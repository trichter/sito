# by TR

from helper import add_doc

from obspy.core.util import gps2DistAzimuth
from scipy.signal import get_window
import logging
import numpy as np
import warnings
from scipy.signal import iirfilter, freqz

log = logging.getLogger(__name__)

r_earth = 6371
def sind(x): return np.sin(x / 180. * np.pi)
def cosd(x): return np.cos(x / 180. * np.pi)
def tand(x): return np.tan(x / 180. * np.pi)
def arcsind(x): return np.arcsin(x) / np.pi * 180
def arccosd(x): return np.arccos(x) / np.pi * 180
def arctand(x): return np.arctan(x) / np.pi * 180
def isint(string):
    string = string.strip()
    if string.startswith('-') or string.startswith('+'):
        string = string[1:]
    return string.isdigit()
def isfloat(string):
    string = string.strip()
    if string.startswith('-') or string.startswith('+'):
        string = string[1:]
    return -1 != string.find('.') == string.rfind('.') and all([part.isdigit() for part in string.split('.')])
def isnumber(string):
    return isint(string) or isfloat(string)

def gps2DistDegree(lat1, lon1, lat2, lon2):
    return arccosd(sind(lat1) * sind(lat2) +
                   cosd(lat1) * cosd(lat2) * cosd(lon1 - lon2))

def gps2dist(lat1, lon1, lat2, lon2):
    """
    Return distance in degree, in km, azimuth 1 and azimuth 2.

    Arguments:
    lat1: Latitude of point A in degrees (positive for northern,
        negative for southern hemisphere)
    lon1: Longitude of point A in degrees (positive for eastern,
        negative for western hemisphere)
    lat2: Latitude of point B in degrees (positive for northern,
        negative for southern hemisphere)
    lon2: Longitude of point B in degrees (positive for eastern,
        negative for western hemisphere)
    return: (Great circle distance in deg, in km, azimuth A->B in degrees,
        azimuth B->A in degrees)
    """
    distm, az1, az2 = gps2DistAzimuth(lat1, lon1, lat2, lon2)
    distdeg = \
        arccosd(sind(lat1) * sind(lat2) + cosd(lat1) * cosd(lat2) * cosd(lon1 - lon2))
    return distdeg, distm / 1000., az1, az2

@add_doc(get_window)
def getWindow(window, N, alpha=0.2):
    """
    Return window of length N

    :param window: 'tukey' for tukey window, for other see below
    :param N: length of window
    :param alpha: alpha parameter in case of tukey window.
    0 -> rectangular window
    1 -> cosine taper

    :return: window (np.array)

    Doc of scipy.signal.get_window:
    """
    if window != 'tukey':
        return get_window(window, N)
    else:
        temp = np.ones(N)
        x = np.linspace(-1., 1., N)
        ind1 = (abs(x) > 1 - alpha) * (x < 0)
        ind2 = (abs(x) > 1 - alpha) * (x > 0)
        temp[ind1] = 0.5 * (1 - np.cos(np.pi * (x[ind1] + 1) / alpha))
        temp[ind2] = 0.5 * (1 - np.cos(np.pi * (x[ind2] - 1) / alpha))
        return temp

def smooth(x, window_len, window='flat'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    
    http://www.scipy.org/Cookbook/SignalSmooth
    """

    #window_len = 2 * window_len + 1
    if window_len % 2 == 0:
        window_len += 1
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len < 3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s = np.r_[x[(window_len - 1) // 2 :0:-1], x, x[-1:-(window_len + 1) // 2:-1]]
    if window == 'flat': #moving average
        w = np.ones(window_len, 'd')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w / w.sum(), s, mode='valid')
    return y

#def smooth2(x, t, window_len, window='flat'):
#    ret = x.copy()
#    N1 = int(window_len / (t[1] - t[0]))
#    N2 = int(window_len / (t[-1] - t[-2]))
#    x = np.r_[x[N1 // 2 :0:-1], x, x[-1:-N2 // 2:-1]]
#    for i in np.arange(x):
#        if i != len(x - 1):
#            N = window_len / (t[i + 1] - t[i])
#        ret[i] = np.mean(x[i:i + window_len])
#    return ret

def timegen(t1, t2, dt=24 * 3600, start=None):
    """
    Generator for every (dt)th (Default: 1) day (UTCDateTime) between t1 and t2
    
    if start is not None starts at this day number 
    """
    if dt == 'day':
        dt = 24 * 3600
    elif dt == 'hour':
        dt = 3600
    t = t1
    if start is not None:
        t += start * dt
    while t1 <= t <= t2:
        yield(t)
        t += dt

def streamtimegen(stream, dt=24 * 3600, start=None, shift=None):
    """
    Generator for streams with time length of dt seconds (Default: 1day).
    
    shift < dt -> Running mean
    """
    if len(stream) == 0:
        return
    elif dt < 0:
        yield stream
        return
    if shift is None:
        shift = dt
    for t in timegen(stream[0].stats.starttime + 0.1, stream[-1].stats.endtime - 0.1, dt=shift, start=start):
        t_next = t + dt
        st_sel = stream.select('%r<=st.starttime<%r' % (t, t_next))
        if len(st_sel) > 0:
            yield st_sel

def daygen(t1, t2, dt=1, start=None):
    """
    Generator for every (dt)th (Default: 1) day (UTCDateTime) between t1 and t2
    
    if start is not None starts at this day number 
    """
    t1 = t1.__class__(t1.date)
    t2 = t2.__class__(t2.date)
    t = t1
    if start is not None:
        t += start * 24 * 3600
    while t1 <= t <= t2:
        yield(t)
        t += dt * 24 * 3600

def streamdaygen(stream, dt=1, start=None):
    """
    Generator for day streams with time length of dt (Default: 1) days.
    """
    if len(stream) == 0:
        return
    elif dt == -1:
        yield stream
        return
    for t_day in daygen(stream[0].stats.starttime + 5, stream[-1].stats.endtime - 5, dt=dt, start=start):
        t_next_day = t_day + dt * 24 * 3600
        ms = stream.select('%r<=st.starttime<%r' % (t_day, t_next_day))
        if len(ms) > 0:
            yield ms

def yeargen(t1, t2):
    """
    Generator for all years (UTCDateTime) between t1 and t2 (includes t1, t2) 
    """
    for year in range(t1.year, t2.year + 1):
        yield(t1.__class__(year, 1, 1))

def streamyeargen(stream):
    """
    Generator for year streams. 
    """
    if len(stream) == 0:
        return
    for t_year in yeargen(stream[0].stats.starttime + 5, stream[-1].stats.endtime - 5):
        t_next_year = t_year.__class__(t_year.year + 1, 1, 1)
        ms = stream.select('%r<=st.starttime<%r' % (t_year, t_next_year))
        if len(ms) > 0:
            yield ms

def streamyeargen2(stream):
    """
    Generator for the tuple (year, year stream). 
    """
    for st_year in streamyeargen(stream):
        t_year = st_year[0].stats.starttime + 5
        t_year = t_year.__class__(t_year.year, 1, 1)
        yield((t_year, st_year))

def calculate(data, operation):
    if np.any(np.isnan(data)):
        log.warning('Found NAN values.')
        data = np.nan_to_num(data)
    if operation in ('sum', 'mean', 'stack', 'psd'):
        ret = np.mean(data, axis=0)
    else:
        ret = operation(data)
        #raise ValueError('Unknown operation.')
    return ret

# http://azitech.wordpress.com/2011/03/15/designing-a-butterworth-low-pass-filter-with-scipy/
def filterResp(freqmin, freqmax, corners=2, zerophase=False, sr=None, N=None, whole=False):
    """
    Butterworth-Bandpass Filter.

    Filter frequency data from freqmin to freqmax using
    corners corners.

    :param data: Data to filter, type numpy.ndarray.
    :param freqmin: Pass band low corner frequency.
    :param freqmax: Pass band high corner frequency.
    :param sr: Sampling rate in Hz.
    :param corners: Filter corners. Note: This is twice the value of PITSA's 
        filter sections
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the number of corners but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    
    Derived from obspy.signal.filter
    """
    df = sr
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high > 1:
        high = 1.0
        msg = "Selected high corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    [b, a] = iirfilter(corners, [low, high], btype='band',
                       ftype='butter', output='ba')
    freqs, values = freqz(b, a, N, whole=whole) #@UnusedVariable
    if zerophase:
        values *= np.conjugate(values)
    return freqs, values
#    if zerophase:
#        firstpass = lfilter(b, a, data)
#        return lfilter(b, a, firstpass[::-1])[::-1]
#    else:
#        return lfilter(b, a, data)

####FIX: only for reading old data
# mapping some trace stats to Q- header fields
SH_OPY_EVENT_IDX = {
    'EVENTNO': 'eventno',
    'EVENT_ID': 'id',
    'DEPTH':'depth',
    'MAGNITUDE':'magnitude',
    'LAT':'latitude',
    'LON':'longitude',
    'ORIGIN':'datetime'}


SH_OPY_IDX = {
    #'LENGTH':'npts',
    #x'SIGN':'I011',
    #'EVENTNO':'event_id',
    'NETWORK': 'network',
    'MARK':'mark',
    #'DELTA':'delta',
    #'CALIB':'calib',
    'DISTANCE':'dist',
    'AZIMUTH':'azi',
    'SLOWNESS':'slowness',
    'INCI':'inci',
    'SIGNOISE':'signoise',
#x 'PWDW':'pwdw',  length of p-wave train in sec
    'DCVREG':'lazi',
    'DCVINCI':'linci',
    'COMMENT':'comment',
    #'STATION':'station',
    #x 'OPINFO':'S002',
    'FILTER':'filter',
    #x 'QUALITY':'quality',
    # 'COMP':'C000',
    # 'CHAN1':'C001',
    # 'CHAN2':'C002',
    # 'BYTEORDER':'C003',
    # 'START':'S021',
    'P-ONSET':'ponset',
    'S-ONSET':'sonset',
    'PLAT':'plat',
    'PLON':'plon',
    'RPIER':'rpier',
    'SUM':'sum'
}


