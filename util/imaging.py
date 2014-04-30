# by TR

from obspy.core import UTCDateTime
try:
    from sito.util import dist2gps
except:
    pass
import logging
import numpy as np
import matplotlib
from matplotlib.colors import Normalize
from matplotlib import cbook
from numpy import ma
log = logging.getLogger(__name__)



def equi(m, lat, lon, radius, indeg=True):
    if indeg:
        radius = radius / 360. * 6371. * 2 * np.pi
    X = []
    Y = []
    for azimuth in range(0, 360):
        glat2, glon2 = dist2gps(radius, azimuth, lat, lon)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
    #~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X, Y = m(X, Y)
    return X, Y

def line(m, lat, lon, azi, start, end, indeg=True):
    if indeg:
        start = start / 360. * 6371. * 2 * np.pi
        end = end / 360. * 6371. * 2 * np.pi
    X = []
    Y = []
    for distance in np.linspace(start, end, 100):
        glat2, glon2 = dist2gps(distance, azi, lat, lon)
        X.append(glon2)
        Y.append(glat2)
    X, Y = m(X, Y)
    return X, Y

def _getUTCListFromSth(stream, time, raisenumber=False):
    N = len(stream)
    if isinstance(time, basestring):
        if 'onset' in time or 'time' in time: #isinstance(relative, basestring):
            ret = stream.getHI(time)
        elif time == 'middle':
            starttime = stream.getHI('starttime')
            endtime = stream.getHI('endtime')
            ret = [starttime[i] + (endtime[i] - starttime[i]) / 2
                   for i in range(N)]
        else:
            raise ValueError('time is string but not expected one.')
    elif isinstance(time, UTCDateTime):
        ret = [time] * N
    elif cbook.iterable(time):
        if np.any([not isinstance(entry, (UTCDateTime, float, int, long)) for entry in time]):
            raise ValueError('time is list, but not of UTCDateTime or float objects.')
        if len(time) != N:
            raise ValueError('time is list, but has not the length of stream.')
        ret = None
    elif isinstance(time, (float, int, long)) and not raisenumber:
        ret = None
    else:
        raise ValueError('time has wrong type.')
    return ret


def getTimeIntervall(stream, start=None, end=None, relative='starttime', ret_rel='utc'):
    """
    Create two lists of UTCDateTimes - start list and end list

    'time' can stand for UTCDateTime, list of UTCDateTimes, header entry out of
    ('ponset', 'sonset', 'startime', 'endtime') or 'middle'
    :param start, end: - None (means start- resp. endtime)
        - time object
        - or seconds relative to param relative
    :param relative: times (if given as seconds=numbers) are taken relative to
        this parameter, is also needed for param ret_rel='relative
        -time object
    :param ret_rel: - 'utc'  output in absolute UTCDateTime
        - 'relative': output in seconds relative to param relative
        - time object: output in seconds relative to time
    :return: start and end list of UTCDateTime or None if stream has length 0
    """
    N = len(stream)
    if N == 0:
        return

    # get list of UTCDateTimes for start_out and end_out
    if start == None:
        start = 'starttime'
    if end == None:
        end = 'endtime'
    start_out = _getUTCListFromSth(stream, start)
    end_out = _getUTCListFromSth(stream, end)

    # get list of UTCDateTimes for relative if needed
    if start_out == None or end_out == None or ret_rel == 'relative':
        relative = _getUTCListFromSth(stream, relative, raisenumber=True)
    # get list of UTCDateTimes for start_out and end_out
    if start_out == None:
        if cbook.iterable(start):
            start_out = [utc + start[i] for i, utc in enumerate(relative)]
        else:
            start_out = [i + start for i in relative]
    if end_out == None:
        if cbook.iterable(start):
            end_out = [utc + end[i] for i, utc in enumerate(relative)]
        else:
            end_out = [i + end for i in relative]

    # convert UTCDateTimes to seconds if ret_rel demands it
    if ret_rel == 'utc':
        return start_out, end_out
    elif ret_rel != 'relative':
        relative = _getUTCListFromSth(stream, ret_rel)
    start_out = [start_out[i] - relative[i] for i in range(N)]
    end_out = [end_out[i] - relative[i] for i in range(N)]
    return start_out, end_out

def getDataWindow(stream, start=None, end=None, relative='starttime'):
    """
    Return array with data in time window (start, end) around relative.

    'time' can stand for UTCDateTime, list of UTCDateTimes, header entry out of
    ('ponset', 'sonset', 'startime', 'endtime') or 'middle'
    :param stream: Stream object with data
    :param start, end: time or float (seconds) relative to param=relative
    :param relative: time, is needed if start or end in seconds (float)
    :return: np.array of shape (N_stream, N_data)
    """
    stream = stream.slice2(start, end, relative=relative)
    N_stream = len(stream)
    if N_stream == 0:
        raise ValueError('Stream has length 0')
    samp = stream.getHI('sampling_rate')
    if min(samp) != max(samp):
        stream.downsample2(min(samp))
        log.warning('Downsampling stream because of differing sampling rate.')
    npts = stream.getHI('npts')
    if min(npts) != max(npts):
        log.warning('Traces in stream have different NPTS. '
                    'Difference: %d samples' % (max(npts) - min(npts)))

    data = np.zeros((N_stream, max(npts)))
    for i, trace in enumerate(stream):
        data[i, :len(trace.data)] = trace.data
    return data



# create colormap Blue -> White -> Red for xcorr plots
cdict = {'red': ((0.0, 0.0, 0.0),
#                 (0.3, 0.5, 0.5),
                  (0.5, 1.0, 1.0),
#                 (0.7, 1.0, 1.0),
                  (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0),
#                 (0.3, 1.0, 1.0),
                  (0.5, 1.0, 1.0),
#                 (0.7, 1.0, 1.0),
                  (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0),
#                 (0.3, 1.0, 1.0),
                  (0.5, 1.0, 1.0),
#                 (0.7, 0.0, 0.0),
                  (1.0, 0.0, 0.0))}
xcorr_cmap = matplotlib.colors.LinearSegmentedColormap('xcorr_cmap', cdict, 256)

class DLogNorm(Normalize):
    """
    Normalize a given positive or negative value to the 0-1 range on a log scale

    negative values are mapped to 0-0.5
    positive values are mapped to 0.5-1

    Derived from:
    matplotlib.colors.LogNorm
    """
    def __init__(self, vmin=None, vmax=None, cmin=1e-5, cmax=1e-5, clip=False):
        """
        If *vmin* or *vmax* is not given, they are taken from the input's
        minimum and maximum value respectively.  If *clip* is *True* and
        the given value falls outside the range, the returned value
        will be 0 or 1, whichever is closer. Returns 0 if::

            vmin==vmax

        Works with scalars or arrays, including masked arrays.  If
        *clip* is *True*, masked values are set to 1; otherwise they
        remain masked.  Clipping silently defeats the purpose of setting
        the over, under, and masked colors in the colormap, so it is
        likely to lead to surprises; therefore the default is
        *clip* = *False*.

        cmin, cmax gives the range of logarithmic plot for positive (cmax)
        and negative (cmin) values. All values  with smaller absolute value
        are mapped to 0.5.
        """
        self.vmin = vmin
        self.vmax = vmax
        self.cmin = cmin
        self.cmax = cmax
        self.clip = clip

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip
        if cbook.iterable(value):
            vtype = 'array'
            val = ma.asarray(value).astype(np.float)
        else:
            vtype = 'scalar'
            val = ma.array([value]).astype(np.float)
        self.autoscale_None(val)
        vmin, vmax = self.vmin, self.vmax
        cmin, cmax = self.cmin * vmin, self.cmax * vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result = 0.0 * val
        else:
            if clip:
                mask = ma.getmask(val)
                val = ma.array(np.clip(val.filled(vmax), vmin, vmax),
                                mask=mask)
            result = 0. * val + 0.5
            result[val > cmax] = (ma.log10(val[val > cmax]) - ma.log10(cmax)) / (np.log10(vmax) - np.log10(cmax)) / 2. + 0.5
            result[val < cmin] = -(ma.log10(-val[val < cmin]) - ma.log10(-cmin)) / (np.log10(-vmin) - np.log10(-cmin)) / 2. + 0.5
        if vtype == 'scalar':
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = self.vmin, self.vmax
        cmin, cmax = self.cmin * vmin, self.cmax * vmax

        if cbook.iterable(value):
            val = np.asarray(value)
            result = 0.0 * val
            result[val > 0.5] = cmax * (vmax / cmax) ** (2. * val[val > 0.5] - 1)
            result[val < 0.5] = cmin * (vmin / cmin) ** (-2. * val[val < 0.5] + 1)
            return result
        else:
            if value == 0.5:
                return 0
            elif value > 0.5:
                return  cmax * (vmax / cmax) ** (2. * value - 1)
            elif value < 0.5:
                return cmin * (vmin / cmin) ** (-2. * value + 1)

    def ticks(self):
        vmin, vmax = self.vmin, self.vmax
        cmin, cmax = self.cmin, self.cmax
        a1 = np.logspace(np.log10(cmax * vmax) + 1, np.log10(vmax), int(-np.log10(cmax)))
        a2 = -np.logspace(np.log10(-cmin * vmin) + 1, np.log10(-vmin), int(-np.log10(cmin)))
        return np.hstack((a1, 0, a2))
