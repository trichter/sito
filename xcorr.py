#!/usr/bin/env python
# by TR

#import logging
from obspy.signal.util import nextpow2
from scipy.fftpack import fft, ifft, fftshift, ifftshift
from sito._xcorr import xcorr as xcorr_c
from sito.util import isnumber, filterResp, smooth, fillArray
import logging
import numpy as np
import obspy.signal
from obspy.core.util.decorator import deprecated



USE_FFTW3 = False

#fft = None
#ifft = None
#def use_fftw3(val=True):
#    global USE_FFTW3, fft, ifft
#    USE_FFTW3 = val
#    print('not ' * (not val) + 'using FFWT3 library.')
#    if val:
#        from sito.util.fftw3_be import fft, ifft
#    else:
#        from scipy.fftpack import fft, ifft
#use_fftw3(USE_FFTW3)

log = logging.getLogger(__name__)

def timeNorm(data, method=None, param=None, recursive=0):
    """
    Calculates normalized data. See Bensen et al.(2007)

    Method is a string. There are the following methods:

    1bit: reduce data to +1 if >0 and -1 if <0
    clip: clip data to the root mean square (rms)
    eventremoval: automatic event detection and removal - if an value is bigger
        than the threshold, the following values are set to zero.
        param: (threshold, number of samples (eg. 30min) to set to zero)
    stalta: automatic event removing with recursive sta/lta trigger
    runningmean: the data is normalized with the running average
        The width of the normalization window determines how much amplitude
        information is retained (N=1 -> 1bit normalization, N very big ->
        rescaled data). Half of the maximum period of the passband filter_ works
        well.
        param: width of window (should be odd)
    runningmean_over_filtered: the data is normalized with the running average
        over the filtered data.
        A band pass filter_ between 20s and 100s period can remove local
        seismicity.
        param: (width of window in seconds, sampling rate, filter_, freq1, freq2)
            filter_: in ('band', 'low', high')
            if filter_ in ('low', 'high') only on frequency is needed
    waterlevel: any amplitude above the waterlevel (multiple of rms) is
        down-weighted by a factor. This procedure is repeated all of the
        waveform data is under the water-level
        param: (water-level factor, reducing factor)

    """
    mask = np.ma.getmask(data)
    if method == '1bit':
        data = np.sign(data)
    elif method == 'clip':
        std = np.std(data)
        data[data > std] = std
        data[data < -std] = -std
    elif method == 'eventremoval':
        if param == None:
            # remove 30 min (at 10Hz) after events if data is bigger than 2000 
            param = (2000, 30 * 60 * 10)
        clip = np.nonzero(abs(data) >= param[0])[0]
        if len(clip) > 0:
            clip = clip[0]
            index = min(clip + param[1], len(data))
            data[clip:index] = 0
            if index < len(data):
                data[index:] = timeNorm(data[index:], method=method, param=param)
    elif method == 'stalta':
        if param is None:
            # STA: 3s at 100Hz, LTA: 10s, trigger on: 1.2, trigger off:1.0 
            param = (100 * 3, 100 * 10, 1.2, 1.0)
        cft = obspy.signal.trigger.recSTALTA(data, param[0], param[1])
        trg = obspy.signal.trigger.triggerOnset(cft, param[2], param[3])
        for on, off in trg:
            data[on:off] = 0
    elif method == 'runningmean':
        if param == None:
            # smooth over 20s at 10Hz
            param = 10 * 10
        smoothed = smooth(np.abs(data), param)
        data /= smoothed
    elif method == 'runningmean_over_filtered':
        if param is None:
            # smooth over 20s at 10Hz over bandpassed data
            param = (10, 10, 'band', 1 / 50., 1 / 15.)
        sr = param[1]
        over = int(param[0] * sr)
        filter_ = param[2]
        if filter_ == 'band':
            data2 = obspy.signal.bandpass(data, param[3], param[4], sr)
        elif filter_ == 'low':
            data2 = obspy.signal.lowpass(data, param[3], sr)
        elif filter_ == 'high':
            data2 = obspy.signal.highpass(data, param[3], sr)
        else:
            raise ValueError("filter_ should be in ('band', 'high', 'low')")
        data /= smooth(np.abs(data2), over)
    elif method == 'waterlevel':
        if param == None:
            # data above 6*rms is recursively reduced by a factor of 10
            param = (6., 10.)
        waterlevel = param[0] * np.std(data)
        indices = np.abs(data) > waterlevel
        if np.any(indices):
            if param[1] == 0:
                data[indices] = 0
            else:
                data[indices] /= param[1]
            data = timeNorm(data, method=method, param=param, recursive=recursive + 1)
    elif method == 'waterlevel_rm':
        if param == None:
            # running mean over 5s at 10Hz data
            # data above 4*rms is recursively reduced by a factor of 10
            param = (5 * 10, 4., 10.)
        running_mean = smooth(np.abs(data), param[0])
        waterlevel = param[1] * np.std(running_mean)
        indices = (running_mean > waterlevel) + (np.abs(data) > waterlevel)
        if np.any(indices):
            param = list(param)
            frac_zeros = 1. * np.count_nonzero(indices) / len(data)
            if param[2] == 0:
                data[indices] = 0
                param[1] *= (1 + frac_zeros)
            else:
                data[indices] /= param[2]
                param[1] *= (1 + frac_zeros * (1 - 1 / param[2]))
            print recursive, frac_zeros, waterlevel
            data = timeNorm(data, method=method, param=param, recursive=recursive + 1)
    elif method == 'waterlevel_env':
        if param == None:
            # data above 4*rms is recursively reduced by a factor of 10
            param = (4., 10.)
        param = list(param)
        if len(param) == 2:
            param.append(0)
            param.append(0)
        env = obspy.signal.cpxtrace.envelope(data)[1][:len(data)]

        # correct std because of zeros
        waterlevel = param[0] * np.std(env) / (1 - param[2])
#        import pylab as plt
#        from imaging import plotTrace
#        from sito import Trace
#        trace = Trace(data=data)
#        trace2 = Trace(data=env)
#        plotTrace(trace)
#        plotTrace(trace2)
#        plt.figure()
#        plt.plot(data)
#        plt.plot(env)
#        plt.hlines(waterlevel, 0, len(data))
#        plt.show()
        indices = env > waterlevel
        frac_zeros = 1. * np.count_nonzero(indices) / len(data)
        if np.any(indices) and frac_zeros > 0.0005 and param[3] < 20:

            if param[1] == 0:
                data[indices] = 0
                #param[0] *= (1 + frac_zeros)
            else:
                data[indices] /= param[2]
                #param[0] *= (1 + frac_zeros * (1 - 1 / param[1]))
            print param[3], frac_zeros, param[2], waterlevel
            param[2] += frac_zeros
            param[3] += 1
            data = timeNorm(data, method=method, param=param)
    elif method == 'waterlevel_env2':
        if param == None:
            # data above 4*rms is recursively reduced by a factor of 10
            param = (4., 10.)
        N = len(data)
        env = obspy.signal.cpxtrace.envelope(data)[1][:N]
        if mask is not False:
            env[mask] = 0.
        num_stds = 96 # 24*4 =^ every 15min
        if N < 86400: # 24*3600
            num_stds = N // 900
        len_parts = N // num_stds # N//96 = N//24//4 =^ 15min
        len_stds = len_parts // 15 # len_parts//15 =^ 1min
        stds = np.array([np.std(env[i:i + len_stds]) for i in np.arange(num_stds) * len_parts])
        if np.min(stds) == 0:
            stds = stds[stds != 0.]
            num_stds = len(stds)
        stds = np.sort(stds)[num_stds // 15:-num_stds // 15]
        stds = stds[stds < np.min(stds) * 2.]
        waterlevel = param[0] * np.mean(stds)
#        import pylab as plt
#        from imaging import plotTrace
#        from sito import Trace
#        trace = Trace(data=data)
#        trace2 = Trace(data=env)
#        plotTrace(trace)
#        plotTrace(trace2)
#        plt.figure()
#        plt.plot(data)
#        plt.plot(env)
#        plt.hlines(waterlevel, 0, len(data))
#        plt.show()
        indices = env > waterlevel
        #frac_zeros = 1. * np.count_nonzero(indices) / len(data)
        if np.any(indices):
            if param[1] == 0:
                # not setting values to zero but masking them
                # -> they will stay zero after spectral whitening
                # and 1bit normalization
                mask = np.ma.mask_or(mask, indices)
                #data[indices] = 0
            else:
                data[indices] /= param[2]
    elif method is not None:
        raise ValueError('The method passed to timeNorm() is not known.')
    return fillArray(data, mask=mask, fill_value=0.)

def spectralWhitening(data, sr=None, smoothi=None, freq_domain=False, apply_filter=None):
    """
    Apply spectral whitening to data.
    
    sr: sampling rate (only needed for smoothing)
    smoothi: None or int
    Data is divided by its smoothed (Default: None) amplitude spectrum.
    """
    if freq_domain:
        mask = False
        spec = data
    else:
        mask = np.ma.getmask(data)
        N = len(data)
        nfft = nextpow2(N)
        spec = fft(data, nfft)
    #df = sr/N
    spec_ampl = np.sqrt(np.abs(np.multiply(spec, np.conjugate(spec))))
    if isinstance(smoothi, basestring) and isnumber(smoothi) and smoothi > 0:
        smoothi = int(smoothi * N / sr)
        spec /= ifftshift(smooth(fftshift(spec_ampl), smoothi))
    else:
        spec /= spec_ampl
    if apply_filter is not None:
        spec *= filterResp(*apply_filter, sr=sr, N=len(spec), whole=True)[1]
    if freq_domain:
        return spec
    else:
        ret = np.real(ifft(spec, nfft)[:N])
        if USE_FFTW3:
            ret = ret.copy()
        return fillArray(ret, mask=mask, fill_value=0.)

#    from pylab import plot, show, subplot
#    freqs = np.fft.fftfreq(nfft, 1. / sr)
#    ax = subplot(211)
##    ax.set_xscale('log')
##    ax.set_yscale('log')
#    plot(freqs, spec_ampl)
#    plot(freqs, ifftshift(smooth(fftshift(spec_ampl), over)), lw=2)
#    ax2 = subplot(212, sharex=ax)
##    ax2.set_xscale('log')
##    ax2.set_yscale('log')
#    plot(freqs, np.abs(spec * np.conjugate(spec)) ** 0.5)
#    plot(freqs, np.ones(nfft), lw=2)
#    show()
#    return np.abs(ifft(spec, nfft)[:N])


def xcorrf(data1, data2, shift=None, shift_zero=0, oneside=False,
           demean=True, window=0, ndat1d=0, ndat2d=0, N1=None, N2=None,
           normalize=True,
           freq_domain=False, transform_back=True,
           stdev1=None, stdev2=None):
    """
    Cross-correlation of numpy arrays data1 and data2 in frequency domain.
    

    We define cross-corelation as:
    xcorr[i] = sum_j (tr1[i+j-shift_zero] * tr2[j])
    The data is demeaned before cross-correlation and the result normalized
    after cross-correlation.

    data1, data2: data
    shift:    maximum samples to shift
              (window for i in the above formula)
    shift_zero: shift tr1 before cross-correlation by this amount of samples to
              the right (this means correlation function is shifted to the
              right or better: the window of what you get of the function
              is shifted to the left)
    oneside:  if True only the right/positive side of the correlation function
              is returned. Overrides parameter shift_zero.
    demean:   if True demean data beforehand
    normalize: if True normalize correlation function
              (1 means perfect correlation)
    window:   Use only data in this window for demeaning and normalizing
              0: window = min(ndat1, ndat2)
              >0: window = this parameter
    ndat1d, ndat2d: If >0 use different values for the length of the arrays when
              calculating the mean (defaults to window parameter)
    return:   numpy array with correlation function of length 2*shift+1 for
              oneside=False and of length shift+1 for oneside=True
    """
    if freq_domain and not transform_back:
        return data1 * np.conjugate(data2)
    elif freq_domain:
        min_size = max(2 * shift + 1 + abs(shift_zero),
                       (N1 + N2) // 2 + shift + abs(shift_zero))
        if len(data1) < min_size:
            raise ValueError('NFFT was not large enough to cover the desired '
                          'xcorr!\nnfft: %d, required minimum: %d' %
                          (len(data1), min_size))
        ret = (ifft(data1 * np.conjugate(data2))).real
    else:
        complex_result = (data1.dtype == np.complex or
                          data2.dtype == np.complex)
        N1 = len(data1)
        N2 = len(data2)
        #if isinstance(data1[0], np.integer) or isinstance(data2[0], np.integer):
        data1 = data1.astype('float64')
        data2 = data2.astype('float64')
        #if (N1-N2)%2==1:
        #    raise ValueError('(N1-N2)%2 has to be 0')
        if window == 0:
            window = min(N1, N2)
        if ndat1d == 0:
            ndat1d = window
        if ndat2d == 0:
            ndat2d = window
        # determine indices for demeaning and normalization
        ind1 = max(0, (N1 - window) // 2)
        ind2 = min(N1, (N1 + window) // 2)
        ind3 = max(0, (N2 - window) // 2)
        ind4 = min(N2, (N2 + window) // 2)

        # demean and normalize data
        if demean:
            data1 -= np.sum(data1[ind1:ind2]) / ndat1d
            data2 -= np.sum(data2[ind3:ind4]) / ndat2d
        if normalize:
            data1 /= np.max(data1[ind1:ind2])
            data2 /= np.max(data2[ind3:ind4])

        # Always use 2**n-sized FFT, perform xcorr
        size = max(2 * shift + 1 + abs(shift_zero),
                   (N1 + N2) // 2 + shift + abs(shift_zero))
        nfft = nextpow2(size)
        IN1 = fft(data1, nfft)
        if USE_FFTW3:
            IN1 = IN1.copy()
        IN1 *= np.conjugate(fft(data2, nfft))
        ret = ifft(IN1)
        if not USE_FFTW3:
            del IN1
        if not complex_result:
            ret = ret.real
    # shift data for time lag 0 to index 'shift'

    ret = np.roll(ret, -(N1 - N2) // 2 + shift + shift_zero)[:2 * shift + 1]
    # normalize xcorr
    if normalize:
        if not freq_domain:
            stdev1 = (np.sum(data1[ind1:ind2] ** 2)) ** 0.5
            stdev2 = (np.sum(data2[ind3:ind4] ** 2)) ** 0.5
#            stdev1 = (np.sum(data1 ** 2)) ** 0.5
#            stdev2 = (np.sum(data2 ** 2)) ** 0.5
        if stdev1 == 0 or stdev2 == 0:
            log.warning('Data is zero!!')
            ret[:] = 0.
        else:
            ret /= stdev1 * stdev2
    if oneside:
        ret = ret[shift:]
    return np.copy(ret)

def acorrf(data, shift, oneside=True, clipdata=False, ** kwargs):
    """
    Auto-correlation of array data in frequency domain.

    clipdata: if True: data2=data[shift:-shift]
              if False: data2=data
    It calls xcorrf.
    See doc for xcorr:
    """
    if clipdata:
        data2 = data[shift:-shift]
    else:
        data2 = data
    return xcorrf(data, data2, shift, oneside=oneside, ** kwargs)

@deprecated
def xcorr_obspy(data1, data2, shift):
    """Cross correlation using ObsPy"""
    return obspy.signal.xcorr(data1, data2, shift, full_xcorr=True)[2]

@deprecated
def xcorrt(data1, data2, shift, shift_zero=0, oneside=False, demean=True, normalize=True, window=0, ndat1d=0, ndat2d=0):
    """
    Cross-correlation of numpy arrays data1 and data2 in time domain.
    """
    if (len(data1) - len(data2)) % 2 == 1:
        raise ValueError('(N1-N2)%2 has to be 0')
    if oneside:
        ret = xcorr_c(data1, data2, (shift + 1) // 2, -((shift + 1) // 2), window, bool(demean), bool(normalize), ndat1d, ndat2d)
        if len(ret) == shift + 2:
            ret = ret[:-1]
    else:
        ret = xcorr_c(data1, data2, shift, shift_zero, window, bool(demean), bool(normalize), ndat1d, ndat2d)
    return ret

@deprecated
def acorrt(data, shift, oneside=True, clipdata=True, ** kwargs):
    """
    Auto-correlation of array data in time domain.

    clipdata: if True: data2=data[shift:-shift]
              if False: data2=data
    It calls xcorrt with parameter oneside=True.
    See doc for xcorr:
    """
    if clipdata:
        data2 = data[shift:-shift]
    else:
        data2 = data
    ret = xcorrt(data, data2, shift, oneside=True, ** kwargs)
    if not oneside:
        ret = np.hstack((ret[::-1], ret[1:]))
    return ret

def getNormFactors(data1, data2, demean=True, num=24):
    """
    return norm factors of xcorr routine of divided data compared to whole set 
    
    The norm factors of the xcorr routine are calculated for the data divided
    into 'num' parts of the data set (each with equal length) and for the whole
    data set and the quotient is returned.
    Only if these are near 1 stacking the xcorr is technically correct.
    """
    N = len(data1)
    if len(data2) != N:
        raise ValueError('Length of data1 has to be the same as length of data2')
    if isinstance(data1[0], np.integer) or isinstance(data2[0], np.integer):
        data1 = 1. * data1
        data2 = 1. * data2
    if demean:
        data1 -= np.mean(data1)
        data2 -= np.mean(data2)
    fac_whole_time = (np.sum(data1 ** 2) * np.sum(data2 ** 2)) ** 0.5
    fac_period = np.zeros(num)

    for i in range(num):
        ind1 = i * N // num
        ind2 = (i + 1) * N // num
        data1p = data1[ind1:ind2] - np.mean(data1[ind1:ind2])
        data2p = data2[ind1:ind2] - np.mean(data2[ind1:ind2])
        fac_period[i] = (np.sum(data1p ** 2) * np.sum(data2p ** 2)) ** 0.5
    fac_period = fac_period * num / fac_whole_time
    return fac_period

