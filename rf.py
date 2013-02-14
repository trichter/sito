#!/usr/bin/env python
# by TR

from numpy import abs, max, maximum, pi
from obspy.signal.util import nextpow2
from scipy.fftpack import fft, ifft
from sito import _toeplitz
from sito.xcorr import xcorrt, acorrt  # ToDo use scipy
import logging
import numpy as np
from sito.util.main import filterResp

log = logging.getLogger(__name__)

def toeplitz(a1, b, a2=None):
    if len(a1.shape) == 1:
        a1 = a1[np.newaxis, :]
    if a2 == None:
        a2 = a1[:, 1:]
    elif len(a2.shape) == 1:
        a2 = a2[np.newaxis, :]
    if len(b.shape) == 1:
        bnew = b[np.newaxis, :]
    ret = _toeplitz.cbto_sl(a1, a2, bnew)  #@UndefinedVariable
    if not a1.dtype == np.complex and not a2.dtype == np.complex:
        ret = np.real(ret)
    if len(b.shape) == 1:
        ret = ret[0, :]
    return ret

def polar(fx, fy, fz):
    """
    Computes the polarization properties of the functions fx, fy, fz
    -> inc,azi,lin

    fx, fy and fz must of course be time aligned, gain corrected and
    have the same type.

    inc is the incidence angle, i.e. the angle which the largest
        eigenvector forms with the the (x,y) plane
    azi is the azimuth of the largest eigenvector, i.e. the angle
        of its projection on the (x,y) plane, counted from the
        direction of fx (0 deg) over the direction of fy (90 deg)
    lin is the rectilinearity
    """
    try:
        from seis.signal import covar3
    except:
        import warnings
        warnings.warn('Could not load seis.signal modul!')
        covmat = np.cov(np.vstack((fx, fy, fz)))
    else:
        covmat = covar3(fx, fy, fz)
    #covmat = np.cov(np.vstack((fx,fy,fz)))
    covmat = covmat.real  # after Bataille & Chiu (1991)
    # Python-only version
    eigval, eigvec = np.linalg.eigh(covmat)
    eigvec = eigvec.transpose()
    eigval = eigval.real
#   # alternatively
#   eigval, eigvec = eig3(covmat)
    def _sortEigenvectors3(eigval, eigvec):
        tmp = [x for x in zip(1 * eigval, 1 * eigvec)]  # 1* to force copy FIXMME
        tmp.sort(reverse=True)
        for i in 0, 1, 2:
            eigval[i], eigvec[i] = tmp[i]
    _sortEigenvectors3(eigval, eigvec)
    # let the largest eigenvector point upwards
    if eigvec[0][2] < 0:
        eigvec = -eigvec
    # compute rectilinearity after Vidale (1986)
    lin = 1. - (eigval[1] + eigval[2]) / eigval[0]
    # "incidence" angle
    inc = np.arctan2(np.sqrt(eigvec[0][0] ** 2 + eigvec[0][1] ** 2), eigvec[0][2]) * 180 / pi
    azi = np.arctan2(eigvec[0][1], eigvec[0][0]) * 180 / pi
    azi = (-90 - azi) % 360  # convert azi from E anticlockwise to azi from N clockwise
    return inc, azi, lin


def gauss(data, sampling_rate, gauss, pad=0, hp=None, return_gauss=False):
    """
    Gauss filter

    gauss       Gauss parameter of Low-pass filter
    pad         multiply number of samples used for fft by 2**pad
    """
    N = len(data)
    nfft = nextpow2(N) * 2 ** pad
    freq = np.fft.fftfreq(nfft, d=1. / sampling_rate)
    gauss = np.exp(maximum(-(0.5 * 2 * pi * freq / gauss) ** 2, -700.))
    if hp:
        hp = filterResp(hp[0], hp[1], corners=hp[2], zerophase=True,
                        sr=sampling_rate, N=nfft, whole=True)[1]
        gauss = hp * gauss
    if return_gauss:
        return gauss
#    import pylab
#    pylab.plot(freq, gauss)
#    pylab.show()
    resp = ifft(gauss * fft(data, nfft), nfft)[:N]
    return resp

def deconvf(rsp_list, src, sampling_rate, water=0.05, gauss=2., tshift=10., pad=0, length=None, normalize=True, normalize_to_src=False, returnall=False):
    """
    Frequency-domain deconvolution using waterlevel method.

    rsp, src    data containing the response and
                source functions, respectively
    water       waterlevel to stabilize the deconvolution
    gauss       Gauss parameter of Low-pass filter
    tshift      shift the resulting function by that amount
    pad         multiply number of samples used for fft by 2**pad
    """
    if length == None:
        length = len(src)
    N = length
    nfft = nextpow2(N) * 2 ** pad
    freq = np.fft.fftfreq(nfft, d=1. / sampling_rate)
    gauss = np.exp(maximum(-(0.5 * 2 * pi * freq / gauss) ** 2, -700.) - 1j * tshift * 2 * pi * freq)

    spec_src = fft(src, nfft)
    spec_src_conj = np.conjugate(spec_src)
    spec_src_water = np.abs(spec_src * spec_src_conj)
    spec_src_water = np.maximum(spec_src_water, max(spec_src_water) * water)

    if normalize_to_src:
        spec_src = gauss * spec_src * spec_src_conj / spec_src_water
        rf_src = ifft(spec_src, nfft)[:N]
        #i1 = int((tshift-1)*sampling_rate)
        #i2 = int((tshift+1)*sampling_rate)
        norm = 1 / max(rf_src)
        rf_src = norm * rf_src

    flag = False
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    rf_list = [ifft(gauss * fft(rsp, nfft) * spec_src_conj / spec_src_water, nfft)[:N] for rsp in rsp_list]
    if normalize:
        if not normalize_to_src:
            norm = 1. / max(rf_list[0])
        for rf in rf_list:
            rf *= norm
    if returnall:
        if not normalize_to_src:
            spec_src = gauss * spec_src * spec_src_conj / spec_src_water
            rf_src = ifft(spec_src, nfft)[:N]
            norm = 1 / max(rf_src)
            rf_src = norm * rf_src
        return rf_list, rf_src, spec_src_conj, spec_src_water, freq, gauss, norm, N, nfft
    elif flag:
        return rf
    else:
        return rf_list

# Gives similar results like a deconvolution with Seismic handler,
# but SH is faster
def deconvt(rsp_list, src, spiking, shift, length=None, normalize=True):
    """
    Time domain deconvolution.

    Calculate Toeplitz auto-correlation matrix of source, invert it, add noise
    and multiply it with cross-correlation vector of response and source.

    In a formula:
    RF = (STS + spiking*I)^-1 * STR
    This function calculates RF.

    N... length
        ( S0   S-1  S-2 ... S-N+1 )
        ( S1   S0   S-1 ... S-N+2 )
    S = ( S2   ...                )
        ( ...                     )
        ( SN-1 ...          S0    )
    R = (R0 R1 ... RN-1)^T
    RF = (RF0 RF1 ... RFN-1)^T
    S... source matrix (shape N*N)
    R... response vector (length N)
    RF... receiver function (deconvolution) vector (length N)
    STS = S^T*S = Toeplitz autocorrelation matrix
    STR = S^T*R = cross-correlation vector
    I... Identity


    :parameters:
    rsp_list    (list of) data containing the response
    src         data of source
    spiking     random noise added to autocorrelation (eg. 1.0, 0.1)
    shift       shift the source by that amount of samples to the left side to get onset in RF at the right time
                (negative -> shift source to the right side)
                shift = (middle of rsp window - middle of src window) + (0 - middle rf window)
    length      number of data points of deconvolution
    normalize   if True normalize all deconvolutions so that the maximum of the
                first deconvolution is 1
    :return:    (list of) deconvolutions (length N)
    """
#    time_rf = winrf[1]-winrf[0] = length / samp
#     shift_zero=int(samp*(
#       (winrsp[1] - winrsp[0] - winsrc[1] + winsrc[0] - time_rf) / 2
#       + winrsp[0] - winsrc[0] - winrf[0]))
# =     shift_zero=int(samp*(
#       (winrsp[1] + winrsp[0] - winsrc[1] - winsrc[0] -winrf[1] - winrf[0]) / 2
#
#

    if length == None:
        length = len(src)
    flag = False
    RF_list = []
    STS = acorrt(src, length, demean=False, clipdata=False)
    STS[0] += spiking
    #print shift, len(src), len(rsp_list), spiking, length
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    for rsp in rsp_list:
        STR = xcorrt(rsp, src, length // 2, shift, demean=False)
        if len(STR) > len(STS):
            STR = np.delete(STR, -1)
        RF = toeplitz(STS, STR)
        RF_list.append(RF)
    if normalize:
        norm = 1 / np.max(np.abs(RF_list[0]))
        for RF in RF_list:
            RF *= norm
    if flag:
        return RF
    else:
        return RF_list

    # comment return commands for testing
    samp = 50.
    from pylab import plot, subplot, show
    subplot(411)
    plot(np.arange(len(src)) / samp, src)
    plot(np.arange(len(rsp)) / samp, rsp)
    subplot(412)
    plot(np.arange(len(STS)) / samp, STS)
    subplot(413)
    plot(np.arange(len(STR)) / samp, STR)
    subplot(414)
    plot(np.arange(len(RF)) / samp, RF_list[0])
    show()


def deconvfAnalyse(rsp, src, sampling_rate, water=0.05, gauss=2., tshift=10., pad=0):
    """
    Frequency-domain deconvolution in pure Python+NumPy. Plote some stuff.

    See deconv.
    """

    #N = len(src)
    #nfft = nextpow2(N) * 2**pad
    #rf_src = fft(src, nfft)
    #rf_src_conj = np.conjugate(rf_src)
    #rf_src=rf_src*rf_src_conj
    #rf_src_water = abs(rf_src)
    #rf_src_water=maximum(rf_src_water, max(rf_src_water)*water)

    #freq= np.fft.fftfreq(nfft, d=1./sampling_rate)
    #gauss = np.exp(maximum(-(0.5*2*pi*freq/gauss)**2,-700.)-1j*tshift*2*pi*freq)
    #rf_src=gauss*rf_src/rf_src_water
    #norm = len(rf_src)/sum(abs(rf_src))
    #rf_src=rf_src*norm

    rf, rf_src, rf_src_conj, rf_src_water, freq, gauss, norm, N, nfft = deconvf(rsp, src, sampling_rate, water, gauss, tshift, pad, normalize_to_src=True, returnall=True)
    rf2 = norm * gauss * fft(rsp, nfft) * rf_src_conj / rf_src_water
    rf_src2 = norm * gauss * fft(src, nfft) * rf_src_conj / rf_src_water

    from pylab import figure
    fig = figure()
    ax1 = fig.add_subplot(311)
    ax1.plot(freq, abs(gauss), label='abs gauss')
    ax1.plot(freq, abs(rf_src_water) / max(abs(rf_src_water)), label='abs fft src water')
    ax1.plot(freq, abs(rf_src2) / max(abs(rf_src2)), label='abs ifft rf_src')
    ax1.plot(freq, abs(rf2) / max(abs(rf2)), label='abs ifft rf')
    ax1.legend()

    #rf_src = ifft(rf_src, nfft)[:N]
    #rf = ifft(rf, nfft)[:N]

    t = np.arange(N) * 1. / sampling_rate
    ax2 = fig.add_subplot(312)
    ax2.plot(t, src, label='src')
    ax2.plot(t, rsp, label='signal')
    ax2.legend()
    ax3 = fig.add_subplot(313)
    ax3.plot(t, rf_src, label='rf_src')
    ax3.plot(t, rf[0], label='rf')
    ax3.legend()
    log.debug(' sum(rsp)/sum(src) in tw = %f |  sum(rf[0])/sum(rf_src) in tw = %f'
              % (sum(rsp) / sum(src), abs(sum(rf[0]) / sum(rf_src))))
    return fig

def mul_synrf(mod, ps=None, **kwargs):
    import sito
    import pylab as plt
    stream = sito.Stream()
    for p in ps:
        stream += synrf(mod, p, **kwargs)
    return stream

def synrf(mod, p=6.4, wave='P', gauss=2., nsamp=1024, fsamp=20., \
          tshft=10., nsv=(3.5, 0.25)):
    """
    Computes a synthetic teleseismic P or SV receiver function.

    Method from seispy.

    Parameters:

    mod             velocity model

    Optional parameters:

    p               horizontal slowness of r.f. is; this value *must*
                    be the angular slowness in seconds/degree and will
                    stored in the extra[10] entry of the trace header
    wave            "P" or "SV".
    gauss           Gauss lowpass parameter used to compute r.f.
    nsamp           number of samples to be used in synthetic computation
    fsamp           sampling frequency in Hz
    tshft           time shift used in the computation of the r.f., i.e.
                    zero delay time corresponds to the sample at the
                    index tshft*fsamp. This value is stored in the 'sec'
                    field of the trace starting time, e.g. a time of
                    0/0/0 0:0:-10 corresponds to tshft of +10 sec.
    nsv             tuple containing the near-surface S velocity and
                    Poisson's ratio
    """
    import seis._rf
    import sito
    from obspy.core import UTCDateTime as UTC

    if type(nsv) is tuple:  nsvs, nspr = nsv
    else:                   nsvs, nspr = nsv, 0.25

    if not wave in ["P", "SV"]:
        raise ValueError, "wave type must be either 'P' or 'SV'"
    if   wave == "P": c = "RFQ"
    else:           c = "RFL"

    # Response of L component, Response of Q component, RF of Q component
    fzz, frr, frf = seis._rf.synrf(mod.z, mod.vp, mod.vs, mod.rh, mod.qp, mod.qs, \
                           p, gauss, nsamp, fsamp, tshft, nsvs, nspr, wave)
    stream = sito.Stream(traces=[sito.Trace(data=fzz),
                                 sito.Trace(data=frr),
                                 sito.Trace(data=frf)])
    for i, tr in enumerate(stream):
        tr.stats.starttime = UTC('2000-01-01') - tshft
        tr.stats.ponset = UTC('2000-01-01')
        tr.stats.sampling_rate = fsamp
        tr.stats.channel = ('LRSP', 'QRSP', c)[i]
        tr.stats.slowness = p
    return stream

def test():
    import pylab as plt
    from sito.stream import read
    from sito.util import cosTaper, get_window
    plt.plot(cosTaper(1000, 0.1))
    plt.plot(get_window(('gaussian', 100), 1000))

    N = 10000
    data = np.concatenate((get_window(('gaussian', 10), N // 4), -2 * get_window(('gaussian', 10), N // 4), get_window(('gaussian', 100), N // 4), np.zeros(N // 4)))
    src = np.concatenate((get_window(('gaussian', 10), N // 10), np.zeros(N * 9 // 10)))

    #deconvfAnalyse(data, src, 100, water=0.1, gauss=10, tshift=10., pad=0) # try gauss=2, 10, 100
    dummy = deconvf(data, src, 100, water=0.01, gauss=2, tshift=10., pad=0)

    ms = read('./tests/data_temp/STREAM.QHD')[0:3]
    deconvfAnalyse(ms[1].data, ms[0].data, 100, water=0.01, gauss=2, tshift=10., pad=0)  # try gauss=2, 10, 100
    plt.show()

if __name__ == '__main__':
    test()

# -*- snip -*-
#def rotate(x, y, z, phi, theta):
#    """Rotation Rz(phi)*Ry(theta).
#
#    Rotate x-axis with angle phi in x-y plane towards y-axis.
#    After that rotate z-axis with angle theta towards x-axis.
#    """
#
#    xn = cosd(theta) * cosd(phi) * x + cosd(theta) * sind(phi) * y - sind(theta) * z
#    yn = -sind(phi) * x + cosd(phi) * y
#    zn = sind(theta) * cosd(phi) * x + sind(theta) * sind(phi) * y + cosd(theta) * z
#
#    return xn, yn, zn

    #defined as the angle measured between the vector pointing from the station
    #to the source and the vector pointing from the station to the north
    #noth n, east e, ba
    #r = e * sin((ba + 180) * 2 * pi / 360) + n * cos((ba + 180) * 2 * pi / 360)
    #t = e * cos((ba + 180) * 2 * pi / 360) - n * sin((ba + 180) * 2 * pi / 360)

#def rotateZNE2LQT(Z, N, E, azi, inci):
#    """Rotate from ZNE to LQT."""
#    L = cosd(inci) * Z - sind(inci) * cosd(azi) * N - \
#        sind(inci) * sind(azi) * E
#    Q = -sind(inci) * Z -cosd(inci) * cosd(azi) * N - \
#        cosd(inci) * sind(azi) * E
#    T = - sind(azi) * N + cosd(azi) * E
#    return L, Q, T
#SH: Q -> -Q , T -> -T
#    L = cosd(inci) * Z - sind(inci) * cosd(azi) * N - \
#        sind(inci) * sind(azi) * E
#    Q = sind(inci) * Z + cosd(inci) * cosd(azi) * N + \
#        cosd(inci) * sind(azi) * E
#    T = sind(azi) * N - cosd(azi) * E


#def rotateLQT2ZNE(L, Q, T, azi, inci):
#    """Rotate from LQT to ZNE."""
#    Z = cosd(inci) * L - sind(inci) * Q
#    N = -sind(inci) * cosd(azi) * L - cosd(inci) * cosd(azi) * Q - \
#        sind(azi) * T
#    E = -sind(inci) * sind(azi) * L - cosd(inci) * sind(azi) * Q + \
#        cosd(azi) * T
#    return Z, N, E
# SH:
#    Z = cosd(inci) * L + sind(inci) * Q
#    N = -sind(inci) * cosd(azi) * L + cosd(inci) * cosd(azi) * Q + \
#        sind(azi) * T
#    E = -sind(inci) * sind(azi) * L + cosd(inci) * sind(azi) * Q - \
#        cosd(azi) * T



#         /  cos<inci>  -cos<azim>*sin<inci>  -sin<azim>*sin<inci>  \
#         |                                                         |
#   R  =  |  sin<inci>   cos<azim>*cos<inci>   sin<azim>*cos<inci>  |
#   <CR> ...
#         |                                                         |
#         \  0           sin<azim>            -cos<azim>
# BackRotation = R^T (R transposed)


#def polar(fx, fy, fz, rltype=1):
    #seis.polar.polarization_py(self[3*i+1].data, self[3*i+2].data, self[3*i].data, rltype)
        # seis.polar.polarization (polarization_py)
        # obspy.signal.eigval
