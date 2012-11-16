# by TR

from matplotlib.mlab import psd
from numpy.fft.helper import fftfreq
from obspy.core import Trace as ObsPyTrace
from obspy.signal.util import nextpow2
from scipy.fftpack import fft, ifft
import scipy.interpolate
from sito import util
from sito.util import filterResp, fillArray
from sito.xcorr import timeNorm, spectralWhitening
import copy
import logging
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.signal
# from mystream import read

log = logging.getLogger(__name__)

# statshf = 'sampling_rate delta calib npts network location station channel starttime'
# statshf_r = statshf + ' endtime'
# shhf_int = 'SIGN EVENTNO MARK'
# shhf_float = 'DISTANCE AZIMUTH SLOWNESS INCI DEPTH MAGNITUDE LAT LON SIGNOISE PWDW DCVREG DCVINCI'
# shhf_str = 'COMMENT OPINFO FILTER QUALITY BYTEORDER P-ONSET S-ONSET ORIGIN'
# shhf = ' '.join([shhf_int, shhf_float, shhf_str])

class Trace(ObsPyTrace):
    """
    Class derieved from obspy.core.Trace with some additional functionality.
    """

#    @classmethod
#    def read(cls, pathname_or_url, format=None, headonly=False, ** kwargs):
#        """Read first trace of waveform file into an Trace object.
#
#        If there are more than one trace, tries to glue them together.
#        See obspy.core.read.
#        """
#
#        mystream = read(pathname_or_url, format=None, headonly=False,
#                        ** kwargs)
#        mytrace = mystream[0]
#        if len(mystream) > 1:
#            log.warning('File ' + pathname_or_url + ' contains ' +
#                        str(len(mystream)) + ' traces.')
#            mystream.merge()
#            mytrace = mystream[0]
#        if len(mystream) > 1:
#            log.error('File ' + pathname_or_url + ' contains ' +
#                      str(len(mystream)) + ' traces of different id.')
#        return cls(mytrace)

    def __init__(self, trace=None, data=np.array([]), header={}):
        """
        Extend Trace.__init__, Trace object can be an argument.
        """
        if trace == None:
            super(Trace, self).__init__(data=data, header=header)
        else:
            self.stats = trace.stats
            # set data without changing npts in stats object
            super(Trace, self).__setattr__('data', trace.data)
        if not 'is_fft' in self.stats:
            self.stats.is_fft = False

    def write(self, filename, format_, **kwargs):
        """
        Saves current trace into a filename.

        Parameters
        ----------
        filename : string
            Name of the output filename.
        format_ : string
            Name of the output format_.
            See :meth:`~obspy.core.stream.Stream.write` method for all possible
            formats.

        Basic Usage
        -----------
        >>> tr = Trace()
        >>> tr.write("out.mseed", format_="MSEED") # doctest: +SKIP
        """
        # we need to import here in order to prevent a circular import of
        # Stream and Trace classes
        from sito.stream import Stream
        Stream([self]).write(filename, format_, **kwargs)

    def print_(self, mod=0):
        """
        Print some header information.

        :param mod: 0-2
        :return: string with information
        """
        dic = copy.deepcopy(self.stats)
        dic['dur'] = self.stats.endtime - self.stats.starttime
        dic['pon'] = self.stats.ponset - self.stats.starttime
        out = '%(network)s.%(station)s.%(location)s.%(channel)s.%(event_id)s ' \
            '| %(dur).1fs %(sampling_rate).1f Hz '
        out2 = ' dist %(dist).1f, pon %(pon).1fs, ' \
            'azi %(azi).1f, inci %(inci).1f'
        if 'event' in dic.keys():
            dic['event_id'] = dic.event.id
        else:
            dic['event_id'] = ''
        if 'mark' in dic.keys():
            dic['mark'] = dic.mark
        else:
            dic['mark'] = 9

        if mod == 0:
            out += '%(npts)d samples | %(filter)s'
            # dic['st'] = self.stats.starttime.isoformat()
        elif mod == 1:
            out += '| %(filter)s | ' + out2
            # dic['st'] = self.stats.starttime.date.isoformat()
        else:
            out += '| ' + out2 + ', lazi %(lazi).1f, linci %(linci).1f, marked:%(mark)d'  # , razi %(razi).1f, azi2 %(azi2).1f'
            # dic['st'] = self.stats.starttime.date.isoformat()
            # dic['razi'] = (-90-dic['lazi']) % 360
            # dic['azi2'] = (dic['lazi']-180) % 360
        return out % (dic)

    def addZeros(self, secs_before, secs_after=None):
        if secs_after is None:
            secs_after = secs_before
        self.data = np.hstack((np.zeros(secs_before * self.stats.sampling_rate),
                               self.data,
                               np.zeros(secs_after * self.stats.sampling_rate)))
        self.stats.npts = len(self.data)
        self.stats.starttime = self.stats.starttime - secs_before

    def getArgMax(self, ret='index', spline=False, spline_enhance=100, func=None):
        data = self.data
        if func is not None:
            data = func(data.copy())
        if not spline:
            argmax = np.argmax(data)
            datmax = data[argmax]
        else:
            n = len(self)
            x2 = np.linspace(0, n, n * spline_enhance + 1)
            f = scipy.interpolate.InterpolatedUnivariateSpline(np.arange(n), data)
            argmax = x2[np.argmax(f(x2))]
            datmax = f(argmax)
        if ret == 'time':
            argmax = argmax * self.stats.delta
        elif ret == 'utc':
            argmax = self.stats.starttime + argmax * self.stats.delta
        return argmax, datmax


    def fft(self, nfft=None):
        if self.stats.is_fft:
            raise ValueError('Trace already ffted.')
        self.stats.npts_data = self.stats.npts
        if nfft is None:
            nfft = nextpow2(self.stats.npts_data)
        self.stats.nfft = nfft
        self.stats.stdev = (np.sum(self.data ** 2)) ** 0.5
        self.data = fft(self.data, nfft, overwrite_x=False)
        self.stats.is_fft = True
        self.stats.filter += 'FFT'

    def ifft(self):
        if not self.stats.is_fft:
            raise ValueError('Trace is no fft.')
        print self.stats.npts_data
        self.data = np.real(ifft(self.data, self.stats.nfft,
                                 overwrite_x=False)[:self.stats.npts_data])
        self.stats.is_fft = False
        self.stats.filter += 'IFFT'

    def fftfreq(self):
        if not self.stats.is_fft:
            raise ValueError('Trace is no fft.')
        return fftfreq(self.stats.npts, 1. / self.stats.sampling_rate)

    def ffttrim(self, min_freq=0, max_freq=100, dec=1):
        if 'freq_min' in self.stats:
            raise ValueError('you can use ffttrim only once.')
        freqs = self.fftfreq()
        freq_bool = (min_freq <= freqs) * (freqs <= max_freq)
        self.data = self.data[freq_bool][::dec]
        self.stats.sampling_rate /= dec
        self.stats.freq_min = min_freq
        self.stats.freq_max = max_freq

    def demean(self):
        """
        Subtract the mean from the trace.
        """
        self.data = self.data - self.data.mean()
        self.stats.filter += 'Dm'

    def detrend(self):
        """
        Detrend trace.
        """
        self.data = scipy.signal.detrend(self.data)  # , axis=-1, type='linear', bp=0)
        self.stats.filter += 'Dt'

    def integrate(self):
        """
        Integrate trace.
        """
        self.data = np.concatenate((self.data[:1], scipy.integrate.cumtrapz(self.data, dx=1. / self.stats.sampling_rate)))
        self.stats.filter += 'Int'


    def filter2(self, freqmin=None, freqmax=None, corners=2, zerophase=False):
        """
        Wrapper for Trace.filter, make entry in self.stats.filter.
        """
        if self.stats.is_fft:
            self.data *= filterResp(freqmin, freqmax, corners=corners, zerophase=zerophase,
                                    sr=self.stats.sampling_rate, N=self.stats.nfft,
                                     whole=True)[1]
            self.stats.filter += 'BP%4.2f,%4.2f,%d,%d' % (freqmin, freqmax, corners, zerophase)
        else:
            mask = np.ma.getmask(self.data)
            if freqmin and freqmax:
                self.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=corners, zerophase=zerophase)
                self.stats.filter += 'BP%4.2f,%4.2f,%d,%d' % (freqmin, freqmax, corners, zerophase)
            elif not freqmin and freqmax:
                self.filter("lowpass", freq=freqmax, corners=corners, zerophase=zerophase)
                self.stats.filter += 'LP%4.2f,%d,%d' % (freqmax, corners, zerophase)
            elif freqmin and not freqmax:
                self.filter("highpass", freq=freqmin, corners=corners, zerophase=zerophase)
                self.stats.filter += 'HP%4.2f,%d,%d' % (freqmin, corners, zerophase)
            self.data = fillArray(self.data, mask=mask, fill_value=0.)

    def downsample2(self, new_sampling_rate):
        """
        Wrapper for Trace.decimate, make entry in trace.stats.filter.
        """
        if self.stats.sampling_rate >= 2 * new_sampling_rate:
            factor = int(self.stats.sampling_rate / new_sampling_rate)
            self.decimate(factor, strict_length=False, no_filter=True)
            self.stats.filter += 'DS%d' % factor

    def taper2(self, zeros=0, taper=0):
        window = self.stats.endtime - self.stats.starttime
        assert window > 2 * zeros + 2 * taper
        w_z = self.stats.starttime + zeros, self.stats.endtime - zeros
        self.slice(None, w_z[0]).data[:] = 0
        self.slice(w_z[1], None).data[:] = 0
        temp = self.slice(*w_z)
        temp.taper(p=2.*taper / (window - 2 * zeros))
        self.slice(*w_z).data[:] = temp.data

    def acorr(self, shift=None, normalize=True, oneside=False):
        if shift is None:
            shift = len(self.data)
        else:
            shift = int(shift * self.stats.sampling_rate)
        size = max(2 * shift + 1, len(self.data) + shift)
        nfft = nextpow2(size)
        IN1 = fft(self.data, nfft)
        IN1 *= np.conjugate(IN1)
        ret = ifft(IN1).real
        # shift data for time lag 0 to index 'shift'
        ret = np.roll(ret, shift)[:2 * shift + 1]
        # normalize xcorr
        if normalize:
            stdev = (np.sum(self.data ** 2)) ** 0.5
            if stdev == 0:
                log.warning('Data is zero!!')
                ret[:] = 0.
            else:
                ret /= stdev ** 2
        if oneside:
            ret = ret[shift:]
        self.data = ret
        self.stats.npts = len(ret)


    def timeNorm(self, method=None, param=None):
        str_param = str(param)
        if param == None:
            str_param = ''
        self.data = timeNorm(self.data, method=method, param=param)
        self.stats.filter += 'TimeNorm%s%s' % (method, str_param)

    def spectralWhitening(self, smoothi=None, apply_filter=None):
        self.data = spectralWhitening(self.data, sr=self.stats.sampling_rate, smoothi=smoothi, freq_domain=self.stats.is_fft, apply_filter=apply_filter)
        self.stats.filter += 'SpecWhite'

    def moveout(self, model='iasp91', phase='Ps', p0=6.4):
        """
        In-place moveout correction.

        :param p0: reference slowness
        """
        if model == 'iasp91':
            model = util.Iasp91()
        phc = util.phase_dict[phase]
        st = self.stats
        i0 = int((st.ponset - st.starttime) * st.sampling_rate)
        self.data[i0:] = util.mocorr(self.data[i0:],
                                         model.z, model.vp, model.vs,
                                         st.slowness, p0, st.sampling_rate, phc).astype(self.data.dtype)
        # trc.data[i0:] = util.mocorr(trc.data[i0:],
        #                model.z, model.vp, model.vs,
        #                p, p0, fs, phc).astype(trc.data.dtype)
        # i0 = -int(float(trc.time)*trc.fsamp)
        st.filter += 'MC%s,%s' % (phase, p0)

    def zmigr(self, model='iasp91', phase='Ps', zmax=750, dz=0.5):
        """
        Z-migration

        It's not working!
        """
        if model == 'iasp91':
            model = util.Iasp91()
        phc = util.phase_dict[phase]
        st = self.stats
        i0 = int((st.ponset - st.starttime) * st.sampling_rate)
        fs = st.sampling_rate
        self.data = util.zmigr(model.z, model.vp, model.vs,
                                    self.data[i0:],
                                    st.slowness, fs, 0., zmax, dz, phc)
        1 / 0
        print('halo2')
        st.npts = len(self.data)
        st.sampling_rate = 1. / dz
        st.filter += 'ZM%s' % phase

    def norm(self, value=1., fak=None):
        """
        Normalize trace.
        """
        if not fak:
            fak = value / np.max(np.abs(self.data))
        self.data = self.data * fak
        self.stats.filter += 'N%s' % fak

    def signoise(self, winsig, winnoise, relative='ponset'):
        """
        Determine signal noise ratio by dividing the maximum in the two windows.
        """
        st = self.stats
        if relative in ('ponset', 'middle'):
            if relative == 'ponset':
                rel_time = getattr(st, relative)
            else:
                rel_time = st.starttime + (st.endtime - st.starttime) / 2
            winsig0 = rel_time - st.starttime + winsig[0]
            winsig1 = rel_time - st.starttime + winsig[1]
            winnoise0 = rel_time - st.starttime + winnoise[0]
            winnoise1 = rel_time - st.starttime + winnoise[1]
        else:
            winsig0 = winsig[0] - st.starttime
            winsig1 = winsig[1] - st.starttime
            winnoise0 = winnoise[0] - st.starttime
            winnoise1 = winnoise[1] - st.starttime
        t = np.arange(self.stats.npts) * 1. / st.sampling_rate
        datasig = self.data[(t >= winsig0) * (t <= winsig1)]
        datanoise = self.data[(t >= winnoise0) * (t <= winnoise1)]
        # ipshell()
        st.signoise = max(abs(datasig)) / max(abs(datanoise))

    def _window(self, start, end, window='tukey', lenslope=10):
        """
        Window between start and end with args passed to util.getWindow function.
        """
        # if relative != 'ok':
        #    start, end = util.getTimeIntervall(Stream([self]), start, end, relative, ttype='secstart')
        #    start = start[0]
        #    end = end[0]
        t = np.linspace(0, self.stats.endtime - self.stats.starttime, self.stats.npts)
        boolarray = (t >= start) * (t <= end)
        lenwindow = len(boolarray[boolarray])
        alpha = 2. * lenslope / (end - start)  # inserted 1- !!!
        if alpha > 1:
            alpha = 1
        elif alpha < 0:
            alpha = 0
        if window == 'tukey':
            self.data[boolarray] *= util.cosTaper(lenwindow, alpha)
        else:
            self.data[boolarray] *= util.get_window(window, lenwindow)
        self.data[boolarray == False] = 0

    def plotTrace(self, *args, **kwargs):
        from sito import imaging
        return imaging.plotTrace(self, *args, **kwargs)

    @util.add_doc(psd)
    def plotPSD(self, ax=None, x_time=True, scale_by_freq=True, Nfft=256 * 16 * 16,
                 pad_to=None,
                 xscale='log', yscale='log', grid=True,
                 xlabel='time (s)', ylabel=None,
                 figtitle='PSD station component date',
                 title_in_axis=False, smooth=True,
                 just_calculate=False, ** kwargs):
        """
        Plot PSD of first trace.

        Doc matplotlib.mlab.psd:
        """
        if self.stats.is_fft:
            pxx = self.data
            if 'freq_min' in self.stats:
                freqs = np.linspace(self.stats.freq_min, self.stats.freq_max, self.stats.npts)
            else:
                freqs = self.fftfreq()
        else:
                pxx, freqs = psd(self.data, NFFT=Nfft, Fs=self.stats.sampling_rate,
                             scale_by_freq=scale_by_freq, pad_to=pad_to)
        if just_calculate:
            return pxx, freqs
        if x_time:
            pxx = pxx[::-1]
            freqs = 1. / freqs[::-1]
        elif 'time' in xlabel:
            xlabel = 'freq (Hz)'
        if smooth:
            pxx = util.smooth(pxx, smooth)

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()

        # ## print title
        if figtitle is not None:
            figtitle = figtitle.replace('station', self.stats.station)
            figtitle = figtitle.replace('component', self.stats.channel[-1])
            try:
                starttime = self.stats.starttime + 0.5
                figtitle = figtitle.replace('time', '%s' % starttime)
                figtitle = figtitle.replace('date', '%s' % starttime.date)
                figtitle = figtitle.replace('year', '%d' % starttime.year)
                figtitle = figtitle.replace('nfft', '%d' % Nfft)
            except:
                pass
            if not title_in_axis:
                fig.suptitle(figtitle, x=0.5,
                             horizontalalignment='center')
                # fig.text(title, 0., 0.95, horizontalalignment = 'left' )
            else:
                ax.text(0.1, 1, figtitle, verticalalignment='top',
                        transform=ax.transAxes)

        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.grid(grid)
        ax.plot(freqs, pxx, **kwargs)

        return ax

