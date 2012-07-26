#!/usr/bin/env python
# by TR

from numpy.random import seed, random
from scipy.signal import correlate
from sito import read
from sito.xcorr import xcorrf, acorrf, xcorrt, acorrt, xcorr_obspy#, use_fftw3
import numpy as np
import os
import unittest
from matplotlib.mlab import psd

class TestCase(unittest.TestCase):
    def setUp(self):
        self.stream = read(os.path.join(os.path.dirname(__file__), 'data', 'PKD_1997_246.mseed'))
        self.stream_down = self.stream.copy()
        self.stream_down.downsample2(1)
        N = self.N = 1001
        self.data1 = np.sin(np.arange(N) / 100.)
        self.data1_dem = self.data1.copy()
        self.data1_dem[:] = self.data1_dem[:] - np.mean(self.data1_dem)
        self.data2 = np.e ** (-(np.arange(N) - 500) ** 2 / 100.) - np.e ** (-(np.arange(N) - 50) ** 2 / 100.) + 5 * np.e ** (-(np.arange(N) - 950) ** 2 / 100.)
        self.data2_dem = self.data2.copy()
        self.data2_dem[:] = self.data2_dem[:] - np.mean(self.data2_dem)


    def test_xcorr(self):
        cor1 = correlate(self.data1, self.data2, 'same')
        cor2 = xcorrf(self.data1.copy(), self.data2.copy(), self.N // 2, demean=False, normalize=False)
#        from pylab import plot, show, subplot
#        subplot(211)
#        plot(data1)
#        plot(data2)
#        subplot(212)
#        plot(cor1)
#        plot(cor2)
#        show()
        np.testing.assert_array_almost_equal(cor1, cor2)

    def test_acorr(self):
        cor1 = xcorrf(self.data1, self.data1, self.N // 4, oneside=True)
        cor2 = acorrf(self.data1, self.N // 4)
#        from pylab import plot, show
#        plot(cor1)
#        plot(cor2)
#        show()
        np.testing.assert_array_almost_equal(cor1, cor2, 3)
        self.assertTrue((cor2[0] - 1) ** 2 < 1e-5)

    def test_xcorr_realdata(self):
        stream = self.stream.copy()
        data = stream.select(component='Z')[0].data
        shift = 100
        cor1 = xcorrf(data, data[shift:-shift], shift)
        cor2 = xcorrt(data, data[shift:-shift], shift)
#        from pylab import plot, show
#        plot(cor1)
#        plot(cor2)
#        show()
        np.testing.assert_array_almost_equal(cor1, cor2)

    def test_time_norm(self):
        stream = self.stream_down.copy()
        stream = stream.select(component='Z')
        stream.timeNorm('1bit')
        self.assertTrue(np.sum(np.abs(stream[0].data)) == stream[0].stats.npts)

    def test_spectral_whitening(self):
        stream = self.stream.copy()
        stream.downsample2(1)
        stream = stream.select(component='Z')
        stream.filter2(0.005, 0.2)
        nfft = stream[0].stats.npts
        pxx1, freqs1 = psd(stream[0].data, Fs=stream[0].stats.sampling_rate, NFFT=nfft,
                           scale_by_freq=True)
        stream0 = stream.copy()
        stream2 = stream.copy()
        stream.spectralWhitening(0.001)
        stream.filter2(0.005, 0.2)
#        use_fftw3()
#        stream2.spectralWhitening(0.001)
#        use_fftw3(False)
        stream2.filter2(0.005, 0.2)

        pxx2, freqs2 = psd(stream[0].data, Fs=stream[0].stats.sampling_rate, NFFT=nfft,
                           scale_by_freq=True)
        stream0, pxx1, pxx2, freqs1, freqs2
#        from pylab import plot, show, subplot, legend
#        subplot(311)
#        plot(stream0[0].data / max(stream0[0].data) / 20)
#        plot(stream[0].data)
#        plot(stream2[0].data)
#        ax = subplot(312)
##        ax.set_xscale('log')
##        ax.set_yscale('log')
#        plot(freqs1, pxx1)
#        ax2 = subplot(313, sharex=ax)
##        ax.set_xscale('log')
##        ax2.set_yscale('log')
#        plot(freqs2, pxx2)
#        show()

#        np.testing.assert_array_almost_equal(stream[0].data, stream2[0].data)

########## these tests can be remove later

    def test_xcorr_obspy(self):
        N = 1001
        data1 = np.sin(np.arange(N) / 100.)
        data2 = np.e ** (-(np.arange(N) - 500) ** 2 / 100.) - np.e ** (-(np.arange(N) - 50) ** 2 / 100.) + 5 * np.e ** (-(np.arange(N) - 950) ** 2 / 100.)
        data1[:] = data1[:] - np.mean(data1)
        data2[:] = data2[:] - np.mean(data2)
        cor1 = correlate(data1, data2, 'same')
        cor2 = xcorr_obspy(data1, data2, N // 2)
        #cor3 = xcorr_obspy(data1, data1, N)
        cor1 *= max(cor2) / max(cor1)
#        from pylab import plot, show, subplot
#        subplot(211)
#        plot(data1)
#        plot(data2)
#        subplot(212)
#        plot(cor1)
#        plot(cor2)
#        plot(cor3)
#        show()
#        np.testing.assert_array_almost_equal(cor1, cor2,1)





    def test_xcorr_xcorrt(self):
        N = 1001
        data1 = np.sin(np.arange(N // 2 + 1) / 100.)
        data2 = np.e ** (-(np.arange(N) - 500) ** 2 / 100.) - np.e ** (-(np.arange(N) - 50) ** 2 / 100.) + 5 * np.e ** (-(np.arange(N) - 950) ** 2 / 100.)
        cor1 = correlate(data1 - np.mean(data1), data2 - np.mean(data2), 'full')
        cor2 = xcorrt(data2, data1, 750, window=N, ndat2d=N // 2 + 1)[::-1]
        cor3 = xcorrt(data1, data2, 750, window=N, ndat1d=N // 2 + 1) #@UnusedVariable
        cor3b = xcorrt(data1, data2, 750, window=0) #@UnusedVariable
        cor3c = xcorrt(data2, data1, 750, window=0)[::-1] #@UnusedVariable
        cor1 *= max(cor2) / max(cor1)

        cor4 = correlate(data2, data1, 'full')
        cor5 = xcorrt(data2, data1, 750, demean=False)
        cor5b = xcorrt(data2, data1, 750, shift_zero= -100, demean=False)
        cor5c = xcorrt(data2, data1, 750, shift_zero=100, demean=False)
        cor4 *= max(cor5) / max(cor4)

        cor7 = correlate(data1, data2, 'full')
        cor8 = xcorrt(data1, data2, 750, demean=False, normalize=False)
        cor10 = xcorrt(data1, data2, 750, oneside=True, demean=False, normalize=False) #@UnusedVariable

#        from pylab import plot, show, subplot, legend
#        subplot(411)
#        plot(data1)
#        plot(data2)
#        subplot(412)
#        plot(cor1, label='scipy.signal all demeaned and normalized')
#        plot(cor2, label='xcorrt ndat1 > ndat2 ndatxd')
#        plot(cor3, label='xcorrt ndat1 < ndat2 ndatxd')
#        plot(cor3b, label='xcorrt ndat1 > ndat2 window = 0')
#        plot(cor3c, label='xcorrt ndat1 < ndat2 window = 0')
#        legend()
#        subplot(413)
#        plot(cor4, label='scipy.signal all normalized')
#        plot(cor5, label='xcorrt')
#        plot(cor5b, label='xcorrt shifted -100')
#        plot(cor5c, label='xcorrt shifted 100')
#        legend()
#        subplot(414)
#        plot(cor7, label='scipy.signal')
#        plot(cor8, label='xcorrt')
#        plot(cor10, label='xcorrt oneside=True')
#        legend()
#        show()
        np.testing.assert_array_almost_equal(cor1, cor2)
        np.testing.assert_array_almost_equal(cor4, cor5)
        np.testing.assert_array_almost_equal(cor7, cor8)
        np.testing.assert_array_almost_equal(cor5[200:300], cor5b[100:200])
        np.testing.assert_array_almost_equal(cor5[200:300], cor5c[300:400])

    def test_xcorr_acorrt(self):
        data1 = np.sin(np.arange(1001) * 2 * np.pi / 500.)
        cor1 = xcorrt(data1, data1, 1000)
        cor2 = acorrt(data1, 1000, oneside=False, clipdata=False)
#        from pylab import plot, show, subplot, legend
#        subplot(211)
#        plot(data1)
#        subplot(212)
#        plot(cor1, label='xcorrt')
#        plot(cor2, label='acorr clip=False')
#        legend()
#        show()
        np.testing.assert_array_almost_equal(cor1, cor2)

    def test_xcorr_xcorrf(self):
        N = 1001
        data1 = np.sin(np.arange(N // 2 + 1) / 100.)
        data2 = np.e ** (-(np.arange(N) - 500) ** 2 / 100.) - np.e ** (-(np.arange(N) - 50) ** 2 / 100.) + 5 * np.e ** (-(np.arange(N) - 950) ** 2 / 100.)
        cor1 = xcorrt(data1, data2, 750)
        cor2 = xcorrf(data1, data2, 750)
        cor3 = xcorrf(data2, data1, 750)[::-1]
        cor4 = xcorrt(data1, data2, 750, demean=False, normalize=False)
        cor5 = xcorrf(data1, data2, 750, demean=False, normalize=False)
        cor6 = xcorrt(data1, data2, 750, shift_zero=100)
        cor7 = xcorrf(data1, data2, 750, shift_zero=100)
        cor8 = xcorrt(data1, data2, 750, window=200)
        cor9 = xcorrf(data1, data2, 750, window=200)
#        from pylab import plot, show, subplot, legend
#        subplot(411)
#        plot(data1)
#        plot(data2)
#        subplot(412)
#        plot(cor1, label='xcorrt')
#        plot(cor2, label='xcorrf ndat1 > ndat2')
#        plot(cor3, label='xcorrf ndat2 > ndat1')
#        legend()
#        subplot(413)
#        plot(cor4, label='xcorrt all not demeaned not normalized')
#        plot(cor5, label='xcorrf')
#        legend()
#        subplot(414)
#        plot(cor6, label='xcorrt shift')
#        plot(cor7, label='xcorrf shift')
#        plot(cor8, label='xcorrt window')
#        plot(cor9, label='xcorrf window')
#        legend()
#        show()
        np.testing.assert_array_almost_equal(cor1, cor2)
        np.testing.assert_array_almost_equal(cor1, cor3)
        np.testing.assert_array_almost_equal(cor4, cor5, 5)
        np.testing.assert_array_almost_equal(cor6, cor7)
        np.testing.assert_array_almost_equal(cor8, cor9)

    def test_xcorr_acorrf(self):
        data1 = np.sin(np.arange(1001) * 2.1 * np.pi / 500.)
        cor1 = acorrt(data1, 1024, oneside=True, clipdata=False)
        cor2 = acorrf(data1, 1024, oneside=True, clipdata=False)
#        from pylab import plot, show, subplot, legend
#        subplot(211)
#        plot(data1)
#        subplot(212)
#        plot(cor1, label='acorrt')
#        plot(cor2, label='acorrf')
#        legend()
#        show()
        print (np.sum((cor1 - cor2) ** 2) / len(cor1)) ** 0.5
        self.assertTrue((np.sum((cor1 - cor2) ** 2) / len(cor1)) ** 0.5 < 0.1)
        #np.testing.assert_array_almost_equal(cor1, cor2, 4)

    def test_xcorr_all_again(self):
        seed(42)
        data1 = random(1200) * np.sin(np.arange(1200)) * np.arange(1200) - 0.5
        data2 = random(1000) * np.sin(np.arange(1000) * 0.4) * 500 - 0.5
        data1 -= np.mean(data1[100:1100])
        data2 -= np.mean(data2)
        data2_pad = np.hstack((np.zeros(100), data2, np.zeros(100)))
        cor1a = correlate(data1, data2, 'valid')
        cor1a2 = xcorr_obspy(data1, data2_pad, 100)
        cor1b = xcorrt(data1, data2, 100)
        cor1c = xcorrf(data1, data2, 100)
        cor1a *= max(cor1b) / max(cor1a)
        cor1a2 *= max(cor1b) / max(cor1a2)
        cor3 = xcorrf(data2, data1, 100)[::-1] #@UnusedVariable
        cor4 = xcorrt(data1, data2, 100, demean=False, normalize=False)
        cor5 = xcorrf(data1, data2, 100, demean=False, normalize=False)
        val = max (cor4)
        cor4 /= val
        cor5 /= val
        cor6 = xcorrt(data1, data2, 100, shift_zero=100)
        cor7 = xcorrf(data1, data2, 100, shift_zero=100)
        cor8 = xcorrt(data1, data2, 500, window=200)
        cor9 = xcorrf(data1, data2, 500, window=200)

#        from pylab import plot, show, subplot, legend
#        subplot(411)
#        plot(data1)
#        plot(data2)
#        subplot(412)
#        plot(cor1a, label='convolve')
#        plot(cor1b, label='xcorrt')
#        plot(cor1c, label='xcorrf')
#        plot(cor1a2, label='xcorr_obspy')
#        plot(cor6, label='xcorrt shift')
#        plot(cor7, label='xcorrf shift')
#        legend()
#        subplot(413)
#        plot(cor4, label='xcorrt all not demeaned not normalized')
#        plot(cor5, label='xcorrf')
#        legend()
#        subplot(414)
#        plot(cor8, label='xcorrt window')
#        plot(cor9, label='xcorrf window')
#        legend()
#        show()
        np.testing.assert_array_almost_equal(cor1a, cor1b)
        np.testing.assert_array_almost_equal(cor1a, cor1a2)
        np.testing.assert_array_almost_equal(cor1b, cor1c)
        np.testing.assert_array_almost_equal(cor1c[20:30], cor7[120:130])
        np.testing.assert_array_almost_equal(cor4, cor5)
        np.testing.assert_array_almost_equal(cor6, cor7)
        np.testing.assert_array_almost_equal(cor8, cor9)




def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
