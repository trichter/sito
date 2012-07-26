#!/usr/bin/env python
# by TR

from numpy.random import random, seed
from sito import util
import sito.util.fftw3_be
import numpy as np
import os
import unittest
from scipy import fftpack
from sito.stream import read
from scipy.fftpack import fft, ifft
from obspy.signal.util import nextpow2
from sito.util.seispy import EarthModel1D, Iasp91

class TestCase(unittest.TestCase):
    """
    Test suite for tspy.util.
    """

    def setUp(self):
        # set specific seed value such that random numbers are reproducible
        seed(42)
        self.Z = random(412) - 0.5
        self.N = random(412) - 0.5
        self.E = random(412) - 0.5
        self.stream = read(os.path.join(os.path.dirname(__file__), 'data', 'PKD_1997_246.mseed'))

### test_main
    def test_util_window(self):
        N = 101
#        from pylab import plot, legend, show
#        plot(util.main.getWindow('tukey', N, 0.8), label='0.8')
#        plot(util.main.getWindow('tukey', N, 0.1), label='0.1')
#        plot(util.cosTaper(N, 0.1), label='cosTaper0.1')
#        plot(util.cosTaper(N, 0.8), label='cosTaper0.8')
#        plot(util.cosTaper(N, 1), label='cosTaper1.0')
#        plot(util.main.getWindow('tukey', N, 1), label='1')
#        plot(util.main.getWindow('hanning', N), label='hanning')
#        legend()
#        show()
        np.testing.assert_array_almost_equal(util.main.getWindow('tukey', N, 0.), util.get_window('boxcar', N))
        np.testing.assert_array_almost_equal(util.main.getWindow('tukey', N, 1.), util.get_window('hanning', N))
        self.assertTrue(np.sum((util.main.getWindow('tukey', N, 0.1) -
                                util.cosTaper(N, 0.1)) ** 2) < 0.2) #cosTaper from Obspy is no normal cosine -> high difference

    def test_util_filterResponse(self):
        sr = self.stream[0].stats.sampling_rate

        # test filterResponse vs filter2
        st = self.stream.copy()
        data = st[0].data.copy()
        N = len(data)
        nfft = nextpow2(len(data))
        values = util.main.filterResp(1, 5, corners=2, sr=20, N=nfft, whole=True)[1]
        st.filter2(1, 5, corners=2)
        data2 = ifft(fft(data, nfft) * values, nfft)
        np.testing.assert_array_almost_equal(st[0].data, data2[:N])

        # test stream interface
        st = self.stream.copy()
        st.fft()
        st.filter2(1, 5, corners=2)
        st.ifft()
        np.testing.assert_array_almost_equal(st[0].data, data2[:N])

        # filtering with filterResponse and zerophase=Trueproduces a peak at
        # the end of data. With nfft=N there is no peak anymore, but still the values are not the same
        st = self.stream.copy()
        freqs, values = util.main.filterResp(1, 5, sr=20, N=nfft, corners=2, whole=True, zerophase=True)
        st.filter2(1, 5, corners=2, zerophase=True)
        data2 = ifft(fft(data, nfft) * values, nfft)
        np.testing.assert_array_almost_equal(st[0].data[:-10 * sr], data2[:N - 10 * sr])

        return

        from numpy.fft.helper import ifftshift, fftfreq
        import matplotlib.pyplot as plt

        real_freqs = (ifftshift(freqs) - np.pi) * 0.5 * sr / np.pi
        freqs2 = fftfreq(nfft, 1 / 20.)
        print real_freqs
        print freqs2

        plt.subplot(411)
        plt.plot(real_freqs, np.abs(values))
        ax = plt.subplot(412)
        plt.plot(data, label='data')
        plt.legend()
        plt.subplot(413, sharex=ax)
        plt.plot(st[0].data, alpha=0.5, label='stream.filter2')
        plt.plot(data2[:N], alpha=0.5, label='filterResponse')
        plt.legend()
        plt.show()

### util.helper
    def test_util_pop(self):
        self.assertEqual(util.pop('echo echo test_pop'), 0)

    def test_util_vectorize_args(self):
        """ Tests the decorator with example pspier function. """
        @util.vectorize_args((0, 1, 2, 3, 4))
        def pspier(depth, slat, slon, slowness, azi, model=util.iasp91):
            tp, ts, xp, xs = model.trace(slowness) #@UnusedVariable
            xs = xs[model.layer(depth)]
            return tuple([xs]) + util.dist2gps(xs, azi, slat, slon)
        # test with scalars
        np.testing.assert_array_almost_equal(np.array(pspier(100, -10, -100, 6.4, 0)), np.array([25.037491884001234, -9.7748325744502864, -100.0]))

        # test with 1D array
        slat = np.array([-10, 0])
        slon = np.array([-160, 10])
        slowness = np.array([6.4, 3.2])
        azi = np.array([13., 170.])
        depth = np.array([200, 200])
        rpier, plat, plon = pspier(depth, slat, slon, slowness, azi)
        #print pspier(200, -10, -160, 6.4, 13)
        #print pspier(200, 0, 10, 3.2, 170)
        #print rpier, plat, plon
        np.testing.assert_array_almost_equal(rpier, np.array([ 53.24810772, 25.94769729]))
        np.testing.assert_array_almost_equal(plat, np.array([-9.53338447, -0.22980792]))
        np.testing.assert_array_almost_equal(plon, np.array([-159.89077001, 10.04052156]))

        # test with 2D array
        slat = np.array([[-10, 0, 75], [30, 40, 50]])
        slon = np.array([[-100, 0, 150], [300, 40, 50]])
        azi = np.array([[0, 10, 20], [30, 40, 50]])
        rpier3, plat3, plon3 = pspier(100, slat, slon, 6.4, azi)
        #print rpier3
        #print plat3
        #print plon3
        #print util.pspier(100, -10, -100, 6.4, 0)
        #print util.pspier(100, 0, 0, 6.4, 10)
        #print util.pspier(100, 75, 150, 6.4, 20)
        #print util.pspier(100, 30, 300, 6.4, 30)
        #print util.pspier(100, 40, 40, 6.4, 40)
        #print util.pspier(100, 50, 50, 6.4, 50)
        np.testing.assert_array_almost_equal(rpier3, np.array([[25.037491884001234, 25.037491884001234, 25.037491884001234], [25.037491884001234, 25.037491884001234, 25.037491884001234]]))
        np.testing.assert_array_almost_equal(plat3, np.array([[-9.7748325744502864, 0.22174660919583036, 75.21139224817631], [30.194936597642542, 40.172334289001896, 50.144424254999045]]))
        np.testing.assert_array_almost_equal(plon3, np.array([[-100.0, 0.039100108343866131, 150.30170770482789], [-59.869743130561837, 40.189416675883059, 50.26915378427686]]))

    def test_util_add_doc(self):
        def bla2():
            """2"""
            pass
        @util.add_doc(bla2)
        def bla12():
            """1"""
            pass
        @util.add_doc('4', bla12)
        def bla34():
            """3"""
            pass
        self.assertEqual(bla12.__doc__, '12')
        self.assertEqual(bla34.__doc__, '3412')

### util.seispy
    def test_util_region(self):
        reg = util.feregion(10, 10)
        self.assertEqual(reg, 'Nigeria')

    def  test_util_ttt(self):
        ttt = util.ttt(50, 50, surface=True)
        self.assertEqual(abs(ttt[0].time - 529.101318359375) < 1e-5, True)

    def  test_util_correct_ttt(self):
        self.assertEqual(util.correct_ttt(), 0.)
        husen99 = EarthModel1D(os.path.join(os.path.dirname(__file__), 'data', 'husen99.dat'))
        cor = util.correct_ttt(newmodel=husen99, oldmodel=util.iasp91)
        cor2 = util.correct_ttt('S', 3.2, newmodel=husen99, oldmodel=util.iasp91)
        self.assertAlmostEqual(cor, -1.36098024643)
        self.assertAlmostEqual(cor2, -3.1553046897)

    def test_util_gps2dist(self):
        """ Tests gps2dist against seispy. """
        try:
            import seis.geo
            from obspy.core.util import gps2DistAzimuth
        except ImportError:
            pass
        else:
            dist_deg = sito.util.gps2DistDegree(10, 20, 30, 40)
            dist_km, az1, az2 = gps2DistAzimuth(10, 20, 30, 40) #@UnusedVariable
            dist_deg_seis, az1_seis, az2_seis = seis.geo.delazi(10, 20, 30, 40)
            #print util.gps2dist(10, 20, 30, 40)
            #print seis.geo.delazi(10, 20, 30, 40)
            self.assertEqual(abs(dist_deg - dist_deg_seis) < 1e-5, True)
            self.assertEqual(abs(az1 - az1_seis) < 0.2, True) # one of both routines is not working exactly
            self.assertEqual(abs(az2 - az2_seis) < 0.2, True)
            #obspy.signal.rotate.gps2DistAzimuth | util.gps2dist vs seis.geo.delazi

    def test_util_pspier(self):
        """ Tests pspier and pspier2. """
        slat = np.array([-10, 0])
        slon = np.array([-160, 10])
        slowness = np.array([6.4, 3.2])
        azi = np.array([13., 170.])
        depth = np.array([200, 200])
        rpier, plat, plon = util.pspier(depth, slat, slon, slowness, azi)
        rpier2, plat2, plon2 = util.pspier2(200, slat, slon, slowness, azi)

        np.testing.assert_array_almost_equal(rpier2, np.array([ 53.25370502, 25.95027721]))
        np.testing.assert_array_almost_equal(plat2, np.array([-9.53337395, -0.22982045]))
        np.testing.assert_array_almost_equal(plon2, np.array([-159.890609, 10.04052357]))
        np.testing.assert_array_almost_equal(rpier, rpier2, decimal=2)
        np.testing.assert_array_almost_equal(plat, plat2, decimal=2)
        np.testing.assert_array_almost_equal(plon, plon2, decimal=2)

    def test_util_trace3(self):
        from obspy.taup.taup import getTravelTimes
        iasp91 = Iasp91(5., 5000)
        slat = np.array([-10, 0])
        slon = np.array([-160, 10])
        slowness = np.array([6.4, 3.2])
        azi = np.array([13., 170.])
        depth = np.array([200, 200])
        rpier, plat, plon = iasp91.pspier(depth, slat, slon, slowness, azi, phase='S')
        rpier2, plat2, plon2 = util.pspier(depth, slat, slon, slowness, azi)
        np.testing.assert_array_almost_equal(rpier, rpier2, 0)
        np.testing.assert_array_almost_equal(plat, plat2, 2)
        np.testing.assert_array_almost_equal(plon, plon2, 2)

        tp, rp, phip = iasp91.trace3(6.4, phase='P', till_turn=False)
        taup = getTravelTimes(phip[-1] * 180 / np.pi, 0)
        self.assertTrue(abs(tp[-1] - taup[0]['time']) < 2)

    def test_util_time2depth(self):
        test = np.array((util.time2depth(20, 'P'), util.time2depth(20, 'Ps'), util.time2depth(20, 'Pppp')))
        np.testing.assert_array_almost_equal(test, np.array((168.7292053, 183.00481858, 76.94140119)))
        np.testing.assert_array_almost_equal(util.time2depth([20, 10], 'P'), np.array((168.7292053, 76.94140119)))
        #print 'depth:'
        #print test
        #print util.time2depth(20, 'P*')
        #print util.time2depth(20, 'S*')

    def test_util_depth2time(self):
        test = np.array((util.depth2time(100, 'S'), util.depth2time(100, 'Sppp'), util.depth2time(100, 'Ssss')))
        np.testing.assert_array_almost_equal(test, np.array((23.75012441, 13.84042822, 47.50024882)))
        #print 'time:'
        #print test
        #print util.depth2time(183, 'P*')
        #print util.depth2time(183, 'S*')

    def test_util_interpolate(self):
        t = np.array([1., 2., 3., 5.])
        z = np.array([0, 1., 3., 8.])
        zneu = util.interpolate(2.25, t, z)
        self.assertEqual(zneu, 1.5)
        zneu = util.interpolate(np.array([-0.2, 2., 4.5, 5, 12.]), t, z)
        np.testing.assert_array_almost_equal(zneu, [0, 1, 6.75, 8, 8.])

    def test_util_model(self):
        model = Iasp91()
        self.assertEqual(model.vp[0] / model.vs[0], 5.8 / 3.36)
        #print(model)

    def test_util_mocorr(self):
        model = Iasp91()
        phc = util.phase_dict['Ps']
        #self.Z = np.sin(np.linspace(0, 10*np.pi, 500))
        #import psmout
        #y=psmout.psmout([self.Z], [8.5], 0., 5.0, 0.01, 1) # much slower becasue of io of modell
        data = util.mocorr(self.Z, model.z, model.vp, model.vs, 8.5, 6.4, 100, phc)
        #print(sum(self.Z), sum(model.z), phc, sum(model.vp), sum(model.vs))
        #print(sum(data[:len(data)*4//5]), len(data))
        self.assertEqual(abs(sum(data[:len(data) * 4 // 5]) + 2.33691887) < 1e-5, True) # very, very strange anomalies

### util.ffwt3
    def test_util_ffwt3(self):
        data = np.arange(10)
        nfft = 16
        fft1 = fftpack.fft(data, nfft)
        fft2 = sito.util.fftw3_be.fft(data, nfft).copy()
        ifft1 = fftpack.ifft(fft1)
        ifft2 = sito.util.fftw3_be.ifft(fft2)
#        from pylab import plot, show, subplot, legend
#        subplot(211)
#        plot(data, label='data')
#        plot(ifft1, label='scipy.fftpack')
#        plot(ifft2, label='fftw3_bg')
#        legend()
#        subplot(212)
#        plot(fft1, label='scipy.fftpack')
#        plot(fft2, label='fftw3_bg')
#        legend()
#        show()
        np.testing.assert_array_almost_equal(fft1, fft2)
        np.testing.assert_array_almost_equal(ifft1, ifft2)

# untested util.main: getWindow
# untested util.helper: parameters, myhash, setRootLogger, runBash
# untested util.imaging: getTimeIntervall, getDataWindow
# tests could be enhanced: nearly all

def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
