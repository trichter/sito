#!/usr/bin/env python
# by TR

from obspy.core import UTCDateTime as UTC
from sito import read, xcorr as xcorr_mod
from sito.data import ParkfieldTest, IPOCTest
from sito.noisexcorr import stack, prepare, noisexcorr, noisexcorrf
import numpy as np
import os
import unittest



class TestCase(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(__file__)
        # Parkfield
        self.pkd = ParkfieldTest(os.path.join(self.path, 'data'),
                                 os.path.join(self.path, 'temp'))
        self.t1_pkd = UTC('1997-09-03')
        self.t2_pkd = UTC('1997-09-03')
        self.ipoc = IPOCTest(os.path.join(self.path, 'data'),
                             os.path.join(self.path, 'temp'))
        self.ipoc2 = IPOCTest(os.path.join(self.path, 'data'),
                             os.path.join(self.path, 'temp/test2'))

        ### Copy test files from location of raw files
        #print os.path.isfile(self.pkd.getDay('PKD', self.t1_pkd) + '.QHD')
        try:
            self.pkd.getRawStream(self.t1_pkd, 'PKD', 'Z', checkfile=True)
        except ValueError:
            print('Copy missing raw files for Parkfield...')
            self.pkd.copyTestDataFiles(self.t1_pkd, self.t2_pkd, component='Z')
        self.t1_ipoc = UTC('2010-05-04')
        self.t1b_ipoc = UTC('2010-05-05')
        self.t2_ipoc = UTC('2010-05-06')
        try:
            self.ipoc.getRawStream(self.t1_ipoc, 'PB01', 'Z', checkfile=True)
        except ValueError:
            print('Copy missing raw files for IPOC...')
            self.ipoc.copyTestDataFiles(self.t1_ipoc, self.t2_ipoc, component='Z')

        #self.t1_ipoc = UTC('2010-05-05')
        #self.t2_ipoc = UTC('2010-05-05')


    # TODO test noisexcorrf vs noisexcorr

    def test_noisexcorr_Parkfield(self):
        correlation = ('PKDZ', 'PKDZ')

        t1 = self.t1_pkd
        t2 = self.t2_pkd
        data = self.pkd

        prepare(data, ('PKD',), t1, t2, filter=(0.1, 1.), component='Z',
                normalize='1bit', param_norm=None, whitening=None)
        noisexcorr(data, [correlation], t1, t2, 200)
        days1 = read(data.getX(('PKDZ', 'PKDZ'), t1) + '.QHD')

        noisexcorr(data, [correlation], t1, t2, 200, period=3600)
        #data.x_day = data.data + '/xcorr/day_test2/%s_day_%d'
        stack(data, [correlation], period=3600)
        days2 = read(data.getX(correlation, t1) + '.QHD')

        stream2 = data.getStream(t1, 'PKD', component='Z')
        norm_factors = xcorr_mod.getNormFactors(stream2[0].data, stream2[0].data)

#        hours1 = read(data.getXHour(correlation[0], correlation[1], t1) + '.QHD') #@UnusedVariable
#        stream1 = data.getRawStream(t1, 'PKD', component='Z')
#        stream1.plotTrace()
#        stream2.plotTrace(scale=0.25)
#        hours1.plot_()
#        days1.plot_()
#        days2.plot_()
#        from pylab import show; show()

#        from sito import ipshell; ipshell()

        # test stacking, 1bit
        np.testing.assert_almost_equal(norm_factors, np.ones(len(norm_factors)), 3)
        np.testing.assert_almost_equal(days1[0].data, days2[0].data, 3)

    def test_noisexcorr_IPOC(self):

        correlation = ('PB01Z', 'PB02Z')

        t1 = self.t1_ipoc
        t2 = self.t2_ipoc
        data = self.ipoc

        prepare(data, ('PB01', 'PB02'), t1, t2, filter=(0.1, 1.), downsample=10, component='Z',
                normalize='1bit', param_norm=None, whitening=True, use_floating_stream=True)
        t2 = t1 = t1 + (t2 - t1) / 2
        noisexcorr(data, [correlation], t1, t2, 200)
        days1 = read(data.getX(correlation, t1) + '.QHD')
        print days1

        noisexcorr(data, [correlation], t1, t2, 200, period=3600)
        #data.x_day = data.data + '/xcorr/day_test2/%s_day_%d'
        stack(data, [correlation], period=3600)
        days2 = read(data.getX(correlation, t1) + '.QHD')

        noisexcorr(data, [correlation[::-1]], t1, t2, 200)
        days3 = read(data.getX(correlation[::-1], t1) + '.QHD')

#        stream2 = data.getStream(t1, 'PB01', component='Z')
#        hours1 = read(data.getXHour(correlation[0], correlation[1], t1) + '.QHD')
#        stream1 = data.getRawStream(t1, 'PB01', component='Z')
#        stream1.plotTrace()
#        stream2.plotTrace(scale=0.25)
#        stream3 = data.getRawStream(t1, 'PB01', component='Z')
#        stream4 = data.getStream(t1, 'PB02', component='Z')
#        stream3.plotTrace()
#        stream4.plotTrace(scale=0.25)
#        hours1.plot_()
#        days = days1 + days2 + days3
#        days.plot_()
#        from pylab import show; show()


        # test stacking, 1bit and if XCORR(PB01,PB02) = XCORR(PB02,PB01)[::-1]
        np.testing.assert_almost_equal(days1[0].data, days2[0].data, 3)
        np.testing.assert_almost_equal(days1[0].data, days3[0].data[::-1], 3)


        # testing noisexcorrf
        t1 = self.t1_ipoc
        t2 = self.t2_ipoc
        data = self.ipoc2
        prepare(data, ('PB01', 'PB02'), t1, t2, filter=(0.1, 1.), downsample=10, component='Z',
                normalize='1bit', param_norm=None, whitening=True, use_floating_stream=True, freq_domain=True)
        t2 = t1 = t1 + (t2 - t1) / 2
        noisexcorrf(data, [correlation], t1, t2, 200)
        days4 = read(data.getX(correlation, t1) + '.QHD')

        # test stacking, 1bit and if XCORR(PB01,PB02) = XCORR(PB02,PB01)[::-1]
#        days = days1 + days2 + days3 + days4
#        days.plot_()
#        from pylab import show; show()
        np.testing.assert_almost_equal(days1[0].data, days4[0].data, 3)


def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
