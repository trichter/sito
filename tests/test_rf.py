#!/usr/bin/env python
# by TR

from numpy.random import random, seed
from sito import rf
from sito.util import cosTaper
import numpy as np
import scipy.linalg
import unittest

class UtilTestCase(unittest.TestCase):

    def setUp(self):
        # set specific seed value such that random numbers are reproducible
        seed(42)
        self.Z = random(412) - 0.5
        self.N = random(412) - 0.5
        self.E = random(412) - 0.5

    def test_rf_toeplitz(self):
        a = np.arange(5) + 5
        r = np.arange(5) + 1
        toep = scipy.linalg.toeplitz(r)
        x = np.dot(scipy.linalg.inv(toep), a)  # compare to scipy.linalg
        x2 = rf.toeplitz(r, a)
        np.testing.assert_array_almost_equal(x, x2)

    def test_rf_polar(self):
        inci, azi, lin = rf.polar(self.E, self.N, self.Z) #@UnusedVariable
        self.assertEqual(abs(azi - 205.718401197) < 1e-5, True)
        self.assertEqual(abs(inci - 42.077846123) < 1e-5, True)
        #print ('ZNE  %s  %s  %s' % (self.Z, self.N, self.E))
        #print ('azi %s  inci %s' % (azi, inci))
        #rms = np.sqrt(np.sum((np.abs(rsp2[Nshift:] - rsp[:-Nshift])) ** 2) /
        #              np.sum((np.abs(rsp[:-Nshift])) ** 2))

    def test_rf_deconvf(self):
        N = len(self.Z)
        data1 = self.N * cosTaper(N, 0.1)
        data2 = self.E * cosTaper(N, 0.1)
        src = np.concatenate((self.Z[:N // 10] * cosTaper(N // 10, 0.1), np.zeros(N - N // 10)))
        rf_data = rf.deconvf([data1, data2], src, 100, water=0.01, gauss=2, tshift=10., pad=0)
        self.assertEqual(np.shape(rf_data), (2, N)) # nonesense test
        #util.deconvfAnalyse(data1, src, 100, water=0.1, gauss=10, tshift=10., pad=0) # try gauss=2, 10, 100
        #import pylab
        #pylab.show



# untested: decovfAnalyse, deconvt

def suite():
    return unittest.makeSuite(UtilTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
