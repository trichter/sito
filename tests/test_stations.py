#!/usr/bin/env python
# by TR

import os.path
import sito
import unittest

class TestCase(unittest.TestCase):

    def setUp(self):
        self.path = os.path.dirname(__file__)

    def test_stations(self):
        """ Test Station class. """
        file1 = os.path.join(self.path, 'data', 'stationdata.txt')
        file2 = os.path.join(self.path, 'temp', 'stationdata_write_test.txt')
        st = sito.Stations.read(file1)
        st.pick('PB01 PB02')
        st.write(file2)
        self.assertEqual(st.getNames() in ['PB02 PB01', 'PB01 PB02'], True)
        self.assertEqual(abs(st['PB01'].latitude + 21.043228) < 1e-5, True)
        self.assertEqual(abs(st['PB02'].longitude + 69.89603) < 1e-5, True)
        #print('str(st) = ' + str(st))

def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
