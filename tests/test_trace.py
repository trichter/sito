#!/usr/bin/env python
# by TR

from sito.stream import read
import os.path
import unittest

class TestCase(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(__file__)

    def test_trace(self):
        file_ = os.path.join(self.path, 'data_temp', 'STREAM.QHD')
        try:
            mt = read(file_)[0]
        except IOError:
            raise IOError('First start the test for data.')
        dummy = 'repr(mt) = ' + str(repr(mt)) + '\nstr(mt) = ' + str(mt)
        dummy += '\nmt.stats = ' + str(mt.stats) + '\nmt.data = ' + str(mt.data)
        dummy += '\nmt.print_(mod=0) = ' + mt.print_(mod=0)
        dummy += '\nmt.print_(mod=1) = ' + mt.print_(mod=1)
        #print (dummy)

def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
