#!/usr/bin/env python
# by TR

from sito import data
import os.path
import unittest

class TestCase(unittest.TestCase):
    def setUp(self):
        print('For tests you need these files:\nstationdata.txt with information for stations PB01, PB02 and PKD and\nhusen99.dat\nPKD_1997_246.mseed\nNEIC_test.txt\n')
        self.path = os.path.dirname(__file__)
        self.eventfile = os.path.join(self.path, 'data_temp', 'event_list.txt')
        self.streamfile = os.path.join(self.path, 'data_temp', 'STREAM')
        #! There has to be at least on test which writes a file from
        #! data.evetPicker to streamfile
        try:
            self.events = data._getEvents(self.eventfile)
        except IOError:
            raise IOError('First start test_events.')
        self.events.pick(latitude= -21.0432, longitude= -69.4874, minval=30, maxval=95, after='2010-1-1', bigger=7.0)

    def test_DataNotImplemented(self):
        self.assertRaises(NotImplementedError, data.Data, './', './')

    def test_data_IPOC(self):
        ipoc = data.IPOC(append='test')
        ipoc.stations.pick('PB01 PB02')
        ipoc.events = self.events
        ms, fail = data.eventPicker(ipoc, window=[-200, 800], write=False)
        if len(fail) == 0:
            ms.write(self.streamfile, 'Q')
            self.assertEqual(len(fail), 0)
#        elif os.path.isfile(self.streamfile + '.QHD'):
#            print (fail[0][-1] + ' -- probably the IPOC archive is not mounted\n')
        else:
            raise Exception(fail[0][-1] + ' -- probably the IPOC archive is not mounted')

#    def test_data_PKD(self):
#        pkd = data.Parkfield(append='test')
#        pkd.stations.pick('PKD')
#        pkd.events = self.events
#        ms2, fail2 = data.eventPicker(pkd, window=[-200, 800], write=False)
#        if len(fail2) > 0:
#            print ms2
#            print fail2
#        self.assertEqual(len(fail2), 0)

def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
