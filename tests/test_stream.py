#!/usr/bin/env python
# by TR

from sito.stream import read
import os.path
import unittest

class TestCase(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(__file__)
        file_ = os.path.join(self.path, 'data_temp', 'STREAM.QHD')
        try:
            self.stream = read(file_)
        except Exception as ex:
            raise IOError('First start the test for data: ' + str(ex))
        self.stream_short = self.stream.copy()
        self.stream.trim2(-50, 200)

    def test_stream(self):
        ms = self.stream_short.copy()
        stationfile = os.path.join(self.path, 'data', 'stationdata.txt')
        gmtfile = os.path.join(self.path, 'temp', 'gmt.txt')
        ms.sort(['event.id', 'station', 'component'])
        dummy1 = ms.getStations(stationfile)
        dummy2 = ms.getEvents()
        #print(dummy1)
        #print(dummy2)
        ms.select(component='Z', logit=True).writeGMTfile(gmtfile, -200, 800)
        #ms.plot2(plotphases='all3')

        #print('repr(ms) = ' + str(repr(ms)))
        #print('str(ms) = ' + str(ms))
        #print('ms.check() = ' + str(ms.check()))

        ms2 = ms.copy()
        ms2.rotateZNE2LQT(-2, 10)
        ms2.trim2(-20, 100)
        #ms2.write('./tests/data_temp/STREAMROT', 'Q')
        #ms2.plot_()
        ms3 = ms2.copy()
        ms3.integrate()
        ms3.filter2(0.03, 0)
        #util.ipshell()
        ms3.receiverf(tshift=20)
        #ms3.moveout()
        ms4 = ms3.select(expr='(st.lazi-st.azi) % 360 < 30 or (st.lazi-st.azi) % 360 > 330', logit=True)
        #ms3.plot_()
        #ms3.select(component='Q').plot_()

        dummy = 'ms4.print_(mod=0) = ' + ms4.print_(mod=0) + '\n'
        dummy += 'ms4.print_(mod=1) = ' + ms4.print_(mod=1) + '\n'
        dummy += 'ms4.print_(mod=2) = ' + ms4.print_(mod=2) + '\n'
        #print (dummy)

        #ms4.write('./tests/data_temp/STREAMDEC', 'Q')
        ms5 = ms4.select(component='Q', logit=True)
        ms5.norm('all')
        ms5.simpleStack(True)
        #ms5.plot_()


        ms5.pspier(50, stationfile)
        #print('you can manipulate ms, ms2, ms3, ms4 now')
        #from sito.util import ipshell
        #ipshell()
    def test_stream_integrate(self):
        ms = self.stream_short.copy()
        ms.integrate()

    def test_stream_window(self):
        ms = self.stream_short.copy()
        #ms.plot_()
        ms.window(0, 20)
        #ms.plot_()
        #from sito import ipshell; ipshell()

def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
