#!/usr/bin/env python
# by TR

from sito import Events
import os.path
import unittest

class TestCase(unittest.TestCase):
    def setUp(self):
        self.path = os.path.dirname(__file__)
        self.eventfile = os.path.join(self.path, 'data_temp', 'event_list.txt')
        self.eventfile2 = os.path.join(self.path, 'data', 'NEIC_test.txt')
        try:
            self.events = Events.read(self.eventfile)
        except IOError:
            print ('Load events froms server...')
            self.events = Events.load(min_datetime='2010-1-1', max_datetime='2011-1-1', min_magnitude=7.0)
            self.events.write(self.eventfile)

    def test_events_read_NEIC(self):
        events3 = Events.read(self.eventfile2, regex=Events.regex_NEIC)
        events4 = Events.read(self.eventfile2, regex=Events.regex_NEIC2)
        self.assertEqual(len(events3), len(events4))

    def test_events_IO(self):
        eventfile3 = os.path.join(self.path, 'temp', 'test_events_IO.txt')
        self.events[:3].write(eventfile3)
        events2 = Events.read(eventfile3)
        self.assertEqual(len(events2), 3)

    def test_events_add(self):
        events_add1 = self.events[:2]
        events_add2 = self.events[2:4]
        self.assertEqual(self.events[:4], events_add1 + events_add2)

    def test_events_pick(self):
        events2 = self.events.pick(latitude= -21., longitude= -69., minval=30, maxval=150, after='2010-05-10 12:00:00', bigger=7.5, replace=False)
        self.assertEqual(len(events2), 1)
        #print (self.events)
        #print ('Some picked events:\n%s' % events2)

    def test_events_plot(self):
        #from pylab import figure, show
        #self.events.plot(-22, -70, show=False)
        #figure()
        self.events.plot(-22, -70, lines=(0, 270), bigmap=True, radius='depth', color='datetime', show=False)
        #show()

def suite():
    return unittest.makeSuite(TestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
