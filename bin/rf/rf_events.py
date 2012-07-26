#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
"""
script for picking events in 30° - 90° range
INPUT: event list from NEIC (Catalog Search, expanded file format with headers)
OUTPUT: event list
"""
import logging
from sito import Events
logging.basicConfig()

def load_events():
    events = Events.load(min_datetime='2000-01-01', min_magnitude=5, max_results=100000)
    events.write('/home/richter/Data/events/events00-11_mag5.txt')
    #    a.getEvents(self, min_datetime=None, max_datetime=None, min_longitude=None, max_longitude=None, min_latitude=None, max_latitude=None, min_depth=None, max_depth=None, min_magnitude=None, max_magnitude=None, magnitude_type=None, author=None,sort_by='datetime', sort_direction='ASC', max_results=100, format='list', **kwargs)
def combine_events():
    events = Events.read('/home/richter/Data/events/NEIC90-10_mag5.5.txt', regex=Events.regex_NEIC2)
    events.pick(before='2000-01-01')
    events += Events.read('/home/richter/Data/events/events00-11_mag5.txt')
    events.pick(bigger=5.5)
    events.write('/home/richter/Data/events/events90-11_mag5.5.txt')
def compare_events():
    events = Events.read('/home/richter/Data/events/NEIC90-10_mag5.5.txt', regex=Events.regex_NEIC2)
    events.pick(after='2009-12-01', before='2010-05-01', bigger=5.5)
    events2 = Events.read('/home/richter/Data/events/events00-11_mag5.txt')
    events2.pick(after='2009-12-01', before='2010-05-01', bigger=5.5)
    print (events)
    print (events2)
    from IPython import embed
    embed()
def pick_events():
    events = Events.read('/home/richter/Data/events/events00-11_mag5.5.txt')
    events.pick(bigger=5.5) #after='2009-12-01', before='2009-12-10',
    events.write('/home/richter/Data/events/events06-11_mag5.5_IPOC.txt')
def pick_rf_events(coords, after, bigger=5.5, label='temp'):
    print 'read events'
    events = Events.read('/home/richter/Data/events/2012_03_02_NEIC90-12-2_mag5.5.text', regex=Events.regex_NEIC2)
    print 'pick events'
    picked1 = events.pick(latitude=coords.latitude, longitude=coords.longitude, minval=27, maxval=93, after=after, bigger=bigger, replace=False)
    print 'write events'
    picked1.write('/home/richter/Data/events/events_27-93_mag5.5_%s.txt' % label)
#    picked2 = events.pick(latitude=coords.latitude, longitude=coords.longitude, minval=90, maxval=160, after=after, bigger=bigger, replace=False)
#    picked2.write('/home/richter/Data/events/events_90-160_mag5.5_%s.txt' % label)

def plot_events(coords):
    events = Events.read('/home/richter/Data/events/events_30-90_mag5.8_IPOC.txt')
    events.plot(coords.latitude, coords.longitude, bigmap=False, circles=(30, 90, 150))



#load_events()
#combine_events()
#compare_events()
from sito import data
#pick_rf_events(data.Parkfield().stations['PKD'], '1996-09-06', 5.5, 'Parkfield')
pick_rf_events(data.IPOC().stations['PB01'], '2006-01-01', 5.5, 'IPOC_new')
