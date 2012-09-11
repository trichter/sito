#!/usr/bin/env python
# by TR
#from glob import glob
#from miic.core.stretch_mod import stretch_mat_creation, velocity_change_estimete, time_windows_creation
#from sito import util
#from sito.stream import Stream, read
#from sito.trace import Trace
#from sito.util.imaging import getDataWindow
#from sito.util.main import daygen, streamyeargen2, streamdaygen, timegen, streamtimegen, \
#    yeargen
#from sito.xcorr import xcorrf, timeNorm
#import logging
#import matplotlib.pyplot as plt
#import numpy as np
#import os.path
#import warnings
#from obspy.core.util.decorator import deprecated
#import itertools
import logging
from sito import Stream, read
from obspy.core.event import readEvents
from sito.data import IPOC
from obspy.core.util.attribdict import AttribDict
import numpy as np
import os.path
import glob
from progressbar import ProgressBar
log = logging.getLogger(__name__)
data = IPOC()

def get_event_id(expr):
    if '/' in expr:
        expr = expr.split('/', 1)[1]
    expr = expr.replace('NLL.', '').replace('Origin#', '')
    return expr


def cut_events(in_, out):
    print 'read events...'
    catalog = readEvents(in_, 'QUAKEML')
    print 'cut events...'
    for event in ProgressBar()(catalog):
        oid = get_event_id(event.origins[0].resource_id.getQuakeMLURI())
        ori = event.origins[0]
        etime = ori.time
        #print 'Select', event
        st = Stream()
        for arrival in ori.arrivals:
            arrival.pick_id.convertIDToQuakeMLURI()
            pick = arrival.pick_id.getReferredObject()
            if not pick:
                print 'FAIL to get pick from arrival'
                continue
            ptime = pick.time
            seed_id = pick.waveform_id.getSEEDString()
            try:
                st1 = Stream(data.client.getWaveform(*(seed_id.split('.') + [ptime - 50, ptime + 250])))
            except Exception as ex:
                print '%s for %s' % (ex, seed_id)
                continue
            st1.merge()
            #print 'load %s %s %.1f' % (seed_id, pick.phase_hint, ptime - etime)
            st1[0].stats['event'] = AttribDict(
                                        id=event.resource_id.resource_id,
                                        origin_id=oid,
                                        etime=etime, ptime=ptime,
                                        lat=ori.latitude, lon=ori.longitude,
                                        depth=ori.depth, rms=ori.quality.standard_error,
                                        mag=event.magnitudes[0].mag)
            st += st1
        st.write(out % oid, 'Q')

def acorr(in_, out, tw, filter_):
    print 'acorr events'
    for fname in ProgressBar()(glob.glob(in_)):
        st1 = read(fname)
        st1.setHI('filter', '')
        st1.filter2(*filter_)
        for tr in st1:
            etime = tr.stats.event.etime
            ptime = tr.stats.event.ptime
            stime = tr.stats.starttime
            start_corr = etime + 2.2 * (ptime - etime)
            end_corr = start_corr + tw
            data_before = tr.slice(stime + 10, ptime - 10).data
            rms_before = np.sqrt(np.sum(data_before ** 2) / len(data_before))
            data_after = tr.slice(start_corr, end_corr).data
            rms_after = np.sqrt(np.sum(data_after ** 2) / len(data_after))
            #if rms_after < rms_before * 1.5:
            #    continue
            tr.stats.event.rms_ratio = rms_after / rms_before
            tr.trim(start_corr, end_corr)
            tr.timeNorm('runningmean', 20)
            tr.taper(p=5. / (end_corr - start_corr))
            tr.addZeros(tw)
            tr.acorr(shift=tw)
        ofname = os.path.splitext(os.path.basename(fname))[0]
        st1.write(out % ofname, 'Q')

def stack(in_, out):
    print 'loading files for stacking...'
    st_sum = Stream()
    st1 = read(in_)
    print 'stack traces...'
    for station in ProgressBar()('PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 PB15 HMBCX MNMCX PATCX PSGCX'.split()):
        #st2 = st1.select(station=station, expr='st.event.rms_ratio > 1.5 and st.sampling_rate>40 and -25<st.event.lat<-18 and -67.8<st.event.lon<-66.0')
        st2 = st1.select(station=station, expr='st.event.rms_ratio > 1.5 and st.sampling_rate>40')
        if len(st2) == 0:
            print 'No traces for station %s' % station
            continue
        tr = st2.simpleStack()
        tr.stats.label = station
        st_sum += tr
    st_sum.write(out, 'Q')


if __name__ == '__main__':
    EVENTS = '/home/richter/Data/picks/seiscomp_quakeml/filter/R?_mag>3.xml'
    OUT = '/home/richter/Data/IPOC/local_events/2011_mag>3/%s'
    #cut_events(EVENTS, OUT)
    IN = '/home/richter/Data/IPOC/local_events/2011_mag>3/*.QHD'
    OUT = '/home/richter/Data/IPOC/local_events/acorr/2011_mag>3_2Hz/%s'
    #acorr(IN, OUT, 50, (2, None))
    IN = '/home/richter/Data/IPOC/local_events/acorr/2011_mag>3_2Hz/*.QHD'
    OUT = '/home/richter/Data/IPOC/local_events/acorr/2011_mag>3_2Hz_stacked'
    stack(IN, OUT)
    #plotting:
    #ms_s.plot_(annotate=True, start=45, figtitle='loacal coda acorr >2Hz', plotinfo=['count'])
