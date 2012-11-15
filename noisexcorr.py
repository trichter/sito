#!/usr/bin/env python
# by TR
from glob import glob
from miic.core.stretch_mod import time_windows_creation, time_stretch_estimate
from sito import util
from sito.stream import Stream, read
from sito.trace import Trace
from sito.util.imaging import getDataWindow
from sito.util.main import daygen, streamyeargen2, streamdaygen, timegen, streamtimegen, \
    yeargen
from sito.xcorr import xcorrf, timeNorm
import logging
import matplotlib.pyplot as plt
import numpy as np
import os.path
import warnings
from obspy.core.util.decorator import deprecated
import itertools
from progressbar import ProgressBar
log = logging.getLogger(__name__)


class FloatingStream(object):
    def __init__(self, data, starttime, station, shift_sec, component, period=24 * 3600, save=0.0, use_get_raw=False):
        """
        Initialize instance.

        :param data: Data object
        :shift_sec: correlation will be calculated by shifting with this maximum amount of seconds
        :component: one component or 'all'
        """
        self.station = station
        self.shift_sec = shift_sec
        self.component = component
        self.period = period
        if use_get_raw:
            self.getRawStream = self.getStream
            self.data_getStream = data.getRawStream
        else:
            self.data_getStream = data.getStream
        #if self.shift_sec != 0:
#        try:
#            self.stream = self.data_getStream(day - shift_sec - reserve, station, component)
#        except ValueError as err:
#            self.stream = Stream()
#            warnings.warn('Day %s  Station %s  Comp %s:  %s' % ((day - shift_sec - reserve).date, self.station, self.component, str(err)))
        self.save = save
        self.stream = Stream()
        self.use_day = True
        try:
            self.getStream(starttime - 48. * 3600)
        except ValueError:
            pass
        try:
            self.getStream(starttime - 24. * 3600)
        except ValueError:
            pass
        self.time = starttime
        self.use_day = False
    def dontGetStream(self):
        self.time += self.period
        self.stream.trim(self.time - self.shift_sec, None)

    def getStream(self, day=None, *args, **kwargs): #@UnusedVariable
        """
        Return a stream for xcorr.

        :param day: UTCDateTime with beginning of the day
        :args, kwargs: not used
        """
        #if self.shift_sec == 0:
        #    return self.data.getStream(day, self.station, component=self.component)
        #else:
        if self.use_day and day:
            self.time = day
        old_time = self.time
        self.time += self.period
#        if len(self.stream) == 0:
#            try:
#                self.stream = self.data_getStream(day, self.station,
#                                                  component=self.component)
#            except ValueError as err:
#                warnings.warn('Day %s  Station %s  Comp %s:  %s' %
#                              (next_day.date, self.station, self.component,
#                               str(err)))
#        elif day - 1 > self.stream[0].stats.endtime:
#            log.warning('Something went wrong getting the stream. Continue...')
#            try:
#                self.stream = self.data_getStream(day, self.station,
#                                                  component=self.component)
#            except ValueError as err:
#                pass


#            st = self.data_getStream(day, self.station,
#                                     component=self.component)
#            if st[0].stats.endtime > self.stream[0].stats.endtime:
#                self.stream += st
#                self.stream.merge()

#        if (len(self.stream) == 0 or
#            self.time + self.shift_sec > self.stream[0].stats.endtime - self.save):

        try:
            self.stream += self.data_getStream(self.time, self.station, component=self.component)
            #print 'to_merge ', self.stream
            self.stream.merge(method=1, interpolation_samples=10)
            if not np.ma.is_masked(self.stream[0].data):
                self.stream[0].data = np.ma.filled(self.stream[0].data, 0.)

        except ValueError as err:
            msg = '%s' % (str(err))
            warnings.warn(msg)
        ret_stream = self.stream.slice(old_time - self.shift_sec, self.time + self.shift_sec)
        self.stream.trim(self.time - self.shift_sec, None)
        #if component == 'all' and len(self.stream) != 3 or component != 'all' and len(self.stream)!=1:
        #    self.stream = Stream()
        if len(ret_stream) != len(self.component):
            raise ValueError('No data for station %s comp %s time %s!' % (
                                self.station, self.component, old_time))
        return ret_stream

def _prepare_stream(stream, output_file,
            freq_domain=True, filter=(None, None), downsample=None, #@ReservedAssignment
            simulate=None,
            eventremoval=None, param_removal=None,
            normalize=None, param_norm=None,
            whitening=None, filter_before_whitening=True,
            use_this_filter_after_whitening=None,
            discard=0.8 * 24 * 3600,
#            edge_smooth=(0, 0),
            trim=(0, 0)):
    for tr in stream:
        tr.data = tr.data.astype('float64')
    # use special parameters for merging, because sometimes there is a
    # gap of one sample in IPOC data
    # method = 1: overlapping data: take the one of second trace
    # fill_value=0 fill gaps with zeros
    stream.merge(method=1, interpolation_samples=10)#, fill_value=0)
    for tr in stream:
        tr.data = util.fillArray(tr.data, fill_value=0.)
        tr.stats.filter = ''
    stream.demean()
    if simulate and tr.stats.station in simulate:
        paz = simulate[tr.stats.station]
        print 'Simulate STS2 for station %s' % tr.stats.station
        stream.simulate(paz_remove=paz['remove'], paz_simulate=paz['simulate'],
                        taper_fraction=0.0003)
    if filter is not None and filter_before_whitening:
        stream.filter2(*filter)
    if downsample:
        stream.downsample2(downsample)
    if eventremoval:
        stream.timeNorm(method=eventremoval, param=param_removal)
    if (freq_domain and discard and
        1. * stream[0].stats.npts / stream[0].stats.sampling_rate < discard):
        return 'Coverage %f station %s day %s'
#    if edge_smooth != (0, 0):
#        for tr in stream:
#            tzeros1 = tr.stats.starttime + edge_smooth[0]
#            tzeros2 = tr.stats.endtime - edge_smooth[0]
#            tsmooth1 = tzeros1 + edge_smooth[1]
#            tsmooth2 = tzeros2 - edge_smooth[1]
#            tr.select(tr.stats.starttime, tzeros1).data[:] = 0
#            tr.select(tzeros2, tr.stats.endtime).data[:] = 0
#            tr_sm = tr.select(tzeros1, tsmooth1)
#            tr_sm.data *= 1 + np.cos(np.linspace(0, np.pi, len(tr_sm)))
#            tr_sm = tr.select(tsmooth2, tzeros2)
#            tr_sm.data *= 1 + np.cos(np.linspace(np.pi, 2 * np.pi, len(tr_sm)))
    if freq_domain and normalize is None and trim == (0, 0):
        stream.fft()
    if whitening:
        if use_this_filter_after_whitening:
            filter = use_this_filter_after_whitening #@ReservedAssignment
        stream.spectralWhitening(smoothi=whitening, apply_filter=filter)
    if normalize:
        stream.timeNorm(method=normalize, param=param_norm)
    if trim != (0, 0):
        for tr in stream:
            tr.trim(tr.stats.starttime + trim[0], tr.stats.endtime - trim[1])
    if freq_domain and (normalize is not None or trim != (0, 0)):
            stream.fft()
    if output_file:
        stream.write(output_file, 'Q')
    else:
        return stream


#def prepare_events(data, stations, events, tw, component='Z', discard=False,
#                    **kwargs):
#    old_year = 0
#    stream = Stream()
#    for event in events:
#        new_year = event.origins[0].time.year
#        if new_year != old_year:
#            new_year = old_year
#            if len(stream) > 0:
#                _prepare_stream(stream, data.getDayEv(event.origins[0].time), discard=discard, **kwargs)
#                stream = Stream()
#        for pick in event.picks:
#            station = pick.WaveformStreamID.station
#            if station not in stations:
#                continue
#            time = pick.time
#            if pick.phase_hint == '??': # TODO ???
#                try:
#                    new_stream = data.getRawStreamFromClient(time + tw[0], time + tw[1], station, component=component)
#                except ValueError:
#                    log.info('Error loading station %s event time %s' % (str(station), time))
#                else:
#                    for tr in new_stream:
#                        tr.stats.event_id = event.resource_id.resource_id
#                    stream += new_stream
#    _prepare_stream(stream, data.getDayEv(event.origins[0].time), discard=discard, **kwargs)

#def eventsxcorrf(data, t1, t2, shift_sec, autocorr=True, pool=None, max_preload=5):
#    """
#    Files have to be in frequency domain!
#    """
#    log.info('Event cross correlation: %s' % util.parameters())
#    xcorr = [] if pool else Stream()
#    for t in yeargen(t1, t2):
#        try:
#            stream = data.getDayEv(t)
#        except ValueError:
#            log.warning('Could not load data for year %s!' % (t))
#            continue
#        while len(stream):
#            event_id = stream[0].stats.event_id
#            event_stream = stream.select("st.event_id == '%s'" % event_id)
#            for tr1, tr2 in (itertools.izip(event_stream, event_stream) if autocorr else
#                             itertools.combinations_with_replacement(event_stream, repeat=2)):
#                if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
#                    raise ValueError('Sampling rate is different.')
#                    # check data                           
#                if tr1.stats.npts != tr2.stats.npts:
#                    log.info('Discard data because of different npts %d vs %d' %
#                             (tr1.stats.npts, tr2.stats.npts))
#                    continue
#                assert tr1.stats.is_fft and tr2.stats.is_fft
#                #log.debug ('Calculating xcorr for %s' % t)
#                correlation = (tr1.stats.station, tr2.stats.station)
#                args = (tr1, tr2, shift_sec, correlation)
#                if pool:
#    #                global counter
#    #                counter += 1
#    #                print 'counter+', counter
#    #                def _callb(res):
#    #                    global counter
#    #                    counter -= 1
#    #                    print 'counter-', counter
#                    raise RuntimeError('POOL not Working')
#                    if len(xcorr) >= max_preload:
#                        xcorr[-max_preload].wait()
#                    xcorr.append(pool.apply_async(_noisexcorrf_traces, args))
#                else:
#                    xcorr.append(_noisexcorrf_traces(*args))
#            for tr in event_stream:
#                stream.remove(tr)
#        if len(xcorr) > 0:
#            data.writeXEv(_get_async_resutls(xcorr) if pool else xcorr,
#                        correlation, t)


@util.add_doc(timeNorm)
def prepare(data, stations, t1, t2, component='all', use_floating_stream=True,
            use_client=True, #arclink_client_for_LVC=None,
            pool=None, max_preload=5, **kwargs):
    """
    Prepare data for cross correlation.
    
    Day files of raw data are loaded (from files data.raw), and the prepared
    data is again written to day files (to files data.raw).

    @param data: sito.data object with attributes (raw, x_prep) set
    @param stations: tuple of stations
    @param t1: UTC object with date of first day
    @param t2: UTC object with date of last day
    @param filter_: tuple of min and max frequencies (or None) for filtering
    @param downsample: downsample to this sampling rate (or None)
    @param component: string out of ('Z', 'N', 'E', 'all')
    @param normalize: method for normalizing the data in time domain (or None)
           see xcorr.timeNorm
    @param param_norm: parameter passed to xcorr.timeNorm()
    @param whitening: bool, apply spectral whitening, default: False
    @type use_floating_stream: bool
    @param use_floating_stream: You should set this to True, if you have day
           files which do not start/end exactly at midnight
           default: False
    @param reserve: parameter passed to FloatingStream()
           should be the time between the end/beginning of your data and
           midnight
           
    Her is the documentation string of xcorr.timeNorm
    """
#            filter=(None, None), downsample=None, #@ReservedAssignment
#            eventremoval=None, param_removal=None,
#            normalize=None, param_norm=None,
#            whitening=None, filter_before_whitening=True, use_this_filter_after_whitening=None,
#            freq_domain=False, discard=0.5):
#    kwargs = dict(filter=filter, downsample=downsample, #@ReservedAssignment
#                  eventremoval=eventremoval, param_removal=param_removal,
#                  normalize=normalize, param_norm=param_norm,
#                  whitening=whitening, filter_before_whitening=filter_before_whitening,
#                  use_this_filter_after_whitening=use_this_filter_after_whitening,
#                  freq_domain=freq_domain, discard=discard)
    log.info('Prepare data for noise correlation: %s' % util.parameters())
    print 'Prepare data for noise correlation...'
    if pool:
        async_results = []
    if use_client:
        kwargs['trim'] = (60, 60)
        for station in ProgressBar()(stations):
            for t_day in daygen(t1, t2):
        #for (station, t_day) in ProgressBar()(itertools.product(stations, daygen(t1, t2))):
                try:
                    stream = data.getRawStreamFromClient(t_day - 60, t_day + 24 * 3600 + 60, station, component=component)
                except ValueError:
                    log.info('Error loading station %s day %s' % (str(station), t_day.date))
                    continue
#                    if station != 'LVC' or not arclink_client_for_LVC:
#                        continue
#                    try:
#                        stream = arclink_client_for_LVC.getWaveform('GE', 'LVC', '10', 'BH' + component,
#                                                                    t_day - 60, t_day + 24 * 3600 + 60)
#                    except:
#                        continue
#                    else:
#                        log.info('Found data on GEOFON.')
                if len(stream) == 0:
                    log.info('Stream length 0 station %s day %s' % (str(station), t_day.date))
                    continue
                ## one day
                def _callback(msg):
                    if msg:
                        log.info(msg % ((stream[0].stats.npts / stream[0].stats.sampling_rate / 24 / 3600),
                                        str(station), t_day.date))
                args = (stream, data.getDay(station, t_day, data.x_prep))
                if pool:
                    if len(async_results) >= max_preload:
                        async_results[-max_preload].wait()
                    async_results.append(
                        pool.apply_async(_prepare_stream, args, kwargs, callback=_callback))
                else:
                    _callback(_prepare_stream(*args, **kwargs))
    else:
        if not use_floating_stream:
            data2 = data
        for station in stations:
            if use_floating_stream:
                data2 = FloatingStream(data, t1, station, 0, component=component, use_get_raw=True)
            for t_day in daygen(t1, t2):
                try:
                    stream = data2.getRawStream(t_day, station, component=component)
                except ValueError:
                    log.info('Error loading station %s day %s' % (str(station), t_day.date))
                    continue
                if len(stream) == 0:
                    log.info('Stream length 0 station %s day %s' % (str(station), t_day.date))
                    continue
                ## one day
                def _callback(msg):
                    if msg:
                        log.info(msg % ((stream[0].stats.npts / stream[0].stats.sampling_rate / 24 / 3600),
                                        str(station), t_day.date))
                args = (stream, data.getDay(station, t_day, data.x_prep))
                if pool:
                    if len(async_results) >= max_preload:
                        async_results[-max_preload].wait()
                    async_results.append(
                        pool.apply_async(_prepare_stream, args, kwargs, callback=_callback))
                else:
                    _callback(_prepare_stream(*args, **kwargs))
    if pool:
        for res in async_results:
            res.wait()

#@deprecated
#def noisexcorr(data, correlations, t1, t2, shift_sec, period=24 * 3600, freq_domain=False, discard=50., test=False):
#    """
#
#    use_method: freq ... frequency domain
#                time ... time domain
#                obspy ... obspy
#    """
#    if freq_domain:
#        return noisexcorrf(data, correlations, t1, t2, shift_sec, period)
#    if period == 'day':
#        period = 24 * 3600
#    elif period == 'hour':
#        period = 3600
#    log.info('Noise cross correlation: %s' % util.parameters())
#    for correlation in correlations:
#        autocorr = correlation[0] == correlation[1]
#        station1 = correlation[0][:-1]
#        station2 = correlation[1][:-1]
#        comp1 = correlation[0][-1]
#        comp2 = correlation[1][-1]
##        data1 = data2 = data
#        data1 = FloatingStream(data, t1, station1, shift_sec, component=comp1, period=period)
#        if not autocorr:
#            data2 = FloatingStream(data, t1, station2, 0, component=comp2, period=period, save=10.)
#        xcorr = Stream()
#        for t in timegen(t1, t2, period):
#            if len(xcorr) > 0 and (t - period).date != t.date and (
#                    (period > 3600 and t.getJulday() == 1) or period <= 3600) and not test:
#                data.writeX(xcorr, correlation, t - period, period=period)
#                xcorr = Stream()
#    # load data
#            try:
#                stream1 = data1.getStream()
#            except ValueError:
#                if not autocorr:
#                    data2.dontGetStream()
#                log.debug('Could not load data for station %s %s' % (correlation[0], t))
#                continue
#            if autocorr:
#                stream2 = stream1.slice(t, t + period)
#            else:
#                try:
#                    stream2 = data2.getStream()
#                except ValueError:
#                    log.debug('Could not load data for station %s %s' % (correlation[1], t))
#                    continue
#    # check data
#            if not (len(stream1) == len(stream2) == 1):
#                log.debug('Stream is not compatible (%d,%d traces)' % (len(stream1), len(stream2)))
#                continue
#            tr1 = stream1[0]
#            tr2 = stream2[0]
#            sr = tr1.stats.sampling_rate
#            # if data for full shift_sec is not availlable slice data further
#            if t - tr1.stats.starttime < shift_sec - 0.1 or tr1.stats.endtime - period - t < shift_sec - 0.1:
#                tr1 = tr1.slice(t, t + period)
#            if sr != tr2.stats.sampling_rate:
#                raise ValueError('Sampling rate is different.')
#            assert not tr1.stats.is_fft and not tr2.stats.is_fft
#                #log.error('Sampling rate differs. This should absolutely not be')
#            if (tr1.stats.npts - tr2.stats.npts < 0 or
#                tr2.stats.npts < discard / 100. * period * sr):
#                log.info('Discard data at time %s because less than %d%% available' % (tr2.stats.starttime + 1, discard))
#                continue
#    # calculate xcorr
#            log.debug('Calculating xcorr for time %s %s' % (t, str(correlation)))
#            cor = xcorrf(tr1.data, tr2.data, int(shift_sec * sr))
#
#            # mirror cor function so that it displays waves from station1 to
#            # station2 for positive lag times
#            tr_cor = Trace(data=cor[::-1], header=tr2.stats.copy())
#            tr_cor.stats.npts = len(cor)
#            tr_cor.stats.sampling_rate = (tr_cor.stats.npts - 1.) / 2 / shift_sec
#            tr_cor.stats.station = '-'.join(correlation)
#            tr_cor.stats.filter += 'Xcorr%s' % shift_sec
#            tr_cor.stats.dist = data.stations.dist(station1, station2, indeg=False)
#            xcorr.append(tr_cor)
#            #util.ipshell()
#        if len(xcorr) > 0 and not test:
#            data.writeX(xcorr, correlation, t, period=period)
#        elif test:
#            return xcorr

def _noisexcorrf_traces(tr1, tr2, shift_sec, correlation):
    sr = tr1.stats.sampling_rate
    cor = xcorrf(tr1.data, tr2.data, int(shift_sec * sr), freq_domain=True,
                 N1=tr1.stats.npts_data, N2=tr2.stats.npts_data,
                 stdev1=tr1.stats.stdev, stdev2=tr2.stats.stdev)

    # mirror cor function so that it displays waves from station1 to
    # station2 for positive lag times
    tr_cor = Trace(data=cor[::-1], header=tr1.stats.copy())
    tr_cor.stats.npts = len(cor)
    tr_cor.stats.sampling_rate = (tr_cor.stats.npts - 1.) / 2 / shift_sec
    tr_cor.stats.station = '-'.join(correlation)
    tr_cor.stats.is_fft = False
    tr_cor.stats.filter += 'IFFTXcorr%s' % shift_sec
    return tr_cor

def _get_async_resutls(async_results):
    xcorr = Stream()
    for res in async_results:
        xcorr.append(res.get())
    return xcorr

#counter = 0
def noisexcorrf(data, correlations, t1, t2, shift_sec, period=24 * 3600,
                 pool=None, max_preload=5):
    """
    Day or period files have to be in frequency domain!
    """
    if period == 'day':
        period = 24 * 3600
    elif period == 'hour':
        period = 3600
    log.info('Noise cross correlation: %s' % util.parameters())
    print 'Noise cross correlation...'
    if period != 24 * 3600:
        raise ValueError('function at the moment only '
                         'working with period=24*3600.')
    for correlation in ProgressBar()(correlations):
        autocorr = correlation[0] == correlation[1]
        station1 = correlation[0][:-1]
        station2 = correlation[1][:-1]
        comp1 = correlation[0][-1]
        comp2 = correlation[1][-1]
#        data1 = FloatingStream(data, t1, station1, shift_sec, component=comp1, period=period)
#        if not autocorr:
#            data2 = FloatingStream(data, t1, station2, 0, component=comp2, period=period)
        xcorr = [] if pool else Stream()
        for t in timegen(t1, t2, period):
            if len(xcorr) > 0 and (t - period).date != t.date and (
                    (period > 3600 and t.julday == 1) or period <= 3600):

                data.writeX(_get_async_resutls(xcorr) if pool else xcorr,
                            correlation, t - period, period=period)
                xcorr = [] if pool else Stream()
            if (not (data.getStream(t, station1, component=comp1, check=True) and
                     (autocorr or data.getStream(t, station2, component=comp2, check=True)))):
                log.debug('No data for %s %s' % (str(correlation), t))
                continue
            try:
                stream1 = data.getStream(t, station1, component=comp1)
                if autocorr:
                    stream2 = stream1
                else:
                    stream2 = data.getStream(t, station2, component=comp2)
            except ValueError:
                log.warning('Could not load data for %s %s' % (str(correlation), t))
                continue
            if not (len(stream1) == len(stream2) == 1):
                log.debug('Stream is not compatible (%d,%d traces)' % (len(stream1), len(stream2)))
                continue
            tr1 = stream1[0]
            tr2 = stream2[0]
            if tr1.stats.sampling_rate != tr2.stats.sampling_rate:
                raise ValueError('Sampling rate is different.')
            # check data                           
            if tr1.stats.npts != tr2.stats.npts:
                log.info('Discard data because of different npts %d vs %d' % (tr1.stats.npts, tr2.stats.npts))
                continue
            assert tr1.stats.is_fft and tr2.stats.is_fft
            log.debug ('Calculating xcorr for %s' % t)
            args = (tr1, tr2, shift_sec, correlation)
            if pool:
#                global counter
#                counter += 1
#                print 'counter+', counter
#                def _callb(res):
#                    global counter
#                    counter -= 1
#                    print 'counter-', counter

                if len(xcorr) >= max_preload:
                    xcorr[-max_preload].wait()
                xcorr.append(pool.apply_async(_noisexcorrf_traces, args))
            else:
                xcorr.append(_noisexcorrf_traces(*args))
        if len(xcorr) > 0:
            data.writeX(_get_async_resutls(xcorr) if pool else xcorr,
                        correlation, t, period=period)

#def xcorr_day(data, correlations, t1, t2, shift_sec, use_floating_stream=True, discard=50, use_method='freq'):
#    """
#
#    use_method: freq ... frequency domain
#                time ... time domain
#                obspy ... obspy
#    """
#    log.info('Cross correlate day wise: %s' % util.parameters())
#    #t1 = t1.__class__(t1.date)
#    #t2 = t2.__class__(t2.date)
#    if not use_floating_stream:
#        data1 = data2 = data
#    for correlation in correlations:
#        autocorr = correlation[0] == correlation[1]
#        station1 = correlation[0][:-1]
#        station2 = correlation[1][:-1]
#        comp1 = correlation[0][-1]
#        comp2 = correlation[1][-1]
#        if use_floating_stream:
#            data1 = FloatingStream(data, t1, station1, shift_sec, component=comp1)
#            if not autocorr:
#                data2 = FloatingStream(data, t1, station2, 0, component=comp2)
#        #t_day = t1 - 24 * 3600
#        xyear = Stream()
#        #while t_day < t2:
#        #    t_day += 24 * 3600
#        for t_day in daygen(t1, t2):
#            if t_day.getJulday() == 1 and len(xyear) > 0:
#                data.writeXDay(xyear, correlation, t_day - 24 * 3600)
#                xyear = Stream()
#            try:
#                stream1 = data1.getStream(t_day, station1, component=comp1)
#                if autocorr:
#                    stream2 = stream1.slice(t_day, t_day + 24 * 3600)
#                else:
#                    stream2 = data2.getStream(t_day, station2, component=comp2)
#            except ValueError:
#                log.debug('Could not load data for %s %s' % (str(correlation), t_day.date))
#                continue
##            if stream1 == None or stream2 == None:
##                log.debug('Could not load data for %s %s' % (str(correlation), t_day.date))
##                continue
#            if not (len(stream1) == len(stream2) == 1):
#                log.debug('Stream is not compatible (%d,%d traces)' % (len(stream1), len(stream2)))
#                continue
#            tr1 = stream1[0]
#            tr2 = stream2[0]
#            sr = tr1.stats.sampling_rate
#            # if data for full shift_sec is not availlable slice data further
#            if t_day - tr1.stats.starttime < shift_sec or tr1.stats.endtime - 24 * 3600 - t_day < shift_sec:
#                tr1 = tr1.slice(t_day, t_day + 24 * 3600)
#            if sr != tr2.stats.sampling_rate:
#                raise ValueError('Sampling rate is different.')
#                #log.error('Sampling rate differs. This should absolutely not be')
#            # check data                           
#            if tr1.stats.npts - tr2.stats.npts < 0 or \
#                tr2.stats.npts < discard / 100. * sr * 24 * 3600:
#                log.info('Discard day %s because less than %d%% available' % ((tr2.stats.starttime + 1).date, discard))
#                continue
#
#            log.debug ('Calculating xcorr for day %s' % t_day.date)
#            if use_method == 'freq':
#                cor = xcorrf(tr1.data, tr2.data, int(shift_sec * sr))
#            elif use_method == 'time':
#                cor = xcorrf(tr1.data, tr2.data, int(shift_sec * sr))
#            elif use_method == 'obspy':
#                tr1 = tr1.slice(t_day, t_day + 24 * 3600)
#                cor = xcorr_obspy(tr1.data, tr2.data, int(shift_sec * sr), True)[2]
#            else:
#                raise ValueError('Choose valid method!')
#
#            # mirror cor function so that it displays waves from station1 to
#            # station2 for positive lag times
#            tr_cor = Trace(data=cor[::-1], header=tr1.stats.copy())
#            tr_cor.stats.npts = len(cor)
#            tr_cor.stats.sampling_rate = (tr_cor.stats.npts - 1.) / 2 / shift_sec
#            tr_cor.stats.station = '-'.join(correlation)
#            tr_cor.stats.filter += 'Xcorr%s' % shift_sec
#            xyear.append(tr_cor)
#            #util.ipshell()
#        if len(xyear) > 0:
#            data.writeXDay(xyear, correlation, t_day)

def getFilters(*args, **kwargs):
    corners = kwargs.get('corners')
    zerophase = kwargs.get('zerophase')
    if len(args) != 1:
        arg = args
    else:
        arg = args[0]
    ret = []
    for i in range(len(arg) - 1):
        if corners is not None and zerophase is not None:
            ret.append((arg[i], arg[i + 1], corners, zerophase))
        elif corners is not None:
            ret.append((arg[i], arg[i + 1], corners))
        elif zerophase is not None:
            raise ValueError('If zerophase is set, corners must also be passed.')
        else:
            ret.append((arg[i], arg[i + 1]))

    return ret


def filter(data, correlations, filters, stack=None, period=24 * 3600): #@ReservedAssignment
    log.info('Filter correlations: %s' % util.parameters())
    for correlation in correlations:
        expr = data.getX(correlation, '????', filter_=None, period=period, stack=stack) + '.QHD'
        files = glob(expr)
        for file_ in files:
            try:
                st = read(file_)
            except Exception as err:
                log.warning('Could not load file, because:/n%s' % str(err))
                continue
            for filter_ in filters:
                st2 = st.copy()
                st2.filter2(*filter_)
                data.writeX(st2, correlation, st[0].stats.endtime, filter=filter_, period=period, stack=stack)

def stack(data, correlations, dt= -1, filters=None, period=24 * 3600, shift=None, onefile=False, yearfiles=False):
    #t1 = t1.__class__(t1.date)
    #t2 = t2.__class__(t2.date)
    log.info('Stack correlations: %s' % util.parameters())
    print 'Stack correlations... '
    if filters is None:
        filters = (None,)
    stack = Stream()
    last_year = None
    for correlation in ProgressBar()(correlations):
        for filter_ in filters:
            try:
                st = read(data.getX(correlation, '*', filter=filter_, period=period) + '.QHD')
            except Exception as err:
                log.warning('Could not load file, because:/n%s' % str(err))
            else:
                for some_traces in streamtimegen(st, dt=dt, start=None, shift=shift):
                    tr = some_traces.calculate('mean')
                    stack.append(tr)
                    this_year = (some_traces[0].stats.starttime).year
                    if last_year is None:
                        last_year = this_year
                    #if yearfiles and (some_traces[0].stats.starttime + period).julday == 1 and len(stack) > 0:
                    if yearfiles and this_year != last_year and len(stack) > 0:
                        data.writeX(stack, correlation, time=some_traces[0].stats.starttime - 365 * 24 * 3600, filter_=filter_, period=period, stack=(dt, shift))
                        last_year = this_year
                        stack = Stream()

                if not onefile:
                    if yearfiles:
                        time = some_traces[0].stats.starttime
                    else:
                        time = None
                    if len(stack) > 0:
                        data.writeX(stack, correlation, time=time, filter=filter_, period=period, stack=(dt, shift))
                    last_year = None
                    stack = Stream()

    if onefile:
        data.writeX(stack, ('all', 'all'), time=None, filter=filters[0], period=period, stack=(dt, shift))


#def stack_hour(data, correlations, t1, t2):
#    #t1 = t1.__class__(t1.date)
#    #t2 = t2.__class__(t2.date)
#    log.info('Stack hour correlations to st: %s' % util.parameters())
#    for correlation in correlations:
#        stream_day = Stream()
#        #t_day = t1 - 24 * 3600
#        #while t_day < t2:
#        #    t_day += 24 * 3600
#        for t_day in daygen(t1, t2):
#            if t_day.getJulday() == 1 and len(stream_day) > 0:
#                data.writeXDay(stream_day, correlation[0], correlation[1], t_day - 24 * 3600)
#                stream_day = Stream()
#            oneday = read(data.getXHour(correlation[0], correlation[1], t_day) + '.QHD')
#            mean_oneday = oneday[0].copy()
#            mean_oneday.data = oneday.calculate('mean')
#            stream_day.append(mean_oneday)
#        data.writeXDay(stream_day, correlation, t_day)

def stack_day(data, correlations, dt= -1, start=None, onefile=False):
    #t1 = t1.__class__(t1.date)
    #t2 = t2.__class__(t2.date)
    log.info('Stack day correlations: %s' % util.parameters())
    if start is not None:
        dt_log = '%s-%s' % (dt, start)
    else:
        dt_log = dt
    stack = Stream()
    for correlation in correlations:
        try:
            days = read(data.getXDay(correlation, '*') + '.QHD')
        except Exception as err:
            log.warning('Could not load file, because:/n%s' % str(err))
        else:
            for somedays in streamdaygen(days, dt=dt, start=start):
                tr = somedays.calculate('mean')
                stack.append(tr)
            if not onefile:
                data.writeXDayStack(stack, correlation, dt_log)
                stack = Stream()
    if onefile:
        data.writeXDayStack(stack, ('all', 'all'), dt_log)


def get_correlations(stations, component='ZNE', stations2=None, only_auto=False, only_cross=False):
    stations = stations.split()
    if stations2 is None:
        stations2 = stations
    else:
        stations2 = stations2.split()
    if component == 'all':
        component = 'ZNE'
    ret = []
    done = []
    for comp in component:
        for st1 in stations:
            for st2 in stations2:
                st1c = st1 + comp
                st2c = st2 + comp
                if only_auto and st1c != st2c:
                    continue
                if only_cross and st1c == st2c:
                    continue
                if not set((st1c, st2c)) in done:
                    ret.append((st1c, st2c))
                    done.append(set((st1c, st2c)))
    return ret

def plotXcorrs(data, correlations, t1, t2, filters=None, filter_now=True, start=None, end=None, select=None, plot_overview=True, plot_years=True,
                        add_to_title='', add_to_file='', show=False, landscape=False, use_dlognorm=True, stack=None, ext='.png', **kwargs):
    figsize = (8.267, 11.693)
    if landscape:
        figsize = figsize[::-1]
    util.checkDir(data.getPlotX(('', ''), t1, stack=stack))
    if filters is None:
        filters = (None,)
    if stack:
        print 'Plot stack xcorrs...'
    else:
        print 'Plot xcorrs...'
    for correlation in ProgressBar()(correlations):
        stations = correlation[0][:-1], correlation[1][:-1]
        try:
            dist = data.stations.dist(*stations)
        except:
            #TaiQ station
            start = -150
            end = 150
        if start is None or end is None:
            if dist >= 120:
                t = (dist // 100) * 50 + 100
            else:
                t = 100
            #print stations, dist, t
            startt = -t
            endt = t
        else:
            startt = start
            endt = end
        if correlation[0] == correlation[1] and (start is None or start < 0):
            startt = 0
        if filter_now:
            try:
                stream_orig = data.readX(correlation, t1, t2, filter=None, stack=stack)
            except IOError as ex:
                print ex
                continue
        for filter_ in filters:
            if filter_now:
                stream = stream_orig.copy()
                if filter_:
                    stream.filter2(*filter_)
            else:
                try:
                    stream = data.readX(correlation, t1, t2, filter=filter_, stack=stack)
                except IOError as ex:
                    print ex
                    continue
            if select:
                stream = stream.select(expr=select)
            #print stream
            if len(stream) > 0:
                if plot_overview:
                    if show:
                        save = None
                        figtitle = add_to_title
                    else:
                        savebase = data.getPlotX(correlation, 'all', filter=filter_, stack=stack) + add_to_file
                        save = savebase + ext
                        savebase = os.path.basename(savebase)
                        figtitle = savebase + ' ' + add_to_title
                    stream.plotXcorr(startt, endt, imshow=True, use_dlognorm=use_dlognorm,
                                     fig=plt.figure(figsize=figsize),
                                     figtitle=figtitle, save=save, show=show,
                                     dateformatter='%Y-%m-%d', #dateformatter='%y %b'
                                     ** kwargs)
                    plt.show()
                if plot_years:
                    for ys in streamyeargen2(stream):
                        t_year, s_year = ys
                        if show:
                            save = None
                            figtitle = add_to_title
                        else:
                            savebase = data.getPlotX(correlation, t_year, filter=filter, stack=stack) + add_to_file
                            save = savebase + ext
                            savebase = os.path.basename(savebase)
                            figtitle = savebase + ' ' + add_to_title
                        s_year.plotXcorr(startt, endt, imshow=True, use_dlognorm=use_dlognorm,
                                     fig=plt.figure(figsize=figsize),
                                     figtitle=figtitle, #'station year ' + add_to_title,
                                     dateformatter='%y-%m-%d', show=False,
                                     save=save,
                                     **kwargs)

#def plotXStacks(data, correlations, start=None, end=None, filters=None, dt= -1, add_to_title='', add_to_file='', **kwargs):
#    if not 'stack' in add_to_title:
#        add_to_title = 'stack ' + add_to_title
#    util.checkDir(data.getPlotXStack(('', ''), dt))
#    if filters is None:
#        filters = (None,)
#    for correlation in correlations:
#        for filter in filters:
#            stations = correlation[0][:-1], correlation[1][:-1]
#            dist = data.stations.dist(*stations)
#            if start is None or end is None:
#                if dist >= 120:
#                    t = (dist // 100) * 50 + 50
#                else:
#                    t = 70
#                startt = -t
#                endt = t
#            else:
#                startt = start
#                endt = end
#            if correlation[0] == correlation[1]:
#                startt = 0
#                #endt *= 2
#            stream = read(data.getXStack(correlation, dt, filter=filter) + '.QHD')
#            if dt == -1:
#                stream.plot_(startt, endt, vmin_rel='vmax',
#                                     #fig=plt.figure(figsize=(8.267, 11.693)),
#                                     figtitle='station ' + add_to_title,
#                                     #dateformatter='%y %b',
#                                     show=False,
#                                     save=data.getPlotXStack(correlation, dt, filter=filter) + add_to_file,
#                                     **kwargs)
#            else:
#                stream.plotXcorr(startt, endt, imshow=True, vmin_rel='vmax',
#                                     fig=plt.figure(figsize=(8.267, 11.693)),
#                                     figtitle='station ' + add_to_title,
#                                     dateformatter='%y-%m-%d', show=False,
#                                     save=data.getPlotXStack(correlation, dt, filter=filter) + add_to_file,
#                                     **kwargs)


def setHIDist(stream, stations, indeg=False, hf='dist', station_splitter='-'):
    for tr in stream:
        stations_str = [st[:-1] for st in
                    tr.stats.station.split(station_splitter)]
        tr.stats[hf] = stations.dist(*stations_str, indeg=indeg)

#def getInformationAbout(data, station, date, station2=None, component='Z'):
#    if station2 == None:
#        station2 = station
#    else:
#        raw_stream2 = data.getRawStream(date, station2, component=component)
#        norm_stream2 = data.getStream(date, station2, component=component)
#        print raw_stream2
#        print norm_stream2
#    correlation = (station + component, station2 + component)
#
#    raw_stream = data.getRawStream(date, station, component=component)
#    norm_stream = data.getStream(date, station, component=component)
#    xcorr = read(data.getXDay(correlation, date))
#    expr = '%r<=st.starttime<%r' % (date - 10, date + 10)
#    xcorr_this_day = xcorr.select(expr=expr)
#    print raw_stream
#    print norm_stream
#    print xcorr_this_day
#    ipshell()

def removeBad(stream, bar_cor=0.8, start=None, end=None, relative='middle'):
    stream2 = stream.copy()
    stream2.trim2(start, end, relative)
    reftr = stream2.calculate('mean')
    st = reftr.stats
    pos_max = np.argmax(np.abs(reftr.data)) / st.sampling_rate - (st.endtime - st.starttime) / 2.
    corr, dt = stretch(stream2, reftr, str_range=0.0, nstr=1, #@UnusedVariable
                       time_windows=((np.abs(pos_max) - 5,), 10),
                       left_side=(pos_max < 0))
    num = 0
    for i, tr in enumerate(stream):
        if corr[0, i] < bar_cor:
            num += 1
            log.info('removing trace no. %d, day %s, because corr=%f<%f' % (i, (tr.stats.starttime + 1).date, corr[0, i], bar_cor))
            stream.remove(tr)
    log.info('%d of %d traces removed' % (num, len(stream) + num))
    #return stream

#def stretch_old(stream, reftr, start=None, end=None, relative='starttime', str_range=0.1, nstr=100, time_windows=None, single_side=False, left_side=False):
#    sr = stream[0].stats.sampling_rate
#    assert sr == reftr.stats.sampling_rate
##    if time_windows is None:
##        tw_mat = np.arange(10 * sr, 20 * sr, dtype=int)[np.newaxis, :] # window from 10s to 20s
#    if time_windows is not None and isinstance(time_windows[1], (float, int)):
#        tw_mat = time_windows_creation(np.array(time_windows[0]) * sr,
#                                       time_windows[1] * sr)
#    else:
##        tw_mat = np.array(time_windows, dtype=int)
##        tw_mat *= sr
##        if tw_mat.shape[1] != 2:
#        raise ValueError('Wrong format for time_window')
#    log.debug('getting data...')
#    data = getDataWindow(stream, start=start, end=end, relative=relative)
#    ref_data = getDataWindow(Stream([reftr]), start=start, end=end, relative=relative)[0, :]
#    if left_side:
#        data = data[:, ::-1]
#        ref_data = ref_data[::-1]
##    if not single_side:
##        data = data[:, len(ref_data) // 2:len(ref_data) // 2 + 10]
##        ref_data = ref_data[len(ref_data) // 2:len(ref_data) // 2 + 10]
##        tw_mat = np.array(((3, 4, 5),))
##        nstr = 3
#    log.debug('create stretching matrices...')
#    str_ref_mat, deltas = stretch_mat_creation(ref_data, str_range=str_range, nstr=nstr, single_side=single_side)
#    log.debug('calculate correlations and time shifts...')
#    corr, dt = velocity_change_estimete(data, tw_mat, str_ref_mat, deltas, single_side=single_side)
##    ipshell()
#    return corr, dt

def stretch(stream, reftr=None, start=None, end=None, relative='starttime', str_range=0.1, nstr=100, time_windows=None, sides='right'):
    sr = stream[0].stats.sampling_rate
    if reftr:
        assert sr == reftr.stats.sampling_rate
    if time_windows is not None and isinstance(time_windows[1], (float, int)):
        tw_mat = time_windows_creation(np.array(time_windows[0]) * int(sr),
                                       time_windows[1] * int(sr))
    else:
        raise ValueError('Wrong format for time_window')
    log.debug('getting data...')
    data = getDataWindow(stream, start=start, end=end, relative=relative)
    ref_data = getDataWindow(Stream([reftr]), start=start, end=end, relative=relative)[0, :] if reftr else None
    log.debug('calculate correlations and time shifts...')
    return time_stretch_estimate(data, ref_trc=ref_data, tw=tw_mat,
                                 stretch_range=str_range, stretch_steps=nstr,
                                 sides=sides)



#def filters(data, station, date, component='Z'):
#    xcorr = read(data.getXDay(station, station2, date))
#    xcorr.

if __name__ == '__main__':
    pass
