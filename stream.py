# by TR

from copy import copy
from obspy.core import AttribDict, UTCDateTime, Stream as ObsPyStream
from obspy.core.util import gps2DistAzimuth
from sito import util, rf, imaging
from sito.trace import Trace
from sito.util import add_doc, cosTaper, calculate
from sito.xcorr import timeNorm, spectralWhitening
import glob
import logging
import numpy as np
import obspy.core  # @UnresolvedImport
import obspy.signal  # @UnresolvedImport
import os.path

log = logging.getLogger(__name__)


def _reprSortedStats(stats):
    """
    Return sorted string representation of Stats object.
    """
    dict_copy = copy(stats)
    priorized_keys = ['network', 'station', 'location', 'channel',
        'starttime', 'endtime', 'sampling_rate', 'delta',
        'npts', 'calib', 'azi', 'inci', 'dist', 'slowness',
        'filter', 'ponset', 'event']
    unwanted_keys = ['_format']
    for key in unwanted_keys:
        if key in dict_copy.keys():
            dict_copy.pop(key)

    # determine longest key name for alignment of all items
    head = ''.join(["%r: %r, " % (key, dict_copy.pop(key)) \
                   for key in priorized_keys if key in dict_copy.keys()])
    head += ''.join(["%r: %r, " % (key, dict_copy.pop(key)) \
                    for key in dict_copy.keys()])
    return 'obspy.core.Stats({%s})' % head[:-2]

# def read(pathname_or_url, format=None, headonly=False, ** kwargs):
#    """Extend read function from stream module. sh entries are moved to stream."""
#    ms = obspy.core.read(pathname_or_url, format, headonly, ** kwargs)
#    if '.' in pathname_or_url:
#        pathname_or_url = pathname_or_url.rsplit('.', 1)[0]
#    if os.path.isfile(pathname_or_url + '.HAT'):
#        with open(pathname_or_url + '.HAT', 'r') as file:
#            content = file.read().split('\n')
#        for i in range(len(ms)):
#            ms[i].stats.update(eval(content[i]))
#    elif ms[0].stats._format == 'Q':
#        for tr in ms:
#            tr.stats['event'] = AttribDict()
#            for key_sh, key_opy in SH_OPY_EVENT_IDX.iteritems():
#                if key_sh in tr.stats.sh.keys():
#                    tr.stats.event[key_opy] = tr.stats.sh[key_sh]
#                    del tr.stats.sh[key_sh]
#            for key_sh, key_opy in SH_OPY_IDX.iteritems():
#                if key_sh in tr.stats.sh.keys():
#                    tr.stats[key_opy] = tr.stats.sh[key_sh]
#                    del tr.stats.sh[key_sh]
#            if len(tr.stats['sh']) == 0:
#                del tr.stats['sh']
#            if 'filter' in tr.stats.keys():
#                #DT;LP:max:25.00;DS by 2;
#                tr.stats.filter = tr.stats.filter.translate(None, ' |;:')
#                tr.stats.filter = tr.stats.filter.replace('DT', 'Dt').replace('max', '').replace('min', '').replace('by', '')
#            if not 'filter' in tr.stats.keys():
#                tr.stats['filter'] = ''
#            if 'comment' in tr.stats.keys():
#                tr.stats.filter += tr.stats.comment
#                del tr.stats['comment']
#            if 'mark' in tr.stats.keys():
#                tr.stats.mark = bool(tr.stats.mark)
#    return Stream(ms)

@add_doc(obspy.core.read)
def read(pathname_or_url, *args, ** kwargs):  # format=None, headonly=False,
    """
    Read waveform files into an Stream object.

    Doc of obspy.core.read:
    """
    pathname_or_url, ext = os.path.splitext(pathname_or_url)
    ms = obspy.core.read(pathname_or_url + ext, *args, ** kwargs)
    ignore_starttime = kwargs.get('ignore_starttime', False)
    content = []
    hatfiles = glob.glob(pathname_or_url + '.HAT')
    hatfiles.sort()
    for hatfile in hatfiles:
        with open(hatfile, 'r') as file_:
            file_content = file_.read().strip('\n')
            file_content = file_content.replace('nan', 'np.nan')
            if file_content != '':
                content.extend(file_content.split('\n'))
    if len(content) > 0 and len(content) == len(ms):
        for i in range(len(ms)):
            st = eval(content[i])  # quiet slow
            if ms[i].stats.get('station') and ms[i].stats.station != st.station or \
               ms[i].stats.get('channel') and ms[i].stats.channel != st.channel or \
               not ignore_starttime and ms[i].stats.get('starttime') and abs(ms[i].stats.starttime - st.starttime) > 1:
                raise ValueError('Error while reading stats from HAT file_')
            ms[i].stats.update(eval(content[i]))
    elif ms[0].stats._format == 'Q':
        from sito.util.main import SH_OPY_EVENT_IDX, SH_OPY_IDX
        for tr in ms:
            tr.stats['event'] = AttribDict()
            for key_sh, key_opy in SH_OPY_EVENT_IDX.iteritems():
                if key_sh in tr.stats.sh.keys():
                    tr.stats.event[key_opy] = tr.stats.sh[key_sh]
                    del tr.stats.sh[key_sh]
            for key_sh, key_opy in SH_OPY_IDX.iteritems():
                if key_sh in tr.stats.sh.keys():
                    tr.stats[key_opy] = tr.stats.sh[key_sh]
                    del tr.stats.sh[key_sh]
            if len(tr.stats['sh']) == 0:
                del tr.stats['sh']
            if 'filter' in tr.stats.keys():
                # DT;LP:max:25.00;DS by 2;
                tr.stats.filter = tr.stats.filter.translate(None, ' |;:')
                tr.stats.filter = tr.stats.filter.replace('DT', 'Dt').replace('max', '').replace('min', '').replace('by', '')
            if not 'filter' in tr.stats.keys():
                tr.stats['filter'] = ''
            if 'comment' in tr.stats.keys():
                tr.stats.filter += tr.stats.comment
                del tr.stats['comment']
            if 'mark' in tr.stats.keys():
                tr.stats.mark = bool(tr.stats.mark)
    ms = Stream(ms)
    # write some event ids if not given
    if len(ms) > 0 and ms[0].stats.get('event') != None and ms[0].stats.event.get('id') == None:
        ms.setEventIDs()
    if len(ms) > 0 and ms[0].stats.get('filter') == None:
        ms.setHI('filter', '')

    if len(ms) > 0 and ms[0].stats.is_fft:
        ms2 = obspy.core.read(pathname_or_url + '_IM' + ext, *args, ** kwargs)
        for i, tr in enumerate(ms):
            tr.data = tr.data + 1j * ms2[i].data
    log.info('Load stream %s with %s traces from %s' % (ms.hash, len(ms), pathname_or_url))
    return ms

class Stream(ObsPyStream):
    """
    Class derieved from obspy.core.Stream with some additional functionality.
    """
    @property
    def hash(self):  # @ReservedAssignment
        return '%03d' % (hash(repr(self)) % 100)

    # @classmethod
    # def read(cls, *args, **kwargs):
    #    return read(*args, **kwargs)

    # @classmethod
    # def stroptions(cls, * args, ** kwargs):
    #    return ' '.join(['%s' % i for i in args] + ['%s=%s' % (i, j) for i, j in kwargs.iteritems()])

    @add_doc(obspy.core.Stream.write)
    def write(self, filename, format_, Qheader=False, ** kwargs):
        """
        Extend obspy.core.Stream.write with logging and writing special header filename.

        Doc of obspy.core.Stream.write:
        """
        util.checkDir(filename)
        if np.any([isinstance(tr.data, np.ma.MaskedArray) for tr in self]):
            log.debug('Filling some masked arrays in stream object before writing.')
            for tr in self:
                tr.data = np.ma.filled(tr.data)
        if Qheader and format_ == 'Q':
            from sito.util.main import SH_OPY_EVENT_IDX, SH_OPY_IDX
            for tr in self:
                if not 'sh' in tr.stats.keys():
                    tr.stats['sh'] = AttribDict({})
                if 'event' in tr.stats.keys():
                    for key_sh, key_opy in SH_OPY_EVENT_IDX.iteritems():
                        if key_opy in tr.stats.event.keys():
                            tr.stats['sh'][key_sh] = tr.stats.event[key_opy]
                for key_sh, key_opy in SH_OPY_IDX.iteritems():
                    if key_opy in tr.stats.keys():
                        tr.stats['sh'][key_sh] = tr.stats[key_opy]
                if 'filter' in tr.stats.keys() and len(tr.stats.filter) > 78:
                    tr.stats.sh.FILTER = tr.stats.filter[:78]
                    tr.stats.sh.COMMENT = tr.stats.filter[78:]
                if 'mark' in tr.stats.keys():
                    tr.stats.sh.MARK = int(tr.stats.mark)
                tr.stats.sh['FILTER'] = ''
        elif format_ == 'SAC':
            for tr in self:
                if not 'sac' in tr.stats.keys():
                    tr.stats['sac'] = AttribDict({})
                # save distance in sac header
                if 'dist' in tr.stats:
                    tr.stats.sac.dist = tr.stats.dist
        if len(self) > 0:
#            from sito import ipshell
#            ipshell()
            if self[0].stats.is_fft:
                st2 = self.copy()
                for tr in st2:
                    tr.data = np.imag(tr.data)
                super(Stream, st2).write(filename + '_IM', format_, ** kwargs)
            super(Stream, self).write(filename, format_, ** kwargs)
        if Qheader and format_ == 'Q':
            for tr in self:
                if 'event' in tr.stats.keys():
                    for key_sh, key_opy in SH_OPY_EVENT_IDX.iteritems():
                        if key_opy in tr.stats.event.keys():
                            del tr.stats['sh'][key_sh]
                for key_sh, key_opy in SH_OPY_IDX.iteritems():
                    if key_opy in tr.stats.keys():
                        del tr.stats['sh'][key_sh]
                if 'filter' in tr.stats.keys() and len(tr.stats.filter) > 78:
                    del tr.stats.sh.COMMENT
                if len(tr.stats['sh']) == 0:
                    del tr.stats['sh']
        if format_ == 'Q':
            towrite = ''
            for tr in self:
                towrite += _reprSortedStats(tr.stats) + '\n'
                # towrite += repr(tr.stats) + '\n'
            with open(filename + '.HAT', 'w') as filename:
                filename.write(towrite)
        log.info('Write stream %s with %s traces to %s' % (self.hash, len(self), filename))

    def writex(self, filename, format_, years=True, stations=True, **kwargs):
        ye = list(set(self.getHI('event.datetime.year')))
        st = self.getStationList()
        if stations and years:
            for s in st:
                for y in ye:
                    self.select(station=s, expr='st.event.datetime.year==%s' % y).write(filename % (s, y), format_, **kwargs)
        elif stations:
            for s in st:
                self.select(station=s).write(filename % s, format_, **kwargs)
        elif not stations and years:
            for y in ye:
                self.select(expr='st.event.datetime.year==%s' % y).write(filename % y, format_, **kwargs)
        else:
            self.write(filename, format_, **kwargs)

    def writey(self, filename, format_, **kwargs):
        self.writex(filename, format_, stations=False, **kwargs)

    def __init__(self, traces=None):
        """
        Extend Stream.__init__. Traces are now converted to Trace objects.
        """
        if not traces == None and len(traces) > 0 and not isinstance(traces[0], Trace):
            traces = [Trace(trace=tr) for tr in traces]
        super(Stream, self).__init__(traces=traces)
        if len(self) > 0 and not self[0].stats.has_key('is_fft'):
            self.setHI('is_fft', False)


    def copy(self, logit=False):
        """
        Return deepcopy of stream.
        """
        newstream = super(Stream, self).copy()
        if logit:
            log.info('Creating new stream %s by copying stream %s' % (newstream.hash, self.hash))
        return newstream

    def print_(self, mod=0):
        """
        Special printing method for header information.

        :param mod: 0-2
        """
        return_string = str(len(self.traces)) + ' Trace(s) in Stream:\n'
        return return_string + "\n".join([tr.print_(mod) for tr in self])

    def getHI(self, hf, operation=None):
        """
        Return list of header information of given header field.

        :param hf: name of header field eg. 'azi', 'event.depth'
        :return: list of header information
        """
        hi = []
#        if '.' in hf:
#            #parts = hf.split('.')
#            for tr in self:
#                hi.append(eval('tr.stats.'+hf))
#        else:
        if hf == 'id':
            hi = [tr.id for tr in self]
        else:
            hi = [eval('tr.stats.' + hf) for tr in self]  # quick and dirty
            # quick but wrong ('.' in hf):
            # hi = [getattr(tr.stats, hf) for tr in self]
            # alternatives: operator.attrgetter
            # http://code.activestate.com/recipes/577346-getattr-with-arbitrary-depth/
        if operation is not None and isinstance(operation, basestring):
            hi = calculate(hi, operation)
        elif operation is not None:
            hi = operation(hi)
        return hi

    def calculate(self, operation, start=None, end=None, relative='starttime',
                   return_trace=True, add_to_stream=False):
        log.info('Calculate %s of stream %s' % (operation, self.hash))
        data = self.getDataWindow(start=start, end=end, relative=relative)
        ret = calculate(data, operation)
        if return_trace:
            ret = Trace(data=ret, header=self[0].stats.copy())
            if add_to_stream:
                self.insert(0, ret)
        return ret

    def getArgMax(self, ret='ind', spline=False, spline_enhance=100, func=None):
        if ret == 'ind' and not spline:
            dtype = int
        elif ret in ('ind', 'time'):
            dtype = float
        elif ret == 'utc':
            dtype = UTCDateTime
        argmax = np.empty(len(self), dtype=dtype)
        dmax = np.empty(len(self))
        for i, tr in enumerate(self):
            argmax[i], dmax[i] = tr.getArgMax(ret=ret, spline=spline,
                                              spline_enhance=spline_enhance,
                                              func=func)
        return argmax, dmax

    def setHI(self, hf, hi):
        """
        Set header information in given header field.

        :param hf: name of header field  eg. 'azi', 'event.depth'
        :param hi: list of entries or entry
        """
        if not isinstance(hi, list) and not isinstance(hi, np.ndarray):
            hi = [hi] * len(self)
        if '.' in hf:
            parts = hf.split('.')
            for i, tr in enumerate(self):
                tr.stats[parts[0]][parts[1]] = hi[i]
        else:
            for i, tr in enumerate(self):
                tr.stats[hf] = hi[i]

    def fft(self, * args, ** kwargs):
        """
        FFT
        """
        log.info('FFT stream %s: %s' % (self.hash, util.parameters()))
        for tr in self:
            tr.fft(*args, ** kwargs)

    def ifft(self, * args, ** kwargs):
        """
        IFFT
        """
        log.info('IFFT stream %s: %s' % (self.hash, util.parameters()))
        for tr in self:
            tr.ifft(*args, ** kwargs)


    def demean(self):
        """
        Subtract the mean from every trace.
        """
        log.info('Demean stream %s' % self.hash)
        for tr in self:
            tr.demean()

    def detrend(self):
        """
        Detrend every trace.
        """
        log.info('Detrend stream %s' % self.hash)
        for tr in self:
            tr.detrend()

    def integrate(self):
        """
        Integrate every trace.
        """
        log.info('Integrate stream %s' % self.hash)
        for tr in self:
            tr.integrate()

    def filter2(self, * args, ** kwargs):
        """
        Wrapper for Stream.filter.
        """
        log.info('Filter stream %s: %s' % (self.hash, util.parameters()))
        for tr in self:
            tr.filter2(*args, ** kwargs)

    def downsample2(self, new_sampling_rate):
        """
        Wrapper for Stream.downsample
        """
        log.info('Downsample stream %s: %s' % (self.hash, new_sampling_rate))
        for tr in self:
            tr.downsample2(new_sampling_rate)

    def taper2(self, *args, **kwargs):
        for tr in self:
            tr.taper2(*args, **kwargs)

    @add_doc(timeNorm)
    def timeNorm(self, *args, **kwargs):
        """
        Time normalization for preparation of cross correlation.

        Doc of sito.xcorr.timeNorm:
        """
        log.info('Time Normalization of stream %s: %s' % (self.hash, util.parameters()))
        for tr in self:
            try:
                tr.timeNorm(*args, **kwargs)
            except Exception as err:
                log.debug('Remove trace at %s because:\n%s' %
                      ((tr.stats.starttime + 0.1).date, str(err)))
                self.remove(tr)

    @add_doc(spectralWhitening)
    def spectralWhitening(self, *args, **kwargs):
        """
        Spectral Whitening for preparation of cross correlation.

        Doc of sito.xcorr.spectralWhite:
        """
        log.info('Spectral Whitening of stream %s: %s' % (self.hash, util.parameters()))
        for tr in self:
            tr.spectralWhitening(*args, **kwargs)


    def trim2(self, start=None, end=None, relative='starttime', pad=False,
              nearest_sample=True):
        """
        Trim time window around ponset.

        :param start, end: seconds relative to around or UTCDateTime
        :param relative: name of header with UTCDateTime or seconds after starttime or
                UTCDateTime or UTCDateTime list
        See util.getTimeIntervall
        See obspy.core.Stream.trim
        """

        log.info('Trim stream %s: %s' % (self.hash, util.parameters()))
        start_utc, end_utc = util.imaging.getTimeIntervall(self, start, end, relative, ret_rel='utc')
        if isinstance(start, UTCDateTime) or type(start) == list:
            start = 'utc'
        if isinstance(end, UTCDateTime) or type(end) == list:
            end = 'utc'
        for i, tr in enumerate(self):
            tr.trim(start_utc[i], end_utc[i], pad, nearest_sample)
            tr.stats.filter += 'TR%s,%s' % (start, end)

    def slice2(self, start=None, end=None, relative='starttime', keep_empty_traces=False):
        """
        Returns new Stream object cut to the given start- and endtime.

        Does not copy the data but only passes a reference. Will by default
        discard any empty traces. Change the keep_empty_traces parameter to
        True to change this behaviour.
        """

        log.info('Select stream %s: %s' % (self.hash, util.parameters()))
        start_utc, end_utc = util.imaging.getTimeIntervall(self, start, end, relative, ret_rel='utc')
        if isinstance(start, UTCDateTime) or type(start) == list:
            start = 'utc'
        if isinstance(end, UTCDateTime) or type(end) == list:
            end = 'utc'
        traces = []
        for i, trace in enumerate(self):
            sliced_trace = trace.slice(start_utc[i], end_utc[i])
            if keep_empty_traces is False and not sliced_trace.stats.npts:
                continue
#            trace.stats.filter += 'SEL%s,%s' % (start, end)
            traces.append(sliced_trace)
        return self.__class__(traces=traces)

    def getDataWindow(self, start, end, relative):
        """
        Return array with data in time window (start, end).

        :param start, end: UTCDateTime or list of UTCDateTimes or floats (seconds)
        :return: numpy.array of shape (N_stream, N_newdata)
        See util.imaging.getDataWindow.
        """

        return util.imaging.getDataWindow(self, start, end, relative)

    def sort(self, keys=('network', 'station', 'location', 'channel',
                         'starttime', 'endtime'), reverse=False, logit=True):
        """
        Override obspy.core.Stream.sort.

        Every sortable header field can be in keys.
        """  # Sort channel reverse (Z,N,E).
        # Check the list and all items.
        if isinstance(keys, basestring):
            keys = [keys]
        if not isinstance(keys, (tuple, list)):
            raise TypeError('keys has to be tuple/list of strings or string')
        # Loop over all keys in reversed order.
        for _i in keys[::-1]:
            if _i == 'component':
                if self[0].stats.channel[-1] in 'ZNE':
                    self.traces.sort(cmp=lambda x, y: cmp(
                                     'ZNE'.index(x.stats.channel[-1]),
                                     'ZNE'.index(y.stats.channel[-1])), reverse=reverse)
                elif self[0].stats.channel[-1] in 'LQT':
                    self.traces.sort(cmp=lambda x, y: cmp(
                                     'LQT'.index(x.stats.channel[-1]),
                                     'LQT'.index(y.stats.channel[-1])), reverse=reverse)
            elif '.' in _i:
                parts = _i.split('.')
                self.traces.sort(key=lambda x: x.stats[parts[0]][parts[1]], reverse=reverse)
            else:
                self.traces.sort(key=lambda x: x.stats[_i], reverse=reverse)
        if logit:
            log.info('Sort stream %s after %s' % (self.hash, util.parameters(only=['keys', 'reverse'])))

    def select(self, expr=None, around=None, indegree=False, time=None, mark=None, event_id=None, logit=False, autocorr=False, xcorr=False, ** kwargs):
        """
        Extend Stream.select

        You cann additionaly choose a list of event_ids or
        a condition expr wherin 'st.' stands for 'trace.stats.'.
        """
        # make given component letter uppercase (if e.g. "z" is given)
        ms = super(Stream, self).select(** kwargs).traces
        if event_id and not isinstance(event_id, list):
            event_id = [event_id]
        if time:
            if isinstance(time, UTCDateTime):
                expr = '%r<=st.starttime' % time
            elif time[1] is None:
                expr = '%r<=st.starttime' % time[0]
            elif time[0] is None:
                expr = 'st.starttime<%r' % time[0]
            else:
                expr = '%r<=st.starttime<%r' % time
        elif autocorr:
            expr = "len(set(st.station.split('-')))==1"
        elif xcorr:
            expr = "len(set(st.station.split('-')))==2"
        if expr:
            expr = expr.replace('st.', 'trace.stats.')
        traces = []
        for trace in ms:
            # skip trace if any given criterion is not matched
            if event_id and trace.stats.event.id not in event_id:
                continue
            if mark != None and mark != bool(trace.stats.mark):
                continue
            if expr and not eval(expr):
                continue
            if around:
                lat, lon, maxdist = around
                if indegree:
                    dist = util.gps2DistDegree(lat, lon,
                                              trace.stats.event.latitude,
                                              trace.stats.event.longitude)
                else:
                    dist = gps2DistAzimuth(lat, lon,
                                           trace.stats.event.latitude,
                                           trace.stats.event.longitude)[0] / 1000.
                if dist > maxdist:
                    continue
            traces.append(trace)
        newstream = self.__class__(traces=traces)
        if logit:
            log.info('Creating new stream %s by selecting traces from stream %s: %s' % (newstream.hash, self.hash, util.parameters()))
        return newstream

    def getStationList(self):
        """
        Return list of stations.
        """
        station_list = []
        for trace in self:
            station = trace.stats['station']
            if station not in station_list: station_list.append(station)
        return station_list

    def getStations(self, file_):
        """
        Return Station instance.
        """
        st_list = self.getStationList()
        if isinstance(file_, basestring):
            import stations
            st = stations.Stations.read(file_)
        else:
            st = file_.copy()
        st.pick(' '.join(st_list))
        return st

    def getEvents(self):
        """
        Return Events instance.
        """
        import events
        return events.Events.readFromStream(self)


    def writeGMTfile(self, filename, start= -20, end=20, reductionFactor=10):
        """
        Write a GMT readable file with data.
        """
        str_list = []
        N_stream = len(self)
        for i, trace in enumerate(self):
            t1 = trace.stats.starttime - trace.stats.ponset  # sh['P-ONSET']
            t2 = trace.stats.endtime - trace.stats.ponset  # sh['P-ONSET']
            N_data = trace.stats.npts
            t = np.linspace(t1, t2, N_data)
            str_list.append('> trace     %d\n' % (i + 1))
            for j, datapoint in enumerate(trace.data):
                if j % reductionFactor == 0 and t[j] > start and t[j] < end:
                    str_list.append('  %f  %f  %f\n' % (t[j], N_stream - i, datapoint))
        with open(filename, 'w') as file_: file_.write(''.join(str_list))

    def check(self, check_ponset=True, check_duplettes=True, delete_incompatible_traces=True):
        """
        Check data for combatibilty for rotation and receiver function methods.

        Data must be in aranged in tuples of 3 with components ZNE.
        Some of the header entries of traces belonging to one station and event
        must be equal.
        Use for example sort(['event.id', 'station', 'component'] before.
        """
        # self.sort(['station', 'eventno', 'component'])
        # if not ( (np.array(self.getHI('npts'))==self[0].getHI('ntps')).all()
        #  and (np.array(self.getHI('sampling_rate'))==self[0].getHI('sampling_rate')).all()):
        #    log.error('Traces must have same sampling_rate and npts.')
        #    return False
        retcode = True
        i = 0
        j = 0
        if check_duplettes:
            while i < len(self) - 1:
                st1 = self[i].stats
                st2 = self[i + 1].stats
                if st1.event.id == st1.event.id and st1.station == st2.station and st1.channel == st2.channel:
                    if st1.npts >= st2.npts:
                        index = i + 1
                        self.pop(i + 1)
                    else:
                        index = i
                    if delete_incompatible_traces:
                        log.debug('Delete incompatible trace %s (duplette)' % (self[index]))
                        self.pop(index)
                        j += 1
                    else:
                        raise ValueError('Stream %s incompatible at least at index %d' % (self.hash, index))
                i += 1
        i = 0
        while i < len(self) - 2:
            st1 = self[i].stats
            st2 = self[i + 1].stats
            st3 = self[i + 2].stats
            if not (
                    (st1.channel.endswith('Z') and st2.channel.endswith('N') and st3.channel.endswith('E') or
                    st1.channel.endswith('L') and st2.channel.endswith('Q') and st3.channel.endswith('T')) and
                    st1.station == st2.station == st3.station and
                    st1.event.id == st2.event.id == st3.event.id and
                    st1.npts == st2.npts == st3.npts and
                    st1.sampling_rate == st2.sampling_rate == st3.sampling_rate and
                    not np.any(np.isnan(self[i].data)) and not np.any(np.isnan(self[i + 1].data)) and not np.any(np.isnan(self[i + 2].data))
                    ) \
                or (check_ponset and not st1.starttime + 10 < st1.ponset < st1.endtime - 30):
                retcode = False
                if delete_incompatible_traces:
                    log.debug('Delete incompatible trace %s' % (self[i]))
                    self.pop(i)
                    if i == len(self) - 3:
                        log.debug('Delete incompatible trace %s' % (self[i]))
                        self.pop(i)
                        log.debug('Delete incompatible trace %s' % (self[i]))
                        self.pop(i)
                else:
                    raise ValueError('Stream %s incompatible at least at index %d' % (self.hash, i))
                j += 1
            else:
                i += 3
        if len(self) % 3 != 0:
            log.warning('Stream must have length 3N for rotation or rf.')
            retcode = False
            if delete_incompatible_traces:
                i = len(self) - 1
                log.debug('Delete incompatible trace %s' % (self[i]))
                self.pop(i)
                j += 1
                if len(self) % 3 != 0:
                    i = len(self) - 1
                    log.debug('Delete incompatible trace %s' % (self[i]))
                    self.pop(i)
                    j += 1
        if not retcode and delete_incompatible_traces:
            log.warning('Delete %s incompatible traces in stream %s, remaining %s traces' % (j, self.hash, len(self)))
        return retcode

    def polar(self, start, end, relative='ponset', retlin=False):
        """
        Calculate polarization.

        Save azimuth (angle between North and largest eigenvalue (L-component)
        in earth surface plane)
        and inclination (angle between vertical and largest eigenvalue)
        in header (lazi and linci) of L-component.
        """
        data = self.getDataWindow(start, end, relative)
        if retlin: lin_list = []
        i = 0
        j = 0
        while j < len(self) - 2:
            try:
                inci, azi, lin = rf.polar(data[i + 2, :],
                                             data[i + 1, :], data[i, :])  # E, N, Z = x, y, z
                if retlin: lin_list.append(lin)
                # wrong documentation in seis.polar.polarization_py ?
                # angle between E (x-axis) clockwise to eigenvalue -> angle between N clockwise to eigenvalue
                self[j].stats.lazi = azi
                self[j + 1].stats.lazi = azi
                self[j + 2].stats.lazi = azi
                # inclination : angle between Z and new L-Component
                self[j].stats.linci = inci
                self[j + 1].stats.linci = inci
                self[j + 2].stats.linci = inci
                i += 3
                j += 3
            except ValueError:
                log.warning('Stream[%d:%d]: error while calculating polarization-> remove traces, remaining %s traces' % (i, i + 3, len(self)))
                log.warning('Stream[%d:%d]: %s' % (i, i + 3, self[j:j + 3]))
                self.pop(j + 2)
                self.pop(j + 1)
                self.pop(j)
                i += 3

        if retlin: return lin_list

    def rotateZNE2LQT(self, start= -2, end=10, relative='ponset', usetheo=False):
        """
        Rotate from ZNE to LQT.

        Call check before.
        See obspy.signal.rotate_ZNE_LQT.
        """
        log.info('Rotate stream %s: %s' % (self.hash, util.parameters()))
        self.polar(start, end, relative)
        if usetheo:
            azi = self.select(component='Z').getHI('azi')
            inci = self.select(component='Z').getHI('inci')
        else:
            azi = self.select(component='Z').getHI('lazi')
            inci = self.select(component='Z').getHI('linci')
        for i in range(len(self) // 3):
            # print azi[i], inci[i]
            Z, N, E = self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data
            L, Q, T = obspy.signal.rotate_ZNE_LQT(Z, N, E, azi[i], inci[i])
            self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data = L, Q, T
        for tr in self:
            tr.stats.filter += 'Rot'
            tr.stats.channel = tr.stats.channel[:-1] + \
                'LQT'['ZNE'.index(tr.stats.channel[-1])]

    def rotateLQT2ZNE(self):
        """
        Backrotate from LQT to ZNE.

        See obspy.signal.rotate_LQT_ZNE.
        """
        log.info('Back rotate stream %s' % (self.hash))
        azi = self.select(component='L').getHI('lazi')
        inci = self.select(component='L').getHI('linci')
        for i in range(len(self) // 3):
            L, Q, T = self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data
            # Z, N, E = util.rotateLQT2ZNE(L, Q, T, azi[i], inci[i])
            Z, N, E = obspy.signal.rotate_LQT_ZNE(L, Q, T, azi[i], inci[i])
            self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data = Z, N, E
        for tr in self:
            tr.stats.filter += 'Tor'
            tr.stats.channel = tr.stats.channel[:-1] + \
                'ZNE'['LQT'.index(tr.stats.channel[-1])]

    @add_doc(rotateZNE2LQT)
    def rotate(self, * args, ** kwargs):
        """
        See Stream.rotateZNE2LQT:
        """
        self.rotateZNE2LQT (*args, ** kwargs)

    @add_doc(rotateLQT2ZNE)
    def brotate(self, * args, ** kwargs):
        """
        See Stream.rotateLQT2ZNE:
        """
        self.rotateLQT2ZNE (*args, ** kwargs)

    def rotateNE2RT(self, start= -2, end=10, relative='ponset'):
        """
        Rotate from NE to RT.

        Call check before.
        See obspy.signal.rotate_NE_RT.
        """
        log.info('Rotate stream from NE to RT %s: %s' % (self.hash, util.parameters()))
        azi = self.select(component='Z').getHI('azi')
        self.polar(start, end, relative)
        for i in range(len(self) // 3):
            # print azi[i], inci[i]
            N, E = self[3 * i + 1].data, self[3 * i + 2].data
            R, T = obspy.signal.rotate_NE_RT(N, E, azi[i])
            # we want to have a right handed system similar to LQT
            # R points away from source
            self[3 * i + 1].data, self[3 * i + 2].data = -R, T
        for tr in self:
            tr.stats.filter += 'RotRT'
            tr.stats.channel = tr.stats.channel[:-1] + \
                'ZRT'['ZNE'.index(tr.stats.channel[-1])]


    def receiverf(self, water=0.01, gauss=2, tshift=10, pad=0, window=None, start= -10, end=30, where='ponset', lenslope=5, return_real=True, normalize=True):
        """
        Apply deconvolution in frequency domain.

        Treat traces with index 3*i as source and traces with index 3*i+1 and
        3*i+2 as response.
        :param water: waterlevel to stabilize the deconvolution (relative to data maximum)
        :param gauss: Gauss parameter of averaging function (std of LP-filter)
        :param tshift: shift the resulting function by that amount
        :param pad: multiply number of samples used for fft by 2**pad.
        :param window, start, end, where, lenslope: use only window (start,end) around
            where of type window (with lenslope seconds of smoothing) of source function
        :param return_real: just use the real part
        See rf.deconv.
        """
        log.info('Deconvolve stream %s: %s' % (self.hash, util.parameters()))
        for i in range(len(self) // 3):
            samp = self[3 * i].stats.sampling_rate
            if window != None:
                src = self[3 * i:3 * i + 1].copy()
                src.window(start, end, where, window, lenslope)
                src = src[0].data
            else:
                src = self[3 * i].data
                start = 0
                end = 0
            rsp = [self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data]
            rf_resp = rf.deconvf(rsp, src, samp, water, gauss, tshift, pad, normalize=normalize)
            # multiply -1 on Q and T component
            if return_real:
                self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data = rf_resp[0].real, -rf_resp[1].real, -rf_resp[2].real
            else:
                self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data = rf_resp[0], -rf_resp[1], -rf_resp[2]
        for tr in self:
            tr.stats.filter += 'RF%s,%s,%swin%s,%s' % (water, gauss, tshift, start, end)
            tr.stats.ponset = tr.stats.starttime + tshift
            tr.stats.tshift = tshift
        if not normalize:
            maxi = np.mean(self[::3].getFunc(np.max))
            for tr in self:
                tr.data /= maxi

    def receivert(self, winsrc=(-20, 80, 5), winrsp=(-20, 80), winrf=(-20, 80), where='ponset', spiking=1, normalize=True):
        """
        Aply deconvolution in time-domain.
        """
        # if not winsrc[1]- winsrc[0]>= winrsp[1]-winrsp[0] >= winrf[1]-winrf[0]:
        #    raise ValueError('Source window has to be bigger than Response window has to be bigger than RF window!')
        log.info('Deconvolve stream %s: %s' % (self.hash, util.parameters()))
        for i in range(len(self) // 3):
            samp = self[3 * i].stats.sampling_rate
            onset = self[3 * i].stats.get(where)
            src = self[3 * i:3 * i + 1].slice(onset + winsrc[0], onset + winsrc[1]).copy()
            src = src[0].data
            # src *= util.main.getWindow('tukey', len(src), 2*winsrc[2]/(winsrc[1]-winsrc[0]))
            src *= cosTaper(len(src), 2 * winsrc[2] / (winsrc[1] - winsrc[0]))
            rsp = self[3 * i:3 * i + 3].slice(onset + winrsp[0], onset + winrsp[1])
            rsp = [tr.data for tr in rsp]
            time_rf = winrf[1] - winrf[0]
            shift = int(samp * ((winrsp[1] - winrsp[0] - winsrc[1] + winsrc[0] - time_rf) / 2 + winrsp[0] - winsrc[0] - winrf[0]))  # shift_zero == 0 for winrsp[1]-winrsp[0] = winsrc[1]-winsrc[0]+time_rf + 2*(-winrsp[0]+winsrc[0]+winrf[0])
            rf_resp = rf.deconvt(rsp, src, spiking, shift, length=int(time_rf * samp), normalize=normalize)
            # multiply -1 on Q and T component
            self[3 * i].data, self[3 * i + 1].data, self[3 * i + 2].data = rf_resp[0], -rf_resp[1], -rf_resp[2]
            for tr in self[3 * i:3 * i + 3]:
                tr.stats.filter += 'RFT%swins%s,%s,%s' % (spiking, winsrc, winrsp, winrf)
                tr.stats.ponset = tr.stats.starttime - winrf[0]
                tr.stats.tshift = -winrf[0]
                tr.stats.npts = len(tr.data)
        if not normalize:
            maxi = np.mean(self[::3].getFunc(np.max))
            for tr in self:
                tr.data /= maxi

    def receiverSH(self, start, end, spiking=1, cut1= -20, cut2=100):
        """
        Deconvolution using Seismic Handler.

        At the moment not working because writing of full header (lazi, linci)
        is not supported anymore.
        """
        log.info('Deconvolve stream %s using SH: %s' % (self.hash, util.parameters()))
        import tempfile
        tempdir = tempfile.gettempdir()
        self.write(tempdir + '/FORSH', 'Q', Qheader=True)
        command = "cd " + tempdir + """ &&
        rm -f ./FORPY.Q?? &&
        rm -f ./FORPY.HAT &&
        mv ./FORSH.HAT ./FORPY.HAT&&
        shc<<END
        read FORSH all
        al all p-onset
        dec3 %s;%s;%s
        cut all %s %s
        write FORPY all
        quit y
END""" % (start, end, spiking, cut1, cut2)
        log.debug(command)
        log.info(util.runBash(command))
        stream = read(tempdir + '/FORPY.QHD', ignore_starttime=True)
        for i, tr in enumerate(stream):
            del tr.stats.sh
            tr.stats.filter += 'RFSH%s,%swin%s,%s' % (-cut1, spiking, start, end)
            tr.stats.ponset = tr.stats.starttime - cut1
            tr.stats.tshift = -cut1
            self[i] = tr

    def moveout(self, model='iasp91', phase='Ps', p0=6.4):
        """
        In-place moveout correction.
        """
        log.info('Moveout correction for stream %s: %s' % (self.hash, util.parameters()))
        if model == 'iasp91':
            model = util.Iasp91()
        for tr in self:
            tr.moveout(model, phase, p0)

    def zmigr(self, model='iasp91', phase='Ps', zmax=750, dz=0.5):
        """
        Z-migration

        Not working!
        """
        log.info('Z-Migration for stream %s: %s' % (self.hash, util.parameters()))
        if model == 'iasp91':
            model = util.Iasp91()
        for tr in self:
            tr.zmigr(model, phase, zmax, dz)

    def writeStationPosition(self, file_):
        """
        Write the station coordinates in stats header fields slat, slon, selev.
        """
        stations = self.getStations(file_)
        for tr in self:
            tr.stats['slat'] = stations[tr.stats.station].latitude
            tr.stats['slon'] = stations[tr.stats.station].longitude
            tr.stats['selev'] = stations[tr.stats.station].elevation

    def pspier(self, depth, filename=None, other_header=None):
        """
        Write piercing point information in header fields rpier, plat, plon

        :param depth: depth of piercing point
        :param file_: file_ with station information, if None the header entries
            must already be set
         """
        log.info('Pspier stream %s: %s' % (self.hash, util.parameters(exclude=['file_'])))
        if filename:
            self.writeStationPosition(filename)
        # if component:
        #    stream = self.select(component=component)
        # else:
        #    stream = self
        slat = np.array(self.getHI('slat'))
        slon = np.array(self.getHI('slon'))
        slowness = np.array(self.getHI('slowness'))
        azi = np.array(self.getHI('azi'))
        rpier, plat, plon = util.pspier(depth, slat, slon, slowness, azi)
        if other_header is None:
            self.setHI('rpier', rpier)
            self.setHI('plat', plat)
            self.setHI('plon', plon)
        else:
            self.setHI('rpier' + other_header, rpier)
            self.setHI('plat' + other_header, plat)
            self.setHI('plon' + other_header, plon)

        for tr in self:
            tr.stats.filter += 'Pspier%d' % depth

    def simpleStack(self, overwrite=True, insert=False, component='all'):
        """
        Calculate sum of traces.

        :param overwrite: if true method overwrites some stats with mean values
        :param insert: if true inserts sum at position 0 in stream
        :return: sumation trace
        """
        stream = self._getStreamFromComponent(component)
        tr = stream[0].copy()
        tr.data = np.mean([i.data for i in stream], axis=0)
        tr.stats.filter += 'Stack%d' % len(stream)
        tr.stats['count'] = len(stream)
        if overwrite:
            for header in 'azi dist lazi ldist event.depth event.latitude event.longitude event.magnitude'.split():
                if (stream[0].stats.has_key(header) or
                    (header.startswith('event.') and
                     stream[0].stats.has_key('event') and
                     stream[0].stats.event.has_key(header[6:]))):
                        setattr(tr.stats, header, np.mean(self.getHI(header)))
                        # exec("tr.stats.%s = np.mean(self.getHI('%s'))" % (header, header))
        if insert:
            self.insert(0, tr)
        return tr

    def phaseStack(self, order=1., component='all'):
        """
        return stack of instantanious phase
        """
        stream = self._getStreamFromComponent(component)
        N = stream[0].stats.npts
        phases = np.zeros((len(stream), N), complex)
        stack_trace = stream[0].copy()
        for i, trace in enumerate(stream):
            analytic_signal, envelope = obspy.signal.cpxtrace.envelope(trace.data)
            phases[i, :] = (analytic_signal / envelope)[:N]
        stack_trace.data = np.abs(np.sum(phases ** order, axis=0)) / len(stream)
        stack_trace.stats.filter += 'PST'
        return stack_trace

    def PWS(self, overwrite=False, order=1, component='Q'):
        """
        return phase weighted stack
        """
        tr1 = self.phaseStack(order=order, component=component)
        tr2 = self.simpleStack(overwrite=overwrite, component=component)
        tr2.data *= tr1.data
        tr2.stats.filter += 'PWS'
        return tr2

    def _getStreamFromComponent(self, component):
        if component in 'QLTZNE':
            return self.select(component=component)
        else:
            return self

    def slantStackOneDir(self, slope, component='Q', p0=7.0, order=1):
        stream = self._getStreamFromComponent(component)
        if slope != 0:
            stream = stream.copy()
            for tr in stream:
                shift = -slope * (tr.stats.slowness ** 2 - p0 ** 2) * tr.stats.sampling_rate
                shiftN = int(round(shift))
                data1 = np.roll(tr.data, shiftN)
                data2 = np.roll(data1, int(shiftN <= shift))
                tr.data = data1 + (data2 - data1) * abs(shift - shiftN)
        ret_tr = stream.PWS(order=order)
        ret_tr.stats.filter += 'slant'
        ret_tr.stats.slope = slope
        ret_tr.stats.slowness = p0
        return ret_tr

    def slantStack(self, slopes=np.linspace(-0.1, 0.1, 21), component='Q', p0=7.0, order=1):
        """
        return slant stack
        """
        stack = self.__class__()
        for slope in slopes:
            tr = self.slantStackOneDir(slope, component, p0, order)
            stack.append(tr)
        return stack

    def getMaxima(self):
        """
        Return np.array of maxima for every trace.
        """
        return np.array([tr.data.argmax() * 1. / tr.stats.sampling_rate for tr in self])

    def getFunc(self, func):
        """
        Return np.array of func for every trace.
        """
        return np.array([func(tr.data) for tr in self])

    def getArg(self, func):
        """
        Return np.array of maxima for every trace.
        """
        return np.array([func(tr.data) * 1. / tr.stats.sampling_rate for tr in self])


    def correlate(self, tr2, shift_sec, start=None, end=None,
                  relative='starttime', demean=False, normalize=True):
        """
        Correlate all traces with given trace and use given maximum shift.
        """
        from sito.xcorr import xcorrf
        data = self.getDataWindow(start, end, relative)
        data2 = Stream([tr2]).getDataWindow(start, end, relative).flatten()
        if shift_sec == 0:
            cors = []
        else:
            cors = Stream()
        for i, tr in enumerate(self):
            sr = tr.stats.sampling_rate
            new_data = xcorrf(data2, data[i, :].flatten(), int(shift_sec * sr),
                              demean=demean, normalize=normalize)
            if shift_sec == 0:
                cors.append(new_data[0])
            else:
                tr_cor = Trace(data=new_data[::-1], header=tr.stats.copy())
                tr_cor.stats.npts = len(new_data)
                tr_cor.stats.sampling_rate = (tr_cor.stats.npts - 1.) / 2 / shift_sec
                tr_cor.stats.is_fft = False
                cors.append(tr_cor)
        return cors

    def correlate_numpy(self, tr2, start, end, start2, end2, relative='starttime',
                        normalize=True, mode='valid'):
        data1 = self.getDataWindow(start, end, relative)
        data2 = Stream([tr2]).getDataWindow(start2, end2, relative).flatten()
        data1[0, :] = Stream([tr2]).getDataWindow(start, end, relative).flatten()
        cors = Stream()
        stdev = 1.
        if normalize:
            cumsum = np.cumsum(data2 ** 2)
            ind = data1.shape[1]
            stdev1 = np.max(cumsum[ind:] - cumsum[:-ind]) ** 0.5
        for i, tr in enumerate(self):
            new_data = np.correlate(data2, data1[i, :].flatten(), mode)
            if normalize:
                stdev = (np.sum(data1[i, :] ** 2)) ** 0.5 * stdev1
            new_data /= stdev
            tr_cor = Trace(data=new_data[::-1], header=tr.stats.copy())
            tr_cor.stats.npts = len(new_data)
            tr_cor.stats.is_fft = False
            cors.append(tr_cor)
        return cors


    def norm(self, method='single', value=1.):
        """
        Normalize all traces in stream to value.

        :param method: 'single' normalize all traces independent.
            'all' normalize all traces with the same factor.
        """
        log.info('Normalize stream %s: %s' % (self.hash, util.parameters()))
        values = np.zeros(len(self))
        if method == 'all':
            for i, tr in enumerate(self):
                values[i] = np.max(np.abs(tr.data))
            fak = value / max(values)
            for tr in self:
                tr.norm(fak=fak)
        else:
            for tr in self:
                tr.norm(value)

    @add_doc(imaging.plot2)
    def plot2(self, * args, ** kwargs):
        """
        Plot always 3 traces in one axis

        Assuming three components in a row. The order is Z,N,E,Z,N,E, etc.
        Plotting around ponset in time window.
        Doc imaging.plot2:
        """
        log.info('Plot2 stream %s: %s' % (self.hash, util.parameters()))
        return imaging.plot2(self, * args, ** kwargs)

    @add_doc(imaging.plotTrace)
    def plotTrace(self, *args, **kwargs):
        log.info('plot first trace of stream %s: %s' % (self.hash, util.parameters()))
        return self[0].plotTrace(*args, **kwargs)

    @add_doc(imaging.Plot.__init__)
    def plot_(self, * args, ** kwargs):
        """
        Plot all traces in one axis around ponset in time window with filling.

        Doc imaging.Plot:
        """
        log.info('Plot stream %s: %s' % (self.hash, util.parameters()))
        return imaging.Plot(self, * args, ** kwargs)

#    @add_doc(imaging.plotComb)
    def plotComb(self, * args, ** kwargs):
        """
        Plot all traces in one axis around ponset in time window with filling.

        Doc imaging.Plot:
        """
        log.info('PlotComb stream %s: %s' % (self.hash, util.parameters()))
        return imaging.plotComb(self, * args, ** kwargs)

#    @add_doc(imaging.plot3)
#    def plot3(self, * args, ** kwargs):
#        """
#        Plot all traces in one axis around ponset in time window with filling.
#
#        Doc imaging.plot3:
#        """
#        log.info('Plot3 stream %s: %s' % (self.hash, util.parameters()))
#        return imaging.plot3(self, * args, ** kwargs)

    @add_doc(imaging.plotRF)
    def plotRF(self, * args, ** kwargs):
        """
        Plot all traces in one axis around ponset in time window with filling.

        Doc imaging.plotRF:
        """
        log.info('PlotRF stream %s: %s' % (self.hash, util.parameters()))
        return imaging.plotRF(self, * args, ** kwargs)

    @add_doc(imaging.plotProfile)
    def plotProfile(self, * args, ** kwargs):
        """
        Plot all traces along profile.

        Doc imaging.plotProfile:
        """
        log.info('PlotProfile stream %s: %s' % (self.hash, util.parameters()))
        return imaging.plotProfile(self, * args, ** kwargs)

    @add_doc(imaging.plotXcorr)
    def plotXcorr(self, * args, ** kwargs):
        """
        Plot all traces along profile.

        Doc imaging.plotXcorr:
        """
        log.info('PlotXcor stream %s: %s' % (self.hash, util.parameters()))
        return imaging.plotXcorr(self, * args, ** kwargs)

    @add_doc(imaging.plotXcorrVsDist)
    def plotXcorrVsDist(self, * args, ** kwargs):
        """
        Plot all traces along profile.

        Doc imaging.plotXcorrVsDist:
        """
        log.info('PlotXcor stream %s: %s' % (self.hash, util.parameters()))
        return imaging.plotXcorrVsDist(self, * args, ** kwargs)

    @add_doc(imaging.plotFilterbands)
    def plotFilterbands(self, * args, ** kwargs):
        """
        Plot all traces along profile.

        Doc imaging.plotXcorr:
        """
        log.info('PlotFilterbands stream %s: %s' % (self.hash, util.parameters()))
        return imaging.plotFilterbands(self, * args, ** kwargs)

    @add_doc(Trace.plotPSD)
    def plotPSD(self, *args, ** kwargs):
        """
        Plot PSD of first trace.

        Doc Trace..plotPSD:
        """
        return self[0].plotPSD(*args, ** kwargs)

    def fftfreq(self):
        return self[0].fftfreq()

    def addXcorrSides(self):
        self.setHI('starttime', UTCDateTime('2000-01-01'))
        for tr in self:
            N = tr.stats.npts
            tr.data = tr.data[N // 2:] + tr.data[N // 2 + N % 2 - 1::-1]

#    @add_doc(imaging.plotStack)
#    def plotStack(self, * args, ** kwargs):
#        """
#        Plot all traces of a stack.
#
#        Doc imaging.plotStack:
#        """
#        log.info('PlotStack stream %s: %s' % (self.hash, util.parameters()))
#        return imaging.plotStack(self, * args, ** kwargs)


    def getBinnedStream(self, bins, bins2=None, header='plon', header2='plat'):
        """
        Return a 'binned' Stream.
        """
        newstream = Stream()
        st = self.copy()
        st.sort(header, logit=False)
        if not bins2:
            for i in range(len(bins) - 1):
                temp = st.select(expr='st.' + header + '>= ' + repr(bins[i]) + ' and st.' + header + '< ' + repr(bins[i + 1]))
                if len(temp) > 0:
                    tr_sum = temp.simpleStack()
                    if 'time' in header:
                        dif_s = 0.5 * (bins[i + 1] - bins[i]) + \
                            (bins[i] - temp[0].stats[header])
                        for header2 in ['starttime', 'endtime', 'ponset', 'sonset']:
                            if header2 in temp[0].stats.keys():
                                tr_sum.stats[header2] = temp[0].stats[header2] + dif_s
                    else:
                        tr_sum.stats[header] = bins[i] + 0.5 * (bins[i + 1] - bins[i])
                    newstream.append(tr_sum)
        else:
            for i in range(len(bins) - 1):
                for j in range(len(bins2) - 1):
                    temp = st.select(expr='st.' + header + '>= ' + repr(bins[i]) + ' and st.' + header + '< ' + repr(bins[i + 1]) +
                                     ' and st.' + header2 + '>= ' + repr(bins2[j]) + ' and st.' + header2 + '< ' + repr(bins2[j + 1]))
                    if len(temp) > 0:
                        tr_sum = temp.simpleStack()
                        tr_sum.stats[header] = bins[i] + 0.5 * (bins[i + 1] - bins[i])
                        tr_sum.stats[header2] = bins2[j] + 0.5 * (bins2[j + 1] - bins2[j])
                        newstream.append(tr_sum)
        if not bins2:
            bins2 = np.array([])
        log.info('Creating new binned stream %s from stream %s: %s' % (newstream.hash, self.hash, util.parameters(add={'bins': 'array%d' % len(bins), 'bins2': 'array%d' % len(bins2)})))
        return newstream

    def setEventIDs(self):
        """
        Set event ids in Stream.
        """
        log.info('Set event id in stream %s' % (self.hash))
        for i, tr in enumerate(self):
            tr.stats.event['id'] = str(i // 3)

    def signoise(self, winsig, winnoise, relative='ponset'):
        """
        Calculate signoise ratio by dividing the maxima of the given time windows.

        Ratio is written in header field 'signoise'.
        """
        log.info('Calculate SigNoise ratio for stream %s: %s' % (self.hash, util.parameters()))
        for tr in self:
            tr.signoise(winsig, winnoise, relative)

    def window(self, start, end, relative='ponset', window='tukey', lenslope=10):
        """
        Window between start and end with args passed to _window function.
        """
        start, end = util.imaging.getTimeIntervall(self, start, end, relative, ret_rel='starttime')
        for i, tr in enumerate(self):
            tr._window(start[i], end[i], window, lenslope)

    def getReasons(self):
        """
        Get a dictionary of reasons for the string to be marked.
        """
        reason_keys = ['SN', 'angle', 'maxL', 'sigQ', 'SN L', 'SN Q', 'broad']
        dict_ = {}
        mark = self.getHI('mark')
        for i, key in enumerate(reason_keys):
            dict_[key] = mark.count(i + 1)
        return dict_

    def setPhase(self, phase):
        """
        Set ponset and dependent information to another phase (eg. PP)
        """
        for trace in self:
            arrival = util.ttt(trace.stats.dist, trace.stats.event.depth).findPhase(phase)
            if arrival != None:
                onset = trace.stats.event.datetime + arrival.time
                stats = AttribDict({'inci': arrival.inci, 'ponset':onset, 'slowness':arrival.slow})
                trace.stats.update(stats)
            else:
                trace.stats['mark'] = True


    def afarm(self, context='normal', signoise=False, signoiseQ=False, maxL=False, sigQ=False, broad=False, angle=360, remove=False):
#   def afarm(self, context='normal', signoise=1.1, signoiseQ=None, maxL=0.5, sigQ=True, broad=True, angle=360, remove=True):

        """
        Mark traces and remove them if they fullfill special conditions.
        marked because ...
        signal:
        mark 1: Signal/Noise ratio on signal is bigger than signoise
        mark 2: difference of theoretical azimuth and azimuth calculated by polarisation is bigger than angle
        RF:
        mark 3: maxL * maximum in (-0.5, 0.5)s is smaller than maximum in (3, 20)
        mark 4: signal on Q
        mark 5: signoise on L
        mark 6: signoise on Q
        mark 7: broad signal

        See source code.
        """
        try:
            self.getHI('mark')
        except KeyError:
            self.setHI('mark', False)
        throwed = self.__class__()
        n_farmed = 0
        if context == 'normal':
            self.signoise([-10, 25], [-50, -20])
            for i in range(len(self) // 3)[::-1]:
                st = self[3 * i].stats
                cond1 = st.signoise <= signoise
                cond2 = abs(st.azi - st.lazi) > angle and abs (st.azi - st.lazi) < 360 - angle
                if cond1 or cond2:
                    if cond1:
                        this_mark = 1
                    elif cond2:
                        this_mark = 2
                    n_farmed += 1
                    self[3 * i].stats['mark'] = self[3 * i + 1].stats['mark'] = self[3 * i + 2].stats['mark'] = this_mark
                    if remove:
                        throwed += self[3 * i:3 * i + 3]
                        self.remove(self[3 * i + 2])
                        self.remove(self[3 * i + 1])
                        self.remove(self[3 * i])
        else:
            self.signoise([0., 5.], [-5, -0.5])
            angle = None
            for i in range(len(self) // 3)[::-1]:
                trL = self[3 * i]
                trQ = self[3 * i + 1]
                pon = trL.stats.ponset
                cond1 = maxL and maxL * abs(trL.slice(pon - 2, pon + 2).max()) <= abs(trL.slice(pon + 3, pon + 20).max())
                cond2 = sigQ and (
                                  abs(trQ.slice(pon, pon + 10).max()) <= abs(trQ.slice(pon + 10, pon + 20).max())
                                  or abs(trQ.slice(pon - 1, pon + 1).max()) >= 2 * abs(trQ.slice(pon + 3, pon + 20).max()))
                cond3 = signoise and trL.stats.signoise <= signoise
                cond4 = signoiseQ and trQ.stats.signoise <= signoiseQ
                cond5 = False
                if broad and not (cond1 or cond2 or cond3 or cond4):
                    tr = trQ.slice(pon - 5, pon + 20)
                    dat = tr.data
                    dat_above = dat > 0.2 * max(dat)
                    dat_under = dat < 0.2 * min(dat)
                    j = 0
                    while j < len(dat):
                        if dat_above[j] == True:
                            length = np.nonzero(dat_above[j:] == False)[0]
                        elif dat_under[j] == True:
                            length = np.nonzero(dat_under[j:] == False)[0]
                        else:
                            j += 1
                            continue
                        if len(length) == 0:
                            length = len(dat) - j
                        else:
                            length = length[0]
                        if length / tr.stats.sampling_rate > broad:
                            cond5 = True
                            break
                        j += length

                if cond1 or cond2 or cond3 or cond4 or cond5:
                    if cond1:
                        this_mark = 3
                    elif cond2:
                        this_mark = 4
                    elif cond3:
                        this_mark = 5
                    elif cond4:
                        this_mark = 6
                    elif cond5:
                        this_mark = 7
                    n_farmed += 1
                    trL.stats['mark'] = trQ.stats['mark'] = self[3 * i + 2].stats['mark'] = this_mark
                    if remove:
                        throwed += self[3 * i:3 * i + 3]
                        self.remove(self[3 * i + 2])
                        self.remove(self[3 * i + 1])
                        self.remove(self[3 * i])
        log.info('Farming - Arguments+Variables: %s' % (util.parameters()))
        log.info('Farming - Farm %s events from stream %s (unmarked %s)' % (n_farmed, self.hash, len(self.select(expr='st.mark==False')) // 3))
        if remove:
            return throwed

    def setHIForHist(self, events, period=24 * 3600):
        start = self[0].stats.starttime + 0.1
        end = self[-1].stats.starttime + 0.1
        hist_list, mag_list = events.produceHist(start, end, period)
        self.setHI('num_events', hist_list)
        self.setHI('max_mag', mag_list)

    def addZeros(self, secs_before, secs_after=None):
        for tr in self:
            tr.addZeros(secs_before, secs_after=secs_after)

    def stretch(self, reftr=None, str_range=0.1, nstr=101, time_windows=None,
                sides='right'):
        from sito.noisexcorr import stretch as stretch_fun
        result = stretch_fun(self, reftr=reftr, str_range=str_range, nstr=nstr,
                             time_windows=time_windows, sides=sides)
        return result



