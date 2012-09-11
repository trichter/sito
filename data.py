from glob import glob
from obspy.core import AttribDict
from obspy.core.util import gps2DistAzimuth
from sito import seismometer
from sito.events import Events
from sito.stations import Stations
from sito.stream import read, Stream
from sito.util import yeargen, timegen
import logging
import os.path
import tempfile
import util
from obspy.seishub import Client

log = logging.getLogger(__name__)

def _getEvents(events, ** kwargs):
    if isinstance(events, basestring):
        events = Events.read(events)
    elif len(kwargs) > 0:
        events = Events.load(events, ** kwargs)
    return events

#def changePaths(data, path):
#    for expr in 'x_hour x_day x_month x_year x_all rf_events rf_results'.split():
#        exec('data.%s = os.path.join(path, os.path.basename(data.%s))' % (expr, expr))
#    return data

FFTW3_WISDOM_FILE = '/home/richter/Data/fftw3_wisdom.txt'

def getCor(st1, st2=None):
    if st2 is None or st1 == st2:
        return st1
    else:
        return st1 + '-' + st2

def getPeriod(period):
    if isinstance(period, basestring):
        return period
    elif period == 24 * 3600:
        return 'day'
    elif period % (24 * 3600) == 0:
        return '%ddays' % (period // (24 * 3600))
    elif period == 3600:
        return 'hour'
    elif period % 3600 == 0:
        return '%dhours' % (period // 3600)
    elif period == -1:
        return 'all'
    elif isinstance(period, (int, float)):
        return '%ds' % period
    else:
        raise TypeError('Argument has invalid type.')

def getTimeFormat(time, period):
    try:
        if period == 'day' or period > 3600:
            return time.year
        else:
            return '%s-%s' % (time.year, time.julday)
    except AttributeError:
        return time

def getFilter(filter_):
    if filter_ is None:
        return ''
    else:
        return ('_filter%s' + '-%s' * (len(filter_) - 1)) % tuple(filter_)

def getStack(time=None, period=None, stack=None):
    ret = ''
    if time is not None:
        ret = '_%s' % getTimeFormat(time, period)
#        if hasattr(time, 'year'):
#            ret += str(time.year)
#        else:
#            ret += time
#        ret += '_'
    if stack is not None:
        try:
            ret += '_%s' % getPeriod(stack)
        except TypeError:
            if stack[1] == None:
                ret += '_%s' % (getPeriod(stack[0]))
            else:
                ret += '_%s_%s' % (getPeriod(stack[0]), getPeriod(stack[1]))
    return ret


class Data(object):
    def __init__(self, data, results, client=None, use_client=False, xcorr_append='', append='', create_symbolic_links=True):
        # raise an Eception if __init__ is not called from a child
        if self.__class__ == Data:
            raise NotImplementedError('This function has to be called from or '
                                      'implemented by the daughter class.')

        self.data = data
        self.results = results
        self.raw = data + '/raw/%s_%d_%03d.mseed'
        self.getstr = self.x_prep = (data + '/xcorr%s' % xcorr_append +
                                     '/prep/%s_%d_%03d')
        self.x_res = xcorr_results = results + '/xcorr%s' % xcorr_append

        self.xcorr = xcorr_results + '/xcorr/%s_%s%s_%s' # period, correlation, filter, time ->1
        #self.x_filter = xcorr_results + '/filter/%s_%s%s_%s' # 1
        self.x_stack = xcorr_results + '/stack/%s_%s%s_stack%s' # period, correlation, filter, number of stacks -> 2
        self.x_plot = xcorr_results + '/plots/%s_%s%s_%s' # 1
        self.x_plot_stack = xcorr_results + '/plots_stack/%s_%s%s_stack%s' # 2

        self.x_ev_prep = self.x_ev_getstr = (data + '/xcorr%s' % xcorr_append +
                                              '/prep/%d')
        self.x_ev_corr = xcorr_results + '/xcorr/%s%s_%s' # correlation, filter, time ->1
        self.x_stack = xcorr_results + '/stack/%s%s_stack%s' # correlation, filter, number of stacks -> 2
        self.x_plot = xcorr_results + '/plots/%s%s_%s' # 1
        self.x_plot_stack = xcorr_results + '/plots_stack/%s%s_stack%s' # 2



#        self.x_day = xcorr_results + '/day/%s_day_%s'
#        self.x_day_stack = xcorr_results + '/stack/%s_stack_%s'
#        self.x_plot_day = xcorr_results + '/plots/%s_day_%s'
#        self.x_plot_day_stack = xcorr_results + '/plots_stack/%s_stack_%s'
#        self.x_hour = xcorr_results + '/hour/%s_hour_%d_%03d'

        self.rf_events = data + '/receiver/events/%s_%s' + append # M5.5_events
        self.rf_results = results + '/receiver/results/%s' + append

        self.client = client
        self.use_client = use_client

        self.stations = self.eventfile = None
        if create_symbolic_links:
            util.checkDir(self.x_res + '/bla')
            try:
                util.checkDir(self.x_prep)
            except OSError:
                import warnings
                warnings.warn('Error with external HD')
            else:
                prepdir = os.path.dirname(self.x_prep)
                if not os.path.islink(prepdir + '/to_results'):
                    os.symlink(self.x_res, prepdir + '/to_results')
                    if not os.path.islink(self.x_res + '/to_prep'):
                        os.symlink(prepdir, self.x_res + '/to_prep')

    def initClient(self, client=None, use_client=True):
        self.use_client = use_client
        if not client:
            with open('/home/richter/.sito/pwd_seishub') as f:
                pwd = f.read().strip('\n')
            client = Client(base_url='http://sec24c74.gfz-potsdam.de:8080',
                            user='richter', password=pwd)
            client = client.waveform
        self.client = client

    def getRawStreamFromClient(self, *args, **kwargs):
        return self.client.getWaveform(*args, **kwargs)


    def setXLogger(self, add=''):
        #logformat = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        from obspy.core import UTCDateTime
        now = UTCDateTime().isoformat()[:19]
        logfile = self.x_res + '/log_%s%s.txt' % (now, add)
        logdfile = self.x_res + '/log_%s%s_debug.txt' % (now, add)
        util.setRootLogger(logfile=logfile, logdebugfile=logdfile, fmode='w')

    def getRawStream(self, date, station, component='Z', endtime=False):
        """
        Return a stream from local files or other database.

        :param date: if endtime==False: The whole day is returned as loaded
                      from day file
                     if endtime!=False: starttime for stream
        :param component: 'Z', 'N', 'E' or 'all'
        :param endtime: read from file only till endtime (if != False)
        :return: one-component or 3-component stream
        """
        raise NotImplementedError('This function has to be implemented by the '
                                  'daughter class.')

    def getStream(self, date, station, component='Z', endtime=False,
                  filename=None, checkfile=False):
        """
        Return a stream from local day files.

        :param date: if endtime==False: The whole day is returned
                     if endtime!=False: starttime for stream
        :param component: 'Z', 'N', 'E' or 'all'
        :param endtime: read from filename only till endtime (if != False)
        :param filename: use this expression as filename (None= use self.raw)
        :return: one-component or 3-component stream
        """

        if filename == None:
            filename = self.getstr + '.QHD'
        filename = filename % (station, date.year, date.julday)
        if not os.path.isfile(filename):
            raise ValueError('No filename for %s %s %s at %s'
                             % (station, component, date.date, filename))
        elif checkfile:
            return True
        try:
            if endtime and date.julday == (endtime - 0.001).julday:
                ms = read(filename, starttime=date, endtime=endtime)
            elif endtime and date.julday != endtime.julday:
                border = date.__class__(date.date) + 24 * 3600
                ms1 = read(filename, starttime=date)#, endtime=border)
                ms2 = self.getStream(border, station, component, endtime)
                ms = ms1 + ms2
                ms.merge()
            else:
                ms = read(filename)
        except:
            raise ValueError('Error reading filename %s' % filename)
        #if component == 'all':
        #    component = 'ZNE'

        NC = len(component)
        if NC == 1:
            ms = ms.select(component=component)
            ms.merge()
            if len(ms) > 1:
                raise ValueError('Gaps in data')
        else:
            if NC == 2:
                ms = (ms.select(component=component[0]) +
                      ms.select(component=component[1]))
            if len(ms) > NC:
                raise ValueError('Gaps in data')
            elif len(ms) < NC:
                raise ValueError('Not enough components')
            Ns = [ms[i].stats.npts for i in range(NC)]
            #N1, N2, N3 = len(ms[0]), len(ms[1]), len(ms[2])
            #Ns = (N1, N2, N3)
            if max(Ns) - min(Ns) > 1:
                raise ValueError('Components have different length')
            elif max(Ns) - min(Ns) == 1:
                for i in range(NC):
                    if Ns[i] > min(Ns):
                        ms[i].data = ms[i].data[:-1]
                        ms[i].stats.ntps -= 1
        if len(ms) == 0:
            raise ValueError('No traces in stream!')
        return ms

    def getDay(self, st, time, filename=None):
        """
        Get filename for One Day File.

        :param st: station
        :param time: UTCDateTime object (year, julday and hour properties)
        :return: string saved in self.prep completed with arguments
        """
        if filename == None:
            filename = self.getstr
        return filename % (st, time.year, time.julday)

    def getDayEv(self, time, filename=None):
        if filename == None:
            filename = self.x_ev_getstr
        return filename % (time.year)


    def getX(self, cor, time=None, filter=None, period=24 * 3600, stack=None): #@ReservedAssignment
        """
        Get filename for xcorr of period.

        :param st1: first station
        :param st2: second station
        :param time: UTCDateTime object (year, julday and hour properties)
        :return: string saved in self.x_hour completed with arguments
        """
        if stack is None:
            return self.xcorr % (getPeriod(period), getCor(*cor), getFilter(filter), getTimeFormat(time, period))
        else:
            return self.x_stack % (getPeriod(period), getCor(*cor), getFilter(filter), getStack(time, period, stack))

    def getXEv(self, cor, time=None, filter=None, stack=None): #@ReservedAssignment
        if stack is None:
            return self.xcorr % (getCor(*cor), getFilter(filter), time.year)
        else:
            return self.x_stack % (getCor(*cor), getFilter(filter), time.year)


#    def getXStack(self, cor, number, filter=None, period=24 * 3600):
#        """
#        Get filename for xcorr of period.
#
#        :param st1: first station
#        :param st2: second station
#        :param time: UTCDateTime object (year, julday and hour properties)
#        :return: string saved in self.x_hour completed with arguments
#        """
#        return self.x_stack % (getPeriod(period), getCor(*cor), getFilter(filter), number)

    def writeXEv(self, stream, * args, **kwargs):
        filename = self.getXEv(*args, **kwargs)
        util.checkDir(filename)
        stream.write(filename, 'Q')

    def writeX(self, stream, * args, **kwargs):
        """
        Write file for xcorr of 1 hour.

        The parameters are passed to getXHour()
        :param stream: stream to write
        :param st1: first station
        :param st2: second station
        :param time: UTCDateTime object (year, julday and hour properties)
        """
        filename = self.getX(*args, **kwargs)
        util.checkDir(filename)
        stream.write(filename, 'Q')

#    def writeXStack(self, stream, * args, **kwargs):
#        """
#        Write file for xcorr of 1 hour.
#
#        The parameters are passed to getXDaySrack()
#        :param stream: stream to write
#        :param st1: first station
#        :param st2: second station
#        :param dt: number of stacked days in one trace (-1 for all)
#        """
#        filename = self.getXStack(* args, **kwargs)
#        util.checkDir(filename)
#        stream.write(filename, 'Q')

    def readX(self, correlation, t1=None, t2=None, period=24 * 3600, select=True, **kwargs):
        st = Stream()
        if t2 is None:
            file_ = self.getX(correlation, t1, period=period, **kwargs) + '.QHD'
            if t1 is None:
                st += read(file_)
            else:
                for file_ in glob(file_):
                    st += read(file_)
        else:
            if period == 'day' or period > 3600:
                iterator = yeargen(t1, t2)
            else:
                iterator = timegen(t1, t2, dt=24 * 3600)
            for t in iterator:
                file_ = self.getX(correlation, t, period=period, **kwargs) + '.QHD'
                try:
                    st += read(file_)
                except (ValueError, IOError):
                    log.warning('An error occured when trying to read file %s' %
                                file_)
            if select:
                st = st.select(expr='%r<st.starttime<%r' % (t1 - 0.1, t2 + 0.1))
        return st

    def getPlotX(self, cor, time=None, filter=None, period=24 * 3600, stack=None): #@ReservedAssignment
        if stack is None:
            return self.x_plot % (getPeriod(period), getCor(*cor), getFilter(filter), getTimeFormat(time, period))
        else:# period, correlation, filter, number of stacks -> 2
            return self.x_plot_stack % (getPeriod(period), getCor(*cor), getFilter(filter), getStack(time, period, stack))

#    def getXDay(self, cor, time):
#        """
#        Get filename for xcorr of 1 day.
#    
#        :param st1: first station
#        :param st2: second station
#        :param time: UTCDateTime object (year, julday and hour properties)
#        :return: string saved in self.x_hour completed with arguments
#        """
#        try:
#            time = time.year
#        except AttributeError:
#            pass
#        return self.x_day % (getCor(*cor), time)
#
#    def getXDayStack(self, cor, dt):
#        """
#        Get filename for xcorr stack.
#    
#        :param st1: first station
#        :param st2: second station
#        :param time: UTCDateTime object (year, julday and hour properties)
#        :return: string saved in self.x_hour completed with arguments
#        """
#        return self.x_day_stack % (getCor(*cor), dt)
#
#    def writeXDay(self, stream, * args, **kwargs):
#        """
#        Write file for xcorr of 1 hour.
#    
#        The parameters are passed to getXHour()
#        :param stream: stream to write
#        :param st1: first station
#        :param st2: second station
#        :param time: UTCDateTime object (year, julday and hour properties)
#        """
#        filename = self.getXDay(* args, **kwargs)
#        util.checkDir(filename)
#        stream.write(filename, 'Q')
#
#    def writeXDayStack(self, stream, * args, **kwargs):
#        """
#        Write file for xcorr of 1 hour.
#    
#        The parameters are passed to getXDaySrack()
#        :param stream: stream to write
#        :param st1: first station
#        :param st2: second station
#        :param dt: number of stacked days in one trace (-1 for all)
#        """
#        filename = self.getXDayStack(* args, **kwargs)
#        util.checkDir(filename)
#        stream.write(filename, 'Q')
#
#    def readDayXcorr(self, correlation, t1, t2, select=True):
#        days = Stream()
#        for t in yeargen(t1, t2):
#            file = self.getXDay(correlation, t) + '.QHD'
#            try:
#                days += read(file)
#            except (ValueError, IOError):
#                log.warning('An error occured when trying to read file %s' %
#                            file)
#        if select:
#            days = days.select(expr='%r<st.starttime<%r' % (t1 + 1, t2 + 1))
#        return days
#
#    def getPlotXDay(self, correlation, time):
#        st1 = correlation[0]
#        st2 = correlation[1]
#        try:
#            time = time.year
#        except AttributeError:
#            pass
#        return self.x_plot_day % (getCor(st1, st2), time)
#
#    def getPlotXDayStack(self, correlation, dt):
##        st1 = correlation[0]
##        st2 = correlation[1]
#        return self.x_plot_day_stack % (getCor(*correlation), dt)
#
#    def getXHour(self, st1, st2, time):
#        """
#        Get filename for xcorr of 1 hour.
#    
#        :param st1: first station
#        :param st2: second station
#        :param time: UTCDateTime object (year, julday and hour properties)
#        :return: string saved in self.x_hour completed with arguments
#        """
#        return self.x_hour % (getCor(st1, st2), time.year, time.julday)
#
#    def writeXHour(self, stream, * args, **kwargs):
#        """
#        Write file for xcorr of 1 hour.
#    
#        The parameters are passed to getXHour()
#        :param stream: stream to write
#        :param st1: first station
#        :param st2: second station
#        :param time: UTCDateTime object (year, julday and hour properties)
#        """
#        filename = self.getXHour(* args, **kwargs)
#        util.checkDir(filename)
#        stream.write(filename, 'Q')


    def writeRFEvents(self, stream, station, time):
        """
        Write file with extracted traces around onsets of events.
        The filename is defined by completing self.rf_events with the arguments
        station and time.year.
        :param stream: stream to write
        :param station: station
        :param time: UTCDateTime object (only year property is used)
        """

        filename = self.rf_events % (station, time.year)
        util.checkDir(filename)
        stream.write(filename, 'Q')

class Parkfield(Data):
    def __init__(self, path=None, project='Parkfield', **kwargs):
        if path is None:
            data = '/home/richter/Data/%s' % project
            results = '/home/richter/Results/%s' % project
        else:
            results = data = '%s/%s' % (path, project)
        super(Parkfield, self).__init__(data=data, results=results, **kwargs)
        if path != None:
            self.raw = ('/home/richter/Data/%s/' % project +
                        'raw/%s_%d_%03d.mseed')
#        self.data = data
#        self.results = results
#        self.raw = data + '/raw/%s_%d_%03d.mseed'
#        self.getstr = self.x_prep = data + '/xcorr/prep/%s_%d_%03d'
#        self.x_hour = results + '/xcorr/hour/%s_hour_%d_%03d'
#        self.x_day = results + '/xcorr/day/%s_day_%d'
#        self.x_month = results + '/xcorr/%s_month'
#        self.x_year = results + '/xcorr/%s_year'
#        self.x_all = results + '/xcorr/%s_all'
#        self.rf_events = data + '/receiver/events/%s_%s' + append # M5.5_events
#        self.rf_results = results + '/receiver/results/%s_%s' + append

        self.events = ('/home/richter/Data/events/events_30-90_mag5.8_'
                       'Parkfield.txt')
        self.stations = Stations.read('/home/richter/Data/stations.txt')
        self.stations.pick('PKD')
        self.paz = seismometer.PAZ_STS2_3

    @util.add_doc(Data.getRawStream)
    def getRawStream(self, date, station, component='Z', endtime=False, checkfile=False):
        return self.getStream(date, station, component=component,
                              endtime=endtime, filename=self.raw, checkfile=checkfile)

class IPOC(Data):
    def __init__(self, path=None, project='IPOC', client='sec24c74',
                  exception=False,
                  use_local_LVC=False, **kwargs):
        if path is None:
            data = '/home/richter/Data/%s' % project
            results = '/home/richter/Results/%s' % project
        else:
            results = data = '%s/%s' % (path, project)
        super(IPOC, self).__init__(data=data, results=results, client=client, **kwargs)

        self.raw_regex = ('/media/PBO/archive/'
                          '{year}/{code}/{station}/{channel}.D/'
                          '{code}.{station}..{channel}.D.{year}.{day:03d}')
        if not use_local_LVC:
            self.raw_regex_LVC = ('/media/PBO/archive/'
                                  '{year}/{code}/{station}/{channel}.D/'
                                  '{code}.{station}.{number}.{channel}.D.{year}'
                                  '.{day:03d}')
        else:
            self.raw_regex_LVC = ('/home/richter/Data/obspyload-data/LVC/'
                                  '{year}/{code}.{station}.{number}.{channel}_'
                                  '{year}_{day:03d}.mseed')

        self.events = '/home/richter/Data/events/2012_03_events_27-93_mag5.5_IPOC.txt'
        self.stations = Stations.read('/home/richter/Data/stations_ipoc.txt')
        self.paz = seismometer.PAZ_STS2_3
        self.exception = exception
        if self.client in (True, 'sec24c74', 'gfz'):
            self.initClient()
            #example usage: client.getWaveform('CX', 'PB01', '', 'WDI', '2009-01-01', '2009-01-02')


    def lookForMseed(self, time, station, channel, code_list='CX'):
        """
        Look for mseed files in remote IPOC repository.
        """
        # /media/PBOarchive/2010/CX/PB07/AE7.D/CX.PB07..AE7.D.2010.001
        year = time.year
        day = time.julday
        # In the case of station LVC, the method is looking for files with
        # channel B?Z/N/E/1/2 and location 10 or 00, network GE
        if station == 'LVC':
            code = 'GE'
            if channel[0] == 'H':
                channel = 'B' + channel[1:]
            if channel[-1] == 'N':
                posDir = ['N', '1']
            elif channel[-1] == 'E':
                posDir = ['E', '2']
            elif channel[-1] == 'Z':
                posDir = ['Z']
            for direction in posDir:
                #for number in ('10', '00'):
                for number in ('10',):
                    channel = channel[:-1] + direction
                    file_ = self.raw_regex_LVC.format(year=year, day=day,
                                                     station=station,
                                                     number=number, code=code,
                                                     channel=channel)
                    if os.path.isfile(file_):
                        return file_
                    log.debug('No file_ for channel %s.%s' % (number, channel))
        # In all other cases he is just looking for the two network codes CX
        # and NC
        else:
            code_list = code_list.split()
            for code in code_list:
                file_ = self.raw_regex.format(year=year, day=day,
                                             station=station, code=code,
                                             channel=channel)
                if os.path.isfile(file_):
                    return file_
                log.debug('No file_ for network %s' % (code))
        log.debug('No file_ for station %s date %s at %s' % (station, time.date, file_))
        return None

    def getChannelFromClient(self, starttime, endtime, station, network='CX',
                             location='', channel='*'):
        ms = self.client.getWaveform(network, station, location, channel,
                                     starttime, endtime)
        if len(ms) == 0:
            raise ValueError('No traces in stream returned by seishub.')
        return ms

    def getRawStreamFromClient(self, starttime, endtime, station, component='Z'):
        if component == 'all':
            component = '?'
        network = 'CX'
        channel = 'HH' + component
        location = ''
        if station == 'LVC':
            #channel 1 2 support
            if component == 'N':
                component = '{N,1}'
            elif component == 'E':
                component = '{E,2}'
            network = 'GE'
            channel = 'BH' + component
            location = '10'
        ms = Stream(self.client.getWaveform(network, station, location, channel,
                                     starttime, endtime))
        if len(ms) == 0:
            raise ValueError('No traces in stream returned by seishub.')
        if station == 'LVC':
            for tr in ms:
                ch = tr.stats.channel
                if ch[-1] in '12':
                    tr.stats.channel = ch[:-1] + ('N' if ch[-1] == '1' else 'E')
        return ms


    @util.add_doc(Data.getRawStream)
    def getRawStream(self, date, station, component='Z', endtime=False, checkfile=False):
        if component == 'all':
            component = 'ZNE'
        NC = len(component)
        if NC > 1:
            stream = Stream()
            if checkfile:
                stream = []
            for comp in component:
                stream.extend(self.getRawStream(date, station, comp,
                                                 endtime))
            if checkfile:
                import numpy as np
                return np.all(stream)
            #if None in stream:
            #    raise ValueError('One or more component is None')
            Ns = [stream[i].stats.npts for i in range(NC)]
            #N1, N2, N3 = len(st_list[0]), len(st_list[1]), len(st_list[2])
            #Ns = (N1, N2, N3)
            if max(Ns) - min(Ns) > 1:
                raise ValueError('Components have different length')
            elif max(Ns) - min(Ns) == 1:
                for i in range(NC):
                    if Ns[i] > min(Ns):
                        stream[i].data = stream[i].data[:-1]
                        stream[i].stats.ntps -= 1
            #return st_list[0] + st_list[1] + st_list[2]
            return stream
        if station == 'LVC':
            log.warning('Using BH channel for LVC')
            file_ = self.lookForMseed(date, station, 'BH' + component)
        else:
            file_ = self.lookForMseed(date, station, 'HH' + component)
        if file_ == None:
            raise ValueError('No IPOC file for %s %s %s' % (
                                station, component, date.date))
        elif checkfile:
            return True
        merge_later = False
        try:
            if endtime and date.julday == endtime.julday:
                ms = read(file_, format='MSEED', starttime=date, endtime=endtime)
            elif endtime and date.julday != endtime.julday:
                border = date.__class__(date.date) + 24 * 3600
                ms1 = read(file_, starttime=date)#, endtime=border)
                ms2 = self.getRawStream(border, station, component, endtime)
                ms = ms1 + ms2
                ms.merge()
                merge_later = True
            else:
                ms = read(file_)
        except (ValueError, TypeError) as ex:
            raise ValueError('Error reading IPOC file %s because:\n%s' % (
                                                            file_, str(ex)))
        if len(ms) == 0:
            raise ValueError('No traces in IPOC stream!')
        if station == 'LVC':
            for tr in ms:
                if tr.stats.channel[-1] == '1':
                    tr.stats.channel = tr.stats.channel[:-1] + 'N'
                elif tr.stats.channel[-1] == '2':
                    tr.stats.channel = tr.stats.channel[:-1] + 'E'
        if any([network == 'NC' for network in ms.getHI('network')]):
            # change network code to CX
            ms.setHI('network', 'CX')
            if merge_later:
                ms.merge()
        return ms

class ParkfieldTest(Parkfield):
    def __init__(self, data_path, temp_path=tempfile.gettempdir(),
                 project='ParkfieldTest'):
        super(ParkfieldTest, self).__init__(path=temp_path, project=project,
                                            append='')
        self.raw = data_path + '/%s' % project + '/raw/%s_%d_%03d.mseed'
        self.stations = Stations.read(data_path + '/stationdata.txt')
        self.stations.pick('PKD')
    def copyTestDataFiles(self, t1, t2, component='Z'):
        pkd = Parkfield()
        t1 = t1.__class__(t1.date)
        t2 = t2.__class__(t2.date)
        for station in self.stations.getNames().split():
            t_day = t1 - 24 * 3600
            while t_day < t2:
                t_day += 24 * 3600
                stream_day = pkd.getRawStream(t_day, station=station,
                                              component=component)
                stream_day.write(self.getDay(station, t_day, self.raw), 'MSEED')

class IPOCTest(Parkfield):
    def __init__(self, data_path, temp_path=tempfile.gettempdir(),
                 project='IPOCTest'):
        super(IPOCTest, self).__init__(path=temp_path, project=project,
                                       append='')
        self.raw = data_path + '/%s' % project + '/raw/%s_%d_%03d.mseed'
        self.stations = Stations.read(data_path + '/stationdata.txt')
        self.stations.pick('PB01 PB02 PB03')
    def copyTestDataFiles(self, t1, t2, component='Z'):
        ipocarchive = IPOC()
        t1 = t1.__class__(t1.date)
        t2 = t2.__class__(t2.date)
        for station in self.stations.getNames().split():
            t_day = t1 - 24 * 3600
            while t_day < t2:
                t_day += 24 * 3600
                stream_day = ipocarchive.getRawStream(t_day, station=station,
                                                      component=component)
                stream_day.write(self.getDay(station, t_day, self.raw), 'MSEED')



def eventPicker(data, component='all', phase='P', window=(-100, 400),
                filter=(None, None), new_sampling_rate=100, write=True, #@ReservedAssignment
                ** kwargs):
    """
    Pick window around onset of events from mseed files.

    The resulting stream is written in seperate files for each station and year.
    :param data: data object with stations property and getRawStream,
                 writeRFEvents methods
    :param events: file with events, Events object or None (in this case kwargs
        have to be defined) - passed to _getEvents
    :param component: 'Z', 'N', 'E' or 'all'
    :param phase: which ponset is used? 'P', 'PP' or 'S' or something else
        consider that events must show this phase for the stations
    :param window: window around pnset in seconds
    :param filter: filter stream between these frequencies
    :param new_sampling_rate: downsample stream to rhis sampling rate
    :param write: if True: everything is written to files
        if False: return stream object
    :kwargs: passed to _getEvents
        - in the end they are passed to events.Events.load function
        if param events == None
    """
    log.info('Start event picker: %s' % util.parameters())
    try:
        log.info('Data used %s' % data.raw)
    except:
        log.info('Data regex used %s' % data.raw_regex)
    log.info('Extraced data for events will be saved in %s' % data.rf_events)
    if data.events == None and len(kwargs) == 0:
        raise Exception('No arguments to determine events!')
    failure_list = []
    if write:
        stream_return = None
    else:
        stream_return = Stream()
    stations = data.stations
    all_events = _getEvents(data.events, **kwargs)
    all_events.sort()
    log.info('Events between %s and %s' % (all_events[0].datetime.date,
                                           all_events[-1].datetime.date))
    first_year = all_events[0].datetime.year
    last_year = all_events[-1].datetime.year
    for station_name, station in stations.items():
        for year in range(first_year, last_year + 1):
            events = all_events.pick(after='%s-1-1' % year, before='%s-1-1' %
                                     (year + 1), replace=False)
            stream_year = Stream()
            for event in events:
                dist = util.gps2DistDegree(station.latitude, station.longitude,
                                           event.latitude, event.longitude)
                baz = gps2DistAzimuth(station.latitude, station.longitude,
                                      event.latitude, event.longitude)[1]
                arrival = util.ttt(dist, event.depth).findPhase(phase)
                if arrival == None:
                    log.warning('Phase %s not present at distance %s depth %s' %
                                (phase, dist, event.depth))
                    arrival = util.ttt(dist, event.depth)[0]
                onset = event.datetime + arrival.time
                t1 = onset + window[0]
                t2 = onset + window[1]
                try:
                    stream = data.getRawStream(t1, station_name, component, t2)
                except Exception as ex:
                    failure_list.append((station_name, event.id, str(ex)))
                    continue
                # Write header entries and basic data processing (filtering, downsampling)
                stats = AttribDict({'event': event, 'station': station_name,
                                    'dist':dist, 'azi':baz,
                                    'inci': arrival.inci,
                                    phase.lower() + 'onset':onset,
                                    'slowness':arrival.slow, 'filter':''})
                for trace in stream:
                    trace.stats.update(stats)
                stream_year.extend(stream)
            if len(stream_year) > 0:
                stream_year.demean()
                stream_year.detrend()
                if filter[0] != None or filter[1] != None:
                    stream_year.filter2(freqmin=filter[0], freqmax=filter[1])
                if new_sampling_rate <= (max(stream_year.getHI('sampling_rate'))
                                         / 2.):
                    stream_year.downsample2(new_sampling_rate)
                if write:
                    data.writeRFEvents(stream_year, station_name,
                                       event.datetime)
                else:
                    stream_return.extend(stream_year)
    if len(failure_list) > 0:
        log.warning('Failed to load the data for:\nstation     event.id     '
                    'reason\n' + '\n'.join([' '.join(entry) for entry in
                                            failure_list]))
    if write:
        return failure_list
    else:
        return stream_return, failure_list


# Parkfield class:
#    @util.add_doc(Data.getStream)
#    def getStream(self, date, station, component='Z', endtime=False):
#        file = self.raw % (date.year, date.julday)
#        if not os.path.isfile(file):
#            raise ValueError('No mseed file for station %s %s' % (station, date.date))
#        try:
#            if endtime and date.julday == (endtime-0.001).julday:
#                ms = read(file, starttime=date, endtime=endtime)
#            elif endtime and date.julday != endtime.julday:
#                border = date.__class__(date.date) + 24 * 3600
#                ms1 = read(file, starttime=date, endtime=border)
#                ms2 = self.getStream(border, station, component, endtime)
#                ms = ms1 + ms2
#                ms.merge()
#            else:
#                ms = read(file)
#        except:
#            raise ValueError('Error reading file %s' % file)
#        if component != 'all':
#            ms = ms.select(component=component)
#            if len(ms)>1:
#                raise ValueError('Gaps in data')
#        else:
#            if len(ms)>3:
#                raise ValueError('Gaps in data')
#            elif len(ms)<3:
#                raise ValueError('Not enough components')
#            N1, N2, N3 = len(ms[0]), len(ms[1]), len(ms[2])
#            Ns = (N1, N2, N3)
#            if max(Ns) - min(Ns) > 1:
#                raise ValueError('Components have different lenght')
#            elif max(Ns) - min(Ns)==1:
#                for i in range(3):
#                    if Ns[i] > min(Ns):
#                        ms[i].data = ms[i].data[:-1]
#                        ms[i].stats.ntps -= 1
#        if len(ms) == 0:
#            raise ValueError('No traces in stream!')
#        return ms
