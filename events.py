# -*- coding: utf-8 -*-
# by TR

from obspy.core import AttribDict, UTCDateTime as UTC
from obspy.core.util import gps2DistAzimuth
from sito import util
from sito.util import add_doc, gps2DistDegree, feregion
import logging
import numpy as np
import obspy.neries  #@UnresolvedImport
import re
from obspy.core.event import Catalog, Event, Origin, Magnitude, EventDescription

log = logging.getLogger(__name__)

#TODO change to obspy.core.event.Catalog class
#### DEPRECIATED
class Events(list):
    """
    Class for holding event data.

    Use class methods event.load or event.read to get event data.
    """
    header = '{:d} Events:  datetime          latitude longitude   depth' + \
        '   mag         event_id     author qual  information   or_id   flynn_region\n'
    #format_ = '{datetime} {datetime_quality} {latitude:-8.3f}  {longitude:-8.3f}  ' + \
    #    '{depth:-6.1f} {depth_quality}  {magnitude:3.1f} {magnitude_type}    ' + \
    #    '{id} {origin_id} {author:4s} {quality} {information} {flynn_region}\n'
    format_ = '{datetime} {datetime_quality:3s} {latitude:-8.3f} {longitude:-8.3f}  ' + \
        '{depth:-6.1f} {depth_quality:2s} {magnitude:3.1f} {magnitude_type:2} ' + \
        '{id} {author:5s} {quality:3} {information:16s} {origin_id:3} {flynn_region}\n'
    #template = Template('$datetime  $latitude  $longitude  $depth  $magnitude  $magnitude_type ' + \
    #    '$id $origin_id $author $flynn_region\n')
    format_GMT = '{longitude:-8.3f}  {latitude:-8.3f}\n'
    format_GMT2 = '{longitude:-8.3f}  {latitude:-8.3f}  {depth:-8.2f}\n'
    header_LATEX = (#'\\multicolumn{{6}}{{c}}{{{:d} Events}} \\\\\n'
                    '%time & latitude \\inu{{\degree}} & '
                    'longitude \\inu{{\\degree}} '
                    '& depth \\inu{{km}} & magnitude & Flynn region \\\\\n'
                    '%\\addlinespace\\midrule\\addlinespace\n')
    header_LATEX2 = (#'\\multicolumn{{6}}{{c}}{{{:d} Events}} \\\\\n'
                    '%time & latitude \\inu{{\degree}} & '
                    'longitude \\inu{{\\degree}} '
                    '& depth \\inu{{km}} & magnitude & RF usage & Flynn region \\\\\n'
                    '%\\addlinespace\\midrule\\addlinespace\n')
    format_LATEX = ('{latex_datetime} & {latitude:-8.3f} & {longitude:-8.3f} & '
                    '{depth:-6.0f} & {magnitude:3.1f} & {flynn_region} \\\\\n')
    format_LATEX2 = ('{latex_datetime} & {latitude:-8.3f} & {longitude:-8.3f} & '
                    '{depth:-6.0f} & {magnitude:3.1f} & \\texttt{{{rf_usage}}} & '
                    '{flynn_region} \\\\\n')

    regex = """
    ^\s*
    (?P<datetime>[-:.\dTZ]+)\s(?P<datetime_quality>...)\s+
    (?P<latitude>[-+\d.]+)\s+
    (?P<longitude>[-+\d.]+)\s+
    (?P<depth>[-+\d.]+)\s(?P<depth_quality>..)\s+
    (?P<magnitude>[\d.]+).
    (?P<magnitude_type>..)\s+
    (?P<id>\w+)\s+
    (?P<author>.{5}).
    (?P<quality>...).
    (?P<information>.{16}).
    (?P<origin_id>...).
    (?P<flynn_region>.+)
    $
    """
    regex_NEIC = """ # NEIC-combatible
    ^[^\d]*
    (?P<datetime>[\d\s.]{20})[\s*?A-Z]+     #(?P<year>\d+)\s+(?P<month>\d+)\s+(?P<day>\d+)\s+(?P<time>[\d.]+)[\s*?A-Z]+
    (?P<latitude>-?[\d.]+)\s+
    (?P<longitude>-?[\d.]+)\s+
    (?P<depth>\d*)[^\d]{3}
    (?P<depthSTD>.{4})\s+       # STD can be empty
    (?P<magnitude>[\d.]+).*$
    """
    regex_NEIC2 = """ # NEIC-combatible http://earthquake.usgs.gov/earthquakes/eqarchives/epic/sr_format.php
    ^.
    (?P<author>.{7}).
    (?P<datetime>.{20})(?P<datetime_quality>..).
    (?P<latitude>.{7}).
    (?P<longitude>.{8}).
    (?P<depth>.{4})(?P<depth_quality>.).
    (?P<STD>.{5})..
    (?P<magnitude>\d.\d).{32}
    (?P<origin_id>\d\d\d).
    (?P<quality>...).
    (?P<information>.{16}).*$
    """
    regex_NEIC_compressed = """ # NEIC combatible
    ^[^\d]*
    (?P<datetime>[\d\s.]{18})[\s*?A-Z]+
    (?P<latitude>-?\d+.\d{3})\s*
    (?P<longitude>-?\d+.\d{3})\s*
    (?P<depth>\d*).{5}
    (?P<depthSTD>.{4}).{11}
    (?P<magnitude>[\d.]{3}).*$         # mag can be empty in seldom cases, in this case, there is no match!!!
    """
    regex_GEOFON = """ # GEOFON-combatible # well apart from N and S in coordinates
    ^(?P<datetime>[\d\s:-]{19})\s+
    (?P<magnitude>\d.\d)\s+
    (?P<latitude>\d+.\d+\s[NS])\s+
    (?P<longitude>\d+.\d+\s[WE])\s+
    (?P<depth>\d+) # STD is empty
    """

    def __init__(self, event_list):
        """
        Init Event instance from event_list of events
        """
        event_list = [AttribDict(event) for event in event_list]
        for event in event_list:
            event.datetime = UTC(event.datetime)
            if feregion is not None:
                event.flynn_region = feregion(event.latitude, event.longitude)
            for item in ['datetime_quality', 'depth_quality', 'magnitude_type', 'author', 'quality', 'information', 'origin_id', 'flynn_region']:
                if event.get(item) == None:
                    event[item] = ''
            #if not event.get('magnitude_type'):
            #    event['magnitude_type'] = 'xx'
            if event.get('id') == None:
                event['id'] = (str(event['datetime']))[:-4].replace('-', '').replace(':', '').replace('.', '')
        super(Events, self).__init__(event_list)

    def __add__(self, events):
        if not isinstance(events, self.__class__):
            raise TypeError
        return self.__class__(super(Events, self).__add__(events))

    def __iadd__(self, events):
        if not isinstance(events, self.__class__):
            raise TypeError
        return self.__class__(super(Events, self).__iadd__(events))

    def __repr__(self):
        return 'Events({0})'.format(super(Events, self).__repr__())

    def __str__(self):
        return self.write(False)

    def __getslice__(self, i, j, k=1):
        return self.__class__(self[i:j:k])

    def get(self, what):
        list_ = []
        for event in self:
            list_.append(event[what])
        if isinstance(list_[0], (float, int)):
            list_ = np.array(list_)
        return list_

    @classmethod
    def read(cls, event_file, regex=regex):
        """
        Read events from event_file with the help of given regular expression.

        :param event_file: location of event_file
        :param regex: regular expression, see Events.regex<Tab>
        """
        with open(event_file, 'r') as f:
            filedata = f.read()
        event_matches = re.finditer(regex, filedata, re.VERBOSE + re.MULTILINE)
        list_ = [i.groupdict() for i in event_matches]
        #util.ipshell()
        for event in list_:  # convert numbers to float and int types
            for key, item in event.iteritems():
                if util.isint(item):
                    event[key] = int(item)
                elif util.isfloat(item):
                    event[key] = float(item)
                else:
                    event[key] = item.strip()
                    #if event[key] == '':
                    #    event[key] = None
                #if key == 'depth' and regex == cls.regex:
                #    event[key] *= 1
        #util.ipshell()
        log.info('Read event information of %d events from events event_file %s' % (len(list_), event_file))
        return cls(list_)

    @classmethod
    @add_doc(obspy.neries.Client.getEvents)
    def load(cls, client=None, ** kwargs):
        """
        Load events from the web with the help of neries client.

        :param client: obspy.neries client, None to create a new
        :kwargs: arguments passed to method getEvents of client

        Doc of obspy.neries.Client.getEvents:
        """
        if client == None:
            client = obspy.neries.Client()
        list_ = client.getEvents(** kwargs)
        events = cls(list_)
        for event in events:
            event['id'] = event['event_id']
            event['depth'] = -event['depth']
            del event['event_id']
        log.info('Load event information of %d events from client %s' % (len(list_), client))
        return events

    def deleteDoublettes(self):
        self.sort('datetime')
        for i in range(len(self) - 1)[::-1]:
            if self[i].id == self[i + 1].id or self[i].latitude == self[i + 1].latitude and self[i].longitude == self[i + 1].longitude:
                log.info('Remove doublette entry with index %s: %s' % (i + 1, self[i + 1]))
                self.remove(self[i + 1])

    def printDoublettes(self):
        self.sort('datetime')
        for i in range(len(self) - 1)[::-1]:
            if self[i].id == self[i + 1].id or self[i].latitude == self[i + 1].latitude and self[i].longitude == self[i + 1].longitude:
                print 'Doublette entry with index %s' % (i + 1)

    @classmethod
    def findDifference(cls, ev1, ev2, dif=10):
        i = j = l = m = 0
        while i < len(ev1) and j < len(ev2):
            if ev1[i].datetime < ev2[j].datetime - dif:
                print 'missing in ev2: ev1[%d]   mag %f' % (i, ev1[i].magnitude)
                i += 1
                l += 1
            elif ev1[i].datetime > ev2[j].datetime + dif:
                print 'missing in ev1: ev2[%d]   mag %f' % (i, ev1[i].magnitude)
                j += 1
                m += 1
            else:
                i += 1
                j += 1
        print i, j
        print len(ev1), len(ev2)
        print 'missing %d in ev1, %d in ev2' % (m, l)

    @classmethod
    def readFromStream(cls, stream):
        """
        Read events from an ObsPy stream.
        """
        list_ = []
        event_id_list = []
        for tr in stream:
            if tr.stats.event.id not in event_id_list:
                list_.append(tr.stats.event)
                event_id_list.append(tr.stats.event.id)

        log.info('Read event information of %d events from stream %s' % (len(list_), stream.hash))
        return cls(list_)

    def write(self, file_, format=format_, header=header):  #@ReservedAssignment
        """
        Write events to file_.

        file_: location of file_
        format: format of entries, see Events.format<Tab>
        header: header, see Events.header, False for no header
        """
        str_list = []
        if header and '{' in header:
            str_list.append(header.format(len(self)))
        if format == self.format_LATEX or format == self.format_LATEX2:
            for event in self:
                event.latex_datetime = str(event.datetime)[:19]
                event.latex_datetime = (event.latex_datetime[:10] + ' ' +
                                        event.latex_datetime[11:])
        for event in self:
            try:
                str_list.append(format.format(** event))
            except KeyError:
                event['origin_id'] = 102
                event['author'] = 'unknown'
                event['flynn_region'] = 'unknown'
                str_list.append(format.format(** event))
                #str_list.append(template.safe_substitute(event))
        output = ''.join(str_list)
        if file_:
            with open(file_, 'w') as f:
                f.write(output)
            log.info('Write events to file_ ' + file_)
        else:
            return output

    def write_LATEX(self, file_, format=format_LATEX, header=header_LATEX):
        return self.write(file_, format=format, header=header)
    def write_LATEX2(self, file_, format=format_LATEX2, header=header_LATEX2):
        return self.write(file_, format=format, header=header)

    def sort(self, keys='datetime', reverse=False, logit=True):
        """
        Sort events after list of keys.
        """
        # Check the list and all items.
        items = ['datetime', 'id', 'latidue', 'longitude', 'depth', 'magnitude']
        msg = "keys must be a list of item strings. Available items to " \
            "sort after: \n" + ' '.join(items)
        if isinstance(keys, basestring):
            keys = [keys]
        if not isinstance(keys, list):
            raise TypeError(msg)
        for _i in keys:
            try:
                items.index(_i)
            except:
                raise TypeError(msg)
        # Loop over all keys in reversed order.
        for _i in keys[::-1]:
            super(Events, self).sort(key=lambda x: x[_i], reverse=reverse)
        if logit:
            log.info('Sort events after %s' % util.parameters(only=['keys', 'reverse']))

    def produceHist(self, start=None, end=None, period=24 * 3600):
        """
        produce Histogramm Information and saves into stream
        """
        self.sort()
        if len(self) == 0:
            return
        if start is None:
            start = UTC(self[0].datetime.date)
        if end is None:
            end = UTC(self[-1].datetime)
        entries = int((end - start) / period) + 1
        hist_list = np.zeros(entries)
        mag_list = np.zeros(entries)
        for event in self:
            #date = UTC(event.datetime.date)
            entry = int((event.datetime - start) / period)
            print entry
            try:
                hist_list[entry] += 1
                mag_list[entry] = max(event.magnitude, mag_list[entry])
            except IndexError:
                pass
        return hist_list, mag_list

    def pick(self, latitude=None, longitude=None, minval=0, maxval=180, indegree=True,
             after='1900-01-01', before='3000-01-01', bigger=0.,
             smaller=10., replace=True):
        """
        Pick events fullfilling the given conditions.

        :param latitude, longitude: coordinates for distance condition
        :param minval, maxval: distance of event has to be between this values
        :param indegree: True if minval and maxval in deg, False if in km
        :param after, before: UTCDateTime objects or strings with time range
        :param bigger, smaller: magnitude range
        :param replace: if True the data in the event list is overwritten
        :return: picked Events instance
        """
        if indegree:
            degorkm = 'deg'
        else:
            degorkm = 'km'
        newdata = []
        dist = 50
        for event in self[::-1]:
            if latitude != None and longitude != None:
                if not indegree:
                    dist = gps2DistAzimuth(event.latitude, event.longitude,
                                           latitude, longitude)[0] / 1000.
                else:
                    dist = gps2DistDegree(event.latitude, event.longitude,
                                          latitude, longitude)
            if bigger <= event.magnitude and smaller >= event.magnitude and \
                dist >= minval and dist <= maxval and \
                UTC(after) <= event.datetime and UTC(before) >= event.datetime:
                    newdata.append(event)
            elif replace:
                self.remove(event)
        if latitude == None:
            latitude = 0
        if longitude == None:
            longitude = 0
        log.info('Pick %d events with distance between %d%s and %d%s from coordinates lat:%5.2f lon:%5.2f; between the dates %s and %s and between the magnitudes %3.1f and %3.1f'
                 % (len(newdata), minval, degorkm, maxval, degorkm, latitude, longitude, after, before, bigger, smaller))
        return self.__class__(newdata[::-1])

    def plot_(self, m, color='datetime', radius='magnitude', alpha=1., show_colorbar=False):
        import pylab as plt
        x, y = m(self.get('longitude'), self.get('latitude'))
        try:
            color_val = self.get(color)
        except KeyError:
            color_val = color
        else:
            if color == 'datetime':
                def UTC2year(utc):
                    import calendar
                    year = utc.year
                    return year + utc.julday / (365. + calendar.isleap(year))
                color_val = np.array([UTC2year(val) for val in color_val])
        try:
            radius_val = self.get(radius)
        except KeyError:
            radius_val = radius
        else:
            if radius == 'magnitude':
                radius_val *= 10  # we will scale the dots by 10 time the magnitude
        m.plot(x, y, ls='', ms=radius_val**0.5, mfc=color_val, marker='o',
               mew=0, alpha=alpha)
        if show_colorbar:
            c = plt.colorbar(orientation='horizontal', shrink=0.4)
            c.set_label(color)

    def plot(self, lat=0., lon=0., bigmap=False, ax=None, show=True, draw_countries=True, **kwargs):
        from mpl_toolkits.basemap import Basemap
        import matplotlib.pyplot as plt
        # lon_0, lat_0 are the center point of the projection.
        # resolution = 'l' means use low resolution coastlines.
        if bigmap:
            m = Basemap(projection='hammer', lon_0=lon, resolution='c', ax=ax)
        else:
            m = Basemap(projection='ortho', lon_0=lon, lat_0=lat, resolution='l', ax=ax)
        #m = Basemap(width=width,height=width,projection='aeqd',
        #    lat_0=lat,lon_0=lon, resolution='l', ax=ax)
        if draw_countries:
            m.drawcountries(linewidth=0.5)
        m.drawcoastlines()
        m.drawmapboundary(fill_color='white')
        #m.drawmapboundary()
        m.fillcontinents(color='#D3D3D3', lake_color='#D3D3D3', zorder=0)  # light gray
        # draw parallels and meridians.
        m.drawparallels(np.arange(-90., 120., 30.))
        m.drawmeridians(np.arange(0., 390., 30.))
        # , lat=0, lon=0, bigmap=False,circles=(30, 90), circles_around=None, lines=None,
        self.plot_(m, lat=lat, lon=lon, bigmap=bigmap, **kwargs)

        #plt.title('Full Disk Orthographic Projection')

        def on_press(event):
            global lat_press, lon_press
            lon_press, lat_press = m(event.xdata, event.ydata, 'inverse')
            print 'position press  lat: %.2f  lon: %.2f' % (lat_press, lon_press)

        def on_release(event):
            lon_release, lat_release = m(event.xdata, event.ydata, 'inverse')
            dist_km = gps2DistAzimuth(lat_press, lon_press, lat_release,
                                      lon_release)[0] / 1000.
            dist_degree = gps2DistDegree(lat_press, lon_press, lat_release,
                                         lon_release)
            if dist_km > 0.1:
                print 'position release lat: %.2f  lon: %.2f' % (lat_release, lon_release)
                print 'Distance between points: %.2f degree or %.2f km' % (dist_degree, dist_km)
        plt.connect('button_press_event', on_press)
        plt.connect('button_release_event', on_release)
        if show:
            plt.show()
#        plt.savefig('ortho_full.png')
#### END DEPRECIATED

import re  #@Reimport
import StringIO
from obspy.core.quakeml import readQuakeML  #@UnresolvedImport
ev_expr = re.compile(r'<event.*>')
mag_expr = re.compile(r'(<magnitude publicID=")(.{0,50})(#.*?">)'
                      '(.*?)(magnitude)(.*?)(magnitude)(.*?</magnitude>)',
                      flags=re.S)
ev_para_expr = re.compile(r'<[eE]ventParameters>')
def readSeisComPEventXML0_6(filename):
    """
    Reads a single SeisComP event XML V0.6 file and returns a ObsPy Catalog object.

    This fixes the following differences to QUAKEML1.1:
    - EventProperties is replaced by eventProperties
    - event start tag is moved behind eventProperties start tag
    - magnitude nodes are moved out of origin nodes
    - an origin_id tag is inserted into every magnitude node
    - the magnitude tag giving the value and uncertainty of the magnitude is
      replaced by a mag tag
    """
    # realy dirty hack
    spans = [0]
    xml = open(filename, 'rt').read()
    m_event = ev_expr.search(xml)
    ins = '<eventParameters>' + m_event.group()
    for m in mag_expr.finditer(xml):
        ins += m.expand(r'\1\2\3<originID>\2</originID>\4mag\6mag\8')
        spans.extend(list(m.span()))
    spans.extend(list(m_event.span()) + [-1])
    xml = ''.join([xml[spans[2 * i]:spans[2 * i + 1]] for i in range(len(spans) // 2)])
    xml = ev_para_expr.sub(ins, xml)
    xml = xml.replace('EventParameters', 'eventParameters')
    temp = StringIO.StringIO(xml)
    return readQuakeML(temp)

regex_GEOFON = r"""
    # regex for text pasted from geofon eartquake bulletin
    ^(?P<time>[\d\s:-]{19,20})\s+
    (?P<magnitude>\d.\d)\s+
    (?P<latitude>\d+.\d+)°(?P<latitude_sign>[NS])\s+
    (?P<longitude>\d+.\d+)°(?P<longitude_sign>[EW])\s+
    (?P<depth>\d+)\s+
    (?P<AM>[AM])\s+
    (?P<flinn>[\w\s]*)$
"""
def read_regex(event_file, regex=regex_GEOFON, creation_info='GEOFON'):
    """
    Read events from event_file with the help of given regular expression.
    """
    with open(event_file, 'r') as f:
        filedata = f.read()
    event_matches = re.finditer(regex, filedata, re.VERBOSE + re.MULTILINE)
    list_ = [i.groupdict() for i in event_matches]
    events = []
    for event in list_:
        # convert numbers to float and int types
        for key, item in event.iteritems():
            if util.isint(item):
                event[key] = int(item)
            elif util.isfloat(item):
                event[key] = float(item)
            else:
                event[key] = item.strip()
        if 'latitude_sign' in event and event['latitude_sign'] == 'S':
            event['latitude'] = -event['latitude']
        if 'longitude_sign' in event and event['longitude_sign'] == 'W':
            event['longitude'] = -event['longitude']
        if 'AM' in event:
            ci = creation_info + (' automatic' if event['AM'] == 'A' else ' manual')
        else:
            ci = creation_info
        ev = Event(event_type='earthquake', creation_info=ci,
                   origins=[Origin(time=UTC(event['time']),
                                   latitude=event['latitude'],
                                   longitude=event['longitude'],
                                   depth=event['depth'])],
                   magnitudes=[Magnitude(mag=event['magnitude'],
                                         magnitude_type='M')],
                   event_descriptions=[EventDescription(event['flinn'],
                                                        'flinn-engdahl region')]
                                       if 'flinn' in event else None
                   )
        events.append(ev)
    events.sort(key=lambda x: x.origins[0].time)
    return Catalog(events)


def getinfo(cat, key):
    key = ('mag' if key == 'magnitude' else
           'standard_error' if key == 'rms' else
           key)
    key1 = 'magnitudes' if key == 'mag' else 'origins'
    key2 = ('quality' if key in ('standard_error', 'azimuthal_gap',
                                 'used_station_count', 'used_phase_count') else
            None)
    if key2:
        return [getattr(getattr(getattr(event, key1)[0], key2), key)
                for event in cat]
    else:
        return [getattr(getattr(event, key1)[0], key) for event in cat]

