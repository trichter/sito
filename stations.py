# by TR

from obspy.core import AttribDict
from obspy.core.util import gps2DistAzimuth
from sito.util import gps2DistDegree
import logging
import re
log = logging.getLogger(__name__)

class Stations(AttribDict):
    """
    Class for holding station data.
    """

    regex = """ # SH-combatible
    \s*(?P<name>\w+)\s*
    lat:\s*(?P<latitude>-?\+?[\d.]+)\s*
    lon:\s*(?P<longitude>-?\+?[\d.]+)\s*
    (?:elevation:\s*(?P<elevation>-?\+?[\d.]+)\s*|)
    (?:info:\s*(?P<info>.*)|).*$"""

    example = """
    PB01  lat: -21.0432  lon: -69.4874  elevation: 900  info: Plate_Boundary_Station_PB01,_Chile
    PB02  lat: -21.3197  lon: -69.8960  elevation: 1015 info: Plate_Boundary_Station_PB02,_Chile
    PB03  lat: -22.0476  lon: -69.7533  elevation: 1460 info: Plate_Boundary_Station_PB03,_Chile
    GRC3 lat:+48.8901739 lon: +11.5858216 elevation: 438.0 array:01 xrel: 26.722525 yrel: -89.032738 name:Graefenberg,_F.R.G.
    GEA3  lat:48.83501053   lon:13.70003414 array:05 xrel:-0.13108 yrel:-1.12033"""

    format_header = 'Stations:\nname      lat       lon     elev    info\n'
    format_file = '{0:5}  {latitude:-8.3f}  {longitude:-8.3f}  {elevation:6.1f}  {info}\n'

    @classmethod
    def read(cls, filename):
        """
        Read station file_ and return instance of Stations.

        Format has to be like in Stations.example.
        """
        with open(filename, 'r') as file_:
            filedata = file_.read()
        # Create an iterator over matches in Stations file_
        st_matches = re.finditer(cls.regex, filedata, re.VERBOSE + re.MULTILINE)
        # Create a list of dictionaries of PDE data
        st_list = [i.groupdict() for i in st_matches]
        st_dic = {}
        for i in st_list:
            st_dic[i['name']] = AttribDict({'latitude':float(i['latitude']), 'longitude':float(i['longitude']), 'info':i['info']})
            try: st_dic[i['name']]['elevation'] = float(i['elevation'])
            except TypeError: pass
        log.info('Read station information of stations %s from file_ %s' % (' '.join(sorted(st_dic.keys())), filename))
        return cls(st_dic)

    def __repr__(self):
        return 'Stations({0})'.format(super(Stations, self).__repr__())
    def __str__(self):
        return self.write(infile=False, header=True)

    def write(self, filename=None, infile=True, header=False):
        """
        Write station information in a filename or return string with information.

        :param filen: filename
        :param infile: write into filename or return string
        :param header: insert a header
        :return: string if infile is False
        """
        str_list = []
        if header: str_list.append(self.format_header)
        for station in self: str_list.append(self.format_file.format(station, ** self[station]))
        if infile:
            with open(filename, 'w') as f:
                f.write(''.join(str_list))
            log.info('Write station data to ' + filename)
        else:
            return ''.join(str_list)

    def pick(self, str_keys, replace=True):
        """
        Pick stations.

        :param str_keys: string with station names eg. 'PB01 PB02'
        :param replace: if True the data in the station list is overwritten
        :return: instance of station data
        """
        newdata = {}
        if not str_keys.startswith('-'):
            log.info('Pick stations ' + str_keys)
            for key in str_keys.split(): newdata[key] = self[key]
        else:
            log.info('Delete stations ' + str_keys[1:] + ' from list')
            for key in self.keys():
                if key not in str_keys[1:].split(): newdata[key] = self[key]
        if replace:
            self.clear()
            for key in newdata.keys(): self[key] = newdata[key]
        return self.__class__(newdata)

    def getNames(self):
        """
        Return station names in one string.
        """
        return ' '.join(self.keys())

    def dist(self, st1, st2, indeg=False):
        dist_deg = gps2DistDegree(self[st1].latitude, self[st1].longitude,
                                 self[st2].latitude, self[st2].longitude)
        dist_km = gps2DistAzimuth(self[st1].latitude, self[st1].longitude, self[st2].latitude, self[st2].longitude)[0] / 1.e3
        if indeg is True:
            return dist_deg
        elif indeg is False:
            return dist_km
        else:
            return dist_km, dist_deg

    def plot(self, basemap, annotate=True, lsize='small', **kwargs_in):
        kwargs = dict(marker='o')
        kwargs.update(kwargs_in)
        for key, val in self.items():
            x, y = basemap(val.longitude, val.latitude)
            basemap.plot((x,), (y,), **kwargs)
            if annotate:
                import matplotlib.pylab as plt
                plt.annotate(key, (x, y), xytext=(3, 3),
                             textcoords='offset points', size=lsize)

def IPOCStations():
    return Stations.read('/home/richter/Data/stations_ipoc.txt')
def ParkfieldStations():
    return Stations.read('/home/richter/Data/stations.txt')

if __name__ == '__main__':
    pass
