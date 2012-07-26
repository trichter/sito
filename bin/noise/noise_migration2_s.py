#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.noisexcorr import setHIDist
from sito.noise_migration import migrate, migrate2, migrate3
from sito.util.imaging import xcorr_cmap, equi
from obspy.signal.util import utlGeoKm
from obspy.signal.filter import envelope

from sito.stream import read
from sito.stations import IPOCStations
from matplotlib import pylab as plt
import numpy as np

def main():
    ipoc = IPOCStations()
    lats = -25, 17
    lons = -74, -66

    lats = np.linspace(-25, -17, 1500)
    lons = np.linspace(-74, -66, 1500)
    velocity = 3
#    print 'x, y=', utlGeoKm(lons[0], lats[0], lons[-1], lats[-1])
#    print 'dx, dy=', utlGeoKm(lons[0], lats[0], lons[-1], lats[-1]) / np.array((len(lons), len(lats)))
#    print 'wave in km = ', (ms[0].stats.endtime - ms[0].stats.starttime) * 3.0
#    print 'dx of wave=', ms[0].stats.delta * 3.0


    from mpl_toolkits.basemap import Basemap
    m = Basemap(llcrnrlon=lons[0], llcrnrlat=lats[0], urcrnrlon=lons[-1], urcrnrlat=lats[-1],
                lat_0= -21, lon_0= -70, projection='stere', resolution='l', ax=None)
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines()
#    m.drawmapboundary(fill_color='aqua')
    m.drawmapboundary()
    m.fillcontinents(color='coral', lake_color='aqua', zorder=0)
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90., 120., 1.))
    m.drawmeridians(np.arange(0., 390., 1.))
    ipoc.plot(m)
    #wiggles = dict(PB02=66, PB04=60, PB07=60)
    wiggles2 = dict(PB02=100, PB03=130, PB04=143, PB07=105, PB09=82, PB10=100, PB14=90, PATCX=125)

#    for st, rad in wiggles.iteritems():
#        m.plot(*equi(m, ipoc[st].latitude, ipoc[st].longitude, 0.5 * rad * velocity, indeg=False))
    for st, rad in wiggles2.iteritems():
        m.plot(*equi(m, ipoc[st].latitude, ipoc[st].longitude, 0.5 * rad * velocity, indeg=False))


    plt.show()
    #plt.savefig('/home/richter/Results/IPOC/xcorr/1bit/migration/migration3_station_%s_2Hz_600x600_interference' % key)
    plt.close()

if __name__ == '__main__':
    main()
