#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.noisexcorr import setHIDist
from sito.noise_migration import migrate, migrate2, migrate3
from sito.util.imaging import xcorr_cmap
from obspy.signal.util import utlGeoKm
from obspy.signal.filter import envelope

from sito.stream import read
from sito.stations import IPOCStations
from matplotlib import pylab as plt
import numpy as np

def main():
    ipoc = IPOCStations()
    ms = read('/home/richter/Results/IPOC/xcorr/1bit/stack/all_stack_-1.QHD')
    setHIDist(ms, ipoc)
#    for tr in ms:
#        st = tr.stats.station
#        if 'PB12' in st or 'PAT' in st or 'HMBC' in st or 'LVC' in st:
#            ms.remove(tr)
    print (list(set(ms[0].stats.station.split('-'))))[0]
    ms = ms.select("(st.station.split('-'))[0][:-1] in 'PB03 PB04 PB07 PB09 PB10 PB14 PATCX'.split()")
    ms = ms.select(autocorr=True)
#    ms.plot_(-200, 200, relative='middle', absolutescale=10, annotate=True)
#    ms.plotXcorrVsDist(-300, 300, absolutescale=200)
#    ms.plotXcorrVsDist(-300, relative='starttime', absolutescale=500)
#    plt.show()
#    ms.filter2(0.05, 0.5)

    #ms = ms.select('st.dist<455')
    ms.addXcorrSides()
#    for tr in ms:
#        tr.data = envelope(tr.data)   
#    ms.normalize()    
    ms.downsample2(5)
#    ms.trim2(100 + np.array(ms.getHI('dist')) / 3., None)
#    ms.plotXcorrVsDist(0, 300, relative='starttime')
#    ms.plot_(absolutescale=1)
#    ms = ms[:1]
#    ms = ms[:1]
#    ms.plot_()
#    plt.show()
#    return
    print ms

    lats = np.linspace(-25, -17, 1500)
    lons = np.linspace(-74, -66, 1500)
#    print 'x, y=', utlGeoKm(lons[0], lats[0], lons[-1], lats[-1])
#    print 'dx, dy=', utlGeoKm(lons[0], lats[0], lons[-1], lats[-1]) / np.array((len(lons), len(lats)))
#    print 'wave in km = ', (ms[0].stats.endtime - ms[0].stats.starttime) * 3.0
#    print 'dx of wave=', ms[0].stats.delta * 3.0
#
##    for tr in ms:
##        tr.data = envelope(tr.data)
##    for key in ipoc:
##        print key
##        ms_sel = ms.select('%r in st.station' % key)
##        print ms_sel
##        if len(ms_sel) > 0:
##    t1 = UTC()
##    data = migrate3(ms, ipoc, lats, lons, 3.0, skip=50.)
##    t2 = UTC()
##    print 'used time:', t2 - t1
##    np.save('/home/richter/Results/IPOC/xcorr/1bit/migration/data_10Hz_3000x3000_interference', data)
##    del data
#
#    t1 = UTC()
#    data = migrate3(ms, ipoc, lats, lons, 3.0, skip=50.)
#    t2 = UTC()
#    print 'used time:', t2 - t1
#    np.save('/home/richter/Results/IPOC/xcorr/1bit/migration/data_some_stations_method3_5Hz_1500x1500_interference.npy', data)


#    data = np.load('/home/richter/Results/IPOC/xcorr/1bit/migration/data_some_stations_method3_5Hz_1500x1500_interference.npy')

#    from mpl_toolkits.basemap import Basemap
#    for key in ipoc:
#        try:
#            data = np.load('/home/richter/Results/IPOC/xcorr/1bit/migration/data_station_%s_2Hz_600x600_interference.npy' % key)
#        except:
#            continue
#        print key
#        m = Basemap(llcrnrlon=lons[0], llcrnrlat=lats[0], urcrnrlon=lons[-1], urcrnrlat=lats[-1],
#                    lat_0= -21, lon_0= -70, projection='stere', resolution='l', ax=None)
#        m.drawcountries(linewidth=0.5)
#        m.drawcoastlines()
#    #    m.drawmapboundary(fill_color='aqua')
#        m.drawmapboundary()
#        m.fillcontinents(color='coral', lake_color='aqua', zorder=0)
#        # draw parallels and meridians.
#        m.drawparallels(np.arange(-90., 120., 1.))
#        m.drawmeridians(np.arange(0., 390., 1.))
#        im = m.imshow(data.transpose(), interpolation="nearest", cmap=xcorr_cmap, alpha=0.8)
#        plt.colorbar(im)
#        ipoc.plot(m)
#        plt.savefig('/home/richter/Results/IPOC/xcorr/1bit/migration/migration3_station_%s_2Hz_600x600_interference' % key)
#        plt.close()

if __name__ == '__main__':
    main()
