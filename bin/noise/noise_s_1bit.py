#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.noisexcorr import (prepare, get_correlations,
                             plotXcorrs, stack_day, setHIDist, noisexcorr,
    noisexcorrf, stack)
from sito.noise_migration import migrate, migrate2, migrate3
from sito.util.imaging import xcorr_cmap
from obspy.signal.util import utlGeoKm
from obspy.signal.filter import envelope

def main():
    method = 'filter0.25-1_1bit_whitening'
    data = IPOC(xcorr_append='/' + method , use_local_LVC=True)
    data.setXLogger(method)
    stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14'
    stations2 = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 MNMCX PATCX PSGCX HMBCX'

    components = 'Z'
    t1 = UTC('2006-01-01')
    t2 = UTC('2011-09-01')
    t2 = UTC('2006-12-01')
#    t1 = UTC('2008-01-01')
#    t2 = UTC('2008-01-10')

    shift = 2000

    # test
    #data.setXLogger('_1bit_test2')
    #stations = 'PB01 PB14 LVC'
    #t1 = UTC('2011-01-01')
    #t2 = UTC('2011-01-20')

    # test2
    #data.setXLogger('_1bit_test2')
    #stations = 'PB01 PB14 LVC'
    #t1 = UTC('2011-01-01')
    #t2 = UTC('2011-01-20')

    period = 'day'
#    prepare(data, stations.split(), t1, t2, filter=(0.25, 1.), downsample=10,
#            component=components, normalize='1bit', norm_param=None, whitening=True,
#            freq_domain=True)
    correlations = get_correlations(stations, components, only_cross=True)
    #print correlations
#    noisexcorr(data, correlations, t1, t2, shift_sec=shift, period=period, freq_domain=True)
#    stack(data, correlations, onefile=False)
    correlations = (('PB03Z', 'PB04Z'),)
    plotXcorrs(data, correlations, t1, t2, start= -50, end=50, plot_overview=True, filter=None, stack_lim=(-0.1, 0.1), plot_years=False, #vmax=1,
               plot_stack=True, plot_psd=True, downsample=2, add_to_title=method + '_50s', add_to_file='_50s.png', show=True, use_dlognorm=False)


#    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
#    stack_day(data, correlations, onefile=True)

#    data.x_plot_day = data.x_res + '/plotsxlim300s/%s_day_%s'
#    plotXcorrs(data, correlations, t1, t2, start= -300, end=300, stack_lim=(-0.05, 0.05),
#               downsample=1, filter=None, plot_years=False)








####### migrate stuff
#    from sito.stream import read
#    from sito.stations import IPOCStations
#    from matplotlib import pylab as plt
#    import numpy as np
#    ipoc = IPOCStations()
#
#    ms = read('/home/richter/Results/IPOC/xcorr/1bit/stack/all_stack_-1.QHD')
#    setHIDist(ms, ipoc)
#    for tr in ms:
#        st = tr.stats.station
#        if 'PB12' in st or 'PAT' in st or 'HMBC' in st or 'LVC' in st:
#            ms.remove(tr)
##    ms.plot_(-200, 200, relative='middle', absolutescale=10, annotate=True)
##    ms.plotXcorrVsDist(-300, 300, absolutescale=200)
#    ms.addXcorrSides()
##    ms.plotXcorrVsDist(-300, relative='starttime', absolutescale=500)
##    plt.show()
#
#    #ms.filter2(0.05, 0.5)
#    ms = ms.select('st.dist<455')
##    for tr in ms:
##        tr.data = envelope(tr.data)
#    ms.downsample2(2)
##    ms.trim2(100 + np.array(ms.getHI('dist')) / 3., None)
##    ms.normalize()
##    ms.plotXcorrVsDist(0, 300, relative='starttime')
##    ms.plot_(absolutescale=1)
##    ms = ms[:1]
##    ms = ms[:1]
##    ms.plot_()
##    plt.show()
##    return
#    print ms
#
#    lats = np.linspace(-25, -17, 600)
#    lons = np.linspace(-74, -66, 600)
#    print 'x, y=', utlGeoKm(lons[0], lats[0], lons[-1], lats[-1])
#    print 'dx, dy=', utlGeoKm(lons[0], lats[0], lons[-1], lats[-1]) / np.array((len(lons), len(lats)))
#    print 'wave in km = ', (ms[0].stats.endtime - ms[0].stats.starttime) * 3.0
#    print 'dx of wave=', ms[0].stats.delta * 3.0
#    t1 = UTC()
#    data = migrate3(ms, ipoc, lats, lons, 3.0, skip=100.)
#    t2 = UTC()
#    print 'used time:', t2 - t1
#
#    np.save('/home/richter/Results/IPOC/xcorr/1bit/migration_data_env', data)
#
#    from mpl_toolkits.basemap import Basemap
#    m = Basemap(llcrnrlon=lons[0], llcrnrlat=lats[0], urcrnrlon=lons[-1], urcrnrlat=lats[-1],
#                lat_0= -21, lon_0= -70, projection='stere', resolution='l', ax=None)
#    m.drawcountries(linewidth=0.5)
#    m.drawcoastlines()
##    m.drawmapboundary(fill_color='aqua')
#    m.drawmapboundary()
#    m.fillcontinents(color='coral', lake_color='aqua', zorder=0)
#    # draw parallels and meridians.
#    m.drawparallels(np.arange(-90., 120., 1.))
#    m.drawmeridians(np.arange(0., 390., 1.))
#    im = m.imshow(data.transpose(), interpolation="nearest", cmap=xcorr_cmap, alpha=0.8)
#    plt.colorbar(im)
#    ipoc.plot(m)
#    plt.show()



if __name__ == '__main__':
    main()
