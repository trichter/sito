#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.noisexcorr import (prepare, get_correlations,
                             plotXcorrs, noisexcorrf, stack, getFilters)
from sito import util
import matplotlib.pyplot as plt
from sito.stream import read
from multiprocessing import Pool
import time

def main():
    stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 HMBCX MNMCX PATCX PSGCX'

    components = 'Z'
    # TOcopilla earthquake: 2007-11-14 15:14
    t1 = UTC('2006-07-01')
    t2 = UTC('2008-12-31')

    shift = 500
    correlations = get_correlations(stations, components)

    method = 'FINAL_filter0.005-10_1bit_Tocopilla'

    data = IPOC(xcorr_append='/' + method, use_local_LVC=False)
    data.setXLogger('_' + method)
    pool = Pool()
    prepare(data, stations.split(), t1, t2, component=components,
            filter=(0.005, 10, 2, True), downsample=20,
            whitening=False,
            normalize='1bit', param_norm=None,
            pool=pool)
    noisexcorrf(data, correlations, t1, t2, shift, pool=pool)
    pool.close()
    pool.join()
    stack(data, correlations, dt=10 * 24 * 3600, shift=5 * 24 * 3600)
    stack(data, correlations, dt= -1)

    filters = None
    #filters = getFilters((0.005, 0.01, 0.1, 1, 5, 10), zerophase=True, corners=2)
#    plotXcorrs(data, correlations, t1, t2, start=None, end=None, filters=filters, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=True, add_to_title='', downsample=None)
    plotXcorrs(data, correlations, t1=None, t2=None, start=None, end=None, filters=filters, plot_overview=True, plot_years=False, use_dlognorm=False,
                      plot_stack=True, plot_psd=True, add_to_title='', downsample=None, stack=('10days', '5days'))


#    ms = read(data.x_day % ('PB03Z', '*') + '.QHD')
#    tr = ms.calculate('mean')
#    tr.plot()
#    ipshell()

#    util.checkDir(data.getPlotX(('', ''), t1))
    #for correlation in correlations:
#        stations = correlation[0][:-1], correlation[1][:-1]
#        dist = data.stations.dist(*stations)
##        if dist >= 120:
##            t = (dist // 100) * 50 + 50
##        else:
##            t = 70
#        t = 200
#        stream = data.readDayXcorr(correlation, t1, t2)
#        if len(stream) > 0:
#            stream.plotXcorr(-t, t, imshow=True, vmax=0.01, vmin_rel='vmax',
#                                 fig=plt.figure(figsize=(8.267, 11.693)),
#                                 figtitle='station ' + method + ' around Tocopilla event',
#                                 dateformatter='%y-%m-%d', show=False,
#                                 save=data.getPlotX(correlation, 'Tocopilla_0.01.png'),
#                                 stack_lim=None)

if __name__ == '__main__':
    main()
