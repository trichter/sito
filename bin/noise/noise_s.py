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
    #stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 PB15 PB16 HMBCX MNMCX PATCX PSGCX LVC TAIQ'
    stations = 'PB01 PB02 PB03 PB04 PB05'
    stations2 = None

    components = 'Z'
    # TOcopilla earthquake:
    #t_Toco=UTC('2007-11-14 15:14:00')
    t1 = UTC('2006-01-01')
    #t2 = UTC('2011-09-01')
    #t1 = UTC('2007-01-01')
    #t2 = UTC('2009-01-01')
    t2 = UTC('2012-01-01')

    shift = 500
    correlations = get_correlations(stations, components, stations2, only_cross=True)

#    method = 'FINAL_filter0.005-5_1bit_whitening_2011+2012'
#    method = 'filter0.01-1_1bit_whitening0.01'
#    method = 'filter0.005_rm20'
#    method = 'filter0.005_1bit'
    method = 'filter0.01-1_water_env2_whitening_1bit_fft'


    data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=False)
    data.setXLogger('_' + method)

#    pool = Pool()
#    prepare(data, stations.split(), t1, t2, component=components,
#            filter=(0.005, 5, 2, True), downsample=20,
##            filter=(1, 10), downsample=None,
##            eventremoval=None, #'waterlevel_env2', param_removal=(10, 0),
#            whitening=True,
#            use_this_filter_after_whitening=(0.005, 5, 2),
#            normalize='1bit', param_norm=None,
#            pool=pool)
#    noisexcorrf(data, correlations, t1, t2, shift, pool=pool)
#    pool.close()
#    pool.join()
#
#    stack(data, correlations, dt=10 * 24 * 3600, shift=5 * 24 * 3600)
#    stack(data, correlations, dt= -1)

    t1p, t2p = t1, t2
#    t1p, t2p = None, None

    filters = None
    filters = getFilters((0.025, 0.05, 0.1, 0.25, 0.5, 1))

    plotXcorrs(data, correlations, t1=t1p, t2=t2p, start=None, end=None, plot_overview=True, plot_years=False, use_dlognorm=False,
                      plot_stack=True, plot_psd=False, add_to_title='', downsample=None, filters=filters, filter_now=False)

    plotXcorrs(data, correlations, t1=t1p, t2=t2p, start=None, end=None, plot_overview=True, plot_years=False, use_dlognorm=False,
                      plot_stack=True, plot_psd=False, add_to_title='', downsample=None, stack=('10days', 'day'), filters=filters, filter_now=False)

    plotXcorrs(data, correlations, t1=t1p, t2=t2p, start=None, end=None, plot_overview=True, plot_years=False, use_dlognorm=False,
                      plot_stack=True, plot_psd=False, add_to_title='', downsample=None, stack=('50days', '5days'), filters=filters, filter_now=False)


#    plotXcorrs(data, correlations, t1, t2, start=0, end=25, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=True, add_to_title=method, downsample=None)
#    plotXcorrs(data, correlations, t1, t2, start=0, end=25, plot_overview=True, plot_years=False, use_dlognorm=True,
#                      plot_stack=True, plot_psd=True, add_to_file='_dlognorm.png', add_to_title=method, downsample=None)


#    plotXcorrs(data, correlations, t1, t2, start=0, end=25, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=True, add_to_title=method, downsample=None, stack=10 * 24 * 3600)
    #plotXStacks(data, correlations, start= -200, end=200, vmax=0.01, dt=10, add_to_file='_scale=0.01_200s.png')
    #plotXStacks(data, correlations, start= -200, end=200, vmax=0.05, dt= -1, add_to_file='_scale=0.05_200s.png')
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

#
#    method = 'rm5_filter0.1-1'
#    data = IPOC(xcorr_append='/tests/' + method, use_local_LVC=True)
#    data.setXLogger('_' + method)
#    prepare(data, stations.split(' '), t1, t2, filter=(0.1, 1.), downsample=10,
#            component=components, normalize='runningmean', norm_param=5 * 10 + 1,
#            use_floating_stream=True)
#    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
#    plotXcorrs(data, correlations, t1, t2, plot_overview=False, plot_stack=True, plot_psd=True, add_to_title=method)
#
#
#    method = 'rm50_filter0.01'
#    data = IPOC(xcorr_append='/tests/' + method, use_local_LVC=True)
#    data.setXLogger('_' + method)
#    prepare(data, stations.split(' '), t1, t2, filter=(0.01, None), downsample=None,
#            component=components, normalize='runningmean', norm_param=50 * 100 + 1,
#            use_floating_stream=True)
#    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
#    plotXcorrs(data, correlations, t1, t2, plot_overview=False, plot_stack=True, plot_psd=True, add_to_title=method)
#
#
#    method = 'rm0.25_filter2'
#    data = IPOC(xcorr_append='/tests/' + method, use_local_LVC=True)
#    data.setXLogger('_' + method)
#    prepare(data, stations.split(' '), t1, t2, filter=(2, None), downsample=None,
#            component=components, normalize='runningmean', norm_param=100 // 4 + 1,
#            use_floating_stream=True)
#    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
#    plotXcorrs(data, correlations, t1, t2, plot_overview=False, plot_stack=True, plot_psd=True, add_to_title=method)


if __name__ == '__main__':
    main()
