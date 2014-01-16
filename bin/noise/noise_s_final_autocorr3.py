#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.noisexcorr import (prepare, get_correlations,
                             plotXcorrs, noisexcorrf, stack)
from sito import util
import matplotlib.pyplot as plt
from sito.stream import read
from multiprocessing import Pool
import time
from sito import seismometer

def main():
    stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 HMBCX MNMCX PATCX PSGCX LVC'
    stations = 'PB09 PB10 PB11 PB12 PB13 PB14 PB15 PB16'
    stations = 'PATCX'
    stations2 = None


    components = 'Z'
    # TOcopilla earthquake: 2007-11-14 15:14
    t1 = UTC('2006-02-01')
    #t2 = UTC('2008-12-31')
    #t2 = UTC('2012-10-01')
    t2 = UTC('2011-12-31')
#    t1 = UTC('2009-05-01')
#    t2 = UTC('2009-05-03')

    shift = 100
    shift = 60
    #correlations = get_correlations(stations, components, stations2, only_auto=True)
    correlations = get_correlations(stations, components, stations2)
    print correlations

    #method = 'FINAL_filter4-6_1bit_auto'
    method = 'FINAL_filter4-6_1bit_auto_3C'
    #method = 'FINAL_filter3-5'

    data = IPOC(xcorr_append='/' + method, use_local_LVC=False)
    data.setXLogger('_' + method)

    pool = Pool()
#    prepare(data, stations.split(), t1, t2, component=components,
#            filter=(2, 4, 2, True), downsample=50,
#            eventremoval='waterlevel_env2', param_removal=(10, 0),
#            whitening=False,
#            normalize='1bit', param_norm=None,
#            pool=pool)
#    noisexcorrf(data, correlations, t1, t2, shift, period=24 * 3600, pool=pool)

#    noisexcorrf(data, correlations, t1, t2, shift, period=5 * 60, pool=pool,
#                max_preload=1000)
    pool.close()
    pool.join()

#    plotXcorrs(data, correlations, t1, t2, start=None, end=None, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=False, add_to_title='', downsample=None)
    #plt.rc('font', size=16)
    plotXcorrs(data, correlations, t1, t2, start=-20, end=20, plot_overview=True, plot_years=False, use_dlognorm=False,
                      plot_stack=True, plot_psd=False, downsample=None, ext='_hg0.02_dis.pdf', vmax=0.02,
                      add_to_title='4-6Hz', ylabel=None)

#    stack(data, correlations, dt= -1)
    #stack(data, correlations, dt=60 * 60, period=5 * 60)
#    stack(data, correlations, dt=24 * 60 * 60, period=5 * 60)
#    plotXcorrs(data, correlations, t1=None, t2=None, start=None, end=None, plot_overview=True, plot_years=False, use_dlognorm=False,
#               plot_stack=True, plot_psd=False, add_to_title='', downsample=None,
#               stack=('10days', '2days'))
#    plotXcorrs(data, correlations, t1, t2, start=0, end=20, plot_overview=True, plot_years=False, use_dlognorm=False,
#               plot_stack=True, plot_psd=False, add_to_title='', downsample=None,
#               period=60 * 5, stack=(60 * 60, None), ext='_hg.png', vmax=0.1)

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
