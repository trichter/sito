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
    stations2 = None


    components = 'Z'
    # TOcopilla earthquake: 2007-11-14 15:14
    t1 = UTC('2006-02-01')
    t2 = UTC('2012-10-01')

    shift = 100
    correlations = get_correlations(stations, components, stations2, only_auto=True)

    method = 'FINAL_filter1-3_1bit_auto'

    data = IPOC(xcorr_append='/' + method, use_local_LVC=False)
    data.setXLogger('_' + method)

    pool = Pool()
    prepare(data, stations.split(), t1, t2, component=components,
            filter=(1, 3, 2, True), downsample=20,
            eventremoval='waterlevel_env2', param_removal=(10, 0),
            whitening=False,
            normalize='1bit', param_norm=None,
            pool=pool)
    noisexcorrf(data, correlations, t1, t2, shift, pool=pool)
    pool.close()
    pool.join()

#    plotXcorrs(data, correlations, t1, t2, start=None, end=None, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=False, add_to_title='', downsample=None)
    plotXcorrs(data, correlations, t1, t2, start=0, end=20, plot_overview=True, plot_years=False, use_dlognorm=False,
                      plot_stack=True, plot_psd=False, add_to_title='', downsample=None, ext='_hg.png', vmax=0.1)
#    stack(data, correlations, dt= -1)
#    stack(data, correlations, dt=10 * 24 * 3600, shift=2 * 24 * 3600)
#    plotXcorrs(data, correlations, t1=None, t2=None, start=None, end=None, plot_overview=True, plot_years=False, use_dlognorm=False,
#               plot_stack=True, plot_psd=False, add_to_title='', downsample=None,
#               stack=('10days', '2days'))
#    plotXcorrs(data, correlations, t1=None, t2=None, start=0, end=20, plot_overview=True, plot_years=False, use_dlognorm=False,
#               plot_stack=True, plot_psd=False, add_to_title='', downsample=None,
#               stack=('10days', '2days'), ext='_hg.png', vmax=0.1)
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
