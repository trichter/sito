#!/usr/bin/env python
# by TR
#@PydevCodeAnalysisIgnore

from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.noisexcorr import (prepare, xcorr_day, get_correlations,
                             plotXcorrs, stack_day, plotXStacks, noisexcorr)
from sito import ipshell, util
import matplotlib.pyplot as plt
from sito.stream import read


def main():
    stations = 'PB03 PB04'

    stations2 = None

    components = 'Z'
    # TOcopilla earthquake: 2007-11-14 15:14
    t1 = UTC('2006-01-01')
    t1 = UTC('2007-11-10')
    t2 = UTC()
    t2 = UTC('2007-11-20')
    shift = 500
    correlations = get_correlations(stations, components, stations2, only_cross=True)


    method = 'filter4-6_water_env2_whitening_1bit'
#    method = 'filter0.01-1_1bit_whitening0.01'
#    method = 'filter0.005_rm20'
#    method = 'filter0.005_1bit'

    data = IPOC(xcorr_append='/Tocopilla/tests/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)
#    prepare(data, stations.split(), t1, t2, component=components,
#            filter=(4, 6), downsample=None,
#            eventremoval='waterlevel_env2', param_removal=(10, 0),
#            #whitening=True,
#            normalize='1bit', param_norm=None)
    correlations = get_correlations(stations, components, stations2, only_auto=True)
#    noisexcorr(data, correlations, t1, t2, shift)
    plotXcorrs(data, correlations, t1, t2, start= -150, end=150, plot_overview=True, plot_years=False, use_dlognorm=True,
                      plot_stack=True, plot_psd=True, add_to_title=method, show=True)

#    stack_day(data, correlations, dt= -1)
#    stack_day(data, correlations, dt=10)

#    plotXStacks(data, correlations, start= -200, end=200, vmax=0.01, dt=10, add_to_file='_scale=0.01_200s.png')
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
