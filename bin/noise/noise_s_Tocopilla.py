#!/usr/bin/env python
# by TR
#@PydevCodeAnalysisIgnore

from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.noisexcorr import (prepare, xcorr_day, get_correlations,
                             plotXcorrs, noisexcorr)
from sito import ipshell, util
import matplotlib.pyplot as plt


def main():

    stations = 'PB01 PB02 PB03 PB04 PB04 PB05 PB06'
    stations2 = 'PB01 PB02 PB03 PB04 PB04 PB05 PB06 PB07 PB08 HMBCX PATCX'

    components = 'Z'
    # TOcopilla earthquake: 2007-11-14 15:14
    t1 = UTC('2007-09-01')
    t2 = UTC('2008-01-31')

    shift = 200
    correlations = get_correlations(stations, components, stations2)


#    method = 'filter0.01-1_1bit'
#    method = 'filter0.01-1_1bit_whitening0.01'
#    method = 'filter2-20_1bit'
#    method = 'filter0.005_1bit'
#    period = 'day'
#
#    data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
#    data.setXLogger('_' + method)
#    prepare(data, stations.split(), t1, t2, filter=(0.005, 1.), downsample=5, whitening=True,
#            component=components, normalize='1bit', norm_param=None)
#    noisexcorr(data, correlations, t1, t2, shift_sec=shift, period=period)

#    correlations = (('PB03Z', 'PB04Z'),)
#    data.x_plot_day = data.x_res + '/plots2/%s_day_%s'
#    plotXcorrs(data, correlations, t1, t2, start=9, end=15, plot_overview=True, filter=(2, None, 2, True), stack_lim=(-0.01, 0.01), downsample=None, plot_years=False,
#                      plot_stack=True, plot_psd=True, add_to_title=method + '_filter2_9-15', add_to_file='_filter2_9-15.png', show=False)




    method = 'filter4-6_1bit'
    period = 'day'
    data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)

#    prepare(data, stations.split(), t1, t2, filter=(4, 6), downsample=None, whitening=None,
#            component=components, normalize='1bit', norm_param=None)
#    noisexcorr(data, correlations, t1, t2, shift_sec=shift, period=period)


    plotXcorrs(data, correlations, t1, t2, start= -50, end=50, plot_overview=True, filter=None, stack_lim=(-0.1, 0.1), plot_years=False,
               plot_stack=True, plot_psd=False, add_to_title=method + '_wodlognorm_50s', add_to_file='_wodlognorm_50s.png', show=True, landscape=True, use_dlognorm=False)
#    correlation = ('PB03Z', 'PB03Z')
#    add_to_title = ''#_againfilter_zoom1'
#    save = False
#    stream = data.readX(correlation, t1, t2, period=period)
#    stream.plotXcorr(0, 100, imshow=True, use_dlognorm=True, filter=(2, 20),
#                                 fig=plt.figure(figsize=figsize),
#                                 figtitle='station ' + add_to_title,
#                                 dateformatter='%y-%m-%d %Hh%M', save=save, show=True #dateformatter=' % y % b'
#                                 )
#    plt.show()




#    correlations = (('PB03Z', 'PB04Z'),)
#    plotXcorrs(data, correlations, UTC('2007-03-01'), UTC('2008-05-01'), start= -200, end=200, vmax=1, plot_overview=True, filter=None, stack_lim=(-0.1, 0.1), plot_years=False,
#                      plot_stack=True, plot_psd=False, downsample=10, add_to_title=method + '_100s', add_to_file='_100s.png', show=False, landscape=True)


#    method = 'filter2-20_rm20'
##    method = 'filter0.005_1bit'
#
#    data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
#    data.setXLogger('_' + method)
#    prepare(data, stations.split(), t1, t2, filter=(2., 20), downsample=None, whitening=None,
#            component=components, normalize='runningmean', norm_param=20 * 100 + 1)
#    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
#    plotXcorrs(data, correlations, t1, t2, plot_overview=True,
#                      plot_stack=True, plot_psd=True, add_to_title=method)


#    util.checkDir(data.getPlotX(('', ''), t1))
#    for correlation in correlations:
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
