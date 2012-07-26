#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.noisexcorr import (prepare, get_correlations,
                             plotXcorrs, stack_day, noisexcorr,
    getFilters, filter, stack)
from sito import util
#import matplotlib.pyplot as plt
from sito.stream import read

def main():
    st1 = 'PB01 PB02 PB05'

    stations = 'PB01 PB02 PB03 PB04 PB05'

    stations2 = None

    components = 'Z'
    # TOcopilla earthquake: 2007-11-14 15:14
    t1 = UTC('2006-01-01')
    t1 = UTC('2007-01-01')
    t1 = UTC('2007-11-01')
    #t2 = UTC('2006-12-31')
#    t1 = UTC('2007-10-07')


    t2 = UTC('2008-12-31')
    t2 = UTC('2007-12-15')
    #t2 = UTC('2007-01-01')
#    t2 = UTC('2007-10-10')
#    t1 = UTC('2008-05-01')
#    t2 = UTC('2008-06-01')


#    t1 = UTC('2010-08-31')
#    t2 = UTC('2010-09-01')
    shift = 500
    stations = 'PB01 PB02 PB03 PB04 PB05'
    correlations = get_correlations(stations, components, stations2, only_cross=True)
    #correlations = (('PB03Z', 'PB04Z'),)

    method = 'filter0.01-1_water_env2_whitening_1bit'
    method = 'filter0.01-1_water_env2_whitening_1bit_fft'
#    method = 'filter0.01-1_1bit_whitening0.01'
#    method = 'filter0.005_rm20'
#    method = 'filter0.005_1bit'

    data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)
#    prepare(data, 'PB01 PB02 PB05'.split(), t1, t2, component=components,
#            filter=(0.01, 1), downsample=20,
#            eventremoval='waterlevel_env2', param_removal=(10, 0),
#            whitening=True,
#            normalize='1bit', param_norm=None,
#            freq_domain=True)
#    noisexcorr(data, correlations, t1, t2, shift, freq_domain=True)
#
#    path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft/xcorr/'
#    import glob
#    from sito.noisexcorr import removeBad
#    import os.path
#    for file in glob.glob(path + '*.QHD'):
#        if not 'filter' in os.path.basename(file):
#            print file
#            ms = read(file)
#            #ms.normalize()
#            removeBad(ms, 0.7)
#            ms.write(os.path.splitext(file)[0], 'Q')

    filters = (None,)
#    plotXcorrs(data, correlations, t1, t2, start= -100, end=100, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=True, add_to_file='', add_to_title='(waterlev whitened 1bit fft)', stack=None, filters=filters, cmax=1e-2)
#    plotXcorrs(data, correlations, t1, t2, start= -100, end=100, plot_overview=True, plot_years=False, use_dlognorm=True,
#                      plot_stack=True, plot_psd=True, add_to_file='_log', add_to_title='(waterlev whitened 1bit fft)', stack=None, filters=filters, cmax=1e-2)

#    stack(data, correlations, dt=50 * 24 * 3600, filters=filters, period=24 * 3600, shift=5 * 24 * 3600, yearfiles=True)
#    stack(data, correlations, dt=10 * 24 * 3600, filters=filters, period=24 * 3600, shift=24 * 3600, yearfiles=True)
#    stack(data, correlations, dt=30 * 24 * 3600, filters=filters, period=24 * 3600, shift=2 * 24 * 3600, yearfiles=True)
#    stack(data, correlations, dt= -1, filters=filters, period=24 * 3600, yearfiles=False)

    filters = getFilters(0.025, 0.05, 0.1, 0.25, 0.5, 1.)
#    filter(data, correlations, filters, stack=None)
#    filter(data, correlations, filters, stack=('50days', '5days'))
#    filter(data, correlations, filters, stack=('30days', '2days'))
#    filter(data, correlations, filters, stack=('10days', 'day'))
    filters = (None,) + tuple(filters)

#    plotXcorrs(data, correlations, t1, t2, start= -200, end=200, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=True, add_to_file='', add_to_title='(waterlev whitened 1bit ftt)', stack=('10days', 'day'), filters=filters, cmax=1e-2)
#    plotXcorrs(data, correlations, t1, t2, start= -200, end=200, plot_overview=True, plot_years=False, use_dlognorm=True,
#                      plot_stack=True, plot_psd=True, add_to_file='_log', add_to_title='(waterlev whitened 1bit fft)', stack=('10days', 'day'), filters=filters, cmax=1e-2)
#    plotXcorrs(data, correlations, t1, t2, start= -200, end=200, plot_overview=True, plot_years=False, use_dlognorm=True,
#                      plot_stack=True, plot_psd=True, add_to_file='_log', add_to_title='(waterlev whitened 1bit fft)', stack=('30days', '2days'), filters=filters, cmax=1e-2)
#    plotXcorrs(data, correlations, t1, t2, start= -200, end=200, plot_overview=True, plot_years=False, use_dlognorm=True,
#                      plot_stack=True, plot_psd=True, add_to_file='_log', add_to_title='(waterlev whitened 1bit fft)', stack=('50days', '5days'), filters=filters, cmax=1e-2)

#    plotXcorrs(data, correlations, t1, t2, start= -100, end=100, plot_overview=True, plot_years=False, use_dlognorm=False,
#                      plot_stack=True, plot_psd=False, add_to_file='', ext='.eps', add_to_title='(waterlev whitened 1bit ftt)', stack=None, filters=(None,), cmax=1e-2, vmax=0.03)

    method = 'filter2-4_water_env2_fft'
    stations = 'PB04'
    # no external HD!!
    data = IPOC(path='/home/richter/Results/', xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)
    prepare(data, stations.split(), t1, t2, component=components,
            filter=(2, 4), downsample=50,
            eventremoval='waterlevel_env2', param_removal=(10, 0),
            whitening=False,
            normalize=None, param_norm=None, freq_domain=True)

    correlations = get_correlations(stations, components, only_auto=True)
#    correlations = correlations + [('PB01Z', 'PB02Z'), ('PB03Z', 'PB04Z'), ('PB03Z', 'PB05Z'), ('PB04Z', 'PB05Z')]
    #correlations = (('PB03Z', 'PB04Z'),)

    noisexcorr(data, correlations, t1, t2, shift, freq_domain=True)
#    stack(data, correlations, dt=10 * 24 * 3600, period=24 * 3600, shift=24 * 3600, yearfiles=True)
    #stack(data, correlations, dt= -1)
    plotXcorrs(data, correlations, t1, t2, start= -12, end=12, plot_overview=True, plot_years=False, use_dlognorm=True, show=False,
                      plot_stack=True, plot_psd=False, add_to_title='(water_1bit)', add_to_file='_log_12s', cmax=1e-4, ext='.png')
    plotXcorrs(data, correlations, t1, t2, start= -12, end=12, plot_overview=True, plot_years=False, use_dlognorm=False, show=False,
                      plot_stack=True, plot_psd=False, add_to_title='(water_1bit)', add_to_file='_12s', cmax=1e-2, ext='.png')

#    plotXcorrs(data, correlations, t1, t2, start= -25, end=25, plot_overview=True, plot_years=False, use_dlognorm=True, show=False,
#                      plot_stack=True, stack=('10days', 'day'), plot_psd=True, add_to_title='(water_1bit)', add_to_file='_log_25s', cmax=1e-2)



    shift = 100
#    correlations = get_correlations(stations, components, stations2)
    method = 'filter4-6_water_env2_1bit_fft'
#    method = 'filter4-6_water_env2_whitening_1bit'
#    method = 'filter0.01-1_1bit_whitening0.01'
#    method = 'filter0.005_rm20'
#    method = 'filter0.005_1bit'
    stations = 'PB01 PB02 PB03 PB04 PB05'

    data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)
#    prepare(data, stations.split(), t1, t2, component=components,
#            filter=(4, 6), downsample=None,
#            eventremoval='waterlevel_env2', param_removal=(10, 0),
#            whitening=False,
#            normalize='1bit', param_norm=None, freq_domain=True)

    correlations = get_correlations(stations, components, only_auto=True)
#    correlations = correlations + [('PB01Z', 'PB02Z'), ('PB03Z', 'PB04Z'), ('PB03Z', 'PB05Z'), ('PB04Z', 'PB05Z')]
    #correlations = (('PB03Z', 'PB04Z'),)

#    noisexcorr(data, correlations, t1, t2, shift, freq_domain=True)
#    stack(data, correlations, dt=10 * 24 * 3600, period=24 * 3600, shift=24 * 3600, yearfiles=True)
    #stack(data, correlations, dt= -1)
#    plotXcorrs(data, correlations, t1, t2, start= -12, end=12, plot_overview=True, plot_years=False, use_dlognorm=True, show=False,
#                      plot_stack=True, plot_psd=False, add_to_title='(water_1bit)', add_to_file='_log_12s', cmax=1e-2, ext='.eps')
#    plotXcorrs(data, correlations, t1, t2, start= -25, end=25, plot_overview=True, plot_years=False, use_dlognorm=True, show=False,
#                      plot_stack=True, stack=('10days', 'day'), plot_psd=True, add_to_title='(water_1bit)', add_to_file='_log_25s', cmax=1e-2)


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
