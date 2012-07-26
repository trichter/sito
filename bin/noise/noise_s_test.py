#!/usr/bin/env python
# by TR
#@PydevCodeAnalysisIgnore

from obspy.core import UTCDateTime as UTC
from sito.data import Parkfield, IPOC
from sito.noisexcorr import (prepare, xcorr_day, get_correlations,
                             plotXcorrs)
from sito.noisexcorr import (prepare, stack_hour, xcorr_day, xcorr_hour,
                             get_correlations, plotXcorrs)
import matplotlib.pyplot as plt

#from sito.util.main import yeargen, streamyeargen, streamyeargen2
#from sito.xcorr import getNormFactors
from sito import ipshell
#from sito.util.helper import exha

def main():
    stations = 'PB01 PB03'
    stations2 = 'PB03'
    components = 'Z'
    t1 = UTC('2010-01-01')
    t2 = UTC('2010-12-31')
    shift = 500
    correlations = get_correlations(stations, components, stations2)


    method = 'filter0.1-1_1bit_whitening0.01'
    data = IPOC(xcorr_append='/tests/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)
    prepare(data, stations.split(), t1, t2, filter=(0.1, 1), downsample=10, whitening=0.01,
            component=components, normalize='1bit', param_norm=None,
            use_floating_stream=True)
    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
    plotXcorrs(data, correlations, t1, t2, plot_overview=False, plot_stack=True, plot_psd=True, add_to_title=method)


    method = 'filter0.1-1_1bit_whitening0.001'
    data = IPOC(xcorr_append='/tests/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)
    prepare(data, stations.split(), t1, t2, filter=(0.1, 1), downsample=10, whitening=0.001,
            component=components, normalize='1bit', param_norm=None,
            use_floating_stream=True)
    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
    plotXcorrs(data, correlations, t1, t2, plot_overview=False, plot_stack=True, plot_psd=True, add_to_title=method)

    method = 'filter0.1-1_1bit_whitening0.1'
    data = IPOC(xcorr_append='/tests/' + method, use_local_LVC=True)
    data.setXLogger('_' + method)
    prepare(data, stations.split(), t1, t2, filter=(0.1, 1), downsample=10, whitening=0.1,
            component=components, normalize='1bit', param_norm=None,
            use_floating_stream=True)
    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
    plotXcorrs(data, correlations, t1, t2, plot_overview=False, plot_stack=True, plot_psd=True, add_to_title=method)


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
