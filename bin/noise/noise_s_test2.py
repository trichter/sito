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
    data = IPOC(xcorr_append='/tests/1bit_filter0.01', use_local_LVC=True)
    data.setXLogger('_1bit0.01Hz')
    stations = 'PB01 PB03'
    stations2 = 'PB03'

    components = 'Z'
    t1 = UTC('2010-01-01')
    t2 = UTC('2010-01-02')
    shift = 500

#    prepare(data, stations.split(), t1, t2, filter=(0.01, None), downsample=None,
#            component=components, normalize='1bit', norm_param=None,
#            use_floating_stream=True)
    correlations = get_correlations(stations, components, stations2)
#    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
    plotXcorrs(data, correlations, t1, t2, plot_overview=False, plot_stack=True, plot_psd=True)


# prepare stream with running mean average
#    prepare(data, stations.split(' '), t1, t2, filter=(0.1, 1.), downsample=10,
#            component=components, normalize='runningmean', norm_param=20 * 10 + 1,
#            use_floating_stream=True)


# - compare: filter = (0.01,0.1); (0.1,1); (1,10) with downsample=None



if __name__ == '__main__':
    main()
