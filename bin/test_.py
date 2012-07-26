#!/usr/bin/env python

from sito.data import IPOC
from sito.noisexcorr import prepare, noisexcorr, plotXcorrs
from obspy.core import UTCDateTime as UTC

t1 = UTC('2008-09-13')
t2 = UTC('2008-09-14')

method = 'filter4-6_water_env2_1bit'
#    method = 'filter0.01-1_1bit_whitening0.01'
#    method = 'filter0.005_rm20'
#    method = 'filter0.005_1bit'
stations = 'PB03 PB04'
correlations = (('PB03Z', 'PB04Z'),)
components = 'Z'
shift = 500
data = IPOC(xcorr_append='/tests/Tocopilla/' + method, use_local_LVC=True)
data.setXLogger('_' + method)
prepare(data, stations.split(), t1, t2, component=components,
        filter=(4, 6), downsample=None,
        eventremoval='waterlevel_env2', param_removal=(10, 0),
        #whitening=True,
        normalize='1bit', param_norm=None)
noisexcorr(data, correlations, t1, t2, shift)
plotXcorrs(data, correlations, t1, t2, start= -50, end=50, plot_overview=True, plot_years=False, use_dlognorm=False,
                  plot_stack=True, plot_psd=True, add_to_title=method, downsample=None)
plotXcorrs(data, correlations, t1, t2, start= -50, end=50, plot_overview=True, plot_years=False, use_dlognorm=True,
                  plot_stack=True, plot_psd=True, add_to_file='_dlognorm.png', add_to_title=method, downsample=None)
