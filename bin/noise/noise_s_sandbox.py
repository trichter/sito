#!/usr/bin/env python
# by TR
#@PydevCodeAnalysisIgnore

from obspy.core import UTCDateTime as UTC
from sito.data import Parkfield, IPOC
from sito.noisexcorr import (prepare, get_correlations,
                             plotXcorrs)
from sito.noisexcorr import (prepare, stack_hour,
                             get_correlations, plotXcorrs)
import matplotlib.pyplot as plt

#from sito.util.main import yeargen, streamyeargen, streamyeargen2
#from sito.xcorr import getNormFactors
#from sito.util.helper import exha


def main():
    data = IPOC(xcorr_append='/tests/1bit_filter0.1-1', use_local_LVC=True)
    data.setXLogger('_1bit')
    stations = 'PB01 PB03'
    stations2 = 'PB03'

    components = 'Z'
    t1 = UTC('2010-01-01')
    t2 = UTC('2010-12-31')
    shift = 500

    prepare(data, stations.split(), t1, t2, filter=(0.1, 1.), downsample=10,
            component=components, normalize='1bit', param_norm=None,
            use_floating_stream=True)
    correlations = get_correlations(stations, components, stations2)
    xcorr_day(data, correlations, t1, t2, shift, use_floating_stream=True)
    plotXcorrs(data, correlations, t1, t2)


# prepare stream with running mean average
#    prepare(data, stations.split(' '), t1, t2, filter=(0.1, 1.), downsample=10,
#            component=components, normalize='runningmean', norm_param=20 * 10 + 1,
#            use_floating_stream=True)


#    prepare(data, stations, t1, t2, filter=(0.1, 1.), downsample=10, component='Z',
#            normalize='1bit', norm_param=None, use_floating_stream=True)

# - compare: filter = (0.01,0.1); (0.1,1); (1,10) with downsample=None


def plot_streams(streams, scales=None):
    fig = plt.figure()
    N = len(streams)
    if scales == None:
        scales = [None] * N
    vlay = 2
    hlay = (N - 1) // vlay + 1
    ax = fig.add_subplot(vlay, hlay, 1)
    for i in range(N - 1):
        fig.add_subplot(vlay, hlay, i + 2, sharex=ax, sharey=ax)
    for i in range(N):
        streams[i].plotTrace(ax=fig.axes[i], absolutescale=scales[i], title_in_axis=True)
    fig.show()

def plotPSD(streams, Nfft=256):
    fig = plt.figure()
    N = len(streams)
    if isinstance(Nfft, int):
        Nfft = [Nfft] * N
    vlay = 2
    hlay = (N - 1) // vlay + 1
    ax = fig.add_subplot(vlay, hlay, 1)
    for i in range(N - 1):
        fig.add_subplot(vlay, hlay, i + 2, sharex=ax)
    for i in range(N):
        streams[i].plotPSD(ax=fig.axes[i], Nfft=Nfft[i], title_in_axis=True,
                           figtitle='station component date  Nfft:nfft')
    fig.show()

def main2():
    data = IPOC(xcorr_append='/1bit', use_local_LVC=True)
    t1 = UTC('2010-01-01')
    stream0_1 = data.getRawStream(t1, 'PB01', component='Z')
    stream0_2 = data.getRawStream(t1, 'PB02', component='Z')

    stream2_1 = data.getStream(t1, 'PB01', component='Z')
    stream2_2 = data.getStream(t1, 'PB02', component='Z')

    plot_streams([stream0_1, stream0_2, stream2_1, stream2_2], [None, None, 0.1, 0.1])
    plotPSD([stream0_1, stream0_2, stream2_1, stream2_2], 4096)

    ipshell()


if __name__ == '__main__':
    main2()
