#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.util import timegen
import matplotlib.pyplot as plt
from sito.debug import embed
import numpy as np
from sito import Trace


def main():
    stations = 'PB01 PB02 PB03 PB04 PB05'

    component = 'Z'
    t1 = UTC('2009-06-01')
    t2 = UTC('2009-07-01')

    data = IPOC('test', use_local_LVC=False)
    data.setXLogger('_test')
    period = 24 * 3600
    ax = None
    plt.ion()
    for station in stations.split():
        pxxs = []
        freqs_old = None
        i = 0
        for t in timegen(t1, t2, period):
            st = data.getRawStreamFromClient(t, t + period, station, component)
            st.merge(method=1, interpolation_samples=10, fill_value='interpolate')
            print st
            pxx, freqs = st.plotPSD(just_calculate=True)
            assert np.all(freqs == freqs_old) or not freqs_old
            freqs_old = freqs
            if max(pxx[4:]) > 1e7:
                print 'discard'
                i += 1
                continue
            pxxs.append(pxx)
        pxx = sum(pxxs) / len(pxxs)
        del pxxs
        tr = Trace(data=pxx,
                   header=dict(is_fft=True, sampling_rate=2 * max(freqs),
                               freq_min=min(freqs), freq_max=max(freqs)))
        ax = tr.plotPSD(ax=ax, label='%s-%d' % (st[0].stats.station, i), figtitle=None)
        plt.draw()
        #embed()
    ax.legend()
    fig = ax.get_figure()
    fig.suptitle('%s  %s  %s to %s' % (stations, component, t1.strftime('%Y-%m-%d'),
                                       t2.strftime('%Y-%m-%d')))
    plt.ioff()
    plt.show()

if __name__ == '__main__':
    main()
