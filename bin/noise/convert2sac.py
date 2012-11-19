#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.noisexcorr import get_correlations
from sito.stream import read
import pylab as plt
from sito.util import checkDir
from obspy.core.util.attribdict import AttribDict

def main():
    stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 PB15 PB16 HMBCX MNMCX PATCX PSGCX LVC'
    #stations = 'PB02 PB10 PB11'
    components = 'Z'
    correlations = get_correlations(stations, components, only_cross=True)

    method = 'FINAL_filter0.005-5_1bit_whitening_2011+2012'

    data = IPOC(xcorr_append='/' + method, use_local_LVC=False)
    #data.setXLogger('_' + method)
    for cor in correlations:
        try:
            st = read(data.getX(cor, time=None, filter=None, period=24 * 3600, stack= -1) + '.QHD')
        except IOError as ex:
            print ex
            continue
        assert len(st) == 1
        sdata = st[0].data
        st[0].data = 0.5 * (sdata + sdata[::-1])
        st.trim2(0, None, 'middle')
        stats = st[0].stats
        stats.station = cor[0][:4] + cor[1][:4]
        stats.starttime = UTC('2012-01-01')
        stats.sac = AttribDict()
        stats.sac.dist = data.stations.dist(cor[0][:-1], cor[1][:-1])
        checkDir(data.x_sac)
        st.write(data.x_sac % cor, 'SAC')
        #st.plot_()
    plt.show()


if __name__ == '__main__':
    main()
