#!/usr/bin/env python
# by TR

from obspy.signal.cpxtrace import envelope
from sito.noisexcorr import setHIDist
import matplotlib.pyplot as plt
from sito.stations import IPOCStations
from sito.stream import read
from operator import itemgetter
from scipy.optimize import curve_fit

from matplotlib.patches import Polygon
import numpy as np

def no_corr_pairs(stream, stations):
    s1 = set()
    s2 = set()
    for tr in stream:
        s1.add(frozenset(tr.stats.station.split('-')))
    for st1 in stations:
        for st2 in stations:
            if st1 != st2:
                s2.add(frozenset({st1 + 'Z', st2 + 'Z'}))
    return s2 - s1

def get_vel(stream):
    ms = stream.copy()
    ms.addXcorrSides()
    for tr in ms:
        tr.data = np.abs(envelope(tr.data))
    dists = ms.getHI('dist')
    maxi = ms.getMaxima()
    v, _ = curve_fit(lambda x, a: x * a, maxi, dists, p0=1)
    return v[0]


#ms = read('/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_whitening/stack/day_PB0[12345]Z-PB0[12345]Z_stack_all.QHD')
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_whitening'
ms = read(path + '/stack/day_*_stack_all.QHD')
output = path + '/xcorr_vs_dist.pdf'

ipoc = IPOCStations()
setHIDist(ms, ipoc)

print 'no correlation for pairs:', no_corr_pairs(ms, ipoc)
v = get_vel(ms)
#v = 3.03093386
print 'velocity:', v
fig = plt.figure(figsize=(10, 12))
plot = ms.plotXcorrVsDist(-300, 300, scale=10, fig=fig,
                          figtitle='%d cross-correlations' % len(ms))
plot.ax.plot((-300, 0, 300), (300 * v, 0, 300 * v), 'r')

d = 30
w = 80
xy = np.array([(d, 0), (300 + d, v * 300), (300 + d + w, v * 300), (d + w, 0)])
polygon1 = Polygon(xy, True, alpha=0.4, color='b', zorder=50)
plot.ax.add_patch(polygon1)
xy[:, 0] *= -1
polygon2 = Polygon(xy, True, alpha=0.4, color='b', zorder=50)
plot.ax.add_patch(polygon2)

plot.ax.set_ylim(0, 720)

plt.savefig(output, bbox_inches='tight')
plt.show()
