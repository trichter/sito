#!/usr/bin/env python
# by TR

from sito.noisexcorr import setHIDist
import matplotlib.pyplot as plt
from sito.stations import IPOCStations
from sito.stream import read

#ms = read('/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_whitening/stack/day_PB0[12345]Z-PB0[12345]Z_stack_all.QHD')
ms = read('/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_whitening/stack/day_*_stack_all.QHD')
ipoc = IPOCStations()
setHIDist(ms, ipoc)
plot = ms.plotXcorrVsDist(-300, 300, scale=10)
plot.ax.set_ylim(0, 700)
print 171 / 58.4
plt.show()
