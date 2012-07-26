#
# Default preamble
#
import matplotlib.pyplot as plt
import obspy, sito
from obspy.core import UTCDateTime as UTC
from sito import read
# -*- snip -*-
#ms = read('/home/richter/Results/IPOC/xcorr/Tocopilla'
#		  '/filter0.005_rm20/stack/PB04Z_stack_*10.QHD')
ms = read('/home/richter/Results/IPOC/xcorr/Tocopilla'
		  '/filter0.005_1bit/day/PB03Z-PB*_day_*.QHD')
#ms.filter2(10., None)
ms.filter2(2., None, corners=6)

#ms.downsample2(10)

ms.plotXcorr(-100, 100, stack_lim=(-0.01, 0.01), vmax=1, use_dlognorm=True, cmax=1e-8)#, vmax=0.001)
plt.show()
