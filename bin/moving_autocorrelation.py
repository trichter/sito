#!/usr/bin/env python
# by TR

from sito.data import IPOC
from obspy.core import UTCDateTime as UTC
import matplotlib.pyplot as plt
from sito.util.main import streamtimegen
from sito.stream import Stream
from progressbar import ProgressBar

data = IPOC()
t_day = UTC('2008-01-01')
station = 'PB01'

stream = data.getRawStreamFromClient(t_day, t_day + 24 * 3600,
                                                     station, component='Z')
stream.setHI('filter', '')
stream.demean()
stream.filter2(0.5, 5)
stream.trim2(0, 5 * 3600)

auto = Stream()
for st in streamtimegen(stream, dt=60, start=None, shift=30, use_slice=True):
    tr = st[0].copy()
    tr.addZeros(60)
    tr.acorr(60)
    auto.append(tr)

print auto
auto.plotXcorr()
stream.plot(type='dayplot')


plt.show()
