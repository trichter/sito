#!/usr/bin/env python
# by TR
"""
Creates a matrix of dayplots (different stations, different dates)
"""

from sito.data import IPOC
from obspy.core import UTCDateTime as UTC
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from progressbar import ProgressBar

redraw = False
# 13 stations,  7 days, 2007, 2008
# days = '2006-12-01 2007-02-01 2007-04-01 2007-06-01 2007-08-01 2007-10-01 2007-12-01'
days = '2008-01-01 2008-02-01 2008-04-01 2008-06-01 2008-08-01 2008-10-01 2008-12-01'
days = '%d-01-01 %d-02-01 %d-04-01 %d-06-01 %d-08-01 %d-10-01 %d-12-01' % ((2012,) * 7)
stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 HMBCX MNMCX PATCX PSGCX LVC'
# 20 stations,  10 days, 2008-2012
# days = '2008-01-01 2008-07-01 2009-01-01 2009-07-01 2010-01-01 2010-07-01 2011-01-01 2011-07-01 2012-01-01 2012-07-01'
# stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 PB15 HMBCX MNMCX PATCX PSGCX LVC'
pic_dir = '/home/richter/Data/IPOC/raw_pictures/'
output = '/home/richter/Results/IPOC/raw_pics/test.pdf'
output = '/home/richter/Results/IPOC/raw_pics/13stations_7days_2012.pdf'
scale_range = 2000
scale_LVC = 20000

days = days.split()
stations = stations.split()
dx = 0.98 / len(days)
dy = 0.98 / len(stations)
data = IPOC()
figsize = (11.69, 16.53)  # A3 # 8.27, 11.69 # A4
fig = plt.figure(figsize=figsize)
for i, station in ProgressBar(len(stations))(enumerate(stations)):
    for j, day in enumerate(days):
        x0 = 0.01 + dx * (j + 0.02)
        y0 = 0.99 - dy * (1 + i)
        if i == 0:
            fig.text(x0 + 0.98 * dx / 2, 0.98, day, ha='center', va='center')
        if j == 0:
            fig.text(0.02, y0 + (0.98 / 2) * dy, station, va='center', ha='center', rotation=90)
        t_day = UTC(day)
        fname = '%s/%s-%s.png' % (pic_dir, station, day)
        try:
            if redraw:
                raise IOError
            img = mpimg.imread(fname)
        except IOError:
            try:
                stream = data.getRawStreamFromClient(t_day, t_day + 24 * 3600,
                                                     station, component='Z')
            except ValueError:
                continue
            fig2 = stream.plot(type='dayplot',
                              vertical_scaling_range=scale_range if station != 'LVC' else scale_LVC,
                              interval=30, handle=True)
            fig2.axes[0].set_position((0, 0, 1, 1))
            fig2.axes[1].set_visible(False)
            fig2.texts[0].set_visible(False)
            fig2.canvas.draw()
            fig2.savefig(fname)
            img = mpimg.imread(fname)
        ax = fig.add_axes([x0, y0, dx * 0.98, dy * 0.98])
        ax.imshow(img)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
fig.savefig(output, dpi=300)
