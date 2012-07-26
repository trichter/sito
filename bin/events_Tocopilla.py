#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito import Events
from sito.data import IPOC
from obspy.core import UTCDateTime as UTC
import matplotlib.pyplot as plt
import logging
logging.basicConfig()

lat_Toc = -22.170
lon_Toc = -69.971

#events = Events.load(min_datetime="2007-01-01", max_datetime="2008-12-31",
#                     min_latitude=lat_Toc - 1, max_latitude=lat_Toc + 1,
#                     min_longitude=lon_Toc - 1., max_longitude=lon_Toc + 1,
#                     max_results=1000000,
#                     min_magnitude=None, max_magnitude=None)
#events.write('/home/richter/Data/events/events_Tocopilla.txt')
events = Events.read('/home/richter/Data/events/events_Tocopilla.txt')
events.pick(latitude=lat_Toc, longitude=lon_Toc, minval=0, maxval=100., indegree=False)
#events.plot(lat_Toc, lon_Toc, circles=(1,))

method = 'filter2-20_1bit'
#method = 'filter0.005_1bit'

data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)

t1 = UTC('2007-11-01')
t2 = UTC('2007-12-01')

period = 1800
correlation = ('PB03Z', 'PB03Z')
stream = data.readX(correlation, t1, t2, period=period)
#stream.filter2(2, 20)
stream.setHIForHist(events, period=period)
figsize = (8.267, 11.693)[::-1]
add_to_title = '_againfilter_zoom1'
#save = data.getPlotXCorr(correlation, 'all') + '_againfilter_zoom1 + events.png'
save = False

stream.plotXcorr(0, 50, imshow=True, use_dlognorm=True, filter=(2, 20),
                                 fig=plt.figure(figsize=figsize),
                                 figtitle='station ' + add_to_title,
                                 dateformatter='%y-%m-%d %Hh%M', save=save, show=True, #dateformatter=' % y % b'
                            plotinfo=('num_events',), #plotinfo_width=0.1, #@UnusedVariable
                 plotlabel=('# events',), #@UnusedVariable
                 plotinfowhere=('right',))
plt.show()

