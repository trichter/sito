#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito.data import IPOC
import numpy as np
from sito.util.main import streamdaygen
import pylab as plt

stations = 'PB01'  # PB04 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB16'
channels = 'WKI WDI WII'
#channels = 'WDI_10'
datafile = '/home/richter/Data/climate/2006-2012_%s_%s.npz'
output = '/home/richter/Results/IPOC/climate/%s.pdf'
calculate = False
show = False

if calculate:
    ipoc = IPOC()
    for station in stations.split():
        for channel in channels.split():
            stream = ipoc.getChannelFromClient('2006-01-01', '2013-01-01',
                                               station=station, channel=channel)
            data = []
            dates = []
            for day in streamdaygen(stream):
                day.merge()
                data.append(np.mean(day[0].data))
                st = day[0].stats.starttime
                et = day[0].stats.endtime
                dates.append(st + (et - st) / 2.)
            np.savez(datafile % (station, channel), dates=dates, data=data)
else:
    #http://stackoverflow.com/questions/7733693/matplotlib-overlay-plots-with-different-scales
    fig, ax = plt.subplots()
    axes = [ax, ax.twinx(), ax.twinx()]
    fig.subplots_adjust(right=0.75)
    axes[-1].spines['right'].set_position(('axes', 1.2))
    axes[-1].set_frame_on(True)
    axes[-1].patch.set_visible(False)
    #colors = ('green', 'red', 'blue')
    colors = 'bgrcmykbgrc'
    for i, channel in enumerate(channels.split()):
        ax = axes[i]
        for j, station in enumerate(stations.split()):
            color = colors[i]
            npzfile = np.load(datafile % (station, channel))
            data = npzfile['data']
            data = data / 100. if channel == 'WKI' else data / 10000. if channel == 'WDI' else data
            dates = [date.toordinal() for date in npzfile['dates']]
            leglabel = station if i == 0 else None
#            print np.where(np.array(dates) < 3000)
#            from IPython import embed
#            embed()
            ax.plot(np.array(dates), data, color=color, label=leglabel)
            if j == 0:
                ax.set_xlim([min(dates), max(dates)])
            #plt.show()
        label = '(temp in degC)' if channel == 'WKI' else '(pressure in bar)' if channel == 'WDI' else '(rel. humidity in %)'
        ax.set_ylabel('%s' % channel + ' ' + label, color=color)
        ax.tick_params(axis='y', colors=color)
    ax = axes[0]
    ax.xaxis_date()
    for label in ax.get_xticklabels():
        label.set_ha('right')
        label.set_rotation(30)
    ax.legend()
    if show and not ' ' in stations:
        plt.show()
    else:
        fig.savefig(output % station)

