#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
"""
script for plotting rfs
"""

from sito import read, imaging
from glob import glob
import matplotlib.pyplot as plt

path = '/home/richter/Results/IPOC/receiver/2012_mag5.5/'
plotdir = path + 'plots/'

def getFig(num=0, ratio=1.5, margin=None, **kwargs):
    axes = [1. - 0.1 * num] + [0.1] * num
    if margin == None:
        margin = [1.7, 0.4, 1.1, 0.3] #left, rigth, bottom, top
    fig = imaging.getFigure(axes, width=15., margin=margin, ratio=ratio,
                            fontsize=12, labelsize='small', **kwargs)
    return fig

def cosmetic(station, fig):
    if station in ('PB05', 'PB09', 'PATCX', 'LVC'):
        if station in ():
            ind = -2
        else:
            ind = -1
        tls = fig.axes[0].get_yticklabels()
        tls[ind].set_visible(False)
    if station == 'PB16':
        fig.axes[0].set_yticks((0, 5, 10))
        #fig.axes[0].set_yticks((), minor=True)
        fig.axes[2].set_yticks((0, 5, 10))
        #fig.axes[2].set_yticks((), minor=True)
        fig.axes[0].set_yticklabels(('0', '5', ''))

def plot(station='*'):
    start = -5
    end = 22
    show = False
    components = 'LQT'
    for file_ in glob(path + station + '_mout.QHD'):
        ms = read(file_)
        station = ms[0].stats.station
        ms = ms.select('not st.mark')
        ms.sort('azi')
        ratio = (len(ms) // 3 * 0.1 + 1.4 + 0.4 * 2.54) / 15
        ratio = min(ratio, 2.)
        fig = getFig(ratio=ratio)
        plot = ms.plotRF(start, end, show=show, fig=fig, scale=1, component=components[1])
        cosmetic(station, plot.fig)
        plot.fig.savefig(plotdir + 'rf_%s_%s.pdf' % (station, components[1]))
        plt.close(plot.fig)
        fig2 = getFig(ratio=ratio)
        plot = ms.plotRF(start, end, show=show, fig=fig2, scale=1, component=components[2])
        cosmetic(station, plot.fig)
        plot.fig.savefig(plotdir + 'rf_%s_%s.pdf' % (station, components[2]))
        plt.close(plot.fig)
        fig3 = getFig(ratio=ratio)
        plot = ms.plotRF(start, end, show=show, fig=fig3, scale=0.9, component=components[0])
        cosmetic(station, plot.fig)
        plot.fig.savefig(plotdir + 'rf_%s_%s.pdf' % (station, components[0]))
        plt.close(plot.fig)

plot()



