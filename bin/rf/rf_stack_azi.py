#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
"""
stack and plot rfs by azimuth

"""
from sito import read, imaging
from glob import glob
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from matplotlib.colorbar import ColorbarBase

path = '/home/richter/Results/IPOC/receiver/2012_mag5.5/'
azi_path = path + '/azi_stack/'
plotdir = azi_path + 'plots/'

def getFig(num=0, ratio=1.5, margin=None, **kwargs):
    axes = [1. - 0.1 * num] + [0.1] * num
    if margin == None:
        margin = [1.7, 1.3, 1.1, 0.3] #left, rigth, bottom, top
    fig = imaging.getFigure(axes, width=15., margin=margin, ratio=ratio,
                            fontsize=12, labelsize='small', **kwargs)
    return fig


def cosmetic(station, fig):
    bins = np.arange(18) * 20 + 10
    if station in ('PB01'):
        fig.axes[0].set_yticks(bins)
        fig.axes[0].set_yticks((), minor=True)
        fig.axes[0].set_yticklabels([str(bi) * ((i + 1) % 2)
                                     for i, bi in enumerate(bins)])
    fig.axes[0].set_ylim((-20, 380))

def plot(station='*'):
    start = -5
    end = 22
    show = False
    for file_ in glob(azi_path + station + '_azi_stack.QHD'):
        ms = read(file_)
        station = ms[0].stats.station
        ratio = (len(ms) * 0.5 + 1.4 + 0.4 * 2.54) / 15
        ratio = min(ratio, 2.)
        fig = getFig(ratio=ratio)
        alpha = None
        num_tr = np.sum(np.array(ms.getHI('count')))
        if num_tr >= 100:
            alpha = ['count', 20, 0, 1., 0.]
        elif num_tr >= 50:
            alpha = ['count', 10, 0, 1., 0.]
        else:
            alpha = ['count', 5, 0, 1., 0.]
        plot = ms.plotRF(start, end, yaxis='azi', ylabel=u'azi (°)', show=show,
                         fig=fig, scale=360 / len(ms), plotinfo=('sum',), plotinfowhere=('top',),
                         alpha=alpha)
        if alpha is not None:
            #http://matplotlib.sourceforge.net/examples/api/colorbar_only.html
            ax2 = plot.fig.add_axes([0.94, 0.2, 0.01, 0.6])
            norm = colors.Normalize(vmin=0, vmax=alpha[1])
            ColorbarBase(ax2, cmap='Greys', norm=norm, extend='max')
        cosmetic(station, plot.fig)
        plot.fig.savefig(plotdir + 'rf_azistack_%s_Q.pdf' % station)
        plt.close(plot.fig)

def calculate(station='*'):
    for file_ in glob(path + station + '_mout.QHD'):
        ms = read(file_).select(component='Q', expr='not st.mark')
        station = ms[0].stats.station
        if len(ms) >= 100:
            bins = np.arange(19) * 20
        elif len(ms) >= 50:
            bins = np.arange(13) * 30
        else:
            bins = np.arange(7) * 60
        ms = ms.getBinnedStream(bins, header='azi')
        ms.write(azi_path + '%s_azi_stack' % station, 'Q')

calculate()
plot()
