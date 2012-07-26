#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
"""
calculate piercing points in specific depth,
create plots (profiles) where rfs are stacked as for longitude (latitude) of
piercing points

functions:
calculate_pspier: calculate piercing points and write to header
plot_pspier: plot piercing point map
calculate_profile: calculate stacked rfs
plot, plot_depth: helper functions
plot2: plot profile(s)
plot3: plot comparison of profiles north and south of -21° latitude
"""

from sito import read, util
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from glob import glob
import os.path
from matplotlib.patches import Rectangle

def calculate_pspier(other_header=False, station='*'):
    from os.path import splitext
    for file_ in glob(path + station + '_mout.QHD'):
        ms = read(file_)
        if other_header:
            ms.pspier(pspier_depth, other_header=str(pspier_depth))
        else:
            ms.pspier(pspier_depth)
        ms.write(splitext(file_)[0], 'Q')

def plot_pspier():
    coordfile = pp_path + 'ppcoords%d.npy' % pspier_depth
    if os.path.isfile(coordfile):
        coords = np.load(coordfile)
    else:
        ms = read(path + '*_mout.QHD').select(component='Q', expr='not st.mark')
        lat = ms.getHI('plat%d' % pspier_depth)
        lon = ms.getHI('plon%d' % pspier_depth)
        coords = np.array([lat, lon])
        np.save(coordfile, coords)
    from sito import map
    m = map.createIPOCMap(show=False, trench=None)
    x, y = m(coords[1, :], coords[0, :])
    m.plot(x, y, 'xr', ms=3)
    plt.gcf().savefig(pp_path + 'map_pp%d.pdf' % pspier_depth)

def calculate_profile():
    ms = read(path + '*_mout.QHD').select(component='Q', expr='not st.mark')
    if header2 is None:
        ms = ms.getBinnedStream(bins1, header=header)
        ms.write(pp_file, 'Q')
    else:
        for i in range(len(bins2) - 1):
            ms2 = ms.select(expr='%f>=st.%s>%f' % (bins2[i], header2, bins2[i + 1]))
            ms3 = ms2.getBinnedStream(bins1, header=header)
            ms3.write(pp_file % (bins2[i], bins2[i + 1]), 'Q')

def plot_depth(ax):
    ax2 = ax.twinx()
    ax2.set_position(ax.get_position(False))
    h = np.array((0, 50, 100, 150, 200))
    h2 = np.arange(20) * 10
    t = util.depth2time(h)
    myLocator = mpl.ticker.FixedLocator(t)
    myMinorLocator = mpl.ticker.FixedLocator(util.depth2time(h2))
    myFormatter = mpl.ticker.FixedFormatter([str(i) for i in h])
    ax2.yaxis.set_major_locator(myLocator)
    ax2.yaxis.set_minor_locator(myMinorLocator)
    ax2.yaxis.set_major_formatter(myFormatter)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_ylabel('depth (km)')


def plot(title=''):
    ms = read(pp_file + '.QHD')
    if stack_lon:
        xlabel = u'longitude of piercing points (°)'
    else:
        xlabel = u'latitude of piercing points (°)'
    plot = ms.plotProfile(start, end, scale=scale, xaxis=header,
                              xlabel=xlabel,
                              figtitle=title, fancy_box=True, box_ax=1)
    ax = plot.fig.axes[0]
    ax.set_xlim(xlim)
    plot_depth(ax)

def plot2():
    global pp_file
    if header2 is None:
        if stack_lon:
            title = u'latitude %d° to %d°' % (bins2[0], bins2[-1])
        else:
            title = None

        plot(title=title)
        plt.savefig(pp_file + '.pdf')
    else:
        pp_file_old = pp_file
        for i in range(len(bins2) - 1):
            pp_file = pp_file_old % (bins2[i], bins2[i + 1])
            title = None
            if header2 and 'plat' in header2:
                title = u'latitude %d° to %d°' % (bins2[i], bins2[i + 1])
            elif header2 and 'plon' in header2:
                title = u'longitude %.1f° to %.1f°' % (bins2[i], bins2[i + 1])
            plot(title=title)
            plt.savefig(pp_file + '.pdf')
            plt.close()
        pp_file = pp_file_old

def plot3():
    assert(header2 is not None)
    #colors = 'bgrcmyk'
    for i in range(len(bins2) - 1):
        ms = read(pp_file % (bins2[i], bins2[i + 1]) + '.QHD')
#        if i == len(bins2) - 2:
#            reverse_y = True
        if i == 0:
            plot = ms.plotProfile(start, end, scale=scale, xaxis=header,
                                  xlabel=u'longitude of piercing points (°)',
                                  plotinfo=(),
                                  figtitle='', topcolor='red', botcolor='w',
                                  show=False, alpha=0.5)
        else:
            ms.plotProfile(start, end, xaxis=header,
                           scale=scale,
                           plotinfo=(),
                           figtitle=None, show=False, topcolor='green', botcolor='w',
                           ax=plot.ax, ax_info=plot.ax_info, alpha=0.3)

    ax = plot.fig.axes[0]
    ax.set_xlim(xlim)
    plot_depth(ax)
    poly1 = Rectangle((0, 0), 1, 1, fc="r", alpha=0.5)
    poly2 = Rectangle((0, 0), 1, 1, fc="g", alpha=0.3)
    label1 = u'%d° to %d°' % (bins2[0], bins2[1])
    label2 = u'%d° to %d°' % (bins2[1], bins2[2])
    leg = ax.legend((poly1, poly2), (label1, label2), loc='lower left', fancybox=True,
                    title='latitude')
    frame = leg.get_frame()
    frame.set_facecolor('wheat')
    frame.set_alpha(0.8)
    # matplotlib.text.Text instances
    for t in leg.get_texts():
        t.set_fontsize('small')
    plot.fig.savefig(pp_file % ('comparison_' + str(bins2[0]), bins2[-1]) + '.pdf')



pspier_depth = 80
stack_lon = True
start = -3
end = 22

path = '/home/richter/Results/IPOC/receiver/2012_mag5.5/'
if stack_lon:
### for longitutional stack
    dif = 0.02
    scale = 2 * dif
    lon1 = -70.8
    lona = -69.5
    lon2 = -68.8
    xlim = (-70.9, -68.7)
    header = 'plon%d' % pspier_depth
    header2 = None
    #header2 = 'plat%d' % pspier_depth
    bins1 = np.linspace(lon1, lon2, int((lon2 - lon1) / dif) + 1)
    bins2 = (-18, -21, -24)
    #bins = np.linspace(lon1, lona, int((lona - lon1) / dif) + 1)
    #difa = 0.05
    #binsa = np.linspace(lona, lon2, int((lon2 - lona) / difa) + 1)[1:]
    #bins_lon = np.hstack((bins, binsa))
    #pp_path = path + 'pspier/
    pp_path = path + 'pspier/stack_lon/'
else:
### for latitudional stack
    dif = 0.5
    scale = 1 * dif
    lat1 = -23. # -25
    lat2 = -18.
    xlim = (-23, -18)
    header = 'plat%d' % pspier_depth
    header2 = 'plon%d' % pspier_depth
    bins1 = np.linspace(lat1, lat2, int((lat2 - lat1) / dif) + 1)
    bins2 = (-69, -69.5, -70, -70.5)
    #pp_path = path + 'pspier/    
    pp_path = path + 'pspier/stack_lat/'
pp_file = pp_path + 'profile_pp_Q_%s' % str(dif)
pp_file = pp_path + '%s_profile_%s%s' % (header, str(dif), '%s')
if header2 is None:
    pp_file = pp_file % ''
else:
    pp_file = pp_file % ('_' + header2 + '_%s_%s')




#calculate_pspier(True)
#plot_pspier()
calculate_profile()
plot2()
#plot3()

#plt.show()
