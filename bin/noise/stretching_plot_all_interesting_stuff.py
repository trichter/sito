#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from sito import read
import numpy as np

import csv

mpl.rcParams.update({'font.size': 15})
#maj_loc = mdates.MonthLocator()
#min_loc = mdates.DayLocator((5, 10, 15, 20, 25))
#maj_fmt = mdates.DateFormatter('%b %Y')
maj_loc = mdates.YearLocator()
min_loc = mdates.MonthLocator()
maj_fmt = mdates.DateFormatter('%Y')
#maj_loc=None


#path = '/home/richter/Results/IPOC/xcorr/PAT_filter%d-%d/stretch/sinus_exp_alt_opt.csv'
path = '/home/richter/Results/IPOC/xcorr/PAT_filter%d-%d/stretch2/sinus_exp_alt2_opt.csv'
#output = '/home/richter/Results/IPOC/xcorr/f6_PAT_filter_analysis_%s_with_vel_perm.pdf'
output = '/home/richter/Results/IPOC/xcorr/f6_PAT_filter_analysis_%s_without_vel_perm.pdf'
path_stack = '/home/richter/Results/IPOC/xcorr/PAT_filter%d-%d/stack/day_PATCXZ_stack_all.QHD'

#    middle of frequency, middle of useful time windows
freq_tws = {3: range(5, 21),  #25
            4: range(3, 21),  #22
            5: range(3, 21),
            6: range(3, 21),  #22
            7: range(2, 21),
            8: [2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15],
            9: range(2, 6),
            10: range(2, 10),
            11: range(2, 10)
            }
#  middle of frequency, values for time windows
values = ('vel_change', 't_phase', 'vel_temp', 't_dec', 'mean_corr', None)
vel_change = {}
t_phase = {}
vel_temp = {}
t_dec = {}
mean_corr = {}
for freq, tws in freq_tws.items():
    vel_change[freq] = []
    t_phase[freq] = []
    vel_temp[freq] = []
    t_dec[freq] = []
    mean_corr[freq] = []
    csv_file = path % (freq - 1, freq + 1)
    for vals in csv.reader(open(csv_file)):
        if not vals[0].startswith('#') and int(vals[1]) + 1 in tws:
            vel_offset, vc, tp, vel_perm, vt, td, mc = [float(v) for v in vals[2:]]
            vel_change[freq].append(100 * abs(vc))
            t_phase[freq].append(365.242 * (abs(tp) - 0.25) - 30.0895)
            #vel_temp[freq].append(100 * abs(vt))
            vel_temp[freq].append(100 * abs(vt + vel_perm))
            t_dec[freq].append(abs(td) / 365.242)
            mean_corr[freq].append(mc)

ylabels = [r'periodic vel change $\epsilon_{\rm{P}}$ (%)',
           r'phase delay $t_{\rm P}-t_{\rm{temp}}$ (days)',
           r'EQ vel change $\epsilon_{\rm{EQ}}$ (%)',
           r'recovery of EQ change $t_{\rm{EQ}}$ (years)',
           r'correlation coefficient $\rm{cc}_{\rm{mean}}$',
           'auto-correlation function']
cmap = mpl.cm.get_cmap('jet')
lss = ['-.', '--', '-', '-', '-', '-', '-', '-', '-']

steps = [(3, 3, 0.1, 3), (3, 3), (), (5, 2), (5, 2, 0.5, 2, 0.1, 2), (2, 1.5, 2, 3), (1, 1.5), (0.1, 2), (2, 2, 2, 2, 4, 2)]

#http://colorbrewer2.org/index.php?type=qualitative&scheme=Set1&n=9
colors = '228, 26, 28; 55, 126, 184; 77, 175, 74; 152, 78, 163; 255, 127, 0; 255, 255, 51; 166, 86, 40; 247, 129, 191; 153, 153, 153'
colors = colors.split(';')
colors = [[float(c2) / 255 for c2 in c.strip().split(',')] for c in colors]

#figure 4 focal colors http://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
colors = ('k', 'r', 'orange', '#D6D600', '#00FF00', 'green', '#00FFFF', '#5F5FFF', '#CC00FF')

# color by colormap
#colors = [cmap(i / 9.) for i in range(9)]

fsize = (15, 8)
fw = 170 / 25.4
fh = fw / 1.61
fsize = (fw, fh)
mpl.rcParams.update({'figure.figsize': fsize, 'font.size': 7, 'lines.linewidth':0.7})
mpl.rcParams.update({'font.size': 7, 'lines.linewidth':0.5})
fig = plt.figure(figsize=fsize)

h = 0.39
w = 0.255
x0 = 0.07
dw = 0.075
y0 = 0.08
dh = 0.11
which_axes = [0, 3, 1, 4, 2, 5]
axes = ([x0, y0 + h + dh, w, h], [x0 + w + dw, y0 + h + dh, w, h], [x0 + 2 * w + 2 * dw, y0 + h + dh, w, h],
        [x0, y0, w, h], [x0 + w + dw, y0, w, h])
axes = axes + ([x0 + 2 * w + 1.65 * dw, y0, w + 0.35 * dw, h],)

[0.18, 0.18, 0.77, 0.77]
for i, key in enumerate(values):
    k = which_axes[i]
    if k is None:
        continue
    ax = fig.add_axes(axes[k])
    ax.set_ylabel(ylabels[i])
    ax.set_xticks([0, 5, 10, 15, 20])
    ax.set_xticks(range(21), minor=True)
    if i < 5:
        ax.set_xlabel('mid of time window (s)')
        for j, (freq, tws) in enumerate(freq_tws.items()):
            l, = ax.plot(tws, locals()[key][freq], label='%d-%dHz' % (freq - 1, freq + 1), lw=1, ls=lss[j], color=colors[j])
            if steps[j]:
                l.set_dashes(steps[j])
            l.set_dash_capstyle('round')
        ax.set_xlim([2, 20])
        if i == 2 or i == 0:
            ax.set_ylim([0, 1.3])
        elif i == 4:
            ax.set_ylim([0, 1])
        elif i == 3:
            ax.set_ylim([0, 4])
            ax.set_yticks(range(5))
            #ax.annotate('phase of temperature', (14, 36), va='top', ha='center')
            #ax.axhline(31, color='k', zorder=-100, lw=1)
            #ax.set_ylim([20, 120])
    else:
        ax.set_xlabel('lag time (s)')
        for j, (freq, tws) in enumerate(freq_tws.items()):
            stream = read(path_stack % (freq - 1, freq + 1))
            stream.trim2(0, 20, relative='middle')
            tr = stream[0]
            lt = np.linspace(0, 20, len(tr))
            data = tr.data
            clip = 0.02
            data[data > clip] = clip
            data[data < -clip] = -clip
            l, = ax.plot(lt, data / clip / 2 + 3 + j, lw=0.5, color=colors[j])
            ddxx = 0.7 * ((freq >= 9) + (freq >= 11))
            ddyy = 0.1 * (freq >= 11)
            ax.annotate('%d-%dHz' % (freq - 1, freq + 1), (19.5, 3 + j - ddyy),
                        va='bottom', ha='right', color=colors[j], size=7)
            print freq, ddxx
            l, = ax.plot((14.6 - ddxx, 16.3 - ddxx), 2 * [3.25 - ddyy + j], lw=1, ls=lss[j], color=colors[j])
            if steps[j]:
                l.set_dashes(steps[j])
            l.set_dash_capstyle('round')
        ax.set_ylim([2.5, 11.5])
        ax.set_yticks([])


#fig.axes[3].legend(loc='center left', bbox_to_anchor=(1.22, 0.5), ncol=2)

for ax in fig.axes:
    ax.xaxis.labelpad = 3
    ax.yaxis.labelpad = 3

fw2 = 85 / 25.4
fh2 = fw2 / 1.6
fsize2 = (fw2, fh2)
fig2 = plt.figure(figsize=fsize2)
ax2 = fig2.add_axes([0.13, 0.18, 0.83, 0.77])
from scipy.optimize import curve_fit
for j, (freq, tws) in enumerate(freq_tws.items()):
    t = np.array(t_phase[freq])
    v = np.array(vel_change[freq])
    ax2.scatter(t, v, s=4, edgecolors='none', label='%d-%dHz' % (freq - 1, freq + 1), c=colors[j])
    f = lambda x, m, y0: m * x + y0
    popt, _ = curve_fit(f, t, v, p0=(0, 0))
    ax2.plot((min(t) - 5, max(t) + 5), (f(min(t) - 5, *popt), f(max(t) + 5, *popt)), color=colors[j], ls='--')
ax2.set_xlim(0, 90)
ax2.set_ylim(0, 0.90)
ax2.set_yticks((0, 0.2, 0.4, 0.6, 0.8))
ax2.set_yticks((0.1, 0.3, 0.5, 0.7), minor=True)
ax2.set_xlabel(r'phase delay $t_{\rm p} - t_{\rm{temp}}$ (days)')
ax2.set_ylabel(r'amplitude $\epsilon_{\rm p}$ (%)')
ax2.legend(fontsize=6, labelspacing=0, frameon=False)


fig3 = plt.figure(figsize=fsize2)
ax3 = fig3.add_subplot(111)
val_eq = []
val_p = []
import numpy as np
for j, (freq, tws) in enumerate(freq_tws.items()):
    index = np.array(tws) >= 10
    val_p.append(np.mean(np.array(vel_change[freq])[index]))
    val_eq.append(np.mean(np.array(vel_temp[freq])[index]))

print freq_tws.keys()
print val_p
print val_eq
ax3.plot(freq_tws.keys(), val_p)
ax3.plot(freq_tws.keys(), val_eq)


fig.savefig(output % 'all')
fig2.savefig(output % 'absvsphase')
fig3.savefig(output % 'absvsfreq')

#plt.show()






