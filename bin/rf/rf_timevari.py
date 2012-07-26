#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
"""
select events from one region and plot rfs
"""


from sito import read, imaging, util
import pylab as plt
from glob import glob
from IPython import embed
from obspy.core import UTCDateTime as UTC
from operator import neg
from sito.stream import Stream
from termcolor import colored
import cPickle
import collections
import matplotlib as mpl
import numpy as np
import scipy.stats
import tempfile
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.transforms
import os
#import sito.debug

def get_region_number(loc='ipoc', dist=0):
    """number of region from distance"""
    if loc == 'ipoc':
        dist = int(dist)
        return ('1' if dist == 800 else
                '2' if dist == 600 else
                '3' if dist == 1400 else
                'X')

def getFig(num=0, ratio=1.5, margin=None, **kwargs):
    """get figure with margin"""
    axes = [1. - 0.1 * num] + [0.1] * num
    margin = margin or [1.7, 0.4, 1.1, 0.3] #left, rigth, bottom, top
    fig = imaging.getFigure(axes, width=15., margin=margin, ratio=ratio,
                            fontsize=12, labelsize='small', **kwargs)
    return fig


def select_rfs_from_region():
    """Select rfs from region."""
    for file_ in glob(PATH + '*_mout.QHD'):
        ms = read(file_).select(expr='not st.mark')
        for region in REGIONS:
            ms2 = ms.select(around=region)
            print len(ms2), 'vs', len(ms), file_
            station = ms2[0].stats.station
            num = get_region_number(dist=region[2])
            ms2.setHI('eventsfrom', region)
            ms2.setHI('eventsfromfe', util.feregion(region[0], region[1]))
            newfile_ = TIMEVARI_PATH + ('rf_%s_R%s' % (station, num))
            ms2.write(newfile_, 'Q')

def cosmetic_to_plot(fig, stream):
    """set y ticklabels as none, plot event as star on the right panel"""
    def UTC2year(utc):
        import calendar
        year = utc.year
        return year + utc.julday / (365. + calendar.isleap(year))
    def getind(utc):
        times = stream.getHI('event.datetime')
        for i, time in enumerate(times):
            if time > utc:
                return i - 0.5
    fig.axes[0].set_yticklabels(' ' * 200)
    # mark Tocopilla event with red line and star
    for ax in fig.axes[:1] + fig.axes[2:4]:
        ax.axhline(getind(T_TOCO), color='r', lw=1, alpha=0.6)
    fig.axes[3].plot((UTC2year(T_TOCO),), (getind(T_TOCO),), marker=(5, 2, 0), #asterix
                      mec='r', ms=3)


def create_symlinks():
    for file_ in glob(TIMEVARI_PATH + 'rf_*_R?.QHD'):
        reg_num = file_[-5]
        station = file_[file_.rindex('_', None, -7) + 1:-7]
        source = TIMEVARI_PATH + ('plots/rf_%s_%s_R%s' %
                                    (station, 'Q', reg_num))
        link_name = source + '_shift.pdf'
        source = source + '.pdf'
        if not os.path.exists(link_name) and os.path.exists(source):
            print('Create link %s' % link_name)
            os.symlink(source, link_name)

def plot_selected_rfs():
    """Plot RF from region """
    start = -5
    end = 22
    show = False
    components = 'LQT'
    for file_ in glob(TIMEVARI_PATH + 'rf_*_R?.QHD'):
        reg_num = file_[-5]
        ms = read(file_)
        station = ms[0].stats.station
        stareg = '%s R%s' % (station, reg_num)
        if PLOT_SHIFTS_IN_AXIS and (not WINDOWS.has_key(stareg) or
                                    not RESULTS.has_key(stareg)):
            continue

        if (station in 'PB12 PB13 PB15 PB16 PATCX'.split() or station == 'PB14'
            and reg_num != '3' or station == 'PATCX' and reg_num != '2'):
            continue
        ratio = (len(ms) // 3 * 0.2 + 1.4 + 0.4 * 2.54) / 15
        ratio = min(ratio, 2.)
        kwargs = dict(show=show, scale=1,
                      plotinfo=('mean', 'azi dist', 'starttime'),
                      plotlabel=('mean', u'azi (°)/dist (°)', 'year'),
                      plotinfowhere=('top', 'right', 'right'),
                      plotinfodicts=[dict(pad=0, size=0.4),
                                     dict(pad=0.1, size=0.8),
                                     dict(pad=0.1, size=0.8)],
                      usehardticks='time',
                      box_trans='fig',
                      box_fs=10)

        fig = getFig(ratio=ratio)
        plot = ms.plotRF(start, end, fig=fig, component=components[1],
                       figtitle='station component R%s' % reg_num, **kwargs)
        cosmetic_to_plot(plot.fig, ms.select(component=components[0]))

        ### highlight phases and plot time shift results
        if PLOT_SHIFTS_IN_AXIS and WINDOWS.has_key(stareg) and RESULTS.has_key(stareg):
            for i, window in enumerate(sorted(zip(*WINDOWS[stareg]))):
                r = RESULTS[stareg][window]
                time_window = r.time
                _unused_sec, dt_window, _unused_cor, _unused_func = window
                plot.fig.axes[0].axvspan(time_window - dt_window,
                                         time_window + dt_window,
                                         facecolor=COLORS_PHASES[i],
                                         alpha=0.3, lw=0.5)
                plot.fig.axes[1].axvspan(time_window - dt_window,
                                         time_window + dt_window,
                                         facecolor=COLORS_PHASES[i],
                                         alpha=0.3, lw=0.5)
            if stareg != 'PB06 R1':
                num_win = len(RESULTS[stareg])
                fac1 = 1.1
                fac2 = 1.2
                if num_win >= 3:
                    Bbox = matplotlib.transforms.Bbox.from_bounds(0. , 0.28, 1.1 * fac1, 1.0 * fac2)
                elif num_win == 2:
                    Bbox = matplotlib.transforms.Bbox.from_bounds(0. , 0.28, 1.1 * fac1, 0.7 * fac2)
                else:
                    Bbox = matplotlib.transforms.Bbox.from_bounds(0. , 0.28, 1.1 * fac1, 0.4 * fac2)
                trans = plot.fig.dpi_scale_trans + plot.fig.transFigure.inverted()
                shift_ax = plot.fig.add_axes(matplotlib.transforms.TransformedBbox(Bbox, trans).bounds)

                    #inset_axes(plot.fig.axes[0], width=0.5, height=0.5, loc=3)
                plot_shifts(shift_ax, stareg)
        elif PLOT_SHIFTS_IN_AXIS:
            plt.close(plot.fig)
            continue

        newfile_ = TIMEVARI_PATH + ('plots/rf_%s_%s_R%s%s.pdf' %
                                    (station, components[1], reg_num,
                                     PLOT_SHIFTS_IN_AXIS * '_shift'))
        try:
            plot.fig.savefig(newfile_)
        except ValueError:
            print('Could not write ' + newfile_)
        print('Ploted file ' + newfile_)
        plt.close(plot.fig)
        if PLOT_SHIFTS_IN_AXIS:
            continue
        fig = getFig(ratio=ratio)
        kwargs.update(dict(plotinfodicts=[dict(pad=0, size=0.4),
                                          dict(pad=0.1, size=0.8),
                                          dict(pad=0.1, size=0.8)]))
        plot = ms.plotRF(start, end, fig=fig, component=components[2],
                       figtitle='station component R%s' % reg_num, **kwargs)
        cosmetic_to_plot(plot.fig, ms.select(component=components[0]))
        newfile_ = TIMEVARI_PATH + ('plots/rf_%s_%s_R%s.pdf' % (station,
                                                          components[2], reg_num))
        try:
            plot.fig.savefig(newfile_)
        except ValueError:
            print('Could not write ' + newfile_)
        print('Ploted file ' + newfile_)
        plt.close(plot.fig)

##### Determine time shifts
def getind(stream, utc):
    """get index of first trace after given time in stream"""
    times = stream.getHI('event.datetime')
    for i, time in enumerate(times):
        if time > utc:
            return i

PH_TYPE = {0:'positive', 1:'negative', 2:'gradient'}
def ph_type(func):
    """Returns phase type, 0 -> positive, 1 -> negative, 2 -> gradient"""
    return 0 if func is None else 1 if func == neg else 2

#statistics named tuple
Stat = collections.namedtuple('Stat', 'mean1 std1 n1 mean2 std2 n2 mean_dif '
                              'std_dif mean_rel std_rel')
def get_stat(time, n, abs_time):
    """get some statistics from time series"""
    time1 = time[:n]
    time2 = time[n:]
    if len(time1) == 0:
        mean1 = std1 = 0
    else:
        mean1 = np.mean(time1)
        std1 = np.std(time1, ddof=1) / len(time1) ** 0.5
    if len(time2) == 0:
        mean2 = std2 = 0
    else:
        mean2 = np.mean(time2)
        std2 = np.std(time2, ddof=1) / len(time2) ** 0.5
    return Stat(mean1, std1, len(time1),
            mean2, std2, len(time2),
            mean2 - mean1, std1 + std2,
            (mean2 - mean1) / abs_time * 100, (std1 + std2) / abs_time * 100)

PERC = 50
Shift = collections.namedtuple('Shift', 'times cors stat')
Results = collections.namedtuple('Results', 'time perc score shift_max '
                                 'shift_old shift_cor shift_cor_score')
def shift_and_correlate(ms, t_event, window):
    """Calculates shifts of phases compared to reference mean trace."""
    n1 = getind(ms, t_event)
    tr_mean = ms.calculate('mean')
    sec, dt, _unused_cor, func = window
    tr_mean_slice = Stream([tr_mean.copy()])
    tr_mean_slice.trim2(sec - dt, sec + dt, relative='ponset')
    time_mean, _ = tr_mean_slice[0].getArgMax(ret='time', spline=True, func=func)
    sec = time_mean = time_mean + sec - dt
    sec1, sec2 = sec - dt, sec + dt

    ms_slice = ms.copy()
    ms_slice.trim2(sec1, sec2, relative='ponset')
    time, maxi = ms_slice.getArgMax(ret='time', spline=True, func=func)
    shift_max = Shift(time - dt, maxi, get_stat(time - dt, n1, sec))

    cor = ms.correlate(tr_mean, dt, start=sec1, end=sec2,
                       relative='ponset')
    time, maxi = cor.getArgMax(ret='time', spline=True)
    shift_old = Shift(time - dt, maxi, get_stat(time - dt, n1, sec))

    cor2 = ms.correlate_numpy(tr_mean, start=sec1, end=sec2,
                              start2=sec1 - dt, end2=sec2 + dt,
                              relative='ponset')
    time, maxi = cor2.getArgMax(ret='time', spline=True)
    shift_cor = Shift(time - dt, maxi, get_stat(time - dt, n1, sec))

    score = scipy.stats.scoreatpercentile(shift_cor.cors, PERC)
    score = round(score, 2)
    score -= score % 0.05
    if score == 1:
        score -= 0.05
    n2 = np.count_nonzero(maxi[:n1] >= score)
    time = time[maxi >= score]
    maxi = maxi[maxi >= score]
    shift_cor_score = Shift(time - dt, maxi, get_stat(time - dt, n2, sec))
    return Results(time_mean, PERC, score, shift_max, shift_old, shift_cor,
                   shift_cor_score)

WINDOWS = {'PB01 R1': ([1.5, 5., 8.9, 15.5], [0.5] * 4, [0.95] * 4, [None] * 4),
           'PB01 R2': ([9], [0.5], [0.95], [neg]),
           'PB01 R3': ([9], [0.5], [0.95], [neg]),
           'PB02 R1': ([4.5, 6.2, 21], [0.5] * 3, [0.95] * 3, [None, neg, None]),
           'PB02 R2': ([4, 7], [0.5] * 2, [0.95] * 2, [None] * 2),
           'PB02 R3': ([4, 7], [0.5] * 2, [0.95] * 2, [None] * 2),
           'PB03 R1': ([8, 19], [0.5] * 2, [0.95] * 2, [None] * 2),
           'PB03 R2': ([7, 11.5, 16.8], [0.5] * 3, [0.95] * 3, [None, neg, None]),
           'PB03 R3': ([4.5, 7, 14], [0.5] * 3, [0.95] * 3, [None] * 3),
           'PB04 R1': ([6, 12.5, 15], [0.5] * 3, [0.95] * 3, [None, neg, None]),
           'PB04 R2': ([16], [0.5] , [0.95] , [neg]),
           'PB04 R3': ([8.5, 16], [0.5] * 2, [0.95] * 2, [None, neg]),
           'PB05 R1': ([6.5, 6.5, 12, 17.5], [0.5] * 4, [0.95] * 4,
                       [np.gradient, None, neg, None]),
           'PB05 R2': ([5, 7, 17], [0.5] * 3 , [0.95] * 3 , [None, neg, None]),
           'PB05 R3': ([5, 12], [0.5] * 2, [0.95] * 2, [None, neg]),
           'PB06 R1': ([8.5, 19], [0.5] * 2, [0.95] * 2, [None, neg]),
           'PB06 R2': ([8, 12], [0.5] * 2, [0.95] * 2, [None, neg]),
           'PB06 R3': ([7.5, 8, 12], [0.5] * 3, [0.95] * 3,
                       [np.gradient, None, neg]),
           'PB07 R1': ([7.5, 11], [0.5] * 2, [0.95] * 2, [None, neg]),
           'PB08 R1': ([7, 20], [0.5] * 2, [0.95] * 2, [None] * 2),
           'PB08 R2': ([10, 16], [0.5] * 2, [0.95] * 2, [None, neg]),
           'PB08 R3': ([10, 16], [0.5] * 2, [0.95] * 2, [None, neg])}
WINDOWS = collections.OrderedDict(sorted(WINDOWS.items(), key=lambda t: t[0]))

RESULTS_FILE = tempfile.gettempdir() + '/' + 'mem_results_rf_timevari.pkl'
RESULTS = {}
def load_shift_results(compute=False):
    global RESULTS
    if not compute:
        try:
            with open(RESULTS_FILE) as f:
                RESULTS = cPickle.load(f)
                return
        except IOError:
            pass
    RESULTS = collections.OrderedDict()
    for stareg in WINDOWS:
        RESULTS[stareg] = {}
        station, region = stareg.split()
        ms = read(TIMEVARI_PATH + ('rf_%s_%s.QHD' % (station, region)))
        ms = ms.select(component='Q')
        for window in zip(*WINDOWS[stareg]):
            RESULTS[stareg][window] = shift_and_correlate(ms, T_TOCO, window)
    with open(RESULTS_FILE, 'w') as f:
        cPickle.dump(RESULTS, f)

str_stat = ('{mean1:>8.3f} +- {std1:.3f} ({n1:d})  vs. '
            '{mean2: .3f} +- {std2:.3f} ({n2:>2d})  =>  '
            'dif = {mean_dif: .3f} +- {std_dif:.3f}, '
            'reldif = ({mean_rel: .3f} +- {std_rel:.3f})%')
def print_shift_results():
    for stareg in RESULTS:
        station, region = stareg.split()
        print
        print(colored('station %s region %s' % (station, region),
                      attrs=['reverse']))
        for window in RESULTS[stareg]:
            r = RESULTS[stareg][window]
            _unused_sec, dt, _unused_cor, func = window
            color = ['blue', 'red', 'green'][ph_type(func)]
            print colored(PH_TYPE[ph_type(func)], color),
            print(' phase at %.3fs, time window %.1fs' % (r.time, 2 * dt))
            print('maxima:')
            print(str_stat.format(**r.shift_max.stat._asdict()))
            np.set_printoptions(precision=4)
            print(r.shift_max.times)
            print(r.shift_max.cors)
            print('old correlations to mean trace:')
            print(str_stat.format(**r.shift_old.stat._asdict()))
            print('correlations to mean trace:')
            print(str_stat.format(**r.shift_cor.stat._asdict()))
            print('{:d}% of correlations are above score {:.2f} to mean trace:'
                  ''.format(r.perc, r.score))
            print(str_stat.format(**r.shift_cor_score.stat._asdict()))

str_stat_latex = ('${mean1:+.3f}\\pm {std1:.3f} \\left({n1:d}\\right)$ &\n'
                  '${mean2:+.3f}\\pm {std2:.3f} \\left({n2:d}\\right)$ &\n'
                  '${mean_dif:+.3f}\\pm {std_dif:.3f}$ &\n'
                  '${mean_rel:+.2f}\\pm {std_rel:.2f}$ \\\\\n ')
def export_shift_results_to_tex_table(just_these=None):
    tex = ''
    for stareg in RESULTS:
        if just_these and stareg not in just_these:
            continue
        tex += stareg + '\n'
        windows = RESULTS[stareg]
        for window in sorted(windows, key=windows.get(0)):
            r = RESULTS[stareg][window]
            _unused_sec, _unused_dt, _unused_cor, func = window
            phase = ph_type(func)
            phase = '*' if phase == 0 else '-' if phase == 1 else '/'
            tex += (' & $\\SI{{{:5.2f}}}s \\left({}\\right)$& maxima & \n'
                    ''.format(r.time, phase) +
                    str_stat_latex.format(**r.shift_max.stat._asdict()) +
                    ' & & correlation & \n' +
                    str_stat_latex.format(**r.shift_cor.stat._asdict()) +
                    ' & & cor. above ${:.2f}$ & \n'.format(r.score) +
                    str_stat_latex.format(**r.shift_cor_score.stat._asdict()) +
                    '\\addlinespace[2pt]\n')
    tex = tex.replace('+', '\phantom{-}').replace('*', '+')
    table_file = TEX_TABLE_FILE_SOME if just_these else TEX_TABLE_FILE
    with open(table_file, 'w') as f:
        f.write(tex)


COLORS = 'blue #808000 red'.split()
COLOR2 = '#800000'
COLORS_PHASES = 'yellow red blue green'.split()
SHIFT_TYPE = 'max cor cor_score'.split()
PLOT_KWARGS = dict(fmt='d', elinewidth=None, capsize=4, ms=4)
PLOT_KWARGS2 = dict(fmt='d', elinewidth=0., capsize=0, ms=2, lw=1)


def plot_shifts(ax, stareg, label=None):
    i = 0.
    #N = 3 * len(RESULTS[stareg]) - 1
    windows = RESULTS[stareg]
    for k, window in enumerate(sorted(windows, key=windows.get(0))):
        r = RESULTS[stareg][window]
        time_window = r.time
        #sec, _unused_dt, _unused_cor, func = window
        #color = ['blue', 'red', 'green'][ph_type(func)]
        y_pos = i #1.*(N - i) / N
        ax.annotate('%.1fs' % time_window, (0.02, y_pos), color='k',
                     xycoords=('axes fraction', 'data'), fontsize=6, va='top')
        ax.axvline(0, lw=1, color='black', zorder= -1)
        ax.axhspan(y_pos - 0.05, y_pos + 0.25, alpha=0.3, facecolor=COLORS_PHASES[k], lw=0.5)
        for j, shift_type in enumerate(SHIFT_TYPE):
            stat = getattr(r, 'shift_' + shift_type).stat
            if stat.n1 <= 1 or stat.n2 <= 1:
                continue
            # blue, cyan, light blue, green
            color = COLORS[j]
            #y_pos = 1.*(N - i) / N
            y_pos = i
            ax.errorbar(stat.mean_dif, y_pos, yerr=None, xerr=stat.std_dif,
                         ecolor=color, mec=color, mfc=color,
                         label=None if not label else shift_type,
                         **PLOT_KWARGS2)
            i += 0.1
    #ax.set_ylim((i + 0.2, -0.1))
    ax.set_ylim((i - 0.05, -0.05))
    ax.set_xticks(())
    ax.set_yticks(())
    xlims = ax.get_xlim()
    xlim_max = max(-xlims[0], xlims[1])
    xlim_max = 0.3 if xlim_max == 0.25 else 0.2 if xlim_max < 0.2 else xlim_max
    ax.set_xlim([-xlim_max, xlim_max])
    xticks = [-xlim_max * 0.5, xlim_max * 0.5]
    ax.set_xticks(xticks)
    #ax.set_xticklabels([str(xtick) + 's' for xtick in xticks])
    ax.set_xlabel('time shift (s)', fontsize=6, labelpad=2)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)
        #tick.set_pad(-10)
    for line in ax.xaxis.get_ticklines():
        line.set_markeredgewidth(1)
        line.set_markersize(4)

def plot_shift_results():
    fig = plt.figure()
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132, sharey=ax1, sharex=ax1)
    ax3 = fig.add_subplot(133, sharey=ax1)
    assert isinstance(ax1, mpl.axes.Axes) #get pydev autocompletion
    assert isinstance(ax2, mpl.axes.Axes)
    assert isinstance(ax3, mpl.axes.Axes)
    ax2.set_yticklabels(())
    ax3.set_yticklabels(())
    ax1.set_xlabel('time shift (s)/ before and after earthquake')
    ax2.set_xlabel('dif of time shift (s)')
    ax3.set_xlabel('dif of time shift (%)')

    ax1.vlines([0], [-5], [100])
    ax2.vlines([0], [-5], [100])
    ax3.vlines([0], [-5], [100])
    #for l in ax2.
    i = k = 0
    for stareg in RESULTS:
        #station, region = stareg.split()
        #print
        #print(colored('station %s region %s' % (station, region), attrs=['reverse']))
        ax1.annotate(stareg, (0.02, i), xycoords=('axes fraction', 'data'))
        for window in RESULTS[stareg]:
            r = RESULTS[stareg][window]
            sec, _unused_dt, _unused_cor, func = window
            color = ['blue', 'red', 'green'][ph_type(func)]
            ax2.annotate('%.1fs' % sec, (0.02, i), color=color,
                         xycoords=('axes fraction', 'data'))
            for j, shift_type in enumerate(SHIFT_TYPE):
                stat = getattr(r, 'shift_' + shift_type).stat
                # blue, cyan, light blue, green
                if j == 0:
                    ax3.annotate('%d/%d' % (stat.n1, stat.n2), (0.02, i),
                                 xycoords=('axes fraction', 'data'))
                color = COLORS[j]
                c2 = COLOR2
                ax1.errorbar(stat.mean1, i, yerr=None, xerr=stat.std1,
                             ecolor=color, mec=color, mfc=color, **PLOT_KWARGS)
                ax1.errorbar(stat.mean2, i + 0.1, yerr=None, xerr=stat.std2,
                             ecolor=c2, mec=c2, mfc=c2,
                             label=None if k > 0 else 'after',
                             **PLOT_KWARGS)
                ax2.errorbar(stat.mean_dif, i, yerr=None, xerr=stat.std_dif,
                             ecolor=color, mec=color, mfc=color,
                             label=None if k >= len(SHIFT_TYPE) else shift_type,
                             **PLOT_KWARGS)
                ax3.errorbar(stat.mean_rel, i, yerr=None, xerr=stat.std_rel,
                             ecolor=color, mec=color, mfc=color, **PLOT_KWARGS)
                i += 0.3
                k += 1
            i += 0.5
        i += 1
    ax1.set_ylim((-5, 100))
    ax1.legend()
    ax2.legend()
    plt.show()


PATH = '/home/richter/Results/IPOC/receiver/2012_mag5.5/'
TIMEVARI_PATH = PATH + 'timevari/'
TEX_TABLE_FILE = '/home/richter/Documents/pics/rf/rf_time_shifts_table.tex'
TEX_TABLE_FILE_SOME = '/home/richter/Documents/pics/rf/rf_time_shifts_some_table.tex'
JUST_THESE = ('PB01 R1', 'PB02 R1', 'PB03 R1')

REGIONS = ((-58.3, -22.0, 800), #1 South Sandwich
           (13.5, -92.1, 1400), #3 Mexico
            (14.1, -91.2, 600)) #2 Guatemala
T_TOCO = UTC('2007-11-14 15:14')
PLOT_SHIFTS_IN_AXIS = True

#select_rfs_from_region()
load_shift_results()

plot_selected_rfs()
PLOT_SHIFTS_IN_AXIS = False
plot_selected_rfs()
#create_symlinks()


#load_shift_results(True)

print_shift_results()
#export_shift_results_to_tex_table()
#export_shift_results_to_tex_table(just_these=JUST_THESE)
#plot_shift_results()
#embed()
#plt.show()


## PKD REGIONS
#REGIONS = ((42.6, 144.8, 200),
#           (-17.4, -174.8, 500),
#           (-30, -178, 500),
#           (34.9, 140.9, 500),
#           (8.3, -83.8, 500),
#           (51.5, -173., 650),
#           (51.3, 145.5, 860),
#           (-23.4, -70.44, 1000))
