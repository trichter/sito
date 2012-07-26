#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
from IPython import embed
from obspy.core import UTCDateTime as UTC
from operator import neg
from sito import read
from sito.stream import Stream
from termcolor import colored
import cPickle
import collections
import matplotlib as mpl
import numpy as np
import pylab as plt
import scipy.stats
import tempfile
#import sito.debug

def getind(stream, utc):
    times = stream.getHI('event.datetime')
    for i, time in enumerate(times):
        if time > utc:
            return i

#def split_stream(stream, utc):
#    i = getind(stream, utc)
#    return stream[:i], stream[i:]

PH_TYPE = {0:'positive', 1:'negative', 2:'gradient'}
def ph_type(func):
    """
    Returns phase type
    
    0 -> positive
    1 -> negative
    2 -> gradient
    """
    return 0 if  func is None else 1 if func == neg else 2


Stat = collections.namedtuple('Stat', 'mean1 std1 n1 mean2 std2 n2 mean_dif '
                              'std_dif mean_rel std_rel')
def get_stat(time, n, abs_time):
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
    """
    Calculates shifts of phases compared to reference mean trace.    
    """
    #ms1, ms2 = split_stream(ms, t_event)
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

#def colored_stat_wrapper(f):
#    def stat2(*args):
#        v = f(*args)
#        if abs(v[0]) > v[1]: print col('1', 'red'),
#        if abs(v[3]) > v[4]: print col('2', 'red'),
#        if abs(v[6]) > v[7]: print col('*', 'green', attrs=['underline', 'dark']),
#        return v
#    return stat2
#def no_wrapper(f): return(f)
##stat_wrapper = no_wrapper
#
#@colored_stat_wrapper
#def stat_description(time1, time2, time):
#    if len(time1) == 0:
#        mean1 = std1 = 0
#    else:
#        mean1 = np.mean(time1)
#        std1 = np.std(time1, ddof=1) / len(time1) ** 0.5
#    if len(time2) == 0:
#        mean2 = std2 = 0
#    else:
#        mean2 = np.mean(time2)
#        std2 = np.std(time2, ddof=1) / len(time2) ** 0.5
#    return (mean1, std1, len(time1),
#            mean2, std2, len(time2),
#            mean2 - mean1, std1 + std2,
#            (mean2 - mean1) / time * 100, (std1 + std2) / time * 100)

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

RESULTS_FILE = tempfile.gettempdir() + '/' + 'mem_results_rf_timevari.pickle'
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
            print('old correlations to mean trace:')
            print(str_stat.format(**r.shift_old.stat._asdict()))
            print('correlations to mean trace:')
            print(str_stat.format(**r.shift_cor.stat._asdict()))
            print('{:d}% of correlations are above score {:.2f} to mean trace:'
                  ''.format(r.perc, r.score))
            print(str_stat.format(**r.shift_cor_score.stat._asdict()))

def export_shift_results_to_tex_table():
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
            print('old correlations to mean trace:')
            print(str_stat.format(**r.shift_old.stat._asdict()))
            print('correlations to mean trace:')
            print(str_stat.format(**r.shift_cor.stat._asdict()))
            print('{:d}% of correlations are above score {:.2f} to mean trace:'
                  ''.format(r.perc, r.score))
            print(str_stat.format(**r.shift_cor_score.stat._asdict()))


COLORS = 'blue #808000 #008000'.split()
COLOR2 = '#800000'
SHIFT_TYPE = 'max cor cor_score'.split()
PLOT_KWARGS = dict(fmt='d', elinewidth=None, capsize=4, ms=4)
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
        station, region = stareg.split()
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


#            pl_ = ms.plotRF()
#            N = len(time)
#            pl_.ax.plot(2 * [sec], (0., N + 1), '-r')
#            pl_.ax.plot(time_m + sec, 0.5 + np.arange(N), 'gx', label='max')
#            pl_.ax.plot(time_old + sec, maxi_old + np.arange(N), 'bx', label='xcorr')
#            pl_.ax.plot(time + sec, maxi + np.arange(N), 'rx', label='correlate')
#            pl_.ax.legend()
#            plt.figure()
#            plt.subplot(211)
#            n = len(maxi)
#            plt.plot(np.arange(n), maxi, 'bo')
#            plt.plot([0, n], (cor_val,) * 2, '-r')
#            plt.subplot(212)
#            plt.plot(np.arange(n), time, 'bo')
#            plt.show()


PATH = '/home/richter/Results/IPOC/receiver/2012_mag5.5/'
TIMEVARI_PATH = PATH + 'timevari/'
T_TOCO = UTC('2007-11-14 15:14')
#load_shift_results(True)
load_shift_results()
#print_shift_results()
plot_shift_results()
#embed()
#plt.show()
