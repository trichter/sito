#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito.data import IPOC, getCor, getFilter, getStack
from obspy.core import UTCDateTime as UTC
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
from sito.noisexcorr import stretch, get_correlations
from sito.util import checkDir
import numpy as np
logging.basicConfig(level=logging.DEBUG)


#method = 'filter4-6_water_env2_1bit'
#method = 'filter4-6_water_env2_1bit_fft'
method = 'filter0.01-1_water_env2_whitening_1bit_fft'
data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit/stretch_t/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit/stretch2/'
path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft/stretch/sw/'


stations = 'PB01 PB02 PB03 PB04 PB05'
#correlations = (('PB01Z', 'PB02Z'),)
#correlations = (('PB03Z', 'PB04Z'), ('PB03Z', 'PB03Z'), ('PB04Z', 'PB04Z'))
correlations = get_correlations(stations, 'Z', only_cross=True)
#correlations = get_correlations(stations, 'Z', only_auto=True)
#correlations = (('PB01Z', 'PB02Z'), ('PB03Z', 'PB04Z'), ('PB03Z', 'PB05Z'), ('PB04Z', 'PB05Z'))
#stack = ('50days', '5days')
stack = ('30days', '2days')
#stack = ('10days', 'day')
#stack = None
filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25), (0.05, 0.1))
filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25))
filters = (None,)

#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #1-2

# after sw phases
#filters = ((0.1, 0.25),)
#filters = (None,)
#tws = (((-133, -63, 139), 4),)*4 #PB03-PB04 left, 1-2,

# after sw
#tws = (((-100, -70, 70), 50),)*5 #1-2 
#tws = (((-100, -50, 50, 100), 50),)*5 #1-3
#tws = (((-150, -100, -70, 70, 100), 50),)*5 #1-4
#tws = (((-150, -100, -80, 80, 100, 150), 50),)*5 #1-5
#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #2-3
#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #2-4
#tws = (((-150, -100, -80, 80, 100, 150), 50),)*5 #2-5
#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #3-4
#tws = (((-150, -100, -70, 70, 100, 150), 50),)*5 #3-5
#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #4-5

# after sw2
#tws = (((-80, 80), 150),)*5 #4-5
#tws = (((-50, 50), 50),) * 10


#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #2-4

#tws = (((-150, -100, -80, 80, 100, 150), 50),)*5 #1-2
#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #1-2
#tws = (((-50, 50), 150),)*4 #PB03-PB04 left

# sw
correlations = (('PB02Z', 'PB03Z'),)
#tws = (((-15, -17, -19), 3), ((-17, -20, -23), 3), ((-15, -17, -19), 3), ((-15, -17, -19), 3)) #PB01-PB02 left
#tws = (((-33, -36), 3), ((-44, -51), 3), ((-33, -36), 3), ((-33, -36), 3)) #PB01-PB03 left
#tws = (((-46, -52), 6), ((-52, -61), 9), ((-46, -52, -58), 6), ((-46, -52), 6)) #PB01-PB04 left
#tws = (((-62, -70), 8), ((-70,), 12), ((-62, -70), 8), ((-62, -70), 8)) #PB01-PB05 left
tws = (((-20, 27, 29, 31), 3), ((-23, 28, 32, 35), 5), ((-22, 27, 29, 31), 3), ((-22, 27, 29, 31), 3)) #PB02-PB03 right!
#tws = (((-33, -38), 5), ((-38, -44, -50), 6), ((-33, -38), 5), ((-33, -38), 5)) #PB02-PB04 left
#tws = (((-55, -60), 5), ((-57, -66), 10), ((-55, -60), 5), ((-50, -55), 5)) #PB02-PB05 left
#tws = (((-15, -17, -19), 3), ((-17, -20, -23), 3), ((-15, -17, -19), 3), ((-15, -17, -19), 3), ((-10, -1, 1), 10)) #PB03-PB04 left
#tws = (((-30, -33, -36), 4), ((-33, -38), 5), ((-30, -33, -36), 4), ((-30, -33), 4)) #PB03-PB05 left
#tws = (((-17, -19, -21, 15), 3), ((-19, -22, -25), 4), ((-17, -19, -21), 3), ((-17, -19), 3)) #PB04-PB05 left

#4-6Hz
#tws = (((-10, -5, -1, 1, 5, 10), 5),)
#tws = (((1, 5), 5),)

#t1 = UTC('2007-10-01')
#t2 = UTC('2008-02-01')
t1 = UTC('2007-01-01')
t2 = UTC('2009-01-01')
#t2 = UTC('2007-09-05')
#t2 = UTC('2011-04-01')
#t2 = UTC('2009-01-01')
period = 24 * 3600
#period = 1800
add_to_file = ''
#add_to_file = '_Tocevent'
add_to_file = '_2007_2008'
#add_to_file = '_end_2007'
ext = '.eps'
nstr = 2001
str_range = 0.02
line = ',-'
#line = '.-'
load_single_file = False
just_plot = False
show = False
#sl_mean = slice(10, 110, None)
sl_mean = slice(None, None, None)
trim_stream = 150

figsize = (8.267, 11.693 * 1.2)[::-1]

thres_corr = 0.35
thres_dt = 0.01

#just_one_side = 'left'
#just_one_side = 'right'
use_left_side = []
use_right_side = []

if not just_plot:
    tws_l = []
    tws_r = []
    for tw in tws:
        values = np.array(tw[0])
        tws_l.append((-values[np.nonzero(values < 0)], tw[1]))
        tws_r.append((values[np.nonzero(values >= 0)], tw[1]))
        if np.any(values < 0):
            use_left_side.append(True)
        else:
            use_left_side.append(False)
        if np.any(values >= 0):
            use_right_side.append(True)
        else:
            use_right_side.append(False)

    for i in range(len(filters)):
        for correlation in correlations:
            if use_right_side[i]:
                if not load_single_file:
                    stream = data.readX(correlation, t1, t2, filter=filters[i], stack=stack)
                else:
                    stream = data.readX(correlation, t1=None, t2=None, filter=filters[i], stack=stack)
                    stream = stream.select(time=(t1, t2))
                if len(stream) == 0:
                    continue
                stream.trim2(-trim_stream, trim_stream, 'middle')
                reftr = stream[sl_mean].calculate('mean')
                corr, dt = stretch(stream, reftr, str_range=str_range, nstr=nstr, time_windows=tws_r[i], single_side=False)
                dates = [time.toordinal() for time in stream.getHI('starttime')]
            else:
                corr = None

            if use_left_side[i]:
                if not load_single_file:
                    stream = data.readX(correlation, t1, t2, filter=filters[i], stack=stack)
                else:
                    stream = data.readX(correlation, t1=None, t2=None, filter=filters[i], stack=stack)
                    stream = stream.select(time=(t1, t2))
                dates = [time.toordinal() for time in stream.getHI('starttime')]
                stream.trim2(-trim_stream, trim_stream, 'middle')
                for tr in stream:
                    tr.data[:] = tr.data[::-1]
#                for tr in stream:
#                    if  UTC('2007-02-18') - tr.stats.starttime > 0:
#                        tr.data = np.hstack(((0,), tr.data[:-1]))
                reftr = stream[sl_mean].calculate('mean')
                corr2, dt2 = stretch(stream, reftr, str_range=str_range, nstr=nstr, time_windows=tws_l[i], single_side=False)
            else:
                corr2 = None
            if corr is None:
                corr_final = corr2
                dt_final = dt2
            elif corr2 is None:
                corr_final = corr
                dt_final = dt
            else:
                corr_final = np.vstack((corr2, corr))
                dt_final = np.vstack((dt2, dt))
            checkDir(path + 'bla')
            np.savez(path + 'data_stretching_%s%s%s%s' % (getCor(*correlation), getFilter(filters[i]), getStack(None, period, stack), add_to_file), tw_start=tws[i][0], tw_width=tws[i][1], dates=dates, corr=corr_final, dt=dt_final)

for i in range(len(filters)):
    for correlation in correlations:
#        labels = ['%s_%d' % (getCor(*correlation), year) for year in range(2007, 2012)]
#        labels.append(st)
        l = '%s%s%s%s' % (getCor(*correlation), getFilter(filters[i]), getStack(None, period, stack), add_to_file)
        print l
#        if len(set(correlation)) == 2:
#            labels.append('%s%s%s%s_left_side' % (getCor(*correlation), getFilter(filters[i]), getStack(None, period, stack), add_to_file))
#            if just_one_side == 'left':
#                labels.remove(labels[0])
#            elif just_one_side == 'right':
#                labels.remove(labels[1])
        try:
            npzfile = np.load(path + 'data_stretching_%s.npz' % l)
#                npzfile2 = np.load(path2 + 'data_stretching_%s.npz' % l)
        except:
            continue
        tw_start, tw_width, corr, dt, dates = None, None, None, None, None
        vars().update(npzfile)
        dt -= 1.
        dt = np.ma.masked_outside(dt, -thres_dt, thres_dt)
        #corr = np.ma.masked_less(corr, thres_corr)
        corr[corr < thres_corr] = thres_corr
        fig = plt.figure(figsize=figsize, subplotpars=mpl.figure.SubplotParams(hspace=0.07))
        plt.suptitle('stretching_%s' % l)

        ax1 = fig.add_subplot(211)
        #ax1.set_position([0.1, 0.5, 0.8, 0.4])
        for j in range(corr.shape[0]):
            ax1.plot_date(dates, corr[j, :], line, label='tw=%d/%ds' % (tw_start[j], tw_start[j] + np.sign(tw_start[j]) * tw_width))
            ax1.set_ylabel('correlation')
            ax1.legend(loc='lower right')
        ax2 = fig.add_subplot(212, sharex=ax1)
        #ax2.set_position([0.1, 0.1, 0.8, 0.4])
        for j in range(corr.shape[0]):
            ax2.plot_date(dates, dt[j, :]*100, line)
#                ax2.plot_date(npzfile2['dates'], npzfile2['dt'][j, :], line)
            ax2.set_ylabel('rel. time shift in %')
            ax2.set_ylim((-thres_dt * 100, thres_dt * 100))
            #plt.legend()

        fig.autofmt_xdate()
        ticks = ax1.get_xticks()
        if len(ticks) > 20:
            ticks = ticks[::len(ticks) // 20]
            ax1.set_xticks(ticks)
        if not show:
            plt.savefig(path + 'stretching_%s%s' % (l, ext))
            plt.close()
if show:
    plt.show()



#
#for t_year, st_year in streamyeargen2(stream):
#    reftr = st_year.calculate('mean') #.select(time=(UTC('2007-11-08'), UTC('2007-11-12')))
#    corr, dt = stretch(st_year, reftr, str_range=0.05, nstr=1000, time_windows=tw, single_side=False)
#    dates = [time.toordinal() for time in st_year.getHI('starttime')]
#
#    np.savez(path + 'data_stretching_PB03-PB04_%d' % t_year.year, dates=dates, corr=corr, dt=dt)




#npzfile = np.load(path + 'data_stretching_PB03_2007.npz')
#corr, dt, dates = None, None, None
#vars().update(npzfile)
#corr[corr < 0.3] = 0.3
#plt.figure()
#ax1 = plt.subplot(211)
#for i in range(corr.shape[0]):
#    plt.plot_date(dates, corr[i, :], '+', label='tw=%ds-%ds' % (tw[0][i], tw[0][i] + tw[1]))
#    plt.ylabel('correlation')
#    plt.legend()
#plt.subplot(212, sharex=ax1)
#for i in range(corr.shape[0]):
#    plt.plot_date(dates, dt[i, :], '+', label='tw=%ds-%ds' % (tw[0][i], tw[0][i] + tw[1]))
#    plt.ylabel('rel. time shift')
#    plt.legend()
#plt.show()


#for st in ('PB03',):# PB04 PB03-PB04'.split():
#    labels = ['%s_%d' % (st, year) for year in range(2007, 2012)]
#    labels.append(st)
#    for l in labels:
#        npzfile = np.load(path + 'data_stretching_%s.npz' % l)
#        corr, dt, dates = None, None, None
#        vars().update(npzfile)
#        corr[corr < 0.3] = 0.3
#        fig = plt.figure(figsize=figsize)
#        plt.suptitle('stretching_%s' % l)
#        ax1 = plt.subplot(211)
#        for i in range(corr.shape[0]):
#            ax1.plot_date(dates, corr[i, :], '+', label='tw=%ds-%ds' % (tw[0][i], tw[0][i] + tw[1]))
#            ax1.set_ylabel('correlation')
#            ax1.legend(loc='lower right')
#        plt.subplot(212, sharex=ax1)
#        for i in range(corr.shape[0]):
#            plt.plot_date(dates, dt[i, :], '+', label='tw=%ds-%ds' % (tw[0][i], tw[0][i] + tw[1]))
#            plt.ylabel('rel. time shift')
#            plt.ylim((0.98, 1.02))
#            #plt.legend()
#        fig.autofmt_xdate()
#        plt.savefig(path + 'stretching_%s%s' % (l, ext))
#        plt.close()

