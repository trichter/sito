#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito.data import IPOC, getCor, getFilter, getStack
from obspy.core import UTCDateTime as UTC
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
from sito.noisexcorr import stretch as stretch_func, get_correlations
from sito.util import checkDir
import numpy as np
from progressbar import ProgressBar
import matplotlib.dates as mdates
logging.basicConfig(level=logging.DEBUG)

def _correct_dates(dates, min_delta=None):
    """Stacking after gap writes wrong date information
    """
    N = len(dates)
    deltas = np.array([dates[i + 1] - dates[i] for i in range(N - 1)])
    if min_delta is None:
        min_delta = np.median(deltas)
    indices = np.nonzero(deltas == 0)
    for index in indices[0][::-1]:
        if index < N:
            dates[index] = dates[index + 1] - min_delta
    return dates

def _insert_zeros(data, dates, min_delta=None):
    N = len(dates)
    deltas = np.array([dates[i + 1] - dates[i] for i in range(N - 1)])
    if min_delta is None:
        min_delta = np.median(deltas)
    indices = np.nonzero(deltas - min_delta >= min_delta - 0.5)
    nums = (np.round(deltas[indices] / min_delta) - 1).astype('int')
    ntws = data.shape[0]
    npts = data.shape[2]
    counter = 0
    for i in range(len(nums)):
        index = indices[0][i]
        num = nums[i]
        data = np.concatenate((data[:, :counter + index + 1, :],
                               np.zeros((ntws, num, npts)),
                               data[:, counter + index + 1:, :]),
                              axis=1)
        counter += num
    return data

mpl.rcParams.update({'font.size': 16})
maj_loc = mdates.MonthLocator()
min_loc = mdates.DayLocator((5, 10, 15, 20, 25))
maj_fmt = mdates.DateFormatter('%b %Y')
#maj_loc=None



#method = 'filter4-6_water_env2_1bit'
#method = 'filter4-6_water_env2_1bit_fft'
#method = 'filter0.01-1_water_env2_whitening_1bit_fft'

#method = '/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft'
method = 'FINAL_filter0.01-0.5_1bit_whitening'
method = 'FINAL_filter0.01-0.5_1bit_auto'
method = 'FINAL_filter1-3_1bit_auto'
method = 'FINAL_filter4-6_1bit_auto'

data = IPOC(xcorr_append='/' + method, use_local_LVC=False)
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit/stretch_t/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit/stretch2/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft/stretch_Toco/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_whitening/stretch_Toco/swcoda/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_auto/stretch_Toco/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter1-3_1bit_auto/stretch_Toco/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch_Toco/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch_Toco_PATCX/5/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch_Toco/'

#stations = 'PB01 PB02 PB03 PB04 PB05'
stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 PB15 PB16 HMBCX MNMCX PATCX PSGCX LVC'
stations = 'PATCX'
#correlations = (('PB01Z', 'PB02Z'),)
#correlations = (('PB03Z', 'PB04Z'), ('PB03Z', 'PB03Z'), ('PB04Z', 'PB04Z'))
#correlations = get_correlations(stations, 'Z', only_cross=True)
correlations = get_correlations(stations, 'Z', only_auto=True)
#correlations = (('PB01Z', 'PB02Z'), ('PB03Z', 'PB04Z'), ('PB03Z', 'PB05Z'), ('PB04Z', 'PB05Z'))

#stack = ('50days', '5days')
#stack = ('30days', '2days')
stack = ('10days', '2days')
#stack = ('10days', 'day')
stack = None

filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25), (0.05, 0.1))
filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25))
filters = (None,)

#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #1-2

# after sw phases
#filters = ((0.1, 0.25),)
#filters = (None,)
#tws = (((-133, -63, 139), 4),)*4 #PB03-PB04 left, 1-2,

# after sw
#tws = (((-150, -120, -90, 70), 50),) * 5 #1-2 
#tws = (((-150, -100, -90, 40, 100), 50),) * 5 #1-3
#tws = (((-200, -150, -120, 70, 100), 50),) * 5 #1-4
#tws = (((-150, -120, 70, 100), 50),) * 5 #1-5
#tws = (((-150, -100, 50, 100), 50),) * 5 #2-3
#tws = (((-150, -100, 50, 100), 50),) * 5 #2-4
#tws = (((-150, -130, 80, 100), 50),) * 5 #2-5
#tws = (((-150, -100, 50, 100), 50),) * 5 #3-4
#tws = (((-150, -100, 50, 100), 50),) * 5 #3-5
#tws = (((-150, -100, 50, 100), 50),) * 5 #4-5

# after sw2
#tws = (((-80, 80), 150),)*5 #4-5
#tws = (((-50, 50), 50),) * 10

#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #2-4

#tws = (((-150, -100, -80, 80, 100, 150), 50),)*5 #1-2
#tws = (((-150, -100, -50, 50, 100, 150), 50),)*5 #1-2
#tws = (((-50, 50), 150),)*4 #PB03-PB04 left

# sw
#correlations = (('PB04Z', 'PB05Z'),)
#tws = (((-18, -20, -22), 3), ((-17, -20, -23), 3), ((-15, -17, -19), 3), ((-15, -17, -19), 3)) #PB01-PB02 left
#tws = (((-36, -39), 3), ((-44, -51), 3), ((-33, -36), 3), ((-33, -36), 3)) #PB01-PB03 left
#tws = (((-52, -58), 6), ((-52, -61), 9), ((-46, -52, -58), 6), ((-46, -52), 6)) #PB01-PB04 left
#tws = (((-70, -78), 8), ((-70,), 12), ((-62, -70), 8), ((-62, -70), 8)) #PB01-PB05 left
#tws = (((-23, 27, 29, 31), 3), ((-23, 28, 32, 35), 5), ((-22, 27, 29, 31), 3), ((-22, 27, 29, 31), 3)) #PB02-PB03 right!
#tws = (((-38, -43), 5), ((-38, -44, -50), 6), ((-33, -38), 5), ((-33, -38), 5)) #PB02-PB04 left
#tws = (((-60, -65), 5), ((-57, -66), 10), ((-55, -60), 5), ((-50, -55), 5)) #PB02-PB05 left
#tws = (((-22, -24, -26), 3), ((-17, -20, -23), 3), ((-15, -17, -19), 3), ((-15, -17, -19), 3), ((-10, -1, 1), 10)) #PB03-PB04 left
#tws = (((-34, -37, -40), 4), ((-33, -38), 5), ((-30, -33, -36), 4), ((-30, -33), 4)) #PB03-PB05 left
#tws = (((-20, -22, -24, 15), 3), ((-19, -22, -25), 4), ((-17, -19, -21), 3), ((-17, -19), 3)) #PB04-PB05 left

#auto 0.05-0.5Hz
#tws = (((60, 90), 50),)

#auto 1-3
#tws = (((5, 10, 15), 5),)
#auto 4-6Hz
#tws = (((-15, -10, -6, 1, 5, 10), 5),)
#tws = (((1, 5, 10, 15), 5),)
#tws = (((1, 5, 10, 15), 5),)
#tws = (((20, 30, 40, 50), 20),)
#tws = (((10, 30, 50), 50),)

USE_TWS_AFTER_SW = False

#t1 = UTC('2007-10-01')
#t2 = UTC('2008-02-01')
t1 = UTC('2007-01-01')
t2 = UTC('2008-12-31')
t_Toco = UTC('2007-11-14 15:14:00')
#t2 = UTC('2007-09-05')
#t2 = UTC('2011-04-01')
#t2 = UTC('2009-01-01')
period = 24 * 3600
#period = 1800
add_to_file = ''
add_to_file = '_Tocevent'
add_to_file = '_2007_2008'
#add_to_file = '_end_2007'
ext = '_line_thres0.png'
MAX_CORR_LINE = True
nstr = 201
str_range = 0.02
line = ',-'
#line = '.-'
load_single_file = False
calculate = False
plot_simple = True
plot_complex = True
show = False
plot_Toco = True
cmap = 'hot_r'
#sl_mean = slice(10, 110, None)
sl_mean = slice(None, None, None)
trim_stream = 200

figsize = (8.267, 11.693 * 1.2)[::-1]
figsize2 = (20, 5)

thres_corr = 0.35
thres_dt = 0.02

thres_cor2 = 0.4
thres_cor2 = 0.
#thres_cor2 = None
lw_Toco = 1.5
#lw_Toco = 3




if calculate:
    for i in range(len(filters)):
        for correlation in ProgressBar()(correlations):
            try:
                if not load_single_file:
                    stream = data.readX(correlation, t1, t2, filter=filters[i], stack=stack)
                else:
                    stream = data.readX(correlation, t1=None, t2=None, filter=filters[i], stack=stack)
                    stream = stream.select(time=(t1, t2))
            except IOError as ex:
                print ex
                print 'Error loading file ... continue with next'
                continue
            print stream
            if len(stream) == 0:
                continue
            if USE_TWS_AFTER_SW:
                dist = data.stations.dist(correlation[0][:-1], correlation[1][:-1])
                if dist > 550:
                    continue
                direct = int(round(dist / 3.))
                tws = (((-direct - 30 - 80, direct + 30), 80),)
            stream.trim2(-trim_stream, trim_stream, 'middle')
            result = stream.stretch(reftr=None, str_range=str_range, nstr=nstr, time_windows=tws[i], sides='right')
            dates = [(time + 12 * 3600).toordinal() for time in stream.getHI('starttime')]
            dates = _correct_dates(dates)
            print np.array([dates[j + 1] - dates[j] for j in range(len(dates) - 1)])
            checkDir(path + 'bla')
            np.savez(path + 'data_stretching_%s%s%s%s' % (getCor(*correlation), getFilter(filters[i]), getStack(None, period, stack), add_to_file), tw_start=tws[i][0], tw_width=tws[i][1], dates=dates, **result)#corr=corr, stretch=stretch)

if not plot_simple and not plot_complex:
    import sys
    sys.exit()
for filter_ in filters:
    for correlation in correlations:
#        labels = ['%s_%d' % (getCor(*correlation), year) for year in range(2007, 2012)]
#        labels.append(st)
        l = '%s%s%s%s' % (getCor(*correlation), getFilter(filter_), getStack(None, period, stack), add_to_file)
        print l
        try:
            npzfile = np.load(path + 'data_stretching_%s.npz' % l)
        except:
            print 'Error for %s-%s' % correlation
            continue
        tw_start, tw_width, corr, stretch, dates, sim_mat, stretch_vec = 7 * [None]
        vars().update(npzfile)

        if plot_simple:
            dt = stretch
            dt = np.ma.masked_outside(dt, -thres_dt, thres_dt)
            #corr = np.ma.masked_less(corr, thres_corr)
            corr[corr < thres_corr] = thres_corr
            fig = plt.figure(figsize=figsize, subplotpars=mpl.figure.SubplotParams(hspace=0.07))
            plt.suptitle('stretching_%s' % l)
            ax1 = fig.add_subplot(211)
            #ax1.set_position([0.1, 0.5, 0.8, 0.4])
            for j in range(corr.shape[1]):
                #print len(dates), len(corr[:, j])
                ax1.plot_date(dates, corr[:, j], line, label='tw=%d/%ds' % (tw_start[j], tw_start[j] + tw_width))
                ax1.set_ylabel('correlation')
                ax1.legend(loc='lower right')
            ax2 = fig.add_subplot(212, sharex=ax1)
            #ax2.set_position([0.1, 0.1, 0.8, 0.4])
            for j in range(corr.shape[1]):
                ax2.plot_date(dates, dt[:, j] * 100, line)
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

        if plot_complex:
            sim_mat = _insert_zeros(sim_mat, dates)
            for i in range(sim_mat.shape[0]):
                l2 = l + ('_%s-%s' % (tw_start[i], tw_start[i] + tw_width))
                fig = plt.figure(figsize=figsize2)
                plt.suptitle('stretching_%s' % l2)
                #ax1 = fig.add_subplot(111)
                ax1 = fig.add_axes([0.06, 0.16, 0.86, 0.74])
                #ax1.set_position([0.1, 0.5, 0.8, 0.4])
                dt = 1.*(dates[-1] - dates[0]) / (len(dates) + 1) / 2
                ds = 1.*(max(stretch_vec) - min(stretch_vec)) / (len(stretch_vec) + 1) / 2
                extent = (dates[0] - dt, dates[-1] + dt, min(stretch_vec) - ds, max(stretch_vec) + ds)
                vmin = thres_cor2 if thres_cor2 is not None else 0.5 * (np.median(sim_mat[i, 0, :]) + np.median(sim_mat[i, -1, :]))
                vmin = max(0, vmin)

                image = ax1.imshow(np.transpose(sim_mat[i, :, :]), cmap=cmap, interpolation='nearest',
                           origin='lower', aspect='auto', vmax=1, vmin=vmin, extent=extent)
                            #vmax=self.vmax, vmin=self.vmin, norm=self.norm)
                if plot_Toco:
                    ax1.axvline(t_Toco.toordinal(), lw=lw_Toco, color='b')
                if MAX_CORR_LINE:
                    ax1.plot(dates, stretch[:, i], color='b', lw=2)
                ax1.xaxis_date()
                ax1.set_xlim(extent[:2])
                ax1.set_ylim(extent[2:])
                ax1.set_ylabel('stretching factor')
                if maj_loc:
                    ax1.xaxis.set_major_locator(maj_loc)
                    ax1.xaxis.set_minor_locator(min_loc)
                    ax1.xaxis.set_major_formatter(maj_fmt)
                fig.autofmt_xdate()
                ticks = ax1.get_xticks()
                if len(ticks) > 20:
                    ticks = ticks[::len(ticks) // 20]
                    ax1.set_xticks(ticks)
                cb = fig.colorbar(image, cax=fig.add_axes([0.94, 0.16, 0.01, 0.74]))
                cb.set_label('correlation')
                if not show:
                    plt.savefig(path + 'stretching_%s%s' % (l2, ext))
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

