#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito import read, Stream
from obspy.core import UTCDateTime as UTC
import numpy as np
import scipy.signal
import obspy.signal
import pylab as plt
import matplotlib as mpl
import os.path
from datetime import date
from sito import debug
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D



############# parameters#####
str_range = 0.01  ###########
#thres_corr = 0.35  ##########
#thres_cor2 = 0.4  ###########
#thres_cor2 = 0.  ############
#thres_dt = (-0.02, 0.02)  ###
#thres_dt = (-0.006, 0.006)  #


def get_ord(time, period=3600):
    day = (time + period / 2)
    day = day.toordinal() + 1.*day.hour / 24 + 1.*day.minute / 24 / 60 + 1.*day.second / 24 / 3600
    return day

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

def trim_stream():
    fname = path + '/xcorr/hour_%s_2011.QHD' % station
    ms = read(fname)
    print ms
    ms.trim(endtime=UTC('2011-05-01'))
    ms = ms[:-1]
    print ms
    ms.write(path + '/xcorr/hour_%s_2011-01', 'Q') % station

def calc_ref():
    fname = path + '/xcorr/hour_%s_%d.QHD'
    for year in range(2007, 2012):
        print year
        ms2 = read(fname % (station, year))
        ms2.trim2(-25, 25, 'middle')
        print 'mean'
        tr = ms2.calculate('mean')
        print 'write'
        tr.write(path + '/xcorr/hour_%s_ref_%d' % (station, year), 'Q')
    ms = read(path + '/xcorr/hour_%s_ref_*.QHD' % station)
    tr = ms.calculate('mean')
    tr.write(path + '/xcorr/hour_%s_ref_2007-2011' % station, 'Q')





def stretching2():
    fname = path + '/xcorr/hour_%s_%d.QHD'
    reftr = None
    #reftr = path + '/xcorr/hour_%s_ref_2008.QHD'
    if reftr:
        reftr = read(reftr % station)[0]
    days1 = 1
    days2 = 1
    year = 2007
    stream = read(fname % (station, year))
    stream.trim2(-25, 25, 'middle')
    tws = (((5, 10, 15), 5),)
    #tws = (((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18), 2),)
    t = stream[0].stats.starttime
    t2 = UTC('%d-%d-%d' % (t.year, t.month, t.day))
    if t2 < t:
        t = t2 + 24 * 3600
    i = 0
    while (t + days2 * 24 * 3600 * 0.5).year != 2012:
        sel = stream.select(time=(t, t + days2 * 24 * 3600))
        if len(sel) > 0.5 * days2 * 48:
            result = sel.stretch(reftr=reftr, str_range=str_range, nstr=101, time_windows=tws[0], sides='right')
            dates = [get_ord(time) for time in sel.getHI('starttime')]
            dates = _correct_dates(dates)
            np.savez(path + '/stretch/1day/data_stretching2_%s_%05d_5swin' % (station, i), tw_start=tws[0][0], tw_width=tws[0][1], dates=dates, **result)  # tw 5s
            #np.savez(path + '/stretch/data_stretching_stretching2_%s_%03d_5swin' % (station, i), tw_start=tws[0][0], tw_width=tws[0][1], dates=dates, **result)
            #np.savez(path + '/stretch/data_stretching_oneref_%s_%03d_5swin' % (station, i), tw_start=tws[0][0], tw_width=tws[0][1], dates=dates, **result)
        i = i + 1
        t = t + days1 * 24 * 3600
        if (t + days2 * 24 * 3600 + 20).year != year:
            year += 1
            stream = stream.select(time=(t, None))
            if year < 2012:
                st2 = read(fname % (station, year))
                st2.trim2(-25, 25, 'middle')
                stream = stream + st2


def stretching():
    fname = path + '/xcorr/hour_%s_2011-01.QHD' % station
    stream = read(fname)
    tws = (((5, 10, 15), 5),)
    #tws = ((range(21), 2),)
    stream.trim2(-25, 25, 'middle')
    #sel = stream.select(time=(None, UTC('%d-05-01' % year)))
    sel = stream
    result = sel.stretch(str_range=str_range, nstr=101, time_windows=tws[0], sides='right')
    dates = [get_ord(time) for time in sel.getHI('starttime')]
    dates = _correct_dates(dates)
    np.savez(path + '/stretch/data_stretching_%s_2011_01-04' % station, tw_start=tws[0][0], tw_width=tws[0][1], dates=dates, **result)  #corr=corr, stretch=stretch)

def calc_temp():
    station = 'PB12'
    channel = 'WKI'
    from sito.data import IPOC
    ipoc = IPOC()
    stream = ipoc.getChannelFromClient('2006-01-01', '2012-01-01',
                                               station=station, channel=channel)
    stream2 = stream.copy()
    day = stream[0].stats.starttime
    day2 = UTC(day.year, day.month, day.day)
    if day2 < day:
        day = day2 + 24 * 3600
    stream2.trim(day, day + 24 * 3600)
    stream2.merge()
    data = []
    while day < stream[-1].stats.endtime:
        st = stream.slice(day, day + 24 * 3600)
        st.merge()
        if len(st) == 0 or len(st[0].data) < 8640:
            print 'skip %s' % day.date
            day = day + 24 * 3600
            continue
        a = st[0].data[:8640]
        if channel == 'WKI':
            a = a / 100.
            a = np.ma.masked_outside(a, 0, 50)
        elif channel == 'WDI':
            a = a / 10000.
            a = np.ma.masked_outside(a, 0, 5)
        else:
            a = np.ma.masked_outside(a, 0, 100)
        a = np.ma.masked_invalid(a)
        data.append(a)
        day = day + 24 * 3600
    #from IPython import embed
    #embed()
    stream2[0].data = np.ma.mean(data, axis=0)
    stream2.write(path + '/climate/%s_%s' % (station, channel), 'Q')
    stream2.plot(method='full')
    plt.show()
#    data = []
#    dates = []
#    for day in streamdaygen(stream):
#        day.merge()
#        data.append(np.mean(day[0].data))
#        st = day[0].stats.starttime
#        et = day[0].stats.endtime
#        dates.append(st + (et - st) / 2.)
#    np.savez(datafile % (station, channel), dates=dates, data=data)

def analyze1():
    npzfile = np.load(path + '/stretch/data_stretching_%s_2010_01-04.npz' % station)
    tw_start, tw_width = npzfile['tw_start'], npzfile['tw_width']
    corr, stretch, dates = npzfile['corr'], npzfile['stretch'], npzfile['dates']
    date0 = dates[0]
    res = []
    res_co = []
    N = corr.shape[1]
    for j in range(N):
        sts = []
        cos = []
        for i in range(int(dates[-1] - date0)):
            select = (date0 + i <= dates) * (dates <= date0 + i + 1)
            st = stretch[select, j]
            st = scipy.signal.detrend(st, type='constant')
            #st = obspy.signal.detrend.simple(st)
            sts.append(st)
            cos.append(corr[select, j])
        res.append(np.array(sts))
        res_co.append(np.array(cos))
    plotdates = dates[dates <= date0 + 1]
    #plotdates = plotdates[:48]
    fig = plt.figure()
#    ax1 = fig.add_subplot(111)
    ax1 = fig.add_axes([0.1, 0.1, 0.67, 0.8])
    for j in range(N):
        color = cmap(1.* j / (N - 1))
        try:
            mm = np.mean(res[j], axis=0)
        except ValueError as ex:
            print j, len(res[j]), ex
            res[j] = [resi for resi in res[j] if len(resi) == 49]
            print len(res[j])
            mm = np.mean(res[j], axis=0)
        ax1.plot(plotdates, -mm, color=color, label='%ds-%ds' % (tw_start[j], tw_start[j] + tw_width))
    ax1.errorbar(plotdates, -np.mean(res[0], axis=0), yerr=np.std(res[0], axis=0, ddof=1) / len(res[0]) ** 0.5)
    #ax1.plot(plotdates, -np.mean(res[1], axis=0), 'g')
    #ax1.errorbar(plotdates, -np.mean(res[1], axis=0), yerr=np.std(res[1], axis=0, ddof=1) / len(res[1]) ** 0.5)
    #ax1.plot(plotdates, -np.mean(res[2], axis=0), 'r')
    #ax1.errorbar(plotdates, -np.mean(res[2], axis=0), yerr=np.std(res[2], axis=0, ddof=1) / len(res[2]) ** 0.5)

    #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#    ax2 = fig.add_subplot(122)
#    try:
#        mm = np.mean(res_co[0], axis=0)
#    except ValueError as ex:
#        print ex
#        res_co[0] = [resi for resi in res_co[0] if len(resi) == 49]
#        mm = np.mean(res_co[0], axis=0)
#    ax2.plot_date(plotdates, mm, 'b')
    #ax2.plot(plotdates, np.mean(res_co[1], axis=0), 'g')
    #ax2.plot(plotdates, np.mean(res_co[2], axis=0), 'r')
    fig.autofmt_xdate()
    ax3 = ax1.twinx()
    clim = read(path + '/climate/PB01_WKI_2010-01.QHD')[0].data
    ax3.plot_date(np.linspace(plotdates[0], plotdates[-1], len(clim)), clim, 'y')
    #clim = read(path + '/PB10_WKI_2011-01.QHD')[0].data
    #ax3.plot_date(np.linspace(plotdates[0], plotdates[-1], len(clim)), clim, 'b')
    #from IPython import embed
    #embed()
    fig.savefig(path + '/plots/hour_%s_2011_01-04_dv_vs_temp_PB01.pdf' % station)
    #plt.show()


def analyze2():
    npzfile = np.load(path + '/stretch/data_stretching_%s_2010.npz' % station)
    year = 2010
    tw_start, tw_width = npzfile['tw_start'], npzfile['tw_width']
    corr, stretch, dates = npzfile['corr'], npzfile['stretch'], npzfile['dates']
    ddays = 30
    ddays2 = 120
    date0 = min(dates[0], dates[48] - 1)
    for yearpart in range(2):
        res = []
        res_co = []
        N = corr.shape[1]
        for j in range(N):
            sts = []
            cos = []
            for i in range(ddays2):
                select = (date0 + i + yearpart * ddays <= dates) * (dates <= date0 + i + 1 + yearpart * ddays)
                st = stretch[select, j]
                st = scipy.signal.detrend(st, type='constant')
                #st = obspy.signal.detrend.simple(st)
                sts.append(st)
                cos.append(corr[select, j])
            res.append(np.array(sts))
            res_co.append(np.array(cos))
        plotdates = dates[dates <= date0 + 1]
        fig = plt.figure()
        #ax1 = fig.add_subplot(111)
        ax1 = fig.add_axes([0.1, 0.1, 0.67, 0.8])
        #import IPython
        #IPython.embed()
        for j in range(N):
            color = cmap(1.* j / (N - 1))
            try:
                m = np.mean(res[j], axis=0)
                co = np.mean(res_co[j], axis=0)
            except ValueError as ex:
                print j, len(res[j]), ex
                res[j] = [resi for resi in res[j] if len(resi) == 49]
                res_co[j] = [resi for resi in res_co[j] if len(resi) == 49]
                print len(res[j])
                m = np.mean(res[j], axis=0)
                co = np.mean(res_co[j], axis=0)
            ax1.plot(plotdates, co, color=color, label='%ds-%ds' % (tw_start[j], tw_start[j] + tw_width))
        #ax1.errorbar(plotdates, -np.mean(res[0], axis=0), yerr=np.std(res[0], axis=0, ddof=1) / len(res[0]) ** 0.5)
        #ax1.plot(plotdates, -np.mean(res[1], axis=0), 'g')
        #ax1.errorbar(plotdates, -np.mean(res[1], axis=0), yerr=np.std(res[1], axis=0, ddof=1) / len(res[1]) ** 0.5)
        #ax1.plot(plotdates, -np.mean(res[2], axis=0), 'r')
        #ax1.errorbar(plotdates, -np.mean(res[2], axis=0), yerr=np.std(res[2], axis=0, ddof=1) / len(res[2]) ** 0.5)

        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #ax2 = fig.add_subplot(122)
        #ax2.plot_date(plotdates, np.mean(res_co[0], axis=0), 'b')
        #ax2.plot(plotdates, np.mean(res_co[1], axis=0), 'g')
        #ax2.plot(plotdates, np.mean(res_co[2], axis=0), 'r')
        fig.autofmt_xdate()
        ax3 = ax1.twinx()
        clim = read(path + '/climate/PB01_WKI_2010-01.QHD')[0].data
        ax3.plot_date(np.linspace(plotdates[0], plotdates[-1], len(clim)), clim, 'y')
        #clim = read(path + '/PB10_WKI_2011-01.QHD')[0].data
        #ax3.plot_date(np.linspace(plotdates[0], plotdates[-1], len(clim)), clim, 'b')
        #from IPython import embed
        #embed()
        fig.savefig(path + '/plots/bla_hour_%s_%d-%02d_correlation.png' % (station, year, yearpart))
        #plt.show()

def analyze3():
    fname = path + '/stretch/data_stretching2_%s_%03d.npz'
    k = 0
    while k <= 60:
        if not os.path.exists(fname % (station, k)):
            k += 1
            continue
        print 'k', k
        npzfile = np.load(fname % (station, k))
        tw_start, tw_width = npzfile['tw_start'], npzfile['tw_width']
        corr, stretch, dates = npzfile['corr'], npzfile['stretch'], npzfile['dates']
        date0 = dates[0]
        res = []
        res_co = []
        N = corr.shape[1]
        for j in range(N):
            sts = []
            cos = []
            for i in range(int(dates[-1] - date0)):
                select = (date0 + i <= dates) * (dates <= date0 + i + 1)
                st = stretch[select, j]
                st = scipy.signal.detrend(st, type='constant')
                #st = obspy.signal.detrend.simple(st)
                sts.append(st)
                cos.append(corr[select, j])
            res.append(np.array(sts))
            res_co.append(np.array(cos))
        plotdates = dates[dates <= date0 + 1]
        #plotdates = plotdates[:48]
        fig = plt.figure()
    #    ax1 = fig.add_subplot(111)
        ax1 = fig.add_axes([0.1, 0.1, 0.67, 0.8])
        for j in range(N):
            color = cmap(1.* j / (N - 1))
            try:
                mm = np.mean(res[j], axis=0)
            except ValueError as ex:
                print k, len(res[j]), ex
                res_new = [resi for resi in res[j] if len(resi) == 49]
                if len(res_new) <= len(res[j]) // 2:
                    plotdates = dates[dates < date0 + 1]
                    res_new = [resi for resi in res[j] if len(resi) == 48]
                res[j] = res_new
                print len(res[j])
                mm = np.mean(res[j], axis=0)
            ax1.plot(plotdates, -mm, color=color, label='%ds-%ds' % (tw_start[j], tw_start[j] + tw_width))
        ax1.errorbar(plotdates, -np.mean(res[0], axis=0), yerr=np.std(res[0], axis=0, ddof=1) / len(res[0]) ** 0.5)
        date1 = date.fromordinal(int(round(dates[0])))
        date2 = date.fromordinal(int(round(dates[-1])))
        ax1.set_title('%s - %s' % (date1, date2))
        ax1.set_ylim([-0.0015, 0.0015])
        #ax1.plot(plotdates, -np.mean(res[1], axis=0), 'g')
        #ax1.errorbar(plotdates, -np.mean(res[1], axis=0), yerr=np.std(res[1], axis=0, ddof=1) / len(res[1]) ** 0.5)
        #ax1.plot(plotdates, -np.mean(res[2], axis=0), 'r')
        #ax1.errorbar(plotdates, -np.mean(res[2], axis=0), yerr=np.std(res[2], axis=0, ddof=1) / len(res[2]) ** 0.5)

        #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #    ax2 = fig.add_subplot(122)
    #    try:
    #        mm = np.mean(res_co[0], axis=0)
    #    except ValueError as ex:
    #        print ex
    #        res_co[0] = [resi for resi in res_co[0] if len(resi) == 49]
    #        mm = np.mean(res_co[0], axis=0)
    #    ax2.plot_date(plotdates, mm, 'b')
        #ax2.plot(plotdates, np.mean(res_co[1], axis=0), 'g')
        #ax2.plot(plotdates, np.mean(res_co[2], axis=0), 'r')
        fig.autofmt_xdate()
        ax3 = ax1.twinx()
        clim = read(path + '/climate/PB01_WKI_2010-01.QHD')[0].data
        ax3.plot_date(np.linspace(plotdates[0], plotdates[-1], len(clim)), clim, 'y')
        #clim = read(path + '/PB10_WKI_2011-01.QHD')[0].data
        #ax3.plot_date(np.linspace(plotdates[0], plotdates[-1], len(clim)), clim, 'b')
        #from IPython import embed
        #embed()
        fig.savefig(path + '/plots/hour_%s2_%03d.png' % (station, k))
        k += 1
        #plt.show()

def fit_temp_to_vel(temp_dates, temp, vel_dates, vel, phase=False):
    temp = np.tile(temp[:-1], 3)
    temp_dates = np.hstack((temp_dates[:-1] - 1, temp_dates[:-1], temp_dates[:-1] + 1))
    f = InterpolatedUnivariateSpline(temp_dates, temp)
    dtemp = max(f(temp_dates)) - min(f(temp_dates))
    print 'dtemp %.2fK' % dtemp
    if phase:
        fopt = lambda x, a, b, c: abs(a) * (f(x - abs(c) / 24) - abs(b))
        popt, _ = curve_fit(fopt, vel_dates, vel, p0=(0.01, 18., 5.))
    else:
        fopt = lambda x, a, b: abs(a) * (f(x) - abs(b)) / dtemp
        popt, _ = curve_fit(fopt, vel_dates, vel, p0=(0.01, 18))
    return popt, fopt

def analyze4():
    fname = path + '/stretch/1day/data_stretching2_%s_%05d_5swin.npz'
    #fname = path + '/stretch/data_stretching2_%s_%03d_5swin.npz'
    smwin = '2s' in fname
    skip_some = False  #13  #False
    k = 0
    temp_color = 'r'
    num_hour = 48  # 48 or 49
    #res = [[], [], []]
    #res_co = [[], [], []]
    #sts = [[], [], []]
    #cos = [[], [], []]
    #res = []
    #res_co = []
    #sts = []
    #cos = []
    num_nans = 0
    while k <= 1715:
        if not os.path.exists(fname % (station, k)):
            k += 1
            continue
        print 'k', k
        with np.load(fname % (station, k)) as npzfile:
            tw_start, tw_width = npzfile['tw_start'], npzfile['tw_width']
            corr, stretch, dates = npzfile['corr'], npzfile['stretch'], npzfile['dates']
        date0 = dates[0]
        N = corr.shape[1]
        for j in range(N):
            for i in range(int(dates[-1] - date0) + 1):
                select = (date0 + i <= dates) * (dates <= date0 + i + 1)
                st = stretch[select, j]
                st = scipy.signal.detrend(st, type='constant')
                if np.any(np.isnan(st)) or np.any(np.isnan(corr[select, j])):
                    num_nans += 1
                    print 'Nan value in corr %d' % num_nans
                    continue
                #if len(sts) = 0:
                #    for l in range(N):
                #        sts.append([])
                try:
                    sts[j].append(st)  #@UndefinedVariable
                    cos[j].append(corr[select, j])  #@UndefinedVariable
                except NameError:
                    sts = [[] for l in range(N)]
                    cos = [[] for l in range(N)]
                    res = [[] for l in range(N)]
                    res_co = [[] for l in range(N)]
                    sts[j].append(st)  #@UndefinedVariable
                    cos[j].append(corr[select, j])
            #res[j].append(np.array(sts))
            #res_co[j].append(np.array(cos))
        print len(sts[0])
        if k == 1:
            plotdates = dates[dates <= date0 + 1]
        k += 1
        #plotdates = plotdates[:48]
    fig = plt.figure()
    #lw = 2
#    ax1 = fig.add_subplot(111)
    plotdates = plotdates + 0.5 / 24  # to get middle of one hour intervals
    ind = np.count_nonzero(plotdates >= round(plotdates[-1]))
    plotdates = np.hstack((plotdates[-ind:] - 1, plotdates[:-ind + 1]))  # to get hours from 0:00 to 0:00
    plotdates = plotdates - 0.00039352  #30s to set date to exact hour
    print 'ind, plotdates', ind, plotdates, len(plotdates)
    ax1 = fig.add_axes([0.1, 0.15, 0.53, 0.79] if not smwin else [0.1, 0.15, 0.4, 0.8])
    if skip_some:
        N = skip_some
    for j in range(N):  # different time windows
        color = cmap(1.* j / (N - 1))
        if len(set(s.shape for s in sts[j])) != 1:
            print len(sts[j]), '- different shapes found!'
            import collections
            print collections.Counter(len(s) for s in sts[j])
            sts[j] = [resi for resi in sts[j] if len(resi) == num_hour]
            print len(sts[j]), 'after deleting all other shapes'
            #if len(res_new) <= len(res[j]) // 2:
            #    plotdates = dates[dates < date0 + 1]
            #    res_new = [resi for resi in res[j] if len(resi) == 48]
        yerr = np.ma.std(sts[j], axis=0, ddof=1) / len(sts[j]) ** 0.5 * 100.
        yerr = np.hstack((yerr[-ind:], yerr[:-ind + 1]))
        sts[j] = -100.*np.ma.mean(sts[j], axis=0)
        sts[j] = np.hstack((sts[j][-ind:], sts[j][:-ind + 1]))
        ax1.plot(plotdates, sts[j], color=color, label='%ds-%ds' % (tw_start[j], tw_start[j] + tw_width))
        if not smwin:
            ax1.errorbar(plotdates, sts[j], color=color, yerr=yerr)
    print 'bla'
    ylim = [-0.12, 0.12] if not smwin else [-0.15, 0.15]
    ax1.set_ylim(ylim)
    handles, labels = ax1.get_legend_handles_labels()
    labels.append('temperature')
    handles.append(Line2D((0,), (0,), color=temp_color))
    labels.append('shifted+scaled')
    handles.append(Line2D((0,), (0,), color=temp_color, dashes=(3, 3)))
    # JGR ax1.legend(handles, labels, bbox_to_anchor=(1.19, 1), loc=2, borderaxespad=0., frameon=True, prop={'size':'small'})
    ax1.legend(handles, labels, bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0., frameon=True)
    #ax1.legend(bbox_to_anchor=(1.14, 1), loc=2, borderaxespad=0)

    fig.autofmt_xdate()


    ax3 = ax1.twinx()
    #clim2 = np.load(open('/home/richter/Data/climate/climate_hour_IQ_mean.npz'))
    clim2 = np.load(open('/home/richter/Data/climate/climate_CHO_hui_mean.npz'))
    clim2_date = clim2['date'] - clim2['date'][0] + round(plotdates[0])
    ax3.plot_date(clim2_date, clim2['temp'], '-', color=temp_color)
    #ax3.plot_date(clim2_date, clim2['temp'] + clim2['temp_std'], '-', color='b')
    #ax3.plot_date(clim2_date, clim2['temp'] - clim2['temp_std'], '-', color='b')
    popt, fopt = fit_temp_to_vel(clim2_date, clim2['temp'], plotdates, sts[0], phase=True)
    popt2, fopt2 = fit_temp_to_vel(clim2_date, clim2['temp'], plotdates, sts[1], phase=True)
    popt3, fopt3 = fit_temp_to_vel(clim2_date, clim2['temp'], plotdates, sts[2], phase=True)
    f = lambda x: x / popt[0] + popt[1]
    print 'popt', popt  #   phase shift  0.3h [ 0.0501  18.05   0.298 ] # CHO temp 1.09h  0.19%
    print 'popt2', popt2  # phase shift  2.7h [ 0.0101  18.06   2.697 ] #          3.5h  0.038%
    print 'popt3', popt3  # phase shift 13.2h [ 0.0074  18.04  13.236 ] #          13.9h  0.028%
    print [f(ylim[0]), f(ylim[1])]
    print fopt(plotdates, *popt)
    ax1.plot(plotdates, fopt(plotdates, *popt), color=temp_color, dashes=(3, 3))
    ax1.plot(plotdates, fopt2(plotdates, *popt2), color=temp_color, dashes=(3, 3))
    ax1.plot(plotdates, fopt3(plotdates, *popt3), color=temp_color, dashes=(3, 3))
    ax1.set_xlim(left=plotdates[0])
    ax3.set_ylim(sorted([f(ylim[0]), f(ylim[1])]))
    xlc = mpl.dates.HourLocator(np.arange(12) * 2)
    ax3.xaxis.set_major_locator(xlc)
    xfmt = mpl.dates.DateFormatter('%H:%M')
    ax3.xaxis.set_major_formatter(xfmt)
    #for l in ax3.get_yticklabels():
    #    l.set_color('r')
    ax3.set_ylabel(u'temperature (°C)')


    ax1.set_title('velocity change (%)')

    # plot correlations
    ax2 = fig.add_axes([0.82, 0.15, 0.15, 0.47])
    for j in range(N):
        color = cmap(1.* j / (N - 1))
        try:
            mm = np.mean(cos[j], axis=0)
        except ValueError as ex:
            cos[j] = [resi for resi in cos[j] if len(resi) == num_hour]
            mm = np.mean(cos[j], axis=0)
        #print mm.shape
        mm = np.hstack((mm[-ind:], mm[:-ind + 1]))
        #print len(plotdates), mm.shape
        ax2.plot(plotdates, mm, color=color)

    ax2.xaxis_date()
    ax2.set_xlim([plotdates[0], plotdates[-1]])
    #ax1.set_xlim([plotdates[0], plotdates[-1]])
    xlc = mpl.dates.HourLocator((4, 12, 20))
    xlc2 = mpl.dates.HourLocator(range(24))
    ax2.xaxis.set_major_locator(xlc)
    ax2.xaxis.set_minor_locator(xlc2)
    #ax2.set_yticks((0.45, 0.5, 0.55, 0.6))
    xfmt = mpl.dates.DateFormatter('%H:%M')
    ax2.xaxis.set_major_formatter(xfmt)
    for label in ax2.get_xticklabels():
        label.set_rotation_mode('anchor')
        label.set_ha('right')
        label.set_rotation(30)
    xlc = mpl.dates.HourLocator(range(0, 24, 3))
    ax1.xaxis.set_major_locator(xlc)
    ax1.xaxis.set_minor_locator(xlc2)

    for label in ax1.get_xticklabels():
        label.set_rotation_mode('anchor')

    ax2.set_title('correlation')
    ax1.set_xlabel('UTC (Chile Time +4h)')
    ax2.set_xlabel('UTC')

    #ax1.xaxis.labelpad = 3
    #ax2.xaxis.labelpad = 3
    #ax3.yaxis.labelpad = 2
    fig.savefig(path + '/plots/hour_%s_TCHO_5swin_1day.pdf' % station)
    #fig.savefig(path + '/plots/hour_%s_TCHO_5swin.pdf' % station)
    print 'num_nans: %s vs len %s' % (num_nans, len(sts[j]))
    #plt.show()

def smooth2(x, window_len=11):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len < 3:
                return x
        s = np.r_[2 * x[0] - x[window_len - 1::-1], x, 2 * x[-1] - x[-1:-window_len:-1]]
        w = np.ones(window_len, 'd')
        y = np.convolve(w / w.sum(), s, mode='same')
        return y[window_len:-window_len + 1]

def smooth(x, y, winlen=11):
    N = len(x) // winlen
    means = np.zeros(N) * 1.0
    ys = np.zeros(N) * 1.0
    for i in range(N):
        part = x[i * winlen:(i + 1) * winlen]
        means[i] = np.ma.mean(part)
        ys[i] = y[int((i + 0.5) * winlen)]
    return means, ys


def analyze_year():
    fname = path + '/stretch/data_stretching_oneref_%s_%03d_5swin.npz'
    k = 0
    corrs = []
    datess = []
    stretchs = []
    while k <= 60:
        if not os.path.exists(fname % (station, k)):
            k += 1
            continue
        print 'k', k
        npzfile = np.load(fname % (station, k))
        corr, stretch, dates = npzfile['corr'], npzfile['stretch'], npzfile['dates']
        corr = corr[:, 1]
        stretch = stretch[:, 1]
        corrs.append(corr)
        datess.append(dates)
        stretchs.append(stretch)
        k += 1
    datess = np.hstack(datess)
    corrs = np.hstack(corrs)
    stretchs = np.hstack(stretchs)
    stretchs = smooth2(stretchs, 48 * 3)
    #stretchs, datess = smooth(stretchs, datess, 48)
    plt.plot(datess, stretchs)
    plt.gcf().autofmt_xdate()
    xfmt = mpl.dates.DateFormatter('%y-%m-%d %H:%M')
    plt.gca().xaxis.set_major_formatter(xfmt)
    plt.savefig(path + '/plots/year_%s_TIQ_5swin.pdf' % station)
    plt.show()

path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto_hour'
station = 'PATCXZ'

cmap = mpl.cm.get_cmap('hot')
cmap = mpl.cm.get_cmap('jet')
cmap = lambda i: 'k' if i < 0.5 else 'g' if i < 0.8 else 'y'

#trim_stream()
#calc_temp()

#stretching()
#analyze2()
#calc_ref()

#stretching2()
#JGR
fw = 85 / 25.4
fh = 0.9 * fw
fw = 1.8 * 85 / 25.4
fh = 0.8 * fw
fsize = (fw, fh)
mpl.rcParams.update({'figure.figsize': fsize, 'font.size': 9, 'lines.linewidth':1})
analyze4()
#analyze_year()








#def analyze_out():
#    fname = path + '/xcorr/hour_PATCXZ_2009-01.QHD'
#    stream = read(fname)
#    tws = (((5, 10, 15), 5),)
#    if stack_before:
#        stream2 = stream.copy()
#        while True:
#            date0 = stream[0].starttime
#            to_stack = stream.select()
#            stream.select('abs(st.starttime-%r) % 24*3600 < 10' % (date0,))
#
#    endd = stream[-1].stats.starttime
#    startd = stream[0].stats.starttime
#    res = [[], [], []]
#    res_co = [[], [], []]
#    for i in range(int((endd - startd) / 24 / 3600)):
#        sstream = stream.select(time=(startd + i * 24 * 3600, startd + (i + 1) * 24 * 3600))
#        result = sstream.stretch(str_range=str_range, nstr=101, time_windows=tws[0], sides='right')
#        stretch = result['stretch']
#        corr = result['corr']
#        for j in range(3):
#            res[j].append(scipy.signal.detrend(stretch[:, j]))
#            res_co[j].append(corr[:, j])
#    plotdates = np.arange(48) / 2.
#    plt.plot(plotdates, np.mean(res[0], axis=0), 'b')
#    plt.plot(plotdates, np.mean(res[1], axis=0), 'g')
#    plt.plot(plotdates, np.mean(res[2], axis=0), 'r')
#
#    plt.show()
#    from IPython import embed
#    embed()
