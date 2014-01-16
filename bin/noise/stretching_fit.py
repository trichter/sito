#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from obspy.core import UTCDateTime as UTC
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.dates as mdates

from scipy.optimize import curve_fit, fmin_l_bfgs_b
from obspy.core.event import readEvents
from sito.noisexcorr import get_correlations
from sito.data import getCor
from scipy.optimize import fmin_tnc
from scipy.optimize.optimize import fmin




mpl.rcParams.update({'font.size': 16})
#maj_loc = mdates.MonthLocator()
#min_loc = mdates.DayLocator((5, 10, 15, 20, 25))
#maj_fmt = mdates.DateFormatter('%b %Y')
maj_loc = mdates.YearLocator()
min_loc = mdates.MonthLocator()
maj_fmt = mdates.DateFormatter('%Y')
#maj_loc=None


path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto_3C/stretch3_10s/'
path = '/home/richter/Results/IPOC/xcorr/PAT_filter9-11/stretch2/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter1-3_1bit_auto/stretch'

FIT = 'sinus_exp_alt'  # 2 = without epsilon perm
stations = 'PATCX'
stations = 'PB05'
combs = get_correlations(stations, 'Z', only_auto=True)

assert FIT in ('sinus', 'exp', 'sinus_exp', 'sinus_exp_alt', 'sinus_exp_alt2', 'events')

####### Fit Sinus
if FIT == 'sinus':
    filename = path + '/data_corrected_stretching_%s.npz'
    output = path + '/data_stretching_%s_without_sinus.npz'
    output_fit = path + '/sinus_%s_%ds.npz'
elif FIT == 'exp':
    event_file = '/home/richter/Data/events/2007-2012_events_in_chile_geofon.xml'
    filename = path + '/data_stretching_%s_without_sinus.npz'
    output = path + '/data_stretching_%s_without_exp.npz'
    output_fit = path + '/exp_%s_%ds.npz'
    events = readEvents(event_file)
    events = events.filter('magnitude > 7.5')
elif FIT == 'sinus_exp':
    event_file = '/home/richter/Data/events/2007-2012_events_in_chile_geofon.xml'
    filename = path + '/data_stretching_%s.npz'  # data_ vs fata_corrected
    output = path + '/data_stretching_%s_without_sinus_exp.npz'
    output_fit = path + '/sinus_exp_%s_%ds.npz'
    events = readEvents(event_file)
    events = events.filter('magnitude > 7.5')
elif FIT.startswith('sinus_exp_alt'):
    event_file = '/home/richter/Data/events/2007-2012_events_in_chile_geofon.xml'
    filename = path + '/data_stretching_%s.npz'
    #filename = path + '/data_stretching_%s_altref.npz'  # nix vs _altref
    output = path + '/data_stretching_%s_without_sinus_exp_alt.npz'
    output_fit = path + '/sinus_exp_alt_%s_%ds.npz'
    events = readEvents(event_file)
    events = events.filter('magnitude > 7.5')
else:
    assert False

popts = {}
for k, comb in enumerate(combs):
    comb = getCor(*comb)
    print
    print comb
    npzfile = np.load(filename % comb)
    tw_start, tw_width, corr, stretch, dates, sim_mat, stretch_vec = 7 * [None]

    vars().update(npzfile)
    stretch = np.ma.masked_where(corr <= 0, stretch)
    stretch = np.ma.masked_where(np.isnan(corr), stretch)
    corr = np.ma.masked_where(corr <= 0, corr)
    corr = np.ma.masked_where(np.isnan(corr), corr)

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharex=ax1)
    ax1.plot(dates, stretch)
    ax1.xaxis_date()
    ax1.xaxis.set_major_locator(maj_loc)
    ax1.xaxis.set_minor_locator(min_loc)
    ax1.xaxis.set_major_formatter(maj_fmt)

    USE_SIGMA = True

    #ax.plot(dates, 0.003 * np.sin((dates + 80) / 365.*2 * np.pi))

    t_Toco = UTC('2007-11-14 15:14:00')

    def exps(evts, t, *args):
        return sum(
            [np.hstack((np.zeros(np.count_nonzero(t < evt)),  # Heaviside function
                        args[3 * i] +  # permanent vel drop
                        abs(args[3 * i + 1]) *  # temporary vel drop
                        np.exp(-1.*(t[t >= evt] - evt) / abs(args[3 * i + 2]) * np.log(10))))  # exponential recovery
            for i, evt in enumerate(evts)])

    def exps2(evts, t, *args):  #without vel perm
        return sum(
            [np.hstack((np.zeros(np.count_nonzero(t < evt)),  # Heaviside function
                        #args[2 * i] +  # permanent vel drop
                        abs(args[2 * i]) *  # temporary vel drop
                        np.exp(-1.*(t[t >= evt] - evt) / abs(args[2 * i + 1]) * np.log(10))))  # exponential recovery
            for i, evt in enumerate(evts)])

    def sinus(t, vel_offset, vel_change, t_phase):
        return (vel_offset + abs(vel_change) * np.sin(2 * np.pi *
                            ((t - UTC('2010-01-01').toordinal()) / t_period -
                             abs(t_phase))))

    def correlation_at_fit(t, sim, stretch_vec, func, args):
        #sim_mat.shape = 3(tws), 1683(days), 201
        #sim = sim_mat[i,:,:]
        #stretch_vec 201
        ds = np.abs(stretch_vec[1] - stretch_vec[0]) / 2
        stretching = func(t, *args)
        index = np.nonzero(np.abs(np.transpose(np.tile(stretching, (len(stretch_vec), 1))) -
                                  (np.tile(stretch_vec, (len(t), 1)))) < ds)
        return np.mean(sim[index])

    if FIT == 'sinus':
        fit_at = (t_Toco.toordinal(), UTC('2009-09-01').toordinal())
        dates_fit = dates[(dates < fit_at[0]) + (dates > fit_at[1])]
        t_period = 365.242
        func = sinus  #(lambda t, vel_offset, vel_change, t_phase:  # Sinus
        p0 = (0., 0.01, .5)
        print 'vel_offset + vel_change * sin((t - t_phase - 2010-01-01) / 1year * 2 * pi)'
        print 'vel_offset   vel_change   t_phase'
    elif FIT == 'exp':
        dates_fit = dates
        N = len(events)
        print 'Number of events: %d' % N
        evts = [ev.origins[0].time.toordinal() for ev in events]
        func = lambda t, *args: exps(evts, t, *args)
        p0 = (0., 0., 30.) * N
        print 'H(t-t_ev) * (vel_perm + vel_temp * exp(-(t-t_ev)/t_dec*ln(10)))'
        print 't_ev   vel_perm   vel_temp   t_dec (to 10%)'
    elif FIT == 'sinus_exp' or FIT == 'sinus_exp_alt':
        dates_fit = dates
        N = len(events)
        print 'Number of events: %d' % N
        t_period = 365.242
        evts = [ev.origins[0].time.toordinal() for ev in events]
        func = lambda t, *args: sinus(t, *(args[:3])) + exps(evts, t, *(args[3:]))
        p0 = (0., 0.01, .5) + (0., 0., 30.) * N
        print 'vel_offset + vel_change * sin((t - t_phase - 2010-01-01) / 1year * 2 * pi) + '
        print 'H(t-t_ev) * (vel_perm + vel_temp * exp(-(t-t_ev)/t_dec*ln(10)))'
        print 'vel_offset   vel_change   t_phase'
        print 't_ev   vel_perm   vel_temp   t_dec (to 10%)'
    elif FIT == 'sinus_exp_alt2':
        dates_fit = dates
        N = len(events)
        print 'Number of events: %d' % N
        t_period = 365.242
        evts = [ev.origins[0].time.toordinal() for ev in events]
        func = lambda t, *args: sinus(t, *(args[:3])) + exps2(evts, t, *(args[3:]))
        func2 = lambda t, *args: sinus(t, *(args[:3])) + exps(evts, t, *(args[3:]))
        p0 = (0., 0.01, .5) + (0., 30.) * N
        print 'vel_offset + vel_change * sin((t - t_phase - 2010-01-01) / 1year * 2 * pi) + '
        print 'H(t-t_ev) * vel_temp * exp(-(t-t_ev)/t_dec*ln(10))'
        print 'vel_offset   vel_change   t_phase'
        print 't_ev   vel_EQ   t_dec (to 10%)'


    for i in range(corr.shape[1]):
        if FIT == 'sinus':
            fit_index = np.where((dates < fit_at[0]) + (dates > fit_at[1]))[0]
            stretch_fit = stretch[fit_index, i]
            corr_fit = corr[fit_index, i]
        else:
            not_masked = np.logical_not(np.ma.getmaskarray(corr[:, i]))
            dates_fit = dates[not_masked]
            stretch_fit = stretch[:, i][not_masked]
            corr_fit = corr[:, i][not_masked]
            sim = sim_mat[i, not_masked, :]

        from sito import debug
        if FIT.startswith('sinus_exp_alt'):
            func(dates_fit, *p0)
            func_min = lambda args:-correlation_at_fit(dates_fit, sim, stretch_vec, func, args)
            # vel_offset, vel_change, t_phase
            # vel_perm vel_temp t_dec
            #p0 = np.array((0., 0.01, .5) + (0., 0., 30.) * N)
            if FIT == 'sinus_exp_alt':
                p0 = np.array((0.0, 0.25 / 100, 150. / t_period,
                               0, 0.6 / 100, 700))
            else:
                p0 = np.array((-0.5 / 100, 0.25 / 100, 150. / t_period,
                               1. / 100, 700))
#            p0 = np.array((-0.7 / 100, 0.6 / 100, 150. / t_period,
#                           0, 0.6 / 100, 700))
#            bnds = ((-0.1 / 100, 0.1 / 100), (0, 1. / 100), (100 / t_period, 200 / t_period),
#                    (0.0, 0.1 / 100), (0.0, 1. / 100), (300, 1200))
            print func_min(p0)
            results = fmin(func_min, p0)
            popt = results
            print func_min(popt), popt
            popt = list(popt) + [-func_min(popt)]
        else:
            sigma = 1. / corr_fit if USE_SIGMA else None
            func(dates_fit, *p0)  # bind lambda, for some reason I dont understand fully this line is needed
    #        import pylab
    #        pylab.plot(dates_fit, stretch_fit)
    #        pylab.plot(dates_fit, func(dates_fit, *p0))
    #        from IPython import embed
    #        embed()

            popt, pcov = curve_fit(func, dates_fit, stretch_fit, p0, sigma=sigma)  #, maxfev=100000)
        if FIT == 'sinus_exp_alt2':
            popt = popt[:3] + [0.] + popt[3:]  # func -> func2
        else:
            func2 = func

        try:
            errors = np.diag(pcov)
        except:
            errors = -1 * np.ones(len(popt))
        if FIT == 'sinus' or 'sinus' in FIT:
            a, b, c = popt[:3]
            ae, be, ce = errors[:3]
            c %= 1
            t1 = ('(%.2f+-%.2f)%%   (%.2f+-%.2f)%%   (%.1f+-%.1f)days' %
                  (100 * a, 100 * ae, 100 * b, 100 * be, c * t_period, ce * t_period))

            t2 = 'max: %.2f -> day %.1f' % ((c + 0.25) % 1, ((c + 0.25) % 1) * t_period)
            t3 = 'min: %.2f -> day %.1f' % ((c + 0.75) % 1, ((c + 0.75) % 1) * t_period)
            print '   '.join((t1, t2, t3))
#        if FIT == 'sinux_exp_alt2':
#            add = 3 * (FIT != 'exp')
#            texts = ['%s   (%.2f+-%.2f)%%   (%.1f+-%.1f)days' %
#                     (str(ev.origins[0].time)[:16],
#                      100 * popt[2 * ii + add], 100 * errors[2 * ii + add],
#                      popt[2 * ii + 1 + add], errors[2 * ii + 1 + add])
#                     for ii, ev in enumerate(events)]
#            print '\n'.join(texts)
#            print
        if FIT == 'exp' or 'exp' in FIT:
            add = 3 * (FIT != 'exp')
            texts = ['%s   (%.2f+-%.2f)%%   (%.2f+-%.2f)%%   (%.1f+-%.1f)days' %
                     (str(ev.origins[0].time)[:16],
                      100 * popt[3 * ii + add], 100 * errors[3 * ii + add],
                      100 * popt[3 * ii + 1 + add], 100 * errors[3 * ii + 1 + add],
                      popt[3 * ii + 2 + add], errors[3 * ii + 2 + add])
                     for ii, ev in enumerate(events)]
            print '\n'.join(texts)
            print

        ax1.plot(dates, func2(dates, *popt), color='k')
        stretch[:, i] = stretch[:, i] - func2(dates, *popt)
        ax2.plot(dates, stretch[:, i])
        #print sim_mat.shape
        for j in range(sim_mat.shape[1]):
            index = -int(round(0.5 * func2(dates[j:j + 1], *popt)[0] / max(stretch_vec) * len(stretch_vec)))
            sim_mat[i, j, :] = np.roll(sim_mat[i, j, :], index)
            if index >= 0:
                sim_mat[i, j, :index] = 0
            else:
                sim_mat[i, j, index:] = 0
        if FIT == 'sinus':
            np.savez(output_fit % (comb, tw_start[i]), dates=dates, sinus=func2(dates, *popt))
        elif FIT == 'exp':
            np.savez(output_fit % (comb, tw_start[i]), dates=dates, exp=func2(dates, *popt))
        else:
            np.savez(output_fit % (comb, tw_start[i]), dates=dates, sinus_exp=func2(dates, *popt))
        popts[(comb, tw_start[i])] = popt


    corr = corr.filled(0)
    stretch = stretch.filled(0)
    if output:
        np.savez(output % comb, tw_start=tw_start, tw_width=tw_width, dates=dates, corr=corr, stretch=stretch, sim_mat=sim_mat, stretch_vec=stretch_vec)  #corr=corr, stretch=stretch)

    fig.autofmt_xdate()

    show_freq = False
    if show_freq:
        pxx, freqs = mpl.mlab.psd(stretch[:, 0], NFFT=256 * 8, Fs=1, scale_by_freq=True)
        N = 256 * 2 * 2 * 2

        plt.figure()
        plt.plot(freqs, pxx)
        fig2 = plt.figure()
        ax = fig2.add_subplot(111)
        ax.plot(np.fft.fftfreq(N), np.abs(np.fft.fft(stretch[:, 0], N)))
        ax.plot(np.fft.fftfreq(N), np.abs(np.fft.fft(stretch[:, 1], N)))
        ax.plot(np.fft.fftfreq(N), np.abs(np.fft.fft(stretch[:, 2], N)))

import csv
f = open(path + FIT + "_opt.csv", "w")
f.write('# stationcomp, tw, vel_offset, vel_change, t_phase, vel_perm, vel_temp, t_dec(to 10perc), mean_corr_at_fit (t_events = %s)\n' % evts[0])
w = csv.writer(f)
for key, val in sorted(popts.items()):
    w.writerow(list(key) + list(val))
plt.show()



####### 5-10s   10-15s    15-20s
####### Results Sinus
#vel_offset + vel_change * sin((t - t_phase - 2010-01-01) / 1year * 2 * pi)
#vel_offset   vel_change   t_phase
# fixed year 365.242days and weighted with correlation as 1/sigma
#(-0.07+-0.00)%   (0.26+-0.00)%   (144.4+-0.0)days   max: 0.65 -> day 235.7   min: 0.15 -> day 53.1
#(-0.07+-0.00)%   (0.18+-0.00)%   (157.5+-0.0)days   max: 0.68 -> day 248.8   min: 0.18 -> day 66.2
#(-0.07+-0.00)%   (0.19+-0.00)%   (159.3+-0.0)days   max: 0.69 -> day 250.6   min: 0.19 -> day 68.0

####### Results Events
## Tocopilla:
#H(t-t_ev) * (vel_perm + vel_temp * exp(-(t-t_ev)/t_dec*ln(10)))
#t_ev   vel_perm   vel_temp   t_dec (to 10%)
#2007-11-14T15:40   (-0.07+-0.00)%   (0.62+-0.00)%   (1003.1+-2958.6)days
#2007-11-14T15:40   (-0.06+-0.00)%   (0.68+-0.00)%   (920.8+-951.6)days
#2007-11-14T15:40   (-0.05+-0.00)%   (0.64+-0.00)%   (879.5+-1053.3)days

####### Results combined fit:
#vel_offset + vel_change * sin((t - t_phase - 2010-01-01) / 1year * 2 * pi) +
#H(t-t_ev) * (vel_perm + vel_temp * exp(-(t-t_ev)/t_dec*ln(10)))
#vel_offset   vel_change   t_phase
#t_ev   vel_perm   vel_temp   t_dec (to 10%)
#(-0.06+-0.00)%   (0.27+-0.00)%   (142.0+-0.0)days   max: 0.64 -> day 233.3   min: 0.14 -> day 50.7
#2007-11-14T15:40   (-0.08+-0.00)%   (0.63+-0.00)%   (982.3+-2720.5)days
#
#(-0.08+-0.00)%   (0.18+-0.00)%   (153.5+-0.0)days   max: 0.67 -> day 244.8   min: 0.17 -> day 62.2
#2007-11-14T15:40   (-0.05+-0.00)%   (0.68+-0.00)%   (908.8+-896.2)days
#
#(-0.08+-0.00)%   (0.19+-0.00)%   (154.2+-0.0)days   max: 0.67 -> day 245.5   min: 0.17 -> day 62.9
#2007-11-14T15:40   (-0.04+-0.00)%   (0.64+-0.00)%   (867.2+-978.2)days
####### Results combined fit without vel_perm:
#vel_offset + vel_change * sin((t - t_phase - 2010-01-01) / 1year * 2 * pi) +
#H(t-t_ev) * vel_temp * exp(-(t-t_ev)/t_dec*ln(10))
#vel_offset   vel_change   t_phase
#
#(-0.14+--100.00)%   (0.28+--100.00)%   (142.2+--365.2)days   max: 0.64 -> day 233.5   min: 0.14 -> day 50.9
#2007-11-14T15:40   (0.00+--100.00)%   (0.55+--100.00)%   (1258.9+--1.0)days
#
#(-0.10+--100.00)%   (0.19+--100.00)%   (152.6+--365.2)days   max: 0.67 -> day 243.9   min: 0.17 -> day 61.3
#2007-11-14T15:40   (0.00+--100.00)%   (0.68+--100.00)%   (770.6+--1.0)days
#
#(-0.10+--100.00)%   (0.19+--100.00)%   (153.6+--365.2)days   max: 0.67 -> day 244.9   min: 0.17 -> day 62.3
#2007-11-14T15:40   (0.00+--100.00)%   (0.65+--100.00)%   (753.0+--1.0)days


#PB05Z
#Number of events: 1
#vel_offset + vel_change * sin((t - t_phase - 2010-01-01) / 1year * 2 * pi) +
#H(t-t_ev) * vel_temp * exp(-(t-t_ev)/t_dec*ln(10))
#vel_offset   vel_change   t_phase
#t_ev   vel_EQ   t_dec (to 10%)
#-0.70170924569
#Optimization terminated successfully.
#         Current function value: -0.721323
#         Iterations: 124
#         Function evaluations: 299
#-0.721322606178 [ -2.28806222e-03   9.38790467e-05   5.17833065e-01   6.29060090e-03
#   1.89128349e+03]
#(-0.23+--100.00)%   (0.01+--100.00)%   (189.1+--365.2)days   max: 0.77 -> day 280.4   min: 0.27 -> day 97.8
#2007-11-14T15:40   (0.00+--100.00)%   (0.63+--100.00)%   (1891.3+--1.0)days
#
#-0.780986465777
#Optimization terminated successfully.
#         Current function value: -0.831725
#         Iterations: 131
#         Function evaluations: 297
#-0.83172532447 [ -1.26184560e-03   5.19063200e-04   4.29030485e-01   3.88431634e-03
#   1.88795755e+03]
#(-0.13+--100.00)%   (0.05+--100.00)%   (156.7+--365.2)days   max: 0.68 -> day 248.0   min: 0.18 -> day 65.4
#2007-11-14T15:40   (0.00+--100.00)%   (0.39+--100.00)%   (1888.0+--1.0)days
#
#-0.645425145451
#Optimization terminated successfully.
#         Current function value: -0.732057
#         Iterations: 113
#         Function evaluations: 279
#-0.732056918799 [ -5.82936473e-04  -1.55551138e-04   6.13380055e-01   1.76243979e-03
#   1.53307741e+03]
#(-0.06+--100.00)%   (-0.02+--100.00)%   (224.0+--365.2)days   max: 0.86 -> day 315.3   min: 0.36 -> day 132.7
#2007-11-14T15:40   (0.00+--100.00)%   (0.18+--100.00)%   (1533.1+--1.0)days




