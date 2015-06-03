#!/usr/bin/env python
from glob import glob
from operator import itemgetter

from matplotlib.dates import date2num, num2date
import numpy as np
import pylab as plt
from sito import read
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
import yaml
import obspy.arclink
import scipy.stats as spstats
import scipy.optimize as spoptimize
import matplotlib as mpl


#def calc_strong_motion_old():
#    t1 = UTC('2007-01-01')
#    t2 = UTC('2011-01-01')
#    data = IPOC()
#    dates1 = []
#    dates2 = []
#    acs = []
#    vels = []
#    while t1 + 60 < t2 - 60:
#        try:
#            ms = data.getChannelFromClient(t1 - 60, t1 + 24 * 3600 + 60, 'PATCX', channel='HLZ')
#        except:
#            pass
#        else:
#            ms.merge(fill_value=0)
#            tr = ms[0]
#            if tr.stats.endtime - tr.stats.starttime < 1000:
#                t1 += 24 * 3600
#                continue
#            tr.stats.filter = ''
#            tr.detrend()
#            tr.filter('highpass', freq=0.1)
#            tr.trim(tr.stats.starttime + 60, tr.stats.endtime - 60)
#            max_acc = np.max(np.abs(tr.data))
#            dates1.append(t1)
#            acs.append(max_acc)
#        try:
#            ms = data.getChannelFromClient(t1 - 60, t1 + 24 * 3600 + 60, 'PATCX', channel='HHZ')
#        except:
#            pass
#        else:
#            ms.merge(fill_value=0)
#            tr = ms[0]
#            if tr.stats.endtime - tr.stats.starttime < 1000:
#                t1 += 24 * 3600
#                continue
#            tr.stats.filter = ''
#            tr.detrend()
#            tr.filter('highpass', freq=0.1)
#            tr.trim(tr.stats.starttime + 60, tr.stats.endtime - 60)
#            max_vel = np.max(np.abs(tr.data))
#            dates2.append(t1)
#            vels.append(max_vel)
#        t1 += 24 * 3600
#
#    np.savez('/home/richter/Results/IPOC/maxima_PATCX.npz', dates_ac=dates1, dates_vel=dates2, vel=vels, ac=acs)

def calc_strong_motion():
    files1 = '/media/PBO/archive/20?[7890]/CX/PATCX/B[HL]Z.D/*.???'
    files2 = '/media/platte/Data/IPOC/raw_PATCX_2011_B-Z/*/*.???'
    acs = []
    vels = []
    for fname in glob(files1) + glob(files2):
        ms = read(fname)
        for tr in ms:
            if tr.stats.endtime - tr.stats.starttime < 1000:
                ms.remove(tr)
                continue
            tr.trim(tr.stats.starttime + 10, tr.stats.endtime - 10)
            tr.stats.filter = ''
            tr.detrend()
            tr.filter('highpass', freq=0.2)
        ms.merge(fill_value=0)
        if len(ms) == 0:
            continue
        tr = ms[0]
        t = tr.stats.starttime + (tr.stats.endtime - tr.stats.starttime) / 2
        maxi = np.max(np.abs(tr.data))
        if 'BH' in fname:
            vels.append((t, maxi))
        elif 'BL' in fname:
            acs.append((t, maxi))
        else:
            print 'wrong filename: ', fname
    dates_vel, vels = zip(*sorted(vels, key=itemgetter(0)))
    dates_ac, acs = zip(*sorted(acs, key=itemgetter(0)))
    np.savez('/home/richter/Results/IPOC/maxima_PATCX_5s.npz',
             dates_ac=date2num(dates_ac), dates_vel=date2num(dates_vel),
             vel=vels, ac=acs)

def calc_strong_motion_Toco():
    t = UTC('2007-11-14')
    data = IPOC()
    channel = 'BHZ'
    acs = {}
    vels = {}

    for station in 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PATCX HMBCX PSGCX MNMCX'.split():
        for channel in ('BLZ', 'BHZ'):
            try:
                ms = data.getChannelFromClient(t - 60, t + 24 * 3600 + 60,
                                               network='GE' if station == 'LVC' else 'CX',
                                               location='10' if station == 'LVC' else '',
                                               station=station, channel=channel)
            except Exception as ex:
                print station, channel, 'sucks'
                continue
            for tr in ms:
                if tr.stats.endtime - tr.stats.starttime < 1000:
                    ms.remove(tr)
                    continue
                tr.trim(tr.stats.starttime + 10, tr.stats.endtime - 10)
                tr.stats.filter = ''
                tr.detrend()
                tr.filter('highpass', freq=0.2)
                tr.trim(tr.stats.starttime + 10, tr.stats.endtime - 10)
            ms.merge(fill_value=0)
            if len(ms) == 0:
                continue
            maxi = float(np.max(np.abs(ms[0].data)))
            if 'BH' in channel:
                vels[station] = maxi / 629145000.0 * 100
            elif 'BL' in channel:
                acs[station] = maxi / 427566.942
    with open(GM_TOCO, 'w') as f:
        yaml.dump({'vels': vels, 'acs':acs}, f, default_flow_style=False)

def plot_strong_motion():
    npzdata = np.load('/home/richter/Results/IPOC/maxima_PATCX.npz')
    #npzdata2 = np.load('/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch/data_corrected_stretching_PATCXZ_altref.npz')
    npzdata2 = np.load('/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch/data_stretching_PATCXZ_without_Tocopilla_and_seasons.npz')

    dates_vel = npzdata['dates_vel']
    dates_ac = npzdata['dates_ac']
    vel = npzdata['vel']
    vel = vel / 629145000.0 * 100
    ac = npzdata['ac']  # correct for sensitivity # cm/s
    ac[ac > 400000] = 0
    ac = ac / 427566.942  # correct for sensitivity # m/s^2
    #print 'dates_vel', [str(d.date()) for d in num2date(dates_vel[vel > 0.2])]
    vel_dict = dict(zip([str(d.date()) for d in num2date(dates_vel[vel > 0.1])], np.round(vel[vel > 0.1], 3)))
    ac_dict = dict(zip([str(d.date()) for d in num2date(dates_ac[ac > 0.005])], np.round(ac[ac > 0.005], 3)))
    from pprint import pprint
    print 'vel'
    pprint(vel_dict)
    print
    print 'ac'
    pprint(ac_dict)


    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.plot(dates_vel, vel, 'r')
    ax1.axhline(y=0.2, color='gray', zorder=-1)
    ax1.set_ylabel('vel')
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax2.plot(dates_ac, ac, 'b')
    ax2.axhline(y=0.02, color='gray', zorder=-1)
    ax2.set_ylabel('ac')
    ax3 = fig.add_subplot(313, sharex=ax1)
    stretch = npzdata2['stretch'][:, 1] * 100
    ax3.plot(npzdata2['dates'], stretch, 'b')
    ax3.xaxis_date()
    fig.autofmt_xdate()
    #plt.show()

def calcr(x, y, f):
    #return (1 - np.sum((y - f(x)) ** 2) / np.sum((y - np.mean(y)) ** 2)) ** 0.5
    return (np.sum((np.mean(y) - f(x)) ** 2) / np.sum((y - np.mean(y)) ** 2)) ** 0.5

def plot_veldrop_vs_sm(combine=False):
    twin = False
    dtypes = (UTC,) + 4 * (float,) + (basestring,)
    data = np.genfromtxt(VELDROP_VS_SM % 'events', dtype=dtypes, names=True)
    date, mag, veldrop, gm, vel, comment = zip(*data)
    gm = 2 * np.array(gm) # correct for wrong sensitivity used when calculating gm
    vel = np.array(vel)
    veldrop = np.array(veldrop)

    dtypes = (basestring,) + 3 * (float,) + (basestring,)
    data = np.genfromtxt(VELDROP_VS_SM % 'stations', dtype=dtypes, names=True)
    st, veldrop2, gm2, vel2, comment = zip(*data)
    gm2 = np.array(gm2)  # here data was corrected in the data file
    vel2 = np.array(vel2)
    veldrop2 = np.array(veldrop2)

    if twin:
        axes1 = [0.12, 0.16, 0.42, 0.68]
        axes2 = [0.55, 0.16, 0.42, 0.68]
    else:
        axes1 = [0.12, 0.16, 0.42, 0.80]
        axes2 = [0.55, 0.16, 0.42, 0.80]
    if combine:
        axes3 = [0.08, 0.18, 0.29, 0.75]
        axes1 = [0.385, 0.18, 0.29, 0.75]
        axes2 = [0.69, 0.18, 0.29, 0.75]
    fig1 = plt.figure()
    ax1 = fig1.add_axes(axes1)
    ax2 = fig1.add_axes(axes2, sharey=ax1)
    print gm, veldrop
    print(len(gm), len(veldrop))
    ax1.scatter(gm, veldrop, s=30, marker='x', color='k')
    ax2.scatter(gm2, veldrop2, s=30, marker='x', color='k')  #facecolor='none')

    fig2 = plt.figure()
    if combine:
        ax3 = fig1.add_axes(axes3, sharey=ax1)
        ax4 = fig2.add_axes(axes2)
    elif twin:
        ax3 = ax1.twiny()
        ax4 = ax2.twiny()
    else:
        ax3 = fig2.add_axes(axes1)
        ax4 = fig2.add_axes(axes2, sharey=ax2)


    ax3.scatter(vel, veldrop, s=30, marker='x', color='gray' if twin else 'k')
    ax4.scatter(vel2, veldrop2, s=30, marker='x', color='gray' if twin else 'k')

    slope, intercept, r_value, p_value, std_err = spstats.linregress(gm, veldrop)
    slope2, _ = spoptimize.curve_fit(lambda x, a: a * x, gm, veldrop, p0=(0,))
    print 'r value, **2', r_value, r_value ** 2
    print 'r_value calculated', calcr(gm, veldrop, lambda x: slope * x + intercept)
    print  'p_value', p_value
    print 'standard deviation', std_err
    print
    print 'r_value intercept=0', calcr(gm, veldrop, lambda x: slope2 * x)
    print
    print
    ax1.annotate(r'$r^2\!=%0.2f$' % r_value ** 2, (0.63, 0.55), None, 'axes fraction')
    x = np.array([0, 0.5])
    ax1.plot(x, slope * x + intercept, 'r-', zorder=-1)
    #ax1.plot(x, slope2 * x, 'b-')
    ax1.set_xlim((0, 0.46))
    ax1.set_ylim((0, 1.2))
    ax2.set_xlim((0, 3))
    if not combine:
        ax1.set_ylabel('velocity drop (%)')
    ax1.set_xlabel(r'peak ground acceleration ($\rm m/\rm s^2\!$)', ha='left')
    #ax1.set_xticklabels(['0.0', '0.1', '0.1', '0.15', '0.2'])
    ax1.annotate('station PATCX\ndifferent events', (0, 1), (5, -5), 'axes fraction', 'offset points', va='top')
    ax2.annotate('Tocopilla event\ndifferent stations', (1, 1), (-5, -5), 'axes fraction', 'offset points', va='top', ha='right')
    for l in ax2.get_yticklabels():
        l.set_visible(False)

    slope, intercept, r_value, p_value, std_err = spstats.linregress(vel, veldrop)
    slope2, _ = spoptimize.curve_fit(lambda x, a: a * x, vel, veldrop, p0=(0,))
    print 'r value, **2', r_value, r_value ** 2
    print  'p_value', p_value
    print 'standard deviation', std_err
    print
    print 'r_value intercept=0', calcr(vel, veldrop, lambda x: slope2 * x)

    ax3.annotate(r'$r^2\!=%0.2f$' % r_value ** 2, (0.72, 0.55), None, 'axes fraction')

    x = np.array([0, 2.5])
    ax3.plot(x, slope * x + intercept, 'r-', zorder=-1)
    #ax2.plot(x, slope2 * x, 'b-', lw=1.5)
    ax3.set_xlim((0, 2.5))
    ax3.set_ylim((0, 1.2))
    ax4.set_xlim((0, 2.5))
    ax3.set_ylabel('velocity drop (%)')
    if combine:
        ax3.set_xlabel(r'dynamic strain ($\mathrm{cm}/\mathrm{s}/\overline{v}_{\mathrm{S}}$)')
    else:
        ax3.set_xlabel('ground velocity (cm/s)', ha='left')
    ax3.set_xticklabels(['0.0', '0.5', '1.0', '1.5', '2.0', ''])
    ax3.annotate('station PATCX\ndifferent events', (0, 1), (5, -5), 'axes fraction', 'offset points', va='top')
    ax4.annotate('Tocopilla event\ndifferent stations', (0, 1), (5, -5), 'axes fraction', 'offset points', va='top')

    for l in ax4.get_yticklabels():
        l.set_visible(False)
    if combine:
        for l in ax1.get_yticklabels():
            l.set_visible(False)
        fig1.savefig(SAVEFIG3)
    else:
        fig1.savefig(SAVEFIG)
        fig2.savefig(SAVEFIG2)

    #plt.show()



VELDROP_VS_SM = '/home/richter/gfz/Results/IPOC/veldrop_vs_groundmotion_%s.txt'
GM_TOCO = '/home/richter/Results/IPOC/maxima_Toco.yaml'
SAVEFIG = '/home/richter/gfz/Results/IPOC/veldrop_vs_groundmotion_corrected.pdf'
SAVEFIG2 = '/home/richter/gfz/Results/IPOC/veldrop_vs_groundvel_corrected.pdf'
SAVEFIG3 = '/home/richter/Results/IPOC/veldrop_vs_gm_ds.pdf'

#calc_strong_motion()
#calc_strong_motion_Toco()

#plot_strong_motion()

# JGR
fw = 85 / 25.4
fh = fw / 1.6
fsize = (fw, fh)
fsize2 = (fw * 1., fh)
mpl.rcParams.update({'figure.figsize': fsize2, 'font.size': 7, 'lines.linewidth':0.7})
plot_veldrop_vs_sm(combine=False)

"""
date     mag dist vel%  gm>1Hz vel>1Hz
2007-11-14 T7.8  0.400  0.207  2.29   +long term effect
2007-12-16  7.1  0.101  0.021  0.41
2008-02-04  6.6  0.153  0.035  0.62
2008-02-16  5.5  0.104  0.026  0.28
2008-03-01  5.7  0.204  0.039  0.77
2008-03-24  5.9  0.378  0.055  0.75
2008-09-10  6.0  0.533  0.101  0.87
2009-04-17  6.0  0.197  0.013  0.31
2009-07-15  5.3  0.185  0.011  0.18
2009-11-13  6.4  0.262  0.028  0.46
2010-02-27 M8.7  0.0    0.008  0.51
2010-03-04  6.3  0.185  0.020  0.34
2010-03-15  4.9?  0.183  0.024  0.24
2010-04-05  5.7  0.088  0.016  0.15
2011-06-20  6.5  0.361  0.070  1.35
"""
