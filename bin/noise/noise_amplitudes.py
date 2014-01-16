#!/usr/bin/env python
# by TR
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.util.main import daygen
import obspy.signal
import numpy as np
from matplotlib.dates import date2num
import pylab as plt

def analyze():
    #stations = 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 PB15 PB16 HMBCX MNMCX PATCX PSGCX LVC TAIQ'
    stations = 'PB01 PB02 PATCX'
    stations = 'PATCX'
    t1 = UTC('2009-01-01')
    t2 = UTC('2010-01-01')
    data = IPOC()
    for station in stations.split():
        hours = [[] for i in range(24)]
        times = []
        levels = []
        for t_day in daygen(t1, t2):
            try:
                stream = data.getRawStreamFromClient(t_day, t_day + 24 * 3600, station, component='Z')
            except ValueError:
                continue
            for tr in stream:
                tr.stats.filter = ''
            stream.demean()
            stream.detrend()
            stream.filter2(4, 6)
            stream.downsample2(5)
            stream.merge()
            tr = stream[0]
            startt = tr.stats.starttime
            endt = tr.stats.endtime
            if endt - startt < 12 * 3600:
                continue
            tr.data = obspy.signal.cpxtrace.envelope(tr.data)[1][:len(tr.data)]
            for hour in range(24):
                tr2 = tr.slice(t_day + hour * 3600, t_day + (hour + 1) * 3600)
                if tr2.stats.endtime - tr2.stats.starttime < 1800:
                    continue
                num_stds = 60  # =^ every minute
                len_parts = len(tr2.data) // 60  # =^ 1min
                len_stds = len_parts // 6  # =^ 10s
                stds = np.array([np.std(tr2.data[i:i + len_stds]) for i in np.arange(num_stds) * len_parts])
                stds = stds[stds != 0.]
                num_stds = len(stds)
                if num_stds < 50:
                    continue
                stds = np.sort(stds)[num_stds // 5:-num_stds // 5]
                stds = stds[stds < np.min(stds) * 2.]
                val = np.mean(stds)
                levels.append(val)
                times.append(date2num(t_day + (0.5 + hour) * 3600))
                hours[hour].append(val)
        errors = np.array([np.std(hours[i], ddof=1) / len(hours[i]) ** 0.5 for i in range(24)])
        hours = np.array([np.mean(hours[i]) for i in range(24)])
        times = np.array(times)
        levels = np.array(levels)
        np.savez('/home/richter/Results/IPOC/xcorr/noise_apmlitudes_%s_4-6Hz.npz' % station,
                 hours=hours, errors=errors, times=times, levels=levels)

def plot():
    station = 'PATCX'
    npz = np.load('/home/richter/Results/IPOC/xcorr/noise_apmlitudes_%s_4-6Hz_corrected.npz' % station)
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    ax1.plot_date(npz['times'], npz['levels'])
    ax2.plot(np.arange(24), npz['hours'])
    ax3.plot(np.arange(24), npz['errors'])
    plt.show()
def correct():
    station = 'PATCX'
    npz = np.load('/home/richter/Results/IPOC/xcorr/noise_apmlitudes_%s_4-6Hz.npz' % station)
    levels = npz['levels']
    ind = levels < 20  # PB01: 6, PB02: 3.5, PATCX: 20
    times = npz['times'][ind]
    levels = levels[ind]
    hours = np.empty((24,))
    errors = np.empty((24,))
    for hour in range(24):
        ind = abs((times - (hour + 0.5) / 24) % 1) < 0.001
        errors[hour] = np.std(levels[ind], ddof=1) / len(levels[ind]) ** 0.5
        hours[hour] = np.mean(levels[ind])
    np.savez('/home/richter/Results/IPOC/xcorr/noise_apmlitudes_%s_4-6Hz_corrected.npz' % station,
            hours=hours, errors=errors,
            times=times, levels=levels)
#analyze()
correct()
plot()
