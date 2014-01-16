#!/usr/bin/env python
# by TR
import numpy as np
import matplotlib.pyplot as plt
from IPython import embed
from obspy.core import UTCDateTime as UTC
from matplotlib.dates import date2num

data_path = '/home/richter/Data/climate/'

def get_dict(data):
    data = dict(data)
    for key, item in data.iteritems():
        if key in ('temp', 'dewp', 'slp', 'stp', 'max', 'min'):
            data[key] = np.ma.masked_equal(item, 9999.9)
        elif key in ('visib', 'wdsp', 'mxspd', 'sndp', 'gust'):
            data[key] = np.ma.masked_equal(item, 999.9)
        elif key == 'prcp':
            data[key] = np.ma.masked_equal(item, 99.99)

        if key in ('temp', 'dewp', 'max', 'min'):  #F -> C
            data[key] = (data[key] - 32.) * 5. / 9.
        elif key in ('wdsp', 'mxspd', 'gust'):  #knots -> km/h
            data[key] *= 1.852
        elif key in ('prcp', 'sndp'):  #inch -> mm
            data[key] *= 25.4
        elif key == 'visb':  # mile -> km
            data[key] *= 1.609344
    return data

def get_dict2(data):
    data = dict(data)
    for key, item in data.iteritems():
        if key == 'temp':
            data[key] = np.ma.masked_where((item < 5) + (item > 500), item)
    return data

def mean_day(dates, temp):
    date0 = dates[0]
    dates_mean = []
    temp_mean = []
    while date0 < dates[-1]:
        select = (date0 <= dates) * (date0 + 1 > dates)
        if len(dates[select]) > 5:
            dates_mean.append(np.ma.mean(dates[select]))
            temp_mean.append(np.ma.mean(temp[select]))
        date0 += 1
    return np.array(dates_mean), np.array(temp_mean)

t_period = 365.242
def sinus(t, vel_offset, vel_change, t_phase):
    return (vel_offset + abs(vel_change) * np.sin(2 * np.pi *
                        ((t - UTC('2010-01-01').toordinal()) / t_period -
                         abs(t_phase))))




fname = data_path + 'climate_hour_IQ'
data = get_dict2(np.load(fname + '.npz'))
data2 = np.load(fname + '_mean.npz')

dates = data['date']
temp = data['temp']
preselect = (UTC('2007-01-01').toordinal() <= dates) * (UTC('2012-01-01').toordinal() >= dates)
dates = dates[preselect]
temp = temp[preselect]

fig = plt.figure()
fig.suptitle('Iquique' if fname[-2:] == 'IQ' else 'Antofagasta')
ax1 = fig.add_subplot(211)
ax1.plot_date(dates, temp, 'k')

dates_mean, temp_mean = mean_day(dates, temp)
ax1.plot_date(dates_mean, temp_mean, 'b')


ax1.set_ylabel('temp in C')
fig.autofmt_xdate()

ax2 = fig.add_subplot(212)
ax2.plot_date(data2['date'], data2['temp'], 'k')
ax2.set_ylabel('temp in C')
for label in ax2.get_xticklabels():
    label.set_rotation(30)
    label.set_ha('right')
fig.savefig(fname + '.pdf')

if True:
    print 'Iquique data'
    from scipy.optimize import curve_fit
    p0 = (18., 3.5, 0.5)
    popt, pcov = curve_fit(sinus, dates_mean, temp_mean, p0)
    print popt[0], abs(popt[1]), abs(popt[2]) * t_period, (abs(popt[2]) - 0.5) * t_period
    # 18.0  3.5  314.3  131.6
    ax1.plot(dates_mean, sinus(dates_mean, *popt), 'r')
if False:
    npzfile = np.load('/home/richter/Data/climate/2006-2012_PB01_WKI.npz')
    ax1.plot(date2num(npzfile['dates']), npzfile['data'] / 100., 'y')


fname = data_path + 'climate_CHO'
data_cho = np.load(fname + '.npz')
data2_cho = np.load(fname + '_145_mean.npz')
data3_cho = np.load(fname + '_sens_145_mean.npz')
ax1.plot_date(data_cho['date'], data_cho['temp'], 'y', zorder= -15)
ax1.plot_date(data_cho['date'], data_cho['temp_sens'], 'b')
dates_mean_cho, temp_mean_cho = mean_day(data_cho['date'], data_cho['temp'])
ax1.plot_date(dates_mean_cho, temp_mean_cho, 'orange')

if True:
    print 'CHO data air'
    from scipy.optimize import curve_fit
    p0 = (14., 4., 0.5)
    popt, pcov = curve_fit(sinus, dates_mean_cho, temp_mean_cho, p0)
    print popt[0], abs(popt[1]), abs(popt[2]) * t_period, (abs(popt[2]) - 0.5) * t_period
    # 18.0  3.5  314.3  131.6
    ax1.plot(dates_mean_cho, sinus(dates_mean_cho, *popt), 'r')
    print 'CHO data at sensor'
    popt, pcov = curve_fit(sinus, data_cho['date'], data_cho['temp_sens'], p0)
    print popt[0], abs(popt[1]), abs(popt[2]) * t_period, (abs(popt[2]) - 0.5) * t_period
    # 18.0  3.5  314.3  131.6
    ax1.plot(dates_mean_cho, sinus(dates_mean_cho, *popt), 'r')

# yearly     Iquique, temp_air  temp_sens     temp_0cm    25cm    50cm   75cm
# tempfluc K    7.0       8.0    8.8
# day phase     131.7     128.8     121.4

# daily:
# tempfluc     3.67         9.08      0.112          38      2.5     0.1    0.05
# phase hour max 19:00     17:35     15:10
#                                   (-2.5h) ->ca. 50cm depth


ax3 = ax2.twinx()
ax3.plot_date(data2_cho['date'] - data2_cho['date'][0] + data2['date'][0], data2_cho['temp'], 'b')
ax4 = ax2.twinx()
ax4.spines["right"].set_position(("axes", 1.05))
ax4.plot_date(data3_cho['date'] - data3_cho['date'][0] + data2['date'][0], data3_cho['temp'], 'r')




#data = get_dict(np.load(data_path + 'climate_CAL.npz'))

#plt.suptitle('Calama')
#ax1 = plt.subplot(411)
#plt.plot_date(data['date'], data['prcp'], 'b')
#plt.ylabel('prcp in mm')
#plt.subplot(412, sharex=ax1)
#plt.plot_date(data['date'], data['min'], 'b')
#plt.plot_date(data['date'], data['temp'], 'k')
#plt.plot_date(data['date'], data['max'], 'r')
#plt.ylabel('temp in C')
#plt.subplot(413, sharex=ax1)
#plt.plot_date(data['date'], data['stp'], 'b')
#plt.ylabel('pressure in mbar')
#plt.subplot(414, sharex=ax1)
#plt.plot_date(data['date'], data['wdsp'], 'b')
#plt.plot_date(data['date'], data['mxspd'], 'r')
#plt.ylabel('wdsp in km/h')

#data = get_dict(np.load('data_path + climate_AF.npz'))

#plt.figure()
#plt.suptitle('Antofagasto')
#ax1 = plt.subplot(411)
#plt.plot_date(data['date'], data['prcp'], 'b')
#plt.ylabel('prcp in mm')
#plt.subplot(412, sharex=ax1)
#plt.plot_date(data['date'], data['min'], 'b')
#plt.plot_date(data['date'], data['temp'], 'k')
#plt.plot_date(data['date'], data['max'], 'r')
#plt.ylabel('temp in C')
#plt.subplot(413, sharex=ax1)
#plt.plot_date(data['date'], data['slp'], 'b')
#plt.ylabel('pressure in mbar')
#plt.subplot(414, sharex=ax1)
#plt.title('wind speed')
#plt.plot_date(data['date'], data['mxspd'], 'b')
#plt.plot_date(data['date'], data['wdsp'], 'r')
#plt.ylabel('wdsp in km/h')

data = get_dict(np.load(data_path + 'climate_IQ.npz'))
plt.figure()
plt.suptitle('Iquique')
ax1 = plt.subplot(411)
plt.plot_date(data['date'], data['prcp'], 'b')
plt.ylabel('prcp in mm')
plt.subplot(412, sharex=ax1)
plt.plot_date(data['date'], data['min'], 'b')
plt.plot_date(data['date'], data['temp'], 'k')
plt.plot_date(data['date'], data['max'], 'r')
plt.ylabel('temp in C')
plt.subplot(413, sharex=ax1)
plt.plot_date(data['date'], data['slp'], 'b')
plt.ylabel('pressure in mbar')
plt.subplot(414, sharex=ax1)
plt.title('wind speed')
plt.plot_date(data['date'], data['mxspd'], 'b')
plt.plot_date(data['date'], data['wdsp'], 'r')
plt.ylabel('wdsp in km/h')

plt.show()

