#!/usr/bin/env python
# by TR
import numpy as np
from sito.util import isfloat, isint
import re
from obspy.core import UTCDateTime as UTC
from matplotlib.dates import date2num
from collections import defaultdict

data_path = '/home/richter/Data/climate/'

def extract_cdo_file():
    file = data_path + 'CDO_CAL.txt'  #AF, IQ

    regex = """
    ^
    (?P<stn>\d+)\s+
    (?P<wban>\d+)\s+
    (?P<date>\d+)\s+
    (?P<temp>[-\d.]+)\s+
    (?P<temp_count>\d+)\s+
    (?P<dewp>[-\d.]+)\s+
    (?P<dewp_count>\d+)\s+
    (?P<slp>[\d.]+)\s+
    (?P<slp_count>\d+)\s+
    (?P<stp>[\d.]+)\s+
    (?P<stp_count>\d+)\s+
    (?P<visib>[\d.]+)\s+
    (?P<visib_count>\d+)\s+
    (?P<wdsp>[\d.]+)\s+
    (?P<wdsp_count>\d+)\s+
    (?P<mxspd>[\d.]+)\s+
    (?P<gust>[\d.]+)\s+
    (?P<max>[-\d.]+)\*?\s+
    (?P<min>[-\d.]+)\*?\s+
    (?P<prcp>[\d.]+)
    (?P<prcp_flag>.)\s+
    (?P<sndp>[\d.]+)\s+
    (?P<frshtt>\d+)
    """

    with open(file, 'r') as f:
        filedata = f.read()
    matches = re.finditer(regex, filedata, re.VERBOSE + re.MULTILINE)
    data_list = [i.groupdict() for i in matches]
    data = {}
    for i, dp in enumerate(data_list):  # convert numbers to float and int types
        for key, item in dp.iteritems():
            if key == 'date':
                dp[key] = UTC(item).toordinal()
            elif item is not None and key != 'frshtt':
                if isint(item):
                    dp[key] = int(item)
                elif isfloat(item):
                    dp[key] = float(item)
        if dp['stn'] == 854420:
            dp['station'] = 'Antofagasta'
        elif dp['stn'] == 854320:
            dp['station'] = 'Calama'
        elif dp['stn'] == 854180:
            dp['station'] = 'Iquique'
        else:
            dp['station'] = 'other'
        if i == 0:
            for key, item in dp.iteritems():
                data[key] = []
        for key, item in dp.iteritems():
            #ipshell()
            data[key].append(item)

    np.savez(data_path + 'climate_CAL', **data)  #AF, IQ


def extract_hour_file():
    fname = data_path + 'HOUR_ALL.txt'
    data = defaultdict(lambda : defaultdict(list))
    with open(fname, 'r') as f:
        for line in f:
            if line[:2] in ('AN', 'IQ'):
                ls = line.split()
                st = ls[0]
                time = date2num(UTC(ls[3] + ls[4]))
                data[st]['date'].append(time)
                data[st]['temp'].append(float(ls[20]))
    np.savez(data_path + 'climate_hour_IQ', **data['IQUIQUE'])
    np.savez(data_path + 'climate_hour_AF', **data['ANTOFAGASTA'])

def get_dict2(data):
    data = dict(data)
    for key, item in data.iteritems():
        if key == 'temp':
            data[key] = np.ma.masked_where((item < 5) + (item > 500), item)
    return data

def mean_hour_file():
    fname = data_path + 'climate_hour_IQ'
    #fname = data_path + 'climate_CHO'
    sens = False
    data = get_dict2(np.load(fname + '.npz'))
    dates = data['date']
    temp = data['temp' + '_sens' * sens]

    preselect = (UTC('2007-01-01').toordinal() <= dates) * (UTC('2012-01-01').toordinal() >= dates)
    dates = dates[preselect]
    temp = temp[preselect]
    date0 = int(dates[0] + 1)
    temps = []
    for i in range(int(dates[-1] - date0)):
        select = (date0 + i <= dates) * (dates <= date0 + i + 1)
        if i == 0 or len(plotdates) != 25:
            plotdates = dates[select]
        temps.append(temp[select])
        if i % 1000 == 0:
            print i
    print len(temps)
    print [len(i) for i in temps]

    temps = [t for t in temps if len(t) == 25]
    print len(temps)
    temps = np.ma.mean(temps, axis=0)
    np.savez(fname + '_sens' * sens + '_mean', date=plotdates, temp=np.array(temps))

def mean_hour_file2():
    #fname = data_path + 'climate_hour_IQ'
    fname = data_path + 'climate_CHO'
    sens = False
    data = get_dict2(np.load(fname + '.npz'))
    dates = data['date']
    temp = data['temp' + '_sens' * sens]

    preselect = (UTC('2007-01-01').toordinal() <= dates) * (UTC('2012-01-01').toordinal() >= dates)
    dates = dates[preselect]
    temp = temp[preselect]
    date0 = int(dates[0] + 1)
    temps = []
    dates2 = []
    for i in range(int(dates[-1] - date0)):
        select = (date0 + i <= dates) * (dates <= date0 + i + 1)
        if i == 0 or len(plotdates) != 1441:
            plotdates = dates[select]
        if len(temp[select]) >= 144:
            temps.append(temp[select])
            dates2.append(dates[select])
        if i % 1000 == 0:
            print i
    print len(temps)
    print [len(i) for i in temps]
    for i in range(len(temps)):
        #print i
        #print plotdates
        #print dates2[i] - dates2[i][0] + plotdates[0]
        #print temps[i]
        temps[i] = np.interp(plotdates, dates2[i] - dates2[i][0] + plotdates[0], temps[i])
    temps = [t for t in temps if len(t) == 1441]
    print len(temps)
    temps_mean = np.ma.mean(temps, axis=0)
    temps2 = np.array(temps) - np.tile(np.ma.mean(temps, axis=1)[:, np.newaxis], len(temps[0]))
    #temps_mean = np.ma.mean(temps2, axis=0) + np.ma.mean(temps2)
    #from IPython import embed
    #embed()
    temps_std = np.ma.std(temps2, axis=0)
    np.savez(fname + '_sens' * sens + '_hui_mean', date=plotdates, temp=np.array(temps_mean), temp_std=np.array(temps_std))


#extract_cdo_file()
#extract_hour_file()
mean_hour_file2()



