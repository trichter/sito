#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito.data import IPOC, getCor, getFilter, getStack
from obspy.core import UTCDateTime as UTC
import logging
from sito.noisexcorr import get_correlations
from sito.util import checkDir
import numpy as np
from progressbar import ProgressBar
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

#method = 'filter4-6_water_env2_1bit'
#method = 'filter4-6_water_env2_1bit_fft'
#method = 'filter0.01-1_water_env2_whitening_1bit_fft'

#method = '/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft'
method = 'FINAL_filter0.01-0.5_1bit_whitening'
method = 'FINAL_filter0.01-0.5_1bit_auto'
method = 'FINAL_filter1-3_1bit_auto'
method = 'FINAL_filter4-6_1bit_auto'
method = 'FINAL_filter4-6_1bit_auto_3C'
method = 'FINAL_filter3-5'
method = 'PAT_filter9-11'
method = 'zerotest_nozero'

data = IPOC(xcorr_append='/' + method, use_local_LVC=False)
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit/stretch_t/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit/stretch2/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft/stretch_Toco/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_whitening/stretch_Toco/swcoda/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_auto/stretch_Toco/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter1-3_1bit_auto/stretch/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto_3C/stretch3_10s/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter3-5/stretch/'
path = '/home/richter/Results/IPOC/xcorr/PAT_filter9-11/stretch3/'
path = '/home/richter/Results/IPOC/xcorr/zerotest_nozero/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch2/'
#path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch_Toco_PATCX/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch_Toco/'

#stations = 'PB01 PB02 PB03 PB04 PB05'
stations = ('PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 '
            'PB14 PB15 PB16 HMBCX MNMCX PATCX PSGCX LVC')
stations = 'PATCX'
#stations = 'PATCX'
#correlations = (('PB01Z', 'PB02Z'),)
#correlations = (('PB03Z', 'PB04Z'), ('PB03Z', 'PB03Z'), ('PB04Z', 'PB04Z'))
#correlations = get_correlations(stations, 'Z', only_cross=True)
correlations = get_correlations(stations, 'Z', only_auto=True)
#correlations = (('PB01Z', 'PB02Z'), ('PB03Z', 'PB04Z'), ('PB03Z', 'PB05Z'), ('PB04Z', 'PB05Z'))

#stack = ('50days', '5days')
#stack = ('30days', '2days')
stack = ('10days', '2days')
stack = ('hour', None)
#stack = ('10days', 'day')
stack = None

filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25), (0.05, 0.1))
filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25))
filters = (None,)
reftr = None
#reftr = 'alternative'
border_time = 25  #7
fit_file = path + '/sinus_exp_alt_%s_%ds.npz'


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
tws = (((5, 10, 15), 5),)
#tws = ((range(21), 2),)
#tws = ((range(11), 2),)
#tws = ((range(0, 41, 2), 10),)
tws1 = tws
tws2 = (((-20, -15, -10, 5, 10, 15), 5),)
#tws = (((1, 5, 10, 15), 5),)
#tws = (((20, 30, 40, 50), 20),)
#tws = (((10, 30, 50), 50),)

############ parameters
USE_TWS_AFTER_SW = False

#t1 = UTC('2007-10-01')
#t2 = UTC('2008-02-01')
t1 = UTC('2007-01-01')
t2 = UTC('2012-10-01')
t1 = UTC('2007-10-01')
t2 = UTC('2007-11-30')
#t1 = UTC('2009-05-01')
#t2 = UTC('2009-10-01')
#t2 = UTC('2007-09-05')
#t2 = UTC('2011-04-01')
#t2 = UTC('2009-01-01')
period = 24 * 3600
#period = 301
#period = 1800
add_to_file = '_Tocevent'
add_to_file = '_2007_2008'
add_to_file = ''
if reftr:
    add_to_file += '_altref'
nstr = 201
str_range = 0.02
load_single_file = False
#sl_mean = slice(10, 110, None)
sl_mean = slice(None, None, None)
trim_stream = 200

############# script starts here

for i in range(len(filters)):
    for correlation in ProgressBar()(correlations):
        if correlation[0] == correlation[1]:
            tws = tws1
        else:
            tws = tws2
        try:
            if not load_single_file:
                stream = data.readX(correlation, t1, t2, period=period, filter=filters[i], stack=stack)
            else:
                stream = data.readX(correlation, t1=None, t2=None, period=period, filter=filters[i], stack=stack)
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
        stretches = None
        if reftr:
            stretches = []
            for tw in tws[i][0]:
                if tw <= border_time:
                    npzfile = np.load(fit_file % (getCor(*correlation), tw))
                stretches.append(-npzfile['sinus_exp'])
        result = stream.stretch(reftr=reftr, stretch=stretches, str_range=str_range, nstr=nstr, time_windows=tws[i], sides='right')
        if period == 24 * 3600:
            dates = [(time + 12 * 3600).toordinal() for time in stream.getHI('starttime')]
        else:
            def get_ord(time):
                day = (time + period / 2)
                day = day.toordinal() + 1.*day.hour / 24 + 1.*day.minute / 24 / 60 + 1.*day.second / 24 / 3600
                return day
            dates = [get_ord(time) for time in stream.getHI('starttime')]
        dates = _correct_dates(dates)
        print np.array([dates[j + 1] - dates[j] for j in range(len(dates) - 1)])
        checkDir(path + 'bla')
        np.savez(path + 'data_stretching_%s%s%s%s' % (getCor(*correlation), getFilter(filters[i]), getStack(None, period, stack), add_to_file), tw_start=tws[i][0], tw_width=tws[i][1], dates=dates, **result)  #corr=corr, stretch=stretch)

