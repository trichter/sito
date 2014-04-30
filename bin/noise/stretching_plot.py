#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito.data import getCor, getFilter, getStack
from obspy.core import UTCDateTime as UTC
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
from sito.noisexcorr import get_correlations
import numpy as np
import matplotlib.dates as mdates
from matplotlib.transforms import offset_copy
from obspy.core.event import readEvents
import itertools
from matplotlib.dates import num2date, date2num
import sys
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon, PathPatch
from matplotlib.path import Path
logging.basicConfig(level=logging.DEBUG)

def _insert_zeros(data, dates, min_delta=None, return_dates=False):
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
        dates = np.concatenate((dates[:index + 1],
                 [dates[index] + min_delta * (j + 1) for j in range(num)],
                 dates[index + 1:]))
        counter += num
    return (data, dates) if return_dates else data


############ some parameters
mpl.rcParams.update({'font.size': 16})
#maj_loc = mdates.MonthLocator()
#min_loc = mdates.DayLocator((5, 10, 15, 20, 25))
#maj_fmt = mdates.DateFormatter('%b %Y')
maj_loc = mdates.YearLocator()
min_loc = mdates.MonthLocator()
maj_fmt = mdates.DateFormatter('%Y')
#maj_loc=None

#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit/stretch_t/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit/stretch2/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft/stretch_Toco/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_whitening/stretch_Toco/swcoda/'
#path = '/home/richter/Results/IPOC/xcorr/FINAL_filter0.01-0.5_1bit_auto/stretch_Toco/'
#path = '/home/richter/Results/IPOC/xcorr/FINAL_filter1-3_1bit_auto/stretch/'
path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch/'
path = '/home/richter/Results/IPOC/xcorr/zerotest_zero/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto_hour/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto_3C/stretch3_10s/'

#path = '/home/richter/Results/IPOC/xcorr/PAT_filter9-11/stretch2/'

#path = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto/stretch_Toco_PATCX/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit_fft/stretch_Toco/'

#stations = 'PB01 PB02 PB03 PB04 PB05'
stations = ('PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 '
            'PB14 PB15 PB16 HMBCX MNMCX PATCX PSGCX LVC')
#stations = 'HMBCX PATCX PSGCX'
stations = 'PATCX'
#stations = 'PB04'

#correlations = get_correlations(stations, 'Z', only_cross=True)
#correlations = get_correlations(stations, 'Z', only_auto=True)
correlations = get_correlations(stations, 'Z', only_auto=True)

#stack = ('50days', '5days')
#stack = ('30days', '2days')
#stack = ('10days', '2days')

#stack = ('10days', 'day')
stack = None
#stack = ('hour', None)

filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25), (0.05, 0.1))
filters = (None, (0.5, 1.), (0.25, 0.5), (0.1, 0.25))
filters = (None,)

period = 24 * 3600

METHOD = 'simple'
#METHOD = 'events'
ALL = True  # all windows for different times in one figure
SIZE = 'A4'  # A4, A3, A41, A42
show = False
alt_ref = False
assert METHOD in 'simple correct woline sinus exp sinus_exp sinus_exp_alt events temp press humi'.split()
for_ = 'DIS'  #'DIS'  #F, C, JGR1, JGR2, JGR3, JGRlegend, False, DIS

# simple
# correct: correct line of highest correlation to be in a specific window
# woline without max corr line
############ fine tuning parameters
plot_all_events = False

alternative_npzfile = (path + '/data_stretching_PATCXZ_without_seasons.npz'
                       if METHOD == 'exp' else
                       path + '/data_stretching_PATCXZ_without_Tocopilla.npz'
                       if METHOD == 'events' else None)
alternative_max_cor = METHOD in ('correct', 'sinus', 'sinus_exp', 'sinus_exp_alt', 'temp', 'press', 'humi')
plot_allinone = ALL
plot_allinone_everysecondticklabel = SIZE[:2] == 'A4'
plot_fit = METHOD in ('sinus', 'exp', 'sinus_exp', 'sinus_exp_alt', 'temp', 'press', 'humi')
plot_fit_file = path + '/' + METHOD + '_%s_%ds.npz'
plot_arrows = for_ in ('C', 'JGR1')
lw_big = 3
lw_small = 1

if for_ in ('JGR1', 'JGR2', 'JGR3', 'JGRlegend', 'DIS'):
    mpl.rcParams.update({'font.size': 7, 'lines.linewidth':1})
    lw_big = 1.5
    lw_small = 0.5
plot_climate_file = None
if METHOD == 'temp':
    plot_climate_file = '/home/richter/Data/climate/2006-2012_PB01_WKI.npz'
if  for_ in ('C', 'JGR1'):
    plot_climate_file = '/home/richter/Data/climate/climate_IQ.npz'
    #plot_climate_file = '/home/richter/Data/climate/climate_CHO.npz'
if METHOD == 'press':
    plot_climate_file = '/home/richter/Data/climate/2006-2012_PB01_WDI.npz'
if METHOD == 'humi':
    plot_climate_file = '/home/richter/Data/climate/2006-2012_PB01_WII.npz'
if for_:
    plot_fit_file2 = path + '/exp_%s_%ds.npz'


plot_events = (METHOD == 'events' or for_ == 'C') and for_ != 'JGR2'
event_file = '/home/richter/Data/events/2007-2012_events_in_chile_geofon.xml'
rotate_event_labels = 0
plot_simple = METHOD == 'simple'

############# parameters#####
str_range = 0.012  ###########
thres_corr = 0.35  ##########
thres_cor2 = 0.4  ###########
thres_cor2 = 0.  ############
thres_dt = (-0.02, 0.02)  ###
thres_dt = (-0.006, 0.006)  #

add_to_file = ''
#add_to_file = '_2007_2008'  # for xcorr
MAX_CORR_LINE = METHOD not in ('woline') and for_ not in ('JGR3',)
temp_color = 'blue'
max_xcorr_color = 'cyan' #'blue'
fit_color = 'orange'
pga_color = '#006000'
line = ',-'
#line = '.-'
plot_complex = True
t_Toco = UTC('2007-11-14 15:14:00')
plot_Toco = METHOD != 'events' and for_ != 'C'
cmap = 'hot_r'
plot_pga = for_ == 'JGR2'

### last hacks
#plot_arrows = False;  #temp_color = None;  #plot_fit = False;  plot_Toco = False;
###

ext = ('.png' if METHOD in ('simple', 'correct') and not MAX_CORR_LINE else
       '_line.png' if METHOD in ('simple', 'correct') else
       '_woline.png' if METHOD == 'woline' else
       '_%s.png' % METHOD if METHOD in ('sinus', 'exp', 'sinus_exp', 'sinus_exp_alt') else
       '_events_%s.png' % SIZE if METHOD == 'events' and SIZE in ('A41', 'A42') else
       '_events.png' if METHOD == 'events' else
       '_temperature_PB01.png' if METHOD == 'temp' else
       '_pressure_PB01.png' if METHOD == 'press' else
       '_humidity_PB01.png' if METHOD == 'humi' else
       None)
#ext = '_woseasons.png'
if ALL:
    ext = '_all' + ext
if for_:
    ext = '_' + for_ + ext
if SIZE[:2] == 'A4' and ALL:
    figsize = (8.267, 11.693 * 1.2)[::-1]
elif SIZE[:2] == 'A4':
    figsize = (8.267 * 0.5, 11.693 * 1.2)[::-1]
else:
    figsize = (20, 10)
if for_ == 'JGR1':
    fw = 170 / 25.4
    fh = 0.3 * fw
    figsize = (fw, fh)
elif for_ in ('JGR2', 'JGR3', 'DIS'):
    fw = 170 / 25.4
    fh = 0.2 * fw
    figsize = (fw, fh)
elif for_ == 'JGRlegend':
    fw = 170 / 25.4
    fh = 0.08 * fw
    figsize = (fw, fh)
if alt_ref:
    add_to_file += '_altref'

if for_ == 'JGRlegend':
    fig = plt.figure(figsize=figsize)
    rect = [0.06, 0.9, 0.93, 0.01]
    #rect = [0.06, 0.8, 0.80, 0.01]
    ax = fig.add_axes(rect, frameon=False)
    labels = ('', '', 'best correlation', 'temperature',
              'fit at correlation matrix', 'exponential part of fit',
              'peak ground acceleration', '', '')
    def myhandler(legend, orig_handle, fontsize, handlebox):
        w, h, x, y = handlebox.width, handlebox.height, handlebox.xdescent, handlebox.ydescent
        xm, ym = x + w / 2, y + h / 2
        hs = (0.1, 0.1, 0.1, 0.3, 3, 1, 0.1, 0.1, 2, 0.1, 0.1, 0.1, 0.1)
        xy = [((xm - 3 + i, ym - 3), (xm - 3 + i, ym + 2 * h - 3)) for i, h in enumerate(hs)]
        codes = [(Path.MOVETO, Path.LINETO) for i in range(len(hs))]
        patch = PathPatch(Path(sum(xy, ()), codes=sum(codes, ())), lw=1, fc='none', color=pga_color)
        handlebox.add_artist(patch)
    mpl.legend.Legend.update_default_handler_map({None: myhandler})
    l = Line2D((0, 0), (1, 1), color='w')
    handles = (l, l,
               Line2D((0, 0), (1, 1), color=max_xcorr_color),
               Line2D((0, 0), (1, 1), color=temp_color),
               Line2D((0, 0), (1, 1), color=fit_color),
               Line2D((0, 0), (1, 1), color=fit_color, dashes=(3, 3)),
               None, l, l)
    ax.legend(handles, labels, bbox_to_anchor=(0., 1., 1., .1), loc=2,
              ncol=5, mode="expand", borderaxespad=0., scatterpoints=1)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.annotate('b)', (0., 0), (1, 5), annotation_clip=False, textcoords='offset points', xycoords='figure fraction', va='bottom', ha='left', size='large', weight='roman')
    plt.savefig(path + 'stretching_%s%s' % ('PATCX', ext.replace('.png', '.pdf')))
    sys.exit()

############# script starts here
print 'Use method %s' % METHOD
if alternative_npzfile:
    print 'Alternative npz file: ON'

fff = 5
for filter_, correlation in itertools.product(filters, correlations):
#        labels = ['%s_%d' % (getCor(*correlation), year) for year in range(2007, 2012)]
#        labels.append(st)
    XCORR = correlation[0] != correlation[1]
    if alternative_npzfile:
        npzfile = np.load(alternative_npzfile)
        l = '%s' % (getCor(*correlation),)
    else:
        l = '%s%s%s%s' % (getCor(*correlation), getFilter(filter_), getStack(None, period, stack), add_to_file)
        try:
            #print path + 'data_stretching_%s.npz' % l
            npzfile = np.load(path + 'data_stretching_%s.npz' % l)
        except Exception as ex:
            print ex
            print 'Error for %s-%s' % correlation
            continue
    tw_start, tw_width, corr, stretch, dates, sim_mat, stretch_vec = 7 * [None]
    vars().update(npzfile)
    sum_to_one_day = False
    if sum_to_one_day:
        sim_mat_new = np.empty((sim_mat.shape[0], 24, sim_mat.shape[2]))
        for i in range(24):
            sim_mat_new[:, i, :] = np.average(sim_mat[:, i::24, :], axis=1)
        sim_mat = sim_mat_new
    if SIZE in ('A41', 'A42'):
        split_date = (dates[-1] - dates[0]) // 2 + dates[0]
        if SIZE == 'A41':
            slice_ = np.where(dates <= split_date)[0]
        else:
            slice_ = np.where(dates >= split_date)[0]
        dates = dates[slice_]
        corr = corr[slice_]
        stretch = stretch[slice_]
        sim_mat = sim_mat[:, slice_, :]
    corr = np.nan_to_num(corr)
    sim_mat = np.nan_to_num(sim_mat)
    if alternative_npzfile:
        stretch = np.ma.masked_where(corr == 0, stretch)
        corr = np.ma.masked_where(corr == 0, corr)
    if alternative_max_cor:
        # insert zeros where there is no data (->white in colorplot)
        sim_mat, dates = _insert_zeros(sim_mat, dates, return_dates=True)
        corr2 = np.zeros(sim_mat.shape[:-1]).transpose()
        stretch2 = np.zeros(sim_mat.shape[:-1]).transpose()
        # calculate corr and stretch again from simmat but change possible range for stretching
        for (ii, ctw) in enumerate(tw_start):
            tmp = sim_mat[ii, :, :]
            # Hard coded thres_dt for some stations
            if 'PATCX' in correlation[0] and 'filter4-6' in path:
                thres_dt = ((-0.01, 0.01) if ctw == 5 else
                            (-0.008, 0.01) if ctw == 10 else
                            (-0.005, 0.006, (t_Toco, UTC('2008-11-01'), 0.0, 0.01))
                              if ctw == 15 else thres_dt)
            # Threshold for line with best correlation
            if thres_dt:
                sel = np.where((stretch_vec > thres_dt[0]) * (stretch_vec < thres_dt[1]))[0]
                tmp_corr_vect = tmp[:, sel].max(axis=1)
                corr2[:, ii] = tmp_corr_vect
                stretch2[:, ii] = stretch_vec[sel][tmp[:, sel].argmax(axis=1)]
                # other threshold in specific time windows (used for PATCX)
                if len(thres_dt) > 2:
                    for min_date, max_date, thres_dt_min, thres_dt_max in thres_dt[2:]:
                        sel = np.where((stretch_vec > thres_dt_min) * (stretch_vec < thres_dt_max))[0]
                        tmp_corr_vect = tmp[:, sel].max(axis=1)
                        #corr2[:, ii] = tmp_corr_vect
                        sel_dates = np.where((dates > min_date.toordinal()) * (dates < max_date.toordinal()))[0]
                        corr2[sel_dates, ii] = tmp_corr_vect[sel_dates]
                        stretch2[sel_dates, ii] = stretch_vec[sel][tmp[sel_dates, :][:, sel].argmax(axis=1)]

            else:
                tmp_corr_vect = tmp.max(axis=1)
                corr2[:, ii] = tmp_corr_vect
                stretch2[:, ii] = stretch_vec[tmp.argmax(axis=1)]
        # mask corrs where there was no data (-> line is not plotted)
        corr = np.ma.masked_where(corr2 == 0, corr2)
        stretch = np.ma.masked_where(corr2 == 0, stretch2)
        # save corrected values in npz file
        if 'PATCX' in correlation[0] and 'filter4-6' in path and not alternative_npzfile:
            np.savez(path + 'data_corrected_stretching_%s%s%s%s' % (getCor(*correlation), getFilter(filter_), getStack(None, period, stack), add_to_file), tw_start=tw_start, tw_width=tw_width, dates=dates, corr=corr2, stretch=stretch2, sim_mat=sim_mat, stretch_vec=stretch_vec)

    elif not alternative_npzfile:
        sim_mat = _insert_zeros(sim_mat, dates)

    if plot_simple and not plot_allinone:
        dt = stretch
        #corr = np.ma.masked_less(corr, thres_corr)
        corr[corr < thres_corr] = thres_corr
        fig = plt.figure(figsize=figsize, subplotpars=
                         mpl.figure.SubplotParams(hspace=0.07))
        plt.suptitle('stretching_%s' % l)
        ax1 = fig.add_subplot(211)
        #ax1.set_position([0.1, 0.5, 0.8, 0.4])
        for j in range(corr.shape[1]):
            #print len(dates), len(corr[:, j])
            ax1.plot_date(dates, corr[:, j], line, label='tw=%d/%ds' %
                          (tw_start[j], tw_start[j] + tw_width))
            ax1.set_ylabel('correlation')
            ax1.legend(loc='lower right')
        ax2 = fig.add_subplot(212, sharex=ax1)
        #ax2.set_position([0.1, 0.1, 0.8, 0.4])
        for j in range(corr.shape[1]):
            ax2.plot_date(dates, dt[:, j] * 100, line)
#                ax2.plot_date(npzfile2['dates'], npzfile2['dt'][j, :], line)
            ax2.set_ylabel('rel. time shift in %')
            ax2.set_ylim((-str_range * 100, str_range * 100))
            #plt.legend()
        fig.autofmt_xdate()
        ticks = ax1.get_xticks()
        if len(ticks) > 20:
            ticks = ticks[::len(ticks) // 20]
            ax1.set_xticks(ticks)
        if not show:
            plt.savefig(path + 'stretching_%s%s' % (l, ext))
            plt.close()

    if plot_events:
        cat = readEvents(event_file)
        if not plot_all_events:
            cat = cat.filter('time < %s' % (t_Toco + 3600)) + cat.filter('time > %s' % (t_Toco + 30 * 24 * 3600), 'time < 2011-10-01')
        fig3, basemap = cat.plot(projection='local', ipoc_stretching=True,
                                 handle=True, colormap='brg', min_mag_size=10,
                                 max_mag_size=20)
        fig3.delaxes(fig3.axes[1])
        import sito.stations
        ipoc = sito.stations.IPOCStations()
        ipoc.pick('PATCX')
        ipoc.plot(basemap, color='k', markersize=8)
        r = rotate_event_labels
        if r == 90:
            axt = fig3.axes[0]
            for child in axt.get_children():
                if isinstance(child, mpl.text.Text):
                    child.set_rotation(r)
                    t = child.get_text()
                    if not 'W' in t and not 'S' in t:
                        child.set_ha('right')
                        child.set_va('bottom')
        plt.title('')
    if not plot_complex:
        continue
    for i in range(sim_mat.shape[0]):
        if plot_allinone:
            l2 = l
        else:
            l2 = l + ('_%s-%s' % (tw_start[i], tw_start[i] + tw_width))
        l3 = '%ss-%ss' % (tw_start[i], tw_start[i] + tw_width)
        if plot_allinone:
            NN = sim_mat.shape[0]
            if plot_climate_file:
                rect = [0.06, 0.16 + 0.74 * (1 - 1. / NN - 1.*i / NN), 0.8, 0.74 / NN]
            else:
                rect = [0.06, 0.16 + 0.74 * (1 - 1. / NN - 1.*i / NN), 0.86, 0.74 / NN]
            if i == 0:
                axes = []
                fig = plt.figure(figsize=figsize)
                ax1 = fig.add_axes(rect)
                axes.append(ax1)
            else:
                ax1 = fig.add_axes(rect, sharex=axes[0], sharey=axes[0])
                axes.append(ax1)
        else:
            fig = plt.figure(figsize=figsize)
            if for_:
                if 'JGR' not in for_ and 'DIS' not in for_:
                    plt.suptitle('stretching factor of ambient noise autocorrelation')
            else:
                plt.suptitle('stretching_%s' % l2)
            if for_ == 'JGR1':
                rect = [0.06, 0.1, 0.80, 0.87]
            elif for_ == 'JGR2':
                rect = [0.06, 0.11, 0.93, 0.86]
            elif for_ == 'JGR3':
                rect = [0.06, 0.11, 0.85, 0.86]
            elif 'DIS' in for_:
                rect = [0.06, 0.095, 0.85, 0.86]
            elif plot_climate_file:
                rect = [0.06, 0.16, 0.8, 0.74]
            else:
                rect = [0.06, 0.16, 0.86, 0.74]
            ax1 = fig.add_axes(rect)
        #ax1.set_position([0.1, 0.5, 0.8, 0.4])
        if not for_:
            ax1.annotate(l3, (0, 1), (10, -10), xycoords='axes fraction',
                         textcoords='offset points', va='top')
        dt = 1.*(dates[-1] - dates[0]) / (len(dates) + 1) / 2
        ds = 1.*(max(stretch_vec) - min(stretch_vec)) / (len(stretch_vec) + 1) / 2
        extent = [dates[0] - dt, dates[-1] + dt, min(stretch_vec) - ds, max(stretch_vec) + ds]
        vmin = thres_cor2 if thres_cor2 is not None else 0.5 * (np.median(sim_mat[i, 0, :]) + np.median(sim_mat[i, -1, :]))
        vmin = max(0, vmin)

        image = ax1.imshow(np.transpose(sim_mat[i, :, :]), cmap=cmap, interpolation='nearest',
                   origin='lower', aspect='auto', vmax=1, vmin=vmin, extent=extent)

                    #vmax=self.vmax, vmin=self.vmin, norm=self.norm)
        if plot_Toco:
            ax1.axvline(t_Toco.toordinal(), color='gray')
        if MAX_CORR_LINE:
            #if thres_dt:
            #    stretch = np.ma.masked_outside(stretch, *thres_dt)
            if for_ == 'JGR3':
                ax1.scatter(dates, stretch[:, i], s=1, marker='o', color=max_xcorr_color, edgecolor='none')
            else:
                #ax1.plot(dates, stretch[:, i], color=max_xcorr_color, lw=lw_small)
                where = np.where(np.ediff1d(dates, to_begin=[0]) >= 4)[0]
                dates_line = np.split(dates, where)
                stretch_line = np.split(stretch[:, i], where)
                for iii in range(len(dates_line)):
                    ax1.plot(dates_line[iii], stretch_line[iii], color=max_xcorr_color, lw=lw_small)

        if plot_fit:
            if METHOD not in ('temp', 'press', 'humi'):
                npzfile = np.load(plot_fit_file % (getCor(*correlation), tw_start[i]))
                ax1.plot(npzfile['dates'], npzfile[METHOD.strip('_alt')], color=fit_color)
            if for_ == 'JGR1':
                npzfile = np.load(plot_fit_file2 % (getCor(*correlation), tw_start[i]))
                ax1.plot(npzfile['dates'], npzfile['exp'], color=fit_color, ls='--')
        if plot_climate_file:
            npzfile = np.load(plot_climate_file % tw_start[i])
            new_ax = ax1.twinx()
            if 'dates' in npzfile:
                dates_ = [date.toordinal() for date in npzfile['dates']]
                temp_data = npzfile['data'] / 100.
            else:
                dates_ = npzfile['date']
                #temp_data = npzfile['temp_sens']
                temp_data = (npzfile['temp'] - 32.) * 5. / 9.
            if METHOD == 'temp' or for_ and temp_color:
                new_ax.plot(dates_, temp_data, color=temp_color, lw=lw_small)
                new_ax.set_ylabel(u'temperature (°C)', color=temp_color)
                ax1.add_artist(new_ax.lines.pop(0))
            elif METHOD == 'press':
                new_ax.plot(dates_, npzfile['data'] / 10000., color=temp_color)
                new_ax.set_ylabel(u'pressure (bar)', color=temp_color)
            elif temp_color:
                new_ax.plot(dates_, npzfile['data'], color=temp_color)
                new_ax.set_ylabel(u'rel. humidity (%)', color=temp_color)
        if plot_pga:
            new_ax = ax1.figure.add_axes(ax1.get_position(True), sharex=ax1, frameon=False)
            new_ax.xaxis.set_visible(False)
            #new_ax.yaxis.set_visible(False)
            npzdata = np.load('/home/richter/Results/IPOC/maxima_PATCX.npz')
            dates_vel = npzdata['dates_vel']
            dates_ac = npzdata['dates_ac']
            vel = npzdata['vel']
            vel = vel / 629145000.0 * 100
            ac = npzdata['ac']  # correct for sensitivity # cm/s
            ac[ac > 400000] = 0
            ac = ac / 427566.942  # correct for sensitivity # m/s^2
            vel_dict = dict(zip([str(d.date()) for d in num2date(dates_vel[vel > 0.1])], np.round(vel[vel > 0.1], 3)))
            ac_dict = dict(zip([str(d.date()) for d in num2date(dates_ac[ac > 0.005])], np.round(ac[ac > 0.005], 3)))
            bot = 0.04
            new_ax.bar(dates_ac - 1, ac, 2, bot, color=pga_color, lw=0)
            new_ax.set_yticks([bot, bot + 0.04, bot + 0.08])
            new_ax.set_yticklabels(['0', '0.04', '0.08'])
            new_ax.set_ylim([0, 0.5])
            #new_ax.axhline(y=0.02, color='gray', zorder= -1)
            #new_ax.set_ylabel('PGA ($m/s^2$)')
        if for_ == 'JGR3':
            # plot rectangle
            t1 = date2num(UTC('2008-06-01'))
            t2 = date2num(UTC('2009-06-01'))
            y1 = -0.001
            y2 = 0.005
            ax1.plot((t1, t2, t2, t1, t1), (y1, y1, y2, y2, y1), 'k', alpha=0.5, lw=2)

        if not plot_allinone or i == sim_mat.shape[0] - 1:
            #ax1.xaxis_date()
            ax1.set_xlim(extent[:2])
            ax1.set_ylim(extent[2:])
            if maj_loc:
                ax1.xaxis.set_major_locator(maj_loc)
                ax1.xaxis.set_minor_locator(min_loc)
                ax1.xaxis.set_major_formatter(maj_fmt)
            ticks = ax1.get_xticks()
            if len(ticks) > 20:
                ticks = ticks[::len(ticks) // 20]
                ax1.set_xticks(ticks)
            rect2 = [0.94, 0.16, 0.01, 0.74]
            if for_ in ('JGR1', 'JGR2', 'JGR3', 'DIS'):
                rect2 = [0.93, 0.1, 0.01, 0.87]
                ax1.yaxis.labelpad = 3
            if for_ == 'JGR3' or for_ == 'DIS':
                rect2 = [0.93, 0.1, 0.01, 0.8]
            if for_ in ('JGR1', 'JGR3', 'DIS'):
                cb = fig.colorbar(image, cax=fig.add_axes(rect2),
                                  extend='min')
                cb.set_label('correlation')
            fig.autofmt_xdate()
            if plot_allinone or True:
#                if plot_fit:
#                    l2 += '_' + METHOD
                #plt.suptitle('stretching PATCX')
                if not plot_allinone:
                    axes = [ax1, ax1, ax1]

                if METHOD in ('temp', 'press', 'humi'):
                    axes[1].set_ylabel('stretching factor (%)', color='b')
                if for_ and for_ != 'JGR2':
                    if NN == 3:
                        axes[1].set_ylabel('velocity decrease (%)')
                    elif NN == 2:
                        axes[1].annotate('velocity decrease (%)', (0, 0.5), (11, 0), textcoords='offset points', xycoords='figure fraction', va='center', ha='center', rotation=90)
                elif for_ == 'JGR2':
                    axes[1].set_ylabel(r'PGA($\rm{m/s^2}$)    vel decr (%)  ')
                else:
                    axes[1].set_ylabel('stretching factor (%)')
                if METHOD == 'exp':
                    axes[0].set_ylim((-0.002, 0.008))
                    axes[0].set_yticks((-0.002, 0, 0.002, 0.004, 0.006, 0.008))
                    axes[0].set_yticklabels(('', '0.0', '0.2', '0.4', '0.6'))
                elif METHOD in ('sinus', 'events'):
                    axes[0].set_ylim((-0.005, 0.005))
                    axes[0].set_yticks((-0.004, -0.002, 0.000, 0.002, 0.004))
                    axes[0].set_yticklabels(('-0.4', '-0.2', '0.0', '0.2', '0.4'))
                    if for_ == 'JGR2':
                        axes[0].set_ylim([-0.007, 0.005])
                        axes[0].set_yticks((-0.002, 0.000, 0.002, 0.004))
                        axes[0].set_yticklabels(('-0.2', '0.0', '0.2', '0.4'))
                elif METHOD in ('sinus_exp',):
                    axes[0].set_ylim((-0.005, 0.007))
                    axes[0].set_yticks((-0.004, -0.002, 0.000, 0.002, 0.004, 0.006))
                    axes[0].set_yticklabels(('-0.4', '-0.2', '0.0', '0.2', '0.4', '0.6'))
                else:
#                    axes[0].set_yticks([-0.005, 0, 0.005])
#                    axes[0].set_yticklabels(['-0.5', '0.0', '0.5'])
#                    axes[0].set_yticks([-0.0075, -0.0025, 0.0025, 0.0075], minor=True)
                     pass
                if for_ == 'JGR3' or for_ == 'DIS':
                    axes[0].set_ylim((-0.021, 0.021))
                    if for_ == 'DIS':
                        axes[0].set_ylim((-0.02, 0.02))
                    axes[0].set_yticks((-0.02, -0.01, 0.0, 0.01, 0.02))
                    axes[0].set_yticklabels(('-2', '-1', '0', '1', '2'))
                elif for_ == 'DIS':
                    axes[0].set_ylim((-0.01, 0.01))
                    axes[0].set_yticks((-0.01, -0.005, 0.0, 0.005, 0.01))
                    axes[0].set_yticklabels(('-1.0', '-0.5', '0.0', '0.5', '1.0'))
                if plot_allinone:
                    for ax in axes[:2]:
                        for label in ax.get_xticklabels():
                            label.set_visible(False)
                    if for_ == 'DIS':
                        axes[0].set_ylim((-0.02, 0.02))
                    axes[0].set_yticks((-0.02, -0.01, 0.0, 0.01, 0.02))
                    axes[0].set_yticklabels(('', '-1', '0', '1', ''))
                if plot_events:
                    cmap2 = plt.get_cmap('brg')
                    times = [ev.origins[0].time for ev in cat]
                    lats = [ev.origins[0].latitude for ev in cat]
                    mags = [ev.magnitudes[0].mag for ev in cat]
                    ax_ev = fig.add_axes(rect, frame_on=False, sharex=axes[0])
                    ax_ev.set_ylim((0, 1))
                    #ax.patch.set_alpha(0.5)
                    #ac.patch.set_facecolor('None')
                    for ii in range(len(times)):
                        time = times[ii]
                        lat = lats[ii]
                        mag = mags[ii]
                        even = ii % 2 == 0

                        if for_ in ('C', 'JGR2'):
                            axes[0].axvline(time.toordinal(), color='grey', zorder=1)
                        else:
                            color = cmap2((time - min(times)) / (max(times) - min(times)))
                            #if not (t_Toco + 3600 < time < t_Toco + 41 * 24 * 3600):
                            ax_ev.axvline(time.toordinal(),
                                        color='k' if lat > -21 else 'grey')
                            s = 50 + 100 * (mag - min(mags)) / (max(mags) - min(mags))
                            ax_ev.scatter(time.toordinal(), 0.05, marker='o',
                                        s=s, c=color, zorder=10)
                            ax_ev.annotate('%.1f' % mag, (time.toordinal(), 0.05), (0, 8 if even else -8),
                                        textcoords='offset points', color=color,
                                        va='bottom' if even else 'top',
                                        ha='center')
                    ax_ev.xaxis.set_visible(False)
                    ax_ev.yaxis.set_visible(False)
                if for_:
                    maj_loc = mdates.YearLocator()
                    min_loc = mdates.MonthLocator()
                    maj_fmt = mdates.DateFormatter('')  #%Y
                elif period >= 24 * 3600:
                    maj_loc = mdates.MonthLocator(range(1, 13, 2))
                    #min_loc = mdates.DayLocator((5, 10, 15, 20, 25))
                    min_loc = mdates.MonthLocator()
                    maj_fmt = mdates.DateFormatter('%b')  #%Y
                else:
                    maj_loc = mdates.DayLocator()
                    min_loc = mdates.DayLocator()
                    maj_fmt = mdates.DateFormatter('%b %D')  #%Y
                ax = axes[-1]
                ax.xaxis.set_major_locator(maj_loc)
                ax.xaxis.set_minor_locator(min_loc)
                ax.xaxis.set_major_formatter(maj_fmt)
                ax.set_xlim(extent[:2])
                fig.canvas.draw()
                labels = [item.get_text() for item in ax.get_xticklabels()]
                for ii, label in enumerate(ax.get_xticklabels()):
                    if plot_allinone_everysecondticklabel and ii % 2 == 1:
                        label.set_visible(False)
                    trans = label.get_transform()
                    offtrans = offset_copy(trans, fig, x=7, y=0, units='points')
                    label.set_transform(offtrans)
                    label.set_x(label.get_position()[0] + 200)
                    label.set_ha('right')
                    label.set_rotation(30)
                from datetime import date
                for year in range(date.fromordinal(int(extent[0])).year,
                                  date.fromordinal(int(extent[1])).year + 1):
                    if year == 2006 and correlation[0] not in 'PB01 PB02 PB03 PB04 PB05 MNMCX PATCX PSGCX LVC HMBCX'.split():
                        continue
                    temp2 = UTC('%s-07-01' % year).toordinal()
                    temp3 = ax.get_ylim()[0]
                    ax.annotate(str(year), (temp2, temp3), (0, -35 + 32 * bool(for_)), annotation_clip=False, textcoords='offset points', xycoords='data', va='top', ha='center')
                ax.set_xticklabels(labels)
            else:
                ax1.set_ylabel('stretching factor')
            if for_ == 'JGR1' and plot_Toco:
                ax.annotate('Mw7.7 Tocopilla EQ', (t_Toco.toordinal(), -0.0033), (30, 0),
                            textcoords='offset points', va='center',
                            arrowprops=dict(arrowstyle="->"))
            if plot_arrows:
                t11 = UTC('2007-03-01').toordinal()
                t22 = t_Toco.toordinal()
                t33 = UTC('2010-01-01').toordinal()
                dttt = 10
                dvvv = 0.0003
                c = 'k'
                aps = dict(arrowstyle='<->', lw=lw_big, color=c)
                dict1 = dict(arrowprops=aps, xycoords='data', zorder=20)
                dict2 = dict(xycoords='data', textcoords='offset points', ha='right', va='center', color=c)
                dict3 = dict(xycoords='data', textcoords='offset points', ha='center', va='bottom', color=c)
                ax.annotate('', xy=(t11, -0.0026 - dvvv), xytext=(t11, 0.0010 + dvvv), **dict1)
                ax.annotate('$2\epsilon_\mathrm{P}$', xy=(t11, 0.001), xytext=(0, 3), **dict2)
                ax.annotate('', xy=(t22, 0. - dvvv), xytext=(t22, 0.006 + dvvv), **dict1)
                ax.annotate('$\epsilon_\mathrm{EQ}$', xy=(t22, 0.006), xytext=(-2, 0), **dict2)
                ypos = -0.0042
                ax.annotate('', xy=(t33 - dttt, ypos), xytext=(t33 + dttt + 68, ypos), **dict1)
                ax.annotate('$t_\mathrm{P}$', xy=(t33, ypos), xytext=(0, 2), **dict3)
                ax.annotate('', xy=(t22 - dttt, ypos), xytext=(t22 + 2 * 365 + dttt, ypos), **dict1)
                ax.annotate('$t_\mathrm{EQ}$', xy=(t22 + 640, ypos), xytext=(0, 2), **dict3)
            for ll in ax.get_xticklines():
                ll.set_markersize(5)
                ll.set_markeredgewidth(lw_big)
            if plot_allinone:
                for ax2 in fig.axes[:3]:
                    for ll in ax2.get_xticklines():
                        ll.set_markersize(5)
                        ll.set_markeredgewidth(lw_big)


                #axes[0].annotate('bla', xy=(0., 0.0), xycoords='axes fraction',
                #            xytext=(1., 0.5), arrowprops=aps)
            if for_ == 'JGR1':
                ax.annotate('a)', (0., 1), (1, -1), annotation_clip=False, textcoords='offset points', xycoords='figure fraction', va='top', ha='left', size='large', weight='roman')
            if for_ == 'DIS':
                HZ = '  1-3Hz' if 'filter1-3' in path else '  4-6Hz' if 'filter4-6' in path else ''
                TW = ' %s-%ss' % (tw_start[i], tw_start[i] + tw_width)
                bbox = dict(boxstyle="round", fc="w", alpha=0.5, ec='none')
                if plot_allinone:
                    for i, ax2 in enumerate(fig.axes[:int(NN)]):
                        if XCORR:
                            tw1, tw2 = tw_start[i], tw_start[i] + tw_width
                            TW = '(%ds,%ds)' % (tw1, tw2)
                        else:
                            TW = ('5-10s', '10-15s', '15-20s')[i]
                        ax2.annotate(TW, (0, 1), (5, -2), textcoords='offset points', xycoords='axes fraction', va='top', ha='left', bbox=bbox, fontsize='small')
                    if XCORR:
                        fig.axes[0].annotate('-'.join(correlation), (0, 1), (15, -2), textcoords='offset points', xycoords='figure fraction', va='top', ha='left', clip_on=False)
                    else:
                        fig.axes[0].annotate(correlation[0] + HZ, (0, 1), (15, -2), textcoords='offset points', xycoords='figure fraction', va='top', ha='left', clip_on=False)
                else:
                    ax.annotate(correlation[0] + HZ + TW, (0, 1), (5, -5), textcoords='offset points', xycoords='axes fraction', va='top', ha='left')
            if not show:
                if for_:
                    plt.savefig(path + 'stretching_%s%s' % (l2, ext.replace('.png', '.pdf')), dpi=600)
                plt.savefig(path + 'stretching_%s%s' % (l2, ext), dpi=300)
                plt.close()

if show:
    plt.show()
