#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
"""
script for calculating receiver functions
INPUT: Event list
       data object

functions:
pick_events -> pick waveforms around ponset of events
calculate_rf -> calculates rf
mout -> move out correction

the rest of the functions are depreciated and moved to other scripts
"""

import sys

from IPython.core import ultratb
sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                                color_scheme='LightBG', call_pdb=1)


import logging
import os.path
import calendar
import numpy as np
from obspy.core import UTCDateTime
from sito import data as mod_data, events, imaging, read, util, Stream
import pylab as plt
import matplotlib as mpl

exha = util.exha
log = logging.getLogger('script_receiver')
top = 'black'
bot = 'gray'
top2 = 'red'
bot2 = 'blue'
top3 = 'white'
bot3 = 'white'

show = False
def main():
    global data, plotdir, mark

##### Parkfield filenames
#    pkd = mod_data.Parkfield()
#    pkd.events = '/home/richter/Data/events/events_30-90_mag5.5_Parkfield.txt'
#    pkd.rf_events = pkd.data + '/receiver/events30-90/%s_%s'
#    pkd.rf_results = pkd.results + '/receiver/results/%s_%s'
#    pkd.rf_results_dir = pkd.results + '/receiver/results/'
#    pkd.rf_events = pkd.data + '/receiver/events30-90_mag5.5/%s_%s'
#    pkd.rf_results = pkd.results + '/receiver/results_mag5.5/%s_%s'
#    pkd.rf_results_dir = pkd.results + '/receiver/results_mag5.5/'
#    pkd.rf_plot_dir = pkd.rf_results_dir + 'plots/'

##### IPOC filenames
    ipoc = mod_data.IPOC()
    ipoc.events = '/home/richter/Data/events/2012_03_events_27-93_mag5.5_IPOC.txt'
    ipoc.rf_events = ipoc.data + '/receiver/2012_events_mag5.5/%s_%s'
    ipoc.rf_results = ipoc.results + '/receiver/2012_mag5.5_RT/%s_%s'
    ipoc.rf_results_dir = ipoc.results + '/receiver/2012_mag5.5_RT/'
    ipoc.rf_plot_dir = ipoc.rf_results_dir + 'plots/'

##### EventPicker
#    data=pkd
#    pick_events()
#    data.events = '/home/richter/Data/events/events_90-160_mag5.5_Parkfield.txt'
#    data.rf_events = data.data + '/receiver/events90-160_mag5.5/%s_%s'
#    pick_events(window = (-100,400), phase = 'PP')

    data = ipoc
#    data.stations.pick('LVC')
#    events2 = events.Events.read(data.events)
#    events2.pick(after='2012-01-01')
#    data.events = events2
    #pick_events()

#    data.events = '/home/richter/Data/events/events_90-160_mag5.5_IPOC.txt'
#    data.rf_events = data.data + '/receiver/events90-160_mag5.5/%s_%s'
#    pick_events(window = (-100,400), phase = 'PP')


##### Constrain events
    #from sito import Events
    #data.events = Events.read(data.events)
    #data.events.pick(after='2010-01-01')


##### Calculating RFs
#    logging.basicConfig()
#    data=pkd
#    plotdir = data.rf_plot_dir
#    util.checkDir(plotdir)
#    calculate_rf()
#    mout()
##    produce_event_files('PKD')
#    create_rf_plots(years=range(1996, 2011))
#
##    create_interval_plot('225<=st.azi<=265.7')
##    time_binned('PKD', 1996.75, 2011.25, 0.5,'225<=st.azi<=265.7')
#
#    # IPOC
    logging.basicConfig()
    data = ipoc
    plotdir = data.rf_plot_dir
    util.checkDir(plotdir)
    calculate_rf(rotateLQT=False)
    mout()
    #mout('LVC')
#    produce_event_files('PB02')
    ###create_rf_plots(components='LQT')  # Do not use!
    create_rf_plots(components='ZRT')

#    create_interval_plot('140<=st.azi<=160')
#    time_binned('PB01', 2006.75, 2011.25, 0.5,'140<=st.azi<=160')

###### Creatibng Event MAP
#    eventfile = '/home/richter/Data/events/events_30-90_mag5.5_IPOC.txt'
#    ev = events.Events.read(eventfile)
#    ev.plot(-22.0, -70.0, bigmap=False, color='k', radius=10, draw_countries=False, circles_around=((-58.3, -22, 800, False, 'b'),))
#    plt.show()


###### Binned profile at PPs
#    data = ipoc
#    lon1 = -70.2
#    lon2 = -68.8
#    dif = 0.02
#    bins = np.linspace(lon1, lon2, int((lon2 - lon1) / dif) + 1)
#    profile_pp_binned('allst', bins, header='plon', xlim=(-70.2, -68.8))


###### Creating plots for different regions
#    ## PKD regions
#    data = pkd
#    station = 'PKD'
#    regions = ((42.6, 144.8, 200),
#               (-17.4, -174.8, 500),
#               (-30, -178, 500),
#               (34.9, 140.9, 500),
#               (8.3, -83.8, 500),
#               (51.5, -173., 650),
#               (51.3, 145.5, 860),
#               (-23.4, -70.44, 1000))
#    ## IPOC regions
#    data = ipoc
#    station = 'PB04'
#    regions = ((-58.3, -22.0, 800),)
#               (13.5, -92.1, 1400),
#                (14.1, -91.2, 600))


#    plotdir = data.rf_results_dir + 'plot_regions2/'
#    util.checkDir(plotdir)
#    stream = read(data.rf_results_dir + '%s_20??_good_nomout.QHD' % station)
#    #stream.filter2(None, 0.2)
#    stream.sort('starttime')
#    for region in regions:
#        str_region = stream.select(around=region)
#        #for i in (26, 25, 4): # PB03 M5.5 South Sandwich
#        #for i in (26, 23, 19, 18, 4): # PB01 M5.5 South Sandwich
#        for i in (23, 30, 0): # PB04 M5.5 South Sandwich
#        #for i in ():
#            str_region.pop(3 * i + 2)
#            str_region.pop(3 * i + 1)
#            str_region.pop(3 * i)
#        ratio = len(str_region) // 3 / 100. * 0.6
#        if len(str_region) // 3 < 150:
#            ratio *= 1.5
#        if len(str_region) // 3 < 50:
#            ratio *= 2. / 1.5
#        if len(str_region) // 3 < 20:
#            ratio *= 5. / 2.
#        margin = [1., 0.2, 1.5, 1.5]  #left, rigth, bottom, top
#        plot = str_region.plotRF(-5, 33, show=False, fig=getFig(ratio=ratio, margin=margin, fontsize=16, labelsize=14),
#                                 absolutescale=3, sumscale=5, topcolor=top, botcolor=bot,
#                                plotinfo=('azi dist', 'starttime'), plotlabel=(u'azi (°)/dist (°)', 'year'), #plotinfo_width=0.1,
#                                usehardticks='time', figtitle='events %dkm around %.2f, %.2f (%s) at %s' % (region[2], region[0], region[1], util.region(region[0], region[1]), station), title_xpos=0.02, title_horalign='left')
##        fig = str_region.plotRF(-10, 60, show=False, fig=getFig(ratio=ratio), scale=6, sumscale=20, topcolor=top, botcolor=bot,
##                                plotinfo=('azi dist', 'starttime'), plotinfo_width=0.1, plotlabel=('time (s)', u'azi/dist', 'year'),
##                                usehardticks='time', figtitle = '%dkm around %.2f, %.2f  filtered 5s-30s' % (region[2],region[0], region[1]))
#        plot.fig.savefig(plotdir + 'rf_5-35_%s_%dkm_around_%.1f_%.1f_Q.eps' % (station, region[2], region[0], region[1]), dpi=300)
#        plt.close(plot.fig)


##### Calculating psds

##### Calculating ***PP*** RFs
#    ipoc = mod_data.IPOC()
#    ipoc.rf_events = ipoc.data + '/receiver/events90-160/%s_%s'
#    ipoc.rf_results = ipoc.results + '/receiver/resultsPP/%s_%s'
#    ipoc.rf_results_dir = ipoc.results + '/receiver/resultsPP/'
#    ipoc.rf_plot_dir = ipoc.rf_results_dir + 'plots/'
#
#    data = ipoc
#    plotdir = data.rf_plot_dir
#    util.checkDir(plotdir)
#    #data.stations.pick('PB03')
#    calculate_rf(pp=True)
#    mout()
#    #produce_event_files('PB03')
#    create_rf_plots()
#    #create_interval_plot('270<=st.azi<=300')
#    #time_binned('PB01', 2006.75, 2011.25, 0.5,'270<=st.azi<=300')


def pick_events(window=(-100, 500), filter=(0.033, 2.), phase='P', new_sampling_rate=20):  #@ReservedAssignment
    logfile = os.path.dirname(data.rf_events) + '/log_pick_events_%s.txt'
    util.checkDir(logfile)
    util.setRootLogger(logfile=logfile % '', logdebugfile=logfile % '_debug')
    log.info('***** Pick events: %s' % util.parameters())
    mod_data.eventPicker(data, component='all', phase=phase, window=window,
                     filter=filter, new_sampling_rate=new_sampling_rate)
@exha()
def calculate_rf(year='*', pp=False, rotateLQT=True, deconvolvef=False):
    logfile = data.rf_results_dir + 'a_log%s.txt'
    util.checkDir(logfile)
    util.setRootLogger(logfile=logfile % '', logdebugfile=logfile % '_debug')
    log.info('***** Calculate RF')
    for station in data.stations.keys():
        stream = read(data.rf_events % (station, year + '.QHD'))
        if pp:
            stream.setPhase('PP')
        stream.pspier(60, data.stations)
        log.info('number of events: %s' % (len(stream) // 3))
        #stream.filter2(0.033, 2.)
        stream.trim2(-20, 100, relative='ponset')
        stream.sort(('event.id', 'station', 'component'))
        stream.check()
        stream.setHI('mark', False)
        if rotateLQT:
            stream.rotateZNE2LQT(-5, 15, usetheo=True)
            stream.afarm(signoise=2.0, remove=True)
        else:
            #stream.rotateZNE2LQT(-5, 15, usetheo=True)
            #stream.afarm(signoise=2.0, remove=True)
            #stream.trim2(-20, 100, relative='ponset')
            #stream.rotateLQT2ZNE(usetheo=True)
            stream.rotateNE2RT()
        #stream.afarm(signoise=2.0, remove=True)
        #stream.trim2(-25, 100, relative='ponset')

        # log.info('number of events after first farm: %s' % (len(stream)//3))
        # util.ipshell()
        #stream.receiverf(water=0.005, gauss=5, tshift=20, pad=0,
        #                 window='tukey', start=-10, end=30, where='ponset', lenslope=5)
        if deconvolvef:
            stream.receiverf()
        else:
            stream.receivert()
        #stream.receiverSH(-10, 80, 1)
        #stream.afarm('rf', signoise=2., signoiseQ=1., maxL=1 / 1.5, sigQ=False, broad=True, remove=False)
        #stream.afarm('rf', signoise=False, signoiseQ=False, maxL=False, sigQ=False, broad=False, remove=False)
        #log.info('number of events after second farm: %s' % (len(stream)//3))

        stream.write(data.rf_results % (station, '') + 'nomout', 'Q')
        #stream.writey(data.rf_results % (station, '%s') + '_nomout', 'Q')
        print stream.getReasons()

@exha()
def mout(station='*', year='*'):
    log.info('***** Move Out correction')
    global rf_stream
    logfile = data.rf_results_dir + 'a_mout_log%s.txt'
    util.checkDir(logfile)
    util.setRootLogger(logfile=logfile % '', logdebugfile=logfile % '_debug')
    rf_stream = read(data.rf_results % (station, year) + 'nomout.QHD')
    if False:
        st2 = rf_stream.copy()
        st2.moveout(phase='Ppps')
        st2.trim2(-20, 100)
        st2.write(data.rf_results % '_Ppps', 'Q')
        st3 = rf_stream.copy()
        st3.moveout(phase='Ppss')
        st3.trim2(-20, 100)
        st3.write(data.rf_results % '_Ppss', 'Q')
    rf_stream.moveout(phase='Ps')
    #rf_stream.trim2(-20, 100)
    rf_stream.writex(data.rf_results % ('%s', '') + 'mout', 'Q', years=False)
    #rf_stream.writex(data.rf_results + '_mout', 'Q')
    #rf_stream = rf_stream.select(expr='st.mark==False')
    #rf_stream.writex(data.rf_results + '_good_mout', 'Q')


@exha()
def produce_event_files(station, year='*'):
    """ load data and throwed data to produce event files with used and not used events """
    log.info('***** Produce event files')
    stream = read(data.rf_results % (station, year) + 'all_nomout.QHD')
    events_bad = events.Events.readFromStream(stream.select(expr='st.mark'))
    events_good = events.Events.readFromStream(stream.select(expr='not st.mark'))
    events_bad.write('%sevents_bad_%s.txt' % (data.rf_results_dir, station), header=False, format_=events_bad.format_GMT)
    events_good.write('%sevents_good_%s.txt' % (data.rf_results_dir, station), header=False, format_=events_good.format_GMT)

def read_rf(force=False):
    global rf_stream
    try:
        rf_stream
    except:
        rf_stream = read(data.rf_results_dir + '*_mout.QHD')
        #rf_stream = read(data.rf_results_dir + 'LVC_mout.QHD') #!!!
    else:
        if force:
            rf_stream = read(data.rf_results_dir + '*_mout.QHD')
    return rf_stream

def getFig(num=0, ratio=1.5, margin=None, **kwargs):
    axes = [1. - 0.1 * num] + [0.1] * num
    if margin == None:
        margin = [1.5, 0.1, 1.5, 0.9]  #left, rigth, bottom, top
    fig = imaging.getFigure(axes, width=15., margin=margin, ratio=ratio,
                            fontsize=12, labelsize='small', **kwargs)
    return fig

def create_rf_plots_years(station, years):
    global rf_stream
    for year in years:
        read(data.rf_results % (station, year))
        create_rf_plots()

@exha()
def create_rf_plots(start= -5, end=22, years=None, components='LQT'):
    """ automatically create pictures for receiver functions """
    log.info('***** Create RF plots')
    util.setRootLogger(logdebugfile=data.rf_results_dir + 'a_log_rf_plots.txt')
    stream = read_rf()
    if len(stream) == 0:
        return False
    if years == None:
        years = ['']
    for station in stream.getStationList():
        for year in years:
            stream_station = stream.select(station=station, expr='not st.mark')
            if year != '':
                stream_station = stream_station.select(expr='st.event.datetime.year==%s' % year)
                year = '_' + str(year)
            stream_station.sort('azi')
            if len(stream_station) < 3:
                continue
            ratio = len(stream_station) // 3 / 100. * 0.6
            if len(stream_station) // 3 < 150:
                ratio *= 1.8
            if len(stream_station) // 3 < 50:
                ratio *= 3 / 1.8
            fig = getFig(ratio=ratio)
            plot = stream_station.plotRF(start, end, show=show, fig=fig, scale=2, sumscale=10, component=components[1], topcolor=top, botcolor=bot, figtitle='station component', fancy_box=True)
            plot.fig.savefig(plotdir + 'rf_%s%s_%s.pdf' % (station, year, components[1]))
            plt.close(plot.fig)
            fig2 = getFig(ratio=ratio)
            plot = stream_station.plotRF(start, end, show=show, fig=fig2, scale=2, sumscale=10, component=components[2], topcolor=top, botcolor=bot, figtitle='station component', fancy_box=True)
            plot.fig.savefig(plotdir + 'rf_%s%s_%s.pdf' % (station, year, components[2]))
            plt.close(plot.fig)
            fig3 = getFig(ratio=ratio)
            plot = stream_station.plotRF(start, end, show=show, fig=fig3, scale=0.9, sumscale=2, component=components[0], topcolor=top, botcolor=bot, figtitle='station component', fancy_box=True)
            plot.fig.savefig(plotdir + 'rf_%s%s_%s.pdf' % (station, year, components[0]))
            plt.close(plot.fig)
    return True

@exha()
def create_interval_plot(expr):
    """ automatically create eps pictures for receiver functions

    use only an azimuth interval and compare sorting after azimuth and time"""
    log.info('***** Create interval plots')
    util.setRootLogger(logdebugfile=data.rf_results_dir + 'a_log_interval_plots.txt')
    stream_all = read_rf()
    for station in data.stations.keys():
        stream = stream_all.select(station=station, expr=expr)
        stream.sort(['azi', 'component'])
        pi = ['dist', 'azi']
        pl = ['time (s)', u'dist (°)', u'azi (°)']
        fig = getFig()
        fig = stream.plotRF(-5, 35, showsum=False, show=show, fig=fig, scale=5, sumscale=50, plotinfo=pi, plotlabel=pl, topcolor=top2, botcolor=bot2)
        fig.savefig(plotdir + 'rf_%s_%s_sorted_azi_Q.eps' % (station, expr))
        plt.close(fig)
        fig2 = getFig()
        fig2 = stream.plotRF(-5, 35, showsum=False, show=show, fig=fig2, scale=5, sumscale=50, component='T', plotinfo=pi, plotlabel=pl, topcolor=top2, botcolor=bot2)
        fig2.savefig(plotdir + 'rf_%s_%s_sorted_azi_T.eps' % (station, expr))
        plt.close(fig2)
        pi = ['starttime', 'azi']
        pl = ['time (s)', 'year', u'azi (°)']

        stream.sort(['starttime', 'component'])
        fig3 = getFig()
        fig3 = stream.plotRF(-5, 35, showsum=False, show=show, fig=fig3, scale=5, sumscale=50, plotinfo=pi, plotlabel=pl, topcolor=top2, botcolor=bot2)
        fig3.savefig(plotdir + 'rf_%s_%s_sorted_time_Q.eps' % (station, expr))
        plt.close(fig3)
        fig4 = getFig()
        fig4 = stream.plotRF(-5, 35, showsum=False, show=show, fig=fig4, scale=5, sumscale=50, component='T', plotinfo=pi, plotlabel=pl, topcolor=top2, botcolor=bot2)
        fig4.savefig(plotdir + 'rf_%s_%s_sorted_time_Q.eps' % (station, expr))
        plt.close(fig4)
@exha()
def profile(stations, expr=None):
    """ plots summation trace for all stations """
    log.info('***** Create profile plots')
    util.setRootLogger(logdebugfile=data.rf_results_dir + 'a_log_profile.txt')
    filename = 'sum_%s_%s' % (stations, expr)
    file_ = data.rf_results_dir + filename
    try:
        sum_stream = read(file_ + '.QHD')
    except:
        stream = read_rf()
        sum_stream = Stream()
        for station in stations.split():
            temp = stream.select(station=station, component='Q', expr=expr)
            sum_stream += temp.simpleStack()
        sum_stream.write(file_, 'Q')
    plot = sum_stream.plotProfile(-2, 21, scale=5)
    plot.fig.savefig(plotdir + filename + '.eps')
    plot.fig.savefig(plotdir + filename + '.png')
    plt.close(plot.fig)

@exha()
def profile_pp_binned(stations, bin, header='plon', expr=None, xlim=None):  #, scale=0.2 @ReservedAssignment
    """ PP profile """
    log.info('***** Create piercing point binned plot')
    util.setRootLogger(logdebugfile=data.rf_results_dir + 'a_log_profile_pp.txt')
    filename = 'profile_pp_%s_%s_dif%s' % (stations, expr, bin[1] - bin[0])
    file_ = data.rf_results_dir + filename
    try:
        stream = read(file_ + '.QHD')
    except:
        stream_all = read_rf()
        stream = stream_all.select(component='Q', expr=expr).getBinnedStream(bin, header=header)
        stream.write(file_, 'Q')
    plot = stream.plotProfile(-2, 21, scale=bin[1] - bin[0], xaxis='plon', xlabel='longitude of piercing points', plotinfo=(), figtitle='profile')
    ax = plot.fig.axes[0]
    if xlim:
        ax.set_xlim(xlim)
    # second ax with depth
    ax2 = ax.twinx()
    h = np.array((0, 50, 100, 150, 200))
    h2 = np.arange(20) * 10
    t = util.depth2time(h)
    myLocator = mpl.ticker.FixedLocator(t)
    myMinorLocator = mpl.ticker.FixedLocator(util.depth2time(h2))
    myFormatter = mpl.ticker.FixedFormatter([str(i) for i in h])
    ax2.yaxis.set_major_locator(myLocator)
    ax2.yaxis.set_minor_locator(myMinorLocator)
    ax2.yaxis.set_major_formatter(myFormatter)
    ax2.set_ylim(ax.get_ylim())
    ax2.set_ylabel('depth (km)')
    plt.show()
    #plot.fig.savefig(plotdir + filename + '.eps')    #st2.plotProfile(-2, 21, scale=0.05, xaxis = 'plon')
    #plot.fig.savefig(plotdir + filename + '.png')
    plt.close(plot.fig)

@exha()
def time_binned(station, start, end, dif=0.5, expr=None):
    log.info('Create time binned plot')
    util.setRootLogger(logdebugfile=data.rf_results_dir + 'a_log_bin_time.txt')
    filename = 'bin_time_%s_%s_dif%s' % (station, expr, dif)
    file_ = data.rf_results_dir + filename
    try:
        binned = read(file_ + '.QHD')
    except:
        binsx = np.linspace(start, end, int((end - start) / dif + 1))
        bins = [UTCDateTime(int(i // 1), 1, 1) + (i % 1) * 24 * 60 * 60 * (365 + calendar.isleap(int(i // 1))) for i in binsx]
        stream = read_rf().select(station=station, component='Q', expr=expr)
        stream.trim2(-10, 50)
        binned = stream.getBinnedStream(bins, header='starttime')
        binned.write(file_, 'Q')
    fig = binned.plotRF(topcolor=top2, fig=getFig(), botcolor=bot2, plotinfo=['sum', 'starttime'], plotlabel=['time', 'count', 'year'], show=show)
    fig.savefig(plotdir + filename + '.eps')
    plt.close(fig)

main()
