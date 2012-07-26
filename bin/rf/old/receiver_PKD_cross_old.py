#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
#@PydevCodeAnalysisIgnore

import logging
from sito import read, Stream, util, imaging
import numpy as np

log = logging.getLogger('receiver_cross')

def pkd1(station, year1, year2, multiples=False, thrown=False, usetheo=True):
    """ used for first processing data """
    resfile = resdir + station
    util.setRootLogger(resdir + 'log_' + station + '.txt')
    stream = Stream()
    for year in range(year1, year2 + 1):
        filename = sourcedir + 'PKD_' + str(year) + '.QHD'
        try:
            stream.extend(read(filename))
        except:
            log.warning('can not read file ' + filename)
    if len(stream) == 0:
        return
    stream.setHI('mark', False)
    stream.writeStationPosition(stationfile)
    stream.pspier(60, stationfile)
    #stream = stream1
    log.info('number of events: %s' % (len(stream) // 3))

    #stream.integrate()
    stream.filter2(0.033, 2.)
    stream.trim2(-50, 300)
    stream.sort(['eventno', 'station', 'component'])
    stream.check()
    stream.rotateZNE2LQT(-5, 15, usetheo=usetheo)
    st_thrown = stream.afarm(signoise=1.5)
    log.info('number of events after first farm: %s' % (len(stream) // 3))
    util.ipshell()
    stream.receiverf(water=0.005, gauss=10, tshift=20, pad=0,
                     window='tukey', start= -10, end=30, where='ponset', lenslope=5)
    stream.sort('azi')
    if thrown:
        st_thrown.receiverf(water=0.005, gauss=10, tshift=20, pad=0,
                    window='tukey', start= -10, end=30, where='ponset', lenslope=5)
    st_thrown = st_thrown + stream.afarm('rf', signoise=4., signoiseQ=1.)
    #st_thrown = stream.afarm('rf', signoise = 3., signoiseQ=1.)
    log.info('number of events after second farm: %s' % (len(stream) // 3))
    if len(stream) == 0:
        log.warning('no events left')
        return
    if multiples:
        st2 = stream.copy()
        st2.moveout(phase='Ppps')
        st2.trim2(-20, 100)
        st2.write(resfile + '_PPPS', 'Q')
        st3 = stream.copy()
        st3.moveout(phase='Ppss')
        st3.trim2(-20, 100)
        st3.write(resfile + '_PPSS', 'Q')
    util.ipshell()
    stream.moveout(phase='Ps')
    util.ipshell()
    stream.trim2(-20, 100)
    stream.write(resfile + '_' + str(year1) + '_MOUT', 'Q')
    if thrown:
        st_thrown.moveout(phase='Ps')
        st_thrown.trim2(-20, 100)
        st_thrown.write(resfile + '_THROWN', 'Q')

def pkd1_2(station, year1, year2, multiples=False, usetheo=True):
    """ used for first processing data """
    resfile = resdir + station
    util.setRootLogger(resdir + 'log_' + station + '.txt')
    stream = Stream()
    for year in range(year1, year2 + 1):
        filename = sourcedir + 'PKD_' + str(year) + '.QHD'
        try:
            stream.extend(read(filename))
        except:
            log.warning('can not read file ' + filename)
    if len(stream) == 0:
        return
    stream.setHI('mark', False)
    stream.writeStationPosition(stationfile)
    stream.pspier(60, stationfile)
    #stream = stream1
    log.info('number of events: %s' % (len(stream) // 3))

    #stream.integrate()
    stream.filter2(0.1, 2.)
    stream.trim2(-50, 300)
    stream.sort(['eventno', 'station', 'component'])
    stream.check()
    stream.rotateZNE2LQT(-5, 15, usetheo=usetheo)
    stream.afarm(signoise=1.2, remove=False)
    stream.receiverf(water=0.005, gauss=10, tshift=20, pad=0,
                     window='tukey', start= -10, end=30, where='ponset', lenslope=5)
    stream.sort('azi')
    stream.afarm('rf', signoise=2., signoiseQ=1., maxL=1 / 1.5,
                           sigQ=False, broad=True, remove=False)
    if multiples:
        st2 = stream.copy()
        st2.moveout(phase='Ppps')
        st2.trim2(-20, 100)
        st2.write(resfile + '_PPPS', 'Q')
        st3 = stream.copy()
        st3.moveout(phase='Ppss')
        st3.trim2(-20, 100)
        st3.write(resfile + '_PPSS', 'Q')
    #stream.moveout(phase='Ps')
    stream.trim2(-20, 100)
    #stream.write(resfile + '_' + str(year1) + '_MOUT', 'Q')
    stream.write(resfile + '_' + str(year1), 'Q')

def pkd_select(station, year1, year2):
    resfile = resdir + station
    util.setRootLogger(resdir + 'log_select_' + station + '.txt')
    stream = Stream()
    for year in range(year1, year2 + 1):
        filename = resfile + '_' + str(year) + '.QHD'
        try:
            stream.extend(read(filename))
        except:
            log.warning('can not read file ' + filename)
    #stream=stream.select(expr='225<st.azi<265.7 and not st.mark')
    stream.sort('ponset')
    stream.downsample2(10)
    stream.write(resfile + '_NOMOUT', 'Q')
    #stream.plotRF(topcolor='red', botcolor='blue', plotinfo=['ponset'], plotlabel=['time', 'year'])
    #util.ipshell()

def ipoc2(station):
    """ load data and throwed data to produce event files with used and not used events """
    resfile = resdir + station
    stream = read(resfile + '_MOUT.QHD')
    thrown = read(resfile + '_THROWN.QHD')

    import pde
    pde_thrown = pde.PDE.readFromStream(thrown)
    pde_selected = pde.PDE.readFromStream(stream)
    pde_thrown.write(resdir + 'events_' + station + '_thrown.txt', format_='GMT')
    pde_selected.write(resdir + 'events_' + station + '_selected.txt', format_='GMT')

def ipoc3(station, suffix_in='_MOUT', suffix_out=''):
    """ automatically create eps pictures for receiver functions """
    try:
        stream = read(resdir + station + suffix_in + '.QHD')
    except:
        return False
    axes = [1. - 0.1 * 2] + [0.1] * 2
    margin = [0.7, 0.1, 0.9, 0.7] #left, rigth, bottom, top
    if len(stream) == 0:
        return False
    fig = imaging.getPublicationFigure(axes, width=17, ratio=0.75, margin=margin, fontsize=10, labelsize=8, usetex=False)
    fig = stream.plotRF(-5, 22, show=False, fig=fig, scale=5, sumscale=50)
    fig.savefig(plotdir + 'rf_%s%s_Q.eps' % (station, suffix_out))
    fig2 = imaging.getPublicationFigure(axes, width=17, ratio=0.75, margin=margin, fontsize=10, labelsize=8, usetex=False)
    fig2 = stream.plotRF(-5, 22, show=False, fig=fig2, scale=5, sumscale=50, component='T')
    fig2.savefig(plotdir + 'rf_%s%s_T.eps' % (station, suffix_out))
    fig3 = imaging.getPublicationFigure(axes, width=17, ratio=0.75, margin=margin, fontsize=10, labelsize=8, usetex=False)
    fig3 = stream.plotRF(-5, 22, show=False, fig=fig3, scale=5, sumscale=50, component='L')
    fig3.savefig(plotdir + 'rf_%s_L.eps' % station)

    return True
    #util.ipshell()

def ipoc4(station):
    """ zmigration """
    try:
        stream = read(resdir + 'M5.5_' + station + '_MOUT.QHD')
    except:
        return False
    stream.zmigr()
    return True

def ipoc5(station, suffix_in='_MOUT', suffix_out='_azi140-160'):
    """ automatically create eps pictures for receiver functions

    use only an azimuth interval and compare sorting after azimuth and time"""
    try:
        stream = read(resdir + station + suffix_in + '.QHD')
    except:
        return False
    axes = [1. - 0.1 * 2] + [0.1] * 2
    margin = [0.7, 0.1, 0.9, 0.7] #left, rigth, bottom, top
    if len(stream) == 0:
        return False
    expr = 'st.azi<=160 and st.azi >= 140'
    stream = stream.select(expr=expr)
    stream.sort(['azi', 'component'])
    pi = ['dist', 'azi']
    pl = ['time (s)', u'dist (°)', u'azi (°)']
    fig = imaging.getPublicationFigure(axes, width=17, ratio=0.75, margin=margin, fontsize=10, labelsize=8, usetex=False)
    fig = stream.plotRF(-1, 11, showsum=False, show=False, fig=fig, scale=5, sumscale=50, plotinfo=pi, plotlabel=pl)
    fig.savefig(plotdir + 'rf_azi140-160_sorted_azi_%s_Q.eps' % (station))

    fig2 = imaging.getPublicationFigure(axes, width=17, ratio=0.75, margin=margin, fontsize=10, labelsize=8, usetex=False)
    fig2 = stream.plotRF(-1, 11, showsum=False, show=False, fig=fig2, scale=5, sumscale=50, component='T', plotinfo=pi, plotlabel=pl)
    fig2.savefig(plotdir + 'rf_azi140-160_sorted_azi_%s_T.eps' % (station))

    pi = ['starttime', 'azi']
    pl = ['time (s)', 'year', u'azi (°)']
    uht = 'time'

    stream.sort(['starttime', 'component'])
    fig3 = imaging.getPublicationFigure(axes, width=17, ratio=0.75, margin=margin, fontsize=10, labelsize=8, usetex=False)
    fig3 = stream.plotRF(-1, 11, showsum=False, show=False, fig=fig3, scale=5, sumscale=50, plotinfo=pi, plotlabel=pl, usehardticks=uht)
    fig3.savefig(plotdir + 'rf_azi140-160_sorted_time_%s_Q.eps' % (station))

    fig4 = imaging.getPublicationFigure(axes, width=17, ratio=0.75, margin=margin, fontsize=10, labelsize=8, usetex=False)
    fig4 = stream.plotRF(-1, 11, showsum=False, show=False, fig=fig4, scale=5, sumscale=50, component='T', plotinfo=pi, plotlabel=pl, usehardticks=uht)
    fig4.savefig(plotdir + 'rf_azi140-160_sorted_time_%s_T.eps' % (station))

def ipoc6(stations, expr=None, label=''):
    """ plots summation trace for all stations """
    st = Stream()
    for station in stations.split():
        temp = read(resdir + station + '_MOUT.QHD').select(component='Q', expr=expr)
        st += temp.sum()
    fname = resdir + 'SUM' + label
    st.write(fname, 'Q')
    st = read(fname + '.QHD')
    fig = st.plotProfile(-2, 21, scale=5)
    fig.savefig(plotdir + 'profile' + label + '.eps')

def ipoc7(stations, expr=None, label=''):
    """ calculate and plot piercing points """
    st = Stream()
    for station in stations.split():
        try:
            st += read(resdir + station + '_MOUT.QHD')
        except:
            pass
    st2 = st.select(component='Q')
    st3 = st2.getBinnedStream(np.arange(-71.025, -67.999, 0.05), header='plon')
    st4 = st2.getBinnedStream(np.arange(-71.025, -67.999, 0.005), header='plon')
    #st5 = st2.getBinnedStream(np.arange(-71.025, -67.999, 0.05), [-25, -21, -18])
    #st5.select(expr='st.plat >= -21 and st.plat <= -18').plotProfile(-2, 21, scale=0.2, xaxis = 'plon', plotcount=True)
    #st5.select(expr='st.plat >= -25 and st.plat <= -21').plotProfile(-2, 21, scale=0.2, xaxis = 'plon', plotcount=True)
    fig1 = st3.plotProfile(-2, 21, scale=0.2, xaxis='plon', plotcount=True)
    fig2 = st4.plotProfile(-2, 21, scale=0.05, xaxis='plon')
    fig1.savefig(plotdir + 'profile' + label + '_deg0.2.eps')
    fig2.savefig(plotdir + 'profile' + label + '_deg0.05.eps')
    #st2.plotProfile(-2, 21, scale=0.05, xaxis = 'plon')



def durchgang():
    for i in 'PB01 PB03'.split():
        ipoc1(i, multiples=True, thrown=True, usetheo=True)
    for i in 'PB02 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 MNMCX PATCX PSGCX HMBCX LVC'.split():
    #for i in 'PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 MNMCX PATCX PSGCX HMBCX LVC'.split():
        ipoc1(i)
    ipoc2('PB01')
    ipoc2('PB03')

    for i in 'PB01 PB03'.split():
        for j, k in (('_MOUT', ''), ('_PPPS', '_ppps'), ('_PPSS', '_ppss')):
            print('%s: %s' % (i, ipoc3(i, j, k)))

    for i in 'PB02 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 MNMCX PATCX PSGCX HMBCX LVC'.split():
        print('%s: %s' % (i, ipoc3(i)))

    ipoc5('PB01')
    ipoc5('PB03')
    ipoc7('PB05 PB04 PB07 PB03 PB06 PB01 PB08', label='')
    ipoc7('PB05 PB04 PB02 PB07 HMBCX PB03 MNMCX PB06 PB01 PB08', label='2')

stationfile = '/home/richter/Data/stations.txt'
sourcedir = '/home/richter/Data/Parkfield/receiver/M5.5_events/'
resdir = ipocresdir = '/home/richter/Data/Parkfield/receiver/M5.5_results_theorot_2/'
plotdir = ipocplotdir = '/home/richter/Results/Parkfield/receiver/M5.5_theorot/plots/'

#binning
def binning():
    import calendar
    from obspy.core import UTCDateTime
    binsx = np.linspace(1996.74, 2011.24, int((2011.24 - 1996.74) / 0.5 + 1))
    bins = [UTCDateTime(int(i // 1), 1, 1) + (i % 1) * 24 * 60 * 60 * (365 + calendar.isleap(int(i // 1))) for i in binsx]
    ms = read('/home/richter/Data/Parkfield/receiver/M5.5_results_theorot_2/PKD_MOUT.QHD')
    ms.downsample2(20)
    #ms2 = ms.select(expr='st.azi>=225 and st.azi<=235.4')
    ms2 = ms.select(expr='st.azi>=225 and st.azi<=265.7 and not st.mark')
    ms4 = ms2.select(component='Q').getBinnedStream(bins, header='starttime')
    #ms5 = ms3.select(component='Q').getBinnedStream(bins, header='starttime')
    #ms4.plotRF(topcolor='red', botcolor='blue', plotinfo=['azi','ponset'], plotlabel=['time', 'azi', 'year'])
    ms4.plotRF(topcolor='red', botcolor='blue', plotinfo=['ponset'], plotlabel=['time', 'year'])
    def xcorr(stream, t1, t2):
        st = stream.copy()
        st.trim2(t1, t2)
        tr_sum = st.sum()
        st.xcorr(tr_sum, 1)
        return st

    list = []
    times = ((5., 7.5), (7., 10), (9.5, 12), (15.5, 18.5), (20.5, 23.5), (27, 29.5), (33.5, 36.5))
    for t1, t2 in times:
        list.append(xcorr(ms4, t1, t2))
    maxima = [st.getMaxima() for st in list]
    for st in list:
        st.plotRF(relative='lagtime0', topcolor='grey', botcolor='white', plotinfo=['sum'], plotlabel=['time', 'count'], showsum=False)
    from pylab import subplot, plot, figure, show, ylabel, legend
    figure()
    a = subplot(611)
    for i in range(4):
        plot(maxima[i], label=str(times[i]))
    ylabel('max of corr')
    legend()
    subplot(612, sharex=a)
    for i in range(3):
        plot(maxima[i + 4], label=str(times[i + 4]))
    legend()

    subplot(613, sharex=a); plot(list[0].getHI('azi')); ylabel('azi')
    subplot(614, sharex=a); plot(list[0].getHI('dist')); ylabel('dist')
    subplot(615, sharex=a); plot(list[0].getHI('event.depth')); ylabel('depth')
    subplot(616, sharex=a); plot(list[0].getHI('sum')); ylabel('count')
    show()
    util.ipshell()
    #xcorr(ms4, 10, 13)
    #xcorr(ms4, 19, 24)
    #xcorr(ms4, 30, 40)
    #xcorr(ms4, 39, 43)


def pkd_cross(): pass


if __name__ == '__main__':
    #for year in range(1996, 2011):
    #    pkd1_2('PKD', year, year)
    # pkd_select('PKD', 1996, 2011)
    binning()
    pkd_cross()
    #util.setRootLogger(resdir + 'log2_' + 'PKD' + '.txt')
    #stream = Stream()
    #for year in range(1996, 2011):
    #    print year
    #    filename = resdir + 'PKD_' + str(year) + '_MOUT.QHD'
    #    print filename
    #    try:
    #        stream.extend(read(filename))
    #    except Exception as exception:
    #        log.warning('can not read file %s:\n%s' %(filename, exception))
    #stream.write(resdir + 'PKD_MOUT', 'Q')

    #ipoc3('PKD', '_MOUT', '')
    #for i in 'PB01 PB03'.split():
    #    ipoc5(i)
    # ipoc6('PB05 PB04 PB07 PB03 PB06 PB01 PB08')
    #ipoc7('PB05 PB04 PB07 PB03 PB06 PB01 PB08', label='')
    #ipoc7('PB05 PB04 PB02 PB07 HMBCX PB03 MNMCX PB06 PB01 PB08', label='2')
    #ipoc2('PB03')
    #ipoc6('PB05 PB04 PB07 PB03 PB06 PB01 PB08', label='')
    #ipoc6('PB05 PB04 PB02 PB07 HMBCX PB03 MNMCX PB06 PB01 PB08', label='2')
    #ipoc6('PB04 PB07 PB03 PB06 PB01', label='3')
    #ipoc6('PB05 PB04 PB07 PB03 PB06 PB01 PB08', 'st.azi<=160 and st.azi >= 140', '_azi140-160')
    #ipoc6('PB05 PB04 PB02 PB07 HMBCX PB03 MNMCX PB06 PB01 PB08', 'st.azi<=160 and st.azi >= 140', '2_azi140-160')
    #ipoc3('PB03')

