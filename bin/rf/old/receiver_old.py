#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
#@PydevCodeAnalysisIgnore

import logging
from mystream import Stream
from sito import util, imaging, read
import numpy as np

log = logging.getLogger('receiver')

def test():
    logging.basicConfig(level=logging.DEBUG)
    stream = read('../tests/data/TEST.QHD')
    # processing = 'trim-20,100 rot-2,10 dec0.005,10,20,0,tukey,-10,50,ponset,5 mout' # rot dec decSH mout
    # process(stream, processing, 'mout', '/home/richter/Dev/tspy/tests/receivertestdata/%s%s')
    stream.filter2(0.05, False)
    stream.trim2(-50, 300)
    stream.rotateZNE2LQT(-2, 20)
    stream.plot_()
    print('stream: ' + str(stream))
    stream.afarm(signoise=1.5)
    print('stream: ' + str(stream))
    stream.receiverf(water=0.005, gauss=10, tshift=20, pad=0,
                    window='tukey', start= -10, end=30, where='ponset', lenslope=5)
    stream.afarm('rf', signoise=4)
    stream.moveout()
    stream.trim2(-20, 100)
    print('stream: ' + str(stream))
    stream.plot_()
    stream.plotRF()
    util.ipshell()

def ipoc1(station, multiples=False, thrown=False, usetheo=True):
    """ used for first processing data """
    resfile = resdir + station
    util.setRootLogger(resdir + 'log_' + station + '.txt')
    try: stream1 = read('IPOC_M5.5_2006_' + station + '.QHD')
    except: stream1 = Stream()
    try: stream2 = read(ipocsourcedir + 'IPOC_M5.5_2007_' + station + '.QHD')
    except: stream2 = Stream()
    try: stream3 = read(ipocsourcedir + 'IPOC_M5.5_2008_' + station + '.QHD')
    except: stream3 = Stream()
    try: stream4 = read(ipocsourcedir + 'IPOC_M5.5_2009_' + station + '.QHD')
    except: stream4 = Stream()
    try: stream5 = read(ipocsourcedir + 'IPOC_M5.5_2010_' + station + '.QHD')
    except: stream5 = Stream()
    stream = stream1 + stream2 + stream3 + stream4 + stream5
    if len(stream) == 0:
        return
    stream.writeStationPosition(stationfile)
    stream.pspier(60, stationfile)
    #stream = stream1
    log.info('number of events: %s' % (len(stream) // 3))

    #stream.integrate()
    stream.filter2(0.033, 2.)
    stream.trim2(-50, 300)
    stream.check()
    stream.rotateZNE2LQT(-5, 15, usetheo=usetheo)
    st_thrown = stream.afarm(signoise=1.5)
    log.info('number of events after first farm: %s' % (len(stream) // 3))
    stream.receiverf(water=0.005, gauss=10, tshift=20, pad=0,
                    window='tukey', start= -10, end=30, where='ponset', lenslope=5)
    stream.sort('azi')
    if thrown:
        st_thrown.receiverf(water=0.005, gauss=10, tshift=20, pad=0,
                    window='tukey', start= -10, end=30, where='ponset', lenslope=5)
    st_thrown = st_thrown + stream.afarm('rf', signoise=4., signoiseQ=1.)
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
    stream.moveout(phase='Ps')
    stream.trim2(-20, 100)
    stream.write(resfile + '_MOUT', 'Q')
    if thrown:
        st_thrown.moveout(phase='Ps')
        st_thrown.trim2(-20, 100)
        st_thrown.write(resfile + '_THROWN', 'Q')

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

stationfile = '/home/richter/Data/stations_ipoc.txt'
sourcedir = ipocsourcedir = '/home/richter/Data/IPOC/receiver/M5.5_events/'
#resdir = ipocresdir = '/home/richter/Data/IPOC/receiver/M5.5_results_theorot/'
#plotdir = ipocplotdir = '/home/richter/Results/IPOC/receiver/M5.5_theorot/plots/'
resdir = ipocresdir = '/home/richter/Data/IPOC/receiver/M5.5_results_polrot/'
plotdir = ipocplotdir = '/home/richter/Results/IPOC/receiver/M5.5_polrot/plots/'
resdir = ipocresdir = '/home/richter/Data/IPOC/receiver/M5.5_results_int/'
plotdir = ipocplotdir = '/home/richter/Results/IPOC/receiver/M5.5_int/plots/'


if __name__ == '__main__':
    durchgang()
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


    util.ipshell()
