#!/usr/bin/env python
# by TR

from sito.data import IPOC
from obspy.core import UTCDateTime as UTC
from obspy.signal import PPSD
import matplotlib as mpl
import numpy as np
from obspy.xseed import Parser

def psd(station, parser):
    data = IPOC()
    ppsd_length = 6 * 3600
    overlap = 0.5
    dt = 3 * 24 * 3600
    t1 = UTC('2006-01-01')
    t2 = UTC('2013-11-01')
    ppsd = None
    print t1, t2
    while t1 < t2:
        try:
            if station != 'LVC':
                stream = data.client.getWaveform('CX', station, '', 'HHZ', t1, t1 + dt + overlap * ppsd_length)
            else:
                stream = data.client.getWaveform('GE', 'LVC', '00', 'BHZ', t1, t1 + dt + overlap * ppsd_length)

        except:
            t1 += dt
            continue
        if ppsd is None:
            ppsd = PPSD(stream[0].stats, parser=parser, skip_on_gaps=True,
                        db_bins=(-200, -50, 0.5),
                        ppsd_length=ppsd_length, overlap=overlap)
        print t1
        ppsd.add(stream)
        t1 += dt
    if ppsd is not None:
        print 'station %s: %d segments' % (station, len(ppsd.times))
        ppsd.save("/home/richter/Results/IPOC/PPSD/ppsd_%s_6h.pkl.bz2" % station, compress=True)
        return True
    else:
        return False

def psd_plot(station):
    ppsd = PPSD.load("/home/richter/Results/IPOC/PPSD/ppsd_%s_6h.pkl.bz2" % station)
    ppsd.plot(show=False, show_coverage=False, show_percentiles=True,
              period_lim=(0.05, 1000), max_percentage=10)
    fig = mpl.pyplot.figure(1)
    ax = fig.axes[0]
    ax.set_ylim((-200, -50))
    ax.set_xticklabels(['%s' % t for t in ax.get_xticks()])
    ax.set_xlabel('period (s)')
    ax.set_ylabel('amplitude (dB)')
    col = ax.collections[0]
    bbox = col.colorbar.ax.get_position()
    fig.delaxes(col.colorbar.ax)
    if COLORBAR:
        fig.colorbar(col, cax=fig.add_axes([bbox.x0 + 0.1, bbox.y0, bbox.width, bbox.height]),
                     extend='max')
        #col.colorbar.set_ticks(range(8))
        col.colorbar.set_label('probability (%)')
    bbox = ax.get_position()
    ax.set_position([bbox.x0, bbox.y0, bbox.width + 0.2 - 0.1 * COLORBAR, bbox.height])
    fig.savefig("/home/richter/Results/IPOC/PPSD/ppsd_%s_6h_cb.png" % station, dpi=150, bbox_inches='tight')
    mpl.pyplot.close()

def correct_sensitivity(station):
    ppsd = PPSD.load("/home/richter/Results/IPOC/PPSD/ppsd_%s_6h_wrong_gain.pkl" % station)
    dB = 10 * np.log10(2516580000.0 / 629145000.0)
    if station == 'PB14':
        dB = 10 * np.log10(588000000.0 / 629145000.0)
    ppsd.spec_bins = ppsd.spec_bins - dB
    ppsd.yedges = ppsd.yedges - dB
    ppsd.save("/home/richter/Results/IPOC/PPSD/ppsd_%s_6h.pkl" % station)




COLORBAR = True
mpl.rcParams.update({'font.size': 14})
stations = 'PATCX PB01 PB02 PB03 PB04'
stations = 'LVC PB09 PB10 PB11 PB12 PB13 PB14 PB15 PB16 HMBCX PATCX PSGCX MNMCX PB02 PB03 PB04 PB05 PB06 PB07 PB08'.split()
stations = 'PSGCX MNMCX PB02 PB03 PB04 PB05 PB06 PB07 PB08'.split()
#stations = ('PB01',)  # PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 PB13 PB14 PB15 PB16 PB17 HMBCX MNMCX PSGCX PATCX LVC'.split()
#stations = 'PB09 LVC'.split()

parser = Parser('/home/richter/Data/paz/ipoc_inventory.dseed')
for st in stations:
    print st
    #correct_sensitivity(st)
    suc = psd(st, parser)
    if suc:
        psd_plot(st)
