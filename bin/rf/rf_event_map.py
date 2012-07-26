#!/usr/bin/env python
# by TR
"""
create event map (selected events vs all events used for calculating rfs)
"""

from sito import events, map, read
import pylab as plt
from glob import glob
import os.path

def create_event_files(station='*'):
    for file_ in glob(path + station + '_mout.QHD'):
        ms = read(file_).select(component='Q', expr='not st.mark')
        station = ms[0].stats.station
        events = ms.getEvents()
        events.write(em_path + 'events_%s.txt' % station)

def create_event_maps(station='*'):
    all_events = events.Events.read(eventfile)
    for file_ in glob(em_path + 'events_%s.txt' % station):
        st_events = events.Events.read(file_)
        station = (os.path.splitext(os.path.basename(file_))[0]).split('_')[1]
        m = map.createRFEventMap(show=False, trench=None)
        all_events.plot_(m, color='b', radius=5, alpha=0.5)
        st_events.plot_(m, color='r', radius=8, alpha=0.5)
        if annotate:
            # fancy box
            from matplotlib import transforms
            ax = m._check_ax()
            offset = transforms.ScaledTranslation(-20. / 72, 20. / 72,
                                                  plt.gcf().dpi_scale_trans)
            trans = ax.transAxes + offset
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
            ax.text(1., 0., station, transform=trans, fontsize=14,
                    verticalalignment='bottom', ha='right', bbox=props)
            plt.gcf().savefig(em_path + 'map_events_%s_anno.pdf' % station)
        else:
            plt.gcf().savefig(em_path + 'map_events_%s.pdf' % station)

def create_region_map():
    all_events = events.Events.read(eventfile)
    import math
    w = math.pi * 4000 * 1000 * 186 / 180.
    dict_basemap_ipoc = dict(lat_0= -21, lon_0= -69.5,
                             projection='aeqd', resolution='c',
                             width=w,
                             height=w)
    m = map.createRFEventMap(show=False, trench=None, circles=(), circles_around=regions, dict_basemap=dict_basemap_ipoc)
    all_events.plot_(m, color='b', radius=5, alpha=0.5)
    ax = m._check_ax()
    ax.text(9.6e6, 1.8e6, 'R1')
    ax.text(4.0e6, 10.5e6, 'R2')
    ax.text(3.2e6, 11.1e6, 'R3')
    plt.gcf().savefig(em_path + 'map_regions.pdf')

regions = ((-58.3, -22.0, 800, False, 'r'), #1 South Sandwich
            (14.1, -91.2, 600, False, 'r'), #2 Guatemala           
           (13.5, -92.1, 1400, False, 'r')) #3 Mexico
annotate = True
eventfile = '/home/richter/Data/events/2012_03_events_27-93_mag5.5_IPOC.txt'
path = '/home/richter/Results/IPOC/receiver/2012_mag5.5/'
em_path = path + 'event_maps/'

#create_event_files()
#create_event_maps()
#annotate = False
#create_event_maps()

create_region_map()
