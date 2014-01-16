#!/usr/bin/env python
# by TR

from sito import map
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.path import Path
import matplotlib.patches as patches
from xml.dom import minidom
from math import cos, pi

def svg2path(svg):
    doc = minidom.parse(svg)  # parseString also exists
    path_strings = [path.getAttribute('d') for path
                    in doc.getElementsByTagName('path')]
    doc.unlink()
    paths = []
    for path_string in path_strings:
        codes = []
        verts = []
        for atom in path_string.split():
            if atom == 'M':
                code = Path.MOVETO
            elif atom == 'Z':
                paths.append(Path(verts + [(0, 0)], codes + [Path.CLOSEPOLY]))
                codes = []
                verts = []
            elif atom == 'C':
                code = Path.CURVE4
            elif atom == 'L':
                code = Path.LINETO
            else:
                coords = tuple(float(c) for c in atom.split(','))
                verts.append(coords)
                codes.append(code)
    return paths


def plot_some_events():
    from obspy.core.event import Catalog, Event, Origin, Magnitude
    from obspy.core import UTCDateTime as UTC

    eqs = """2008-09-10T16:12:03    6.0    -20.40    -69.40     40
    2008-03-24T20:39:06    5.9    -20.10    -69.20     85
    2008-03-01T19:51:59    5.7    -20.10    -69.60     15
    2008-02-15T16:54:04    5.5    -23.00    -70.20     32
    2008-02-04T17:01:30    6.6    -20.20    -70.00     36
    2007-12-16T08:09:16    7.1    -22.80    -70.00     14
    2007-11-14T15:40:51    7.8    -22.34    -70.06     37"""  #GEOFON:-22.30    -69.80
    events = []
    for eq in eqs.split('\n'):
        time, mag, lat, lon, depth = eq.split()
        ev = Event(event_type='earthquake', creation_info='GEOFON',
                    origins=[Origin(time=UTC(time), latitude=float(lat),
                                    longitude=float(lon), depth=float(depth))],
                    magnitudes=[Magnitude(mag=float(mag), magnitude_type='M')])
        events.append(ev)
    cat = Catalog(events[::-1])
    #print cat
    #cat.plot(projection='local')
    lons = [ev.origins[0].longitude for ev in cat]
    lats = [ev.origins[0].latitude for ev in cat]
    dates = [ev.origins[0].time for ev in cat]
    mags = [ev.magnitudes[0].mag for ev in cat]
    #map.createIPOCMap(show=False, ll=(-25, -71), ur=(-18, -68.5), earthquake=(lons, lats, mags), cities='ipoc',
    #                  save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_some_earthquakes_with_mag.png')
    #map.createIPOCMap(show=False, ll=(-25, -71), ur=(-18, -68.5), earthquake=(lons, lats), cities='ipoc',
    #                  save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_some_earthquakes.png')
    #map.createIPOCMap(show=False, ll=(-25, -71), ur=(-18, -68.5), earthquake=(lons[:1], lats[:1]), cities='ipoc',
    #                  save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_just_Toco.png')

#plot_some_events()

def get_slipmodel():
    fname = '/home/richter/Data/bernd_slip/jgrb17176-sup-0002-ds01.txt'
    import numpy
    data = numpy.loadtxt(fname, skiprows=1)
    lat = data[:, 0]
    lon = data[:, 1]
    z = data[:, 3]
    #return (lon, lat, z) + ((0.5, 1, 1.5, 2, 2.5, 3), 'gray')
    return (lon, lat, z) + ((0.5, 1, 1.5, 2, 2.5, 3), 'gray')

def plot_results1(map_dic=None, scale=1):
    from sito.stations import IPOCStations
    from matplotlib.lines import Line2D
    # see /home/richter/Results/IPOC/veldrop_vs_groundmotion_stations.txt
    map_dic2 = dict(show=False, ll=(-25, -71.25), ur=(-18, -68.25),
                    slip=get_slipmodel(), cities='ipoc', stations=None)
    if map_dic is not None:
        map_dic2.update(map_dic)
    map_dic = map_dic2
    m = map.createIPOCMap(**map_dic)
    marker = 'o'
    scale2 = 2.
    ipoc = IPOCStations()
    ipoc_no_data = ipoc.pick('PB09 PB10 PB11 PB12 PB13 PB14 PB15', replace=False)
    ipoc_no_data.plot(m, marker=marker, mfc='w', ms=8 * scale, zorder=20)
    #1-3Hz
    ipoc_bad_data = ipoc.pick('MNMCX ', replace=False)
    ipoc_no_signal = ipoc.pick('LVC PB02 PB06 PB08 PSGCX', replace=False)  #<0.2%
    ipoc_small_signal = ipoc.pick('HMBCX PATCX PB01 PB03 PB04 PB05 PB07', replace=False)  #0.2-0.5%
    ipoc_big_signal = ipoc.pick('PB04 PB05', replace=False)  #>0.5%

    ipoc_bad_data.plot(m, marker=marker, mfc='#00FFFF', ms=8 * scale * scale2, zorder=20)
    ipoc_no_signal.plot(m, marker=marker, mfc='y', ms=8 * scale * scale2, zorder=20)
    ipoc_small_signal.plot(m, marker=marker, mfc='orange', ms=8 * scale * scale2, zorder=20)
    ipoc_big_signal.plot(m, marker=marker, mfc='r', ms=8 * scale * scale2, zorder=20)

    #4-6Hz
    ipoc_bad_data = ipoc.pick('', replace=False)
    ipoc_no_signal = ipoc.pick('HMBCX MNMCX PB01 PB06 LVC', replace=False)  #<0.2%
    ipoc_small_signal = ipoc.pick('PB04 PB05 PB07 PB08 PSGCX ', replace=False)  #0.2-0.5%
    ipoc_big_signal = ipoc.pick('PB02 PB03 PATCX', replace=False)  #>0.5%

    ipoc_bad_data.plot(m, marker=marker, mfc='#00FFFF', ms=8 * scale, zorder=21, annotate=False)
    ipoc_no_signal.plot(m, marker=marker, mfc='y', ms=8 * scale, zorder=21, annotate=False)
    ipoc_small_signal.plot(m, marker=marker, mfc='orange', ms=8 * scale, zorder=21, annotate=False)
    ipoc_big_signal.plot(m, marker=marker, mfc='r', ms=8 * scale, zorder=21, annotate=False)

    colors = 'w #00FFFF y orange r'.split()
    labels = 'no data,bad data,<0.2%,0.2-0.5%,>0.5%,1Hz-3Hz\n4Hz-6Hz'.split(',')
    legend_kwargs = dict(loc='center left', bbox_to_anchor=(1.22, 0.5))

    def myhandler(legend, orig_handle, fontsize, handlebox):
        w, h, x, y = handlebox.width, handlebox.height, handlebox.xdescent, handlebox.ydescent
        xm, ym = x + w / 2, y + h / 2
        s1, s2 = 4 * scale, 4 * scale * scale2
        a_list = [Line2D((xm, xm + 0.7 * w), (ym, 0.8 * ym + h), color='k'),
                  Line2D((xm,), (ym,), marker=marker, ms=2 * s2, color='w'),
                  Line2D((xm, xm + 0.7 * w), (ym, ym - 0.5 * h), color='k'),
                  Line2D((xm,), (ym,), marker=marker, ms=2 * s1, color='w')]
        for a in a_list:
            handlebox.add_artist(a)
    mpl.legend.Legend.update_default_handler_map({None: myhandler})
    plt.gca().legend([Line2D((0,), (1,), marker=marker, ls='', mfc=c, ms=8 * scale) for c in colors] +
              [None], labels, numpoints=1, **legend_kwargs)
    mpl.rcParams.update({'lines.linewidth':1.})
    m.drawmapscale(-68.7, -24., -70, -21.5, 50, fontsize=7, yoffset=0.005 * (m.ymax - m.ymin))
    plt.gcf().savefig('/home/richter/Documents/pics/maps/ipoc/ipoc_map_results3.pdf')

def plot_small_map_for_F():
    from sito.stations import IPOCStations
    from matplotlib.lines import Line2D

    ipoc = IPOCStations()
    ipoc.pick('PATCX')
    m = map.createIPOCMap(show=False, ll=(-24, -70.9), ur=(-18, -68.5),
                          earthquake='Tocopilla_position', stations=ipoc,
                          cities='Tocopilla')

def plot_salar_map():
    map_dic = dict(show=False, ll=(-21.4, -70.3), ur=(-20.7, -69.8),
                   figsize=(fw, 1.61 * fw * 0.7),
                   margin=(0.05, 0.05, 0.9, 0.9), lw=2,
                   station_markersize=4,
                   grid=0.2, grid_labels=True, grid_lw=0.2, slip=None, earthquake=None,
                   countries=None, coastlines=None,
                   elevation_args=(1592 * 2, None, False), elevation_offset=1000,
                   shaded=True, shaded_args=(90, 45, 0.7),
                   colormap=cm.binary,
                   elevation='/home/richter/Data/map/salar_90m.tif',
                   elev_oceans=False,
                   stations=None,
                   spines_lw=2,
                   loffset=1000)
    from sito.data import IPOC
    m = map.createIPOCMap(**map_dic)
    chos = [('CHO1', -21.094050, -70.102000, 653),
            ('CHO2', -21.105933, -70.096900, 620),
            ('CHO3', -21.106233, -70.097517, 625)]
    kw = dict(bbox=dict(boxstyle="round", fc="w", alpha=0.5, ec='none'))
    IPOC().stations.plot(m, mfc='w', ms=4, zorder=10, lsize='small', kwargs_an=kw)
    #ASTER GDEM is a product of METI and NASA.
    for station, lat, lon, height in chos[1:2]:
        x, y = m(lon, lat)
        m.plot((x,), (y,), marker='o', mfc='w', ms=4, zorder=10)
        plt.annotate(station, (x, y), xytext=(3, 3), textcoords='offset points', size='small', **kw)
    #plt.annotate('ASTER GDEM is a product of METI and NASA.', (1, 0), xycoords='axes fraction', ha='right', va='bottom', size='xx-small')
    plt.annotate('Salar Grande', (0.62, 0.56), rotation=-80, xycoords='axes fraction',
                 ha='center', va='center', size='small', color='k')
    path = svg2path('/home/richter/Documents/pics/maps/ipoc/salar/salar.svg')[0]
    lat21, lon70, px2deg = 862 - 15, 557 + 15, (1.177425 - 0.869511) / 1055
    for i in range(len(path)):
        x, y = path.vertices[i, :]
        lon = -70 + (x - lon70) * px2deg / cos(21. / 180 * pi)
        lat = -21 - (y - lat21) * px2deg
        path.vertices[i, :] = m(lon, lat)
    patch = patches.PathPatch(path, facecolor='none', lw=0.8, ec='r', alpha=0.5, zorder=50)
    plt.gca().add_patch(patch)
    mpl.rcParams.update({'lines.linewidth':1.})
    m.drawmapscale(-70.2, -21.35, -70, -21, 10, fontsize=7, yoffset=0.005 * (m.ymax - m.ymin))
    plt.gcf().savefig('/home/richter/Documents/pics/maps/ipoc/salar_map.pdf', dpi=1200)

def plot_PATCX_map():
    map_dic = dict(show=False, ll=(-20.91, -70.22), ur=(-20.69, -70.08),
                   figsize=(fw, 1.61 * fw * 0.7),
                   margin=(0.05, 0.05, 0.9, 0.9), lw=2,
                   station_markersize=4,
                   grid=0.1, grid_labels=True, grid_lw=0.2, slip=None, earthquake=None,
                   countries=None, coastlines=None,
                   elevation_args=(1592 * 2, None, False), elevation_offset=1000,
                   shaded=True, shaded_args=(90, 45, 0.7),
                   colormap=cm.binary,
                   elevation='/home/richter/Data/map/PATCX_90m.tif',
                   elev_oceans=False,
                   stations='ipoc', spines_lw=2)
    m = map.createIPOCMap(**map_dic)
    plt.annotate('ASTER GDEM is a product of METI and NASA.', (1, 0), xycoords='axes fraction', ha='right', va='bottom', size='xx-small')
    plt.gcf().savefig('/home/richter/Documents/pics/maps/ipoc/PATCX_map_90m.pdf', dpi=600)

#map.createIPOCMap(show=False, save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_bare.pdf')
#map.createFancyIPOCMap(show=False, shaded=False, save='/home/richter/Documents/pics/maps/ipoc/ipoc_map.pdf')
#map.createFancyIPOCMap(show=False, save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_shaded.pdf')

# JGR
fw = 85 / 25.4
# dis
fwd = 140 / 25.4
#font.size : 8.0
mpl.rcParams.update({'font.size': 6, 'lines.linewidth':0.5})
#mpl.rcParams.update({'font.size': 9, 'lines.linewidth':1})
margin = (0.1, 0.075, 0.65, 0.85)
#map.createFancyIPOCMap(show=False, elevation='/home/richter/Data/map/nchile_228m.grd',
#                       save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_shaded.pdf',
#                       figsize=(fwd, 1.61 * fwd * 0.9),
#                       margin=(0.1, 0.05, 0.62, 0.9),
#                       shaded_args=(315, 45, 0.4),
#                       lw=1, spines_lw=2, dpi=200, station_markersize=6,
#                       station_labelsize=9, loffset=10000)

# for JGR
map_dic = dict(figsize=(fw, 1.61 * fw * 0.8), margin=(0.05, 0.05, 0.6, 0.9),
               lw=0.5, station_markersize=3, spines_lw=2, loffset=10000)
# for dis
map_dic = dict(figsize=(fw, 1.61 * fw * 0.8), margin=(0.05, 0.05, 0.6, 0.9),
               lw=0.5, station_markersize=3, spines_lw=1, loffset=10000)


#map.createFancyIPOCMap(show=False, save=None)

#plot_results1()
plot_results1(map_dic, scale=0.5)
#plot_salar_map()
#plot_PATCX_map()


#plot_small_map_for_F()
#plt.show()
