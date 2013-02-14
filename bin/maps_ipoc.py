#!/usr/bin/env python
# by TR

from sito import map
import matplotlib.pyplot as plt

#map.createIPOCMap(show=False, save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_bare.pdf')
#map.createFancyIPOCMap(show=False, shaded=False, save='/home/richter/Documents/pics/maps/ipoc/ipoc_map.pdf')
#map.createFancyIPOCMap(show=False, save='/home/richter/Documents/pics/maps/ipoc/ipoc_map_shaded.pdf')


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

def plot_results1():
    from sito.stations import IPOCStations
    from matplotlib.lines import Line2D

    m = map.createIPOCMap(show=False, ll=(-25, -71), ur=(-18, -68.5),
                          earthquake='Tocopilla_position', cities='ipoc',
                          stations=None)
    ipoc = IPOCStations()
    ipoc_no_data = ipoc.pick('PB09 PB10 PB11 PB12 PB13 PB14 PB15', replace=False)
    #1-3Hz
    ipoc_bad_data = ipoc.pick('LVC MNMCX ', replace=False)
    ipoc_no_signal = ipoc.pick('PB06 PB08 PSGCX', replace=False)  #<0.2%
    ipoc_small_signal = ipoc.pick('HMBCX PATCX PB01 PB02 PB03 PB07', replace=False)  #0.2-0.5%
    ipoc_big_signal = ipoc.pick('PB04 PB05', replace=False)  #>0.5%

    ipoc_no_data.plot(m, mfc='w', ms=8, zorder=20)
    ipoc_bad_data.plot(m, mfc='#00FFFF', ms=8, zorder=20)
    ipoc_no_signal.plot(m, mfc='y', ms=8, zorder=20)
    ipoc_small_signal.plot(m, mfc='orange', ms=8, zorder=20)
    ipoc_big_signal.plot(m, mfc='r', ms=8, zorder=20)

    #4-6Hz
    for st in ipoc:
        ipoc[st].longitude += 0.15
    ipoc_bad_data = ipoc.pick('PB04', replace=False)
    ipoc_no_signal = ipoc.pick('HMBCX MNMCX PB01 PB06 PSGCX LVC', replace=False)  #<0.2%
    ipoc_small_signal = ipoc.pick('PB02 PB03 PB05 PB07 PB08', replace=False)  #0.2-0.5%
    ipoc_big_signal = ipoc.pick('PATCX', replace=False)  #>0.5%

    ipoc_bad_data.plot(m, marker='s', mfc='#00FFFF', ms=8, zorder=20, annotate=False)
    ipoc_no_signal.plot(m, marker='s', mfc='y', ms=8, zorder=20, annotate=False)
    ipoc_small_signal.plot(m, marker='s', mfc='orange', ms=8, zorder=20, annotate=False)
    ipoc_big_signal.plot(m, marker='s', mfc='r', ms=8, zorder=20, annotate=False)

    colors = 'w #00FFFF y orange r'.split()
    labels = 'no data,bad data,<0.2%,0.2-0.5%,>0.5%,1Hz-3Hz,4Hz-6Hz'.split(',')
    plt.gca().legend([Line2D((0,), (1,), marker='o', ls='', mfc=c, ms=8) for c in colors] +
              [Line2D((0,), (1,), marker=m, ls='', mfc='w', ms=8) for m in 'os'],
              labels, numpoints=1, loc='lower right')
    plt.gcf().savefig('/home/richter/Documents/pics/maps/ipoc/ipoc_map_results.png')

def plot_small_map_for_F():
    from sito.stations import IPOCStations
    from matplotlib.lines import Line2D

    ipoc = IPOCStations()
    ipoc.pick('PATCX')
    m = map.createIPOCMap(show=False, ll=(-24, -70.9), ur=(-18, -68.5),
                          earthquake='Tocopilla_position', stations=ipoc,
                          cities='Tocopilla')



#plot_results1()
plot_small_map_for_F()
plt.show()
