#!/usr/bin/python
# by TR

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from obspy.imaging.beachball import Beach
from matplotlib import patches, collections, colors, colorbar
from osgeo import gdal
import scipy as sp
import scipy.interpolate
import os.path

def createColormapFromGPF(file_):
    data = sp.loadtxt(file_)
    cdict = {'red': np.take(data, (0, 1, 1), axis=1),
             'green': np.take(data, (0, 2, 2), axis=1),
             'blue': np.take(data, (0, 3, 3), axis=1)}
    name = os.path.splitext(os.path.basename(file_))[0]
    return colors.LinearSegmentedColormap(name, cdict)

def plotTrench(map_, lons, lats, size=10000, sep=50000, side=1, **kwargs):
    x, y = map_(lons, lats)
    N = len(x)
    N_new = int(((y[-1] - y[0]) ** 2 + (x[-1] - x[0]) ** 2) ** 0.5 / size * 10)
    param = np.linspace(0, 1, N)
    param_new = np.linspace(0, 1, N_new)
    # interpolate with splines   
    #s=None for perfect fit, 4.2
    tck = (scipy.interpolate.splprep([x, y], u=param, s=N ** 4.3))[0]
    x, y = scipy.interpolate.splev(param_new, tck)
    # interpolate with polyfit        
    #polyx = np.polyfit(param, x, 30)
    #polyy = np.polyfit(param, y, 30)
    #x = np.polyval(polyx, param_new)
    #y = np.polyval(polyy, param_new)                    
    map_.plot(x, y, 'k')
    polys = []
    i2 = 0
    while True:
        try:
            i1 = np.nonzero(np.logical_and((x - x[i2]) ** 2 + (y - y[i2]) ** 2 >
                                           sep ** 2, np.arange(N_new) > i2)
                            )[0][0]
            i2 = np.nonzero(np.logical_and((x - x[i1]) ** 2 + (y - y[i1]) ** 2 >
                                           size ** 2, np.arange(N_new) > i1)
                            )[0][0]
        except IndexError:
            break
        angle = np.arctan2(y[i2] - y[i1], x[i2] - x[i1]) + side * np.pi / 3.
        x2 = x[i1] + size * np.cos(angle)
        y2 = y[i1] + size * np.sin(angle)
        xy = np.array(((x[i1], y[i1]), (x[i2], y[i2]), (x2, y2)))
        polys.append(patches.Polygon(xy))
    col = collections.PatchCollection(polys, match_original=False, zorder=20,
                                      **kwargs)
    plt.gca().add_collection(col)


def createFancyIPOCMap(**kwargs_in):
    kwargs = dict(figsize=(19 / 2.54, 27 / 2.54),
                  margin=(0.1, 0.075, 0.65, 0.85),
                  cities='ipoc', trench='ipoc',
                  earthquake='Tocopilla',
                  #elevation='/home/richter/Data/map/SRTM_nchile_250m.tif',
                  elevation='/home/richter/Data/map/nchile.grd',
                  #elevation='/home/richter/Data/map/CleanTopo2_nchile.tif',
                  #colormap='/home/richter/Data/cm/wiki-2.0.gpf',
                  colormap='/home/richter/Data/cm/bath_111_tpglarm.gpf',
                  #colormap=cm.GMT_globe
                  )
    kwargs.update(kwargs_in)
    return createIPOCMap(**kwargs)

def createIPOCMap(**kwargs_in):
    dict_basemap_ipoc = dict(lat_0= -21, lon_0= -69.5,
                             projection='stere', resolution='h')
    kwargs = dict(ll=(-25, -72), ur=(-17, -68),
                  figsize=(15 / 2.54, 18 * 1.5 / 2.54),
                  margin=(0.15, 0.075, 0.7, 0.85),
                  dict_basemap=dict_basemap_ipoc,
                  grid=True,
                  watercolor=None,
                  fillcontinents=None,
                  stations='ipoc')
    kwargs.update(kwargs_in)
    if kwargs['stations'] == 'ipoc':
        from sito.stations import IPOCStations
        stations = IPOCStations()
        stations['PB13/PB16'] = stations['PB13']
        del stations['PB13']
        del stations['PB16']
        kwargs['stations'] = stations
    return createMap(**kwargs)

def createRFEventMap(circles=(27, 93), circles_around=None, **kwargs_in):
    #http://matplotlib.sourceforge.net/basemap/doc/html/api/basemap_api.html
    w = np.pi * 6371 * 1000 * 186 / 180.
    dict_basemap_ipoc = dict(lat_0= -21, lon_0= -69.5,
                             projection='aeqd', resolution='c',
                             width=w,
                             height=w)
    kwargs = dict(figsize=(15 / 2.54, 15 / 2.54),
                  margin=(0.01, 0.01, 0.98, 0.98),
                  dict_basemap=dict_basemap_ipoc,
                  grid=10, grid_labels=False,
                  mapboundary=False,
                  countries=False,
                  coastlines=False,
                  lw=0.5,
                  use_lsmask=True,
                  watercolor='blue',
                  fillcontinents='lightgrey',
                  stations=None
                  #ll=(-21 - 50, -70 - 50), ur=(-21 + 50, -70 + 50)
                  )
    kwargs.update(kwargs_in)
    m = createMap(**kwargs)
    if circles or circles_around:
        from sito import util
    if circles:
        for radius in circles:
            x, y = util.imaging.equi(m, -21, -69.5, radius)
            plt.plot(x, y, 'k')
    if circles_around:
    #tuple of circles(tuple len 5) (lat, lon, radius, indeg, color/line style)
        for circle in circles_around:
            X, Y = util.imaging.equi(m, *circle[:4])
            plt.plot(X, Y, circle[4])
#    if lines != None:
#        for azi in lines:
#            X, Y = util.imaging.line(m, -21, -70, azi, 0, 180)
#            plt.plot(X, Y, 'k')
    return m


def createMap(ll=None, ur=None, figsize=None, margin=None,
               ax=None,
               dict_basemap=None,
               countries=True, lw=1, coastlines=True,
               mapboundary=True, grid=False, grid_labels=True,
               watercolor='#ADD8E6', fillcontinents='coral', use_lsmask=False,
               stations=None, cities=None, trench=None, earthquake=None,
               elevation=None, elevation_args=(7000, 500, True),
               shaded=True, shaded_args=(315, 45, 0.2),
               colormap=None,
               title=None,
               show=True, save=None):
    if figsize:
        fig = plt.figure(figsize=figsize)
        if margin:
            ax = fig.add_axes(margin)
    else:
        fig = None
    if dict_basemap is None:
        dict_basemap = dict(llcrnrlon= -180, llcrnrlat= -60, urcrnrlon=180,
                            urcrnrlat=60,
                            lat_0=0, lon_0=0,
                            projection='hammer', resolution='l')
    elif ll is not None:
        print dict_basemap
        up_dict = dict(llcrnrlon=ll[1], llcrnrlat=ll[0],
                       urcrnrlon=ur[1], urcrnrlat=ur[0])
        dict_basemap.update(up_dict)
    m = Basemap(ax=ax, **dict_basemap)
    if fillcontinents and use_lsmask:
        m.drawlsmask(land_color=fillcontinents, ocean_color='white')
    elif fillcontinents:
        m.fillcontinents(color=fillcontinents, lake_color=watercolor, zorder=0)
    if countries:
        m.drawcountries(linewidth=lw)
    if coastlines:
        m.drawcoastlines(linewidth=lw)
    if mapboundary:
        m.drawmapboundary(fill_color=watercolor, zorder= -10)
    if grid:
        if grid is True:
            grid = 1.
        if grid_labels:
            m.drawparallels(np.arange(-90., 90., grid),
                            labels=[True, True, False, False], linewidth=lw)
            m.drawmeridians(np.arange(0., 390., grid),
                            labels=[False, False, True, True], linewidth=lw)
        else:
            m.drawparallels(np.arange(-90., 90., grid),
                            labels=[False, False, False, False], linewidth=lw)
            m.drawmeridians(np.arange(0., 390., grid),
                            labels=[False, False, False, False], linewidth=lw)
    if stations:
        stations.plot(m, mfc='b', ms=5, zorder=10)
    if earthquake == 'Tocopilla':
        x, y = m(-70.06, -22.34)
        x2, y2 = m(-69.5, -24.2) #x2, y2 = m(-70.4, -21.5)
        m.plot((x, x2), (y, y2), 'k-', lw=1) #line
        m.plot((x), (y), 'r*', ms=15, zorder=10) #epicenter
        b = Beach([358, 26, 109], xy=(x2, y2), width=50000, linewidth=1)
        b.set_zorder(10)
        plt.gca().add_collection(b)
        plt.annotate('14 Nov. 2007\n M 7.6', (x2, y2), xytext=(22, 0),
                     textcoords='offset points', va='center')
    elif earthquake == 'Tocopilla_position':
        x, y = m(-70.06, -22.34)
        m.plot((x), (y), 'r*', ms=15, zorder=10) #epicenter
    elif earthquake:
        xs, ys = zip(*[m(lon, lat) for lon, lat in zip(*earthquake[:2])])
        m.plot(xs, ys, 'r*', ms=15, zorder=10) #epicenters
        if len(earthquake) == 3:
            for i in range(len(xs)):
                x, y, mag = xs[i], ys[i], earthquake[2][i]
                plt.annotate('M%.1f' % mag, xy=(x, y), xytext=(-5, -15),
                           textcoords='offset points', va='center')
        elif len(earthquake) > 3:
            for i in range(len(xs)):
                x, y, mag, date = xs[i], ys[i], earthquake[2][i], earthquake[3][i]
                plt.annotate('M%.1f\n%s' % (mag, date.date), xy=(x, y), xytext=(10, 0),
                           textcoords='offset points', va='center')
    if elevation:
        vmax, resol, bar = elevation_args
        geo = gdal.Open(elevation)
        topoin = geo.ReadAsArray()
        coords = geo.GetGeoTransform()
        nlons = topoin.shape[1]; nlats = topoin.shape[0]
        delon = coords[1]
        delat = coords[5]
        lons = coords[0] + delon * np.arange(nlons)
        lats = coords[3] + delat * np.arange(nlats)[::-1] # reverse lats
        # create masked array, reversing data in latitude direction
        # (so that data is oriented in increasing latitude,
        # as transform_scalar requires).
        topoin = np.ma.masked_less(topoin[::-1, :], -11000)
        #topoin = array[::-1, :]
        # transform DEM data to a 250m native projection grid
        nx = int((m.xmax - m.xmin) / resol) + 1
        ny = int((m.ymax - m.ymin) / resol) + 1
        topodat = m.transform_scalar(topoin, lons, lats, nx, ny, masked=True)
        if elevation == 'CleanTopo2_nchile.tif':
            # -10,701m(0) to 8,248m(18948)
            topodat = topodat - 10701
        if (isinstance(colormap, basestring) and os.path.isfile(colormap) and
            colormap.endswith('.gpf')):
            colormap = createColormapFromGPF(colormap)
        if shaded:
            azdeg, altdeg, fraction = shaded_args
            ls = colors.LightSource(azdeg=azdeg, altdeg=altdeg)
            # convert data to rgb array including shading from light source.
            rgb = colormap(0.5 * topodat / vmax + 0.5)
            rgb1 = ls.shade_rgb(rgb, elevation=topodat, fraction=fraction)
            rgb[:, :, 0:3] = rgb1
            m.imshow(rgb, zorder=1)
        else:
            m.imshow(topodat, cmap=colormap, zorder=1, vmin= -vmax, vmax=vmax)
        if bar:
            ax = plt.gca()
            fig = plt.gcf()
            cb_ax = fig.add_axes([0.8, 0.3, 0.02, 0.4])
            cb = colorbar.ColorbarBase(cb_ax, cmap=colormap,
                                       norm=colors.Normalize(vmin= -vmax,
                                                             vmax=vmax))
            cb.set_label('elevation and bathymetry')
            ticks = (np.arange(20) - 10) * 1000
            ticks = ticks[np.abs(ticks) <= vmax]
            cb.set_ticks(ticks)
            tls = ['%dm' % i if i % 2000 == 0 else '' for i in ticks]
            cb.set_ticklabels(tls)
            plt.axes(ax)  # make the original axes current again
    if cities == 'ipoc':
        cities = 'Antofagasta Tocopilla Iquique Calama Pisagua Arica'.split()
        lons = [-70.4, -70.2, -70.1525, -68.933333, -70.216667, -70.333333]
        lats = [-23.65, -22.096389, -20.213889, -22.466667, -19.6, -18.483333]
        x, y = m(lons, lats)
        m.plot(x[:len(cities)], y[:len(cities)], 'o', ms=5, mfc='w', mec='k',
               zorder=10)
        for i in range(len(cities)):
            if 'Anto' in cities[i]:
                plt.annotate(cities[i], (x[i], y[i]), xytext=(5, 0),
                                     textcoords='offset points', va='top')
            else:
                plt.annotate(cities[i], (x[i], y[i]), xytext=(-5, 0),
                                     textcoords='offset points', ha='right')
    if trench == 'ipoc':
        coords = np.transpose(np.loadtxt('/home/richter/Data/map/'
                                         'nchile_trench.txt'))
        plotTrench(m, coords[0, :], coords[1, :], facecolors='k',
                   edgecolors='k')
    if title:
        plt.title(title, y=1.05)
    if save:
        plt.savefig(save)
    if show:
        plt.show()
    return m



def plotEvents(catalog, client=None, **kwargs):
    """
    Plot events on map.
    
    When client is given it shows the waveform when the event is clicked.
    IPython may not be started in pylab mode for this feature working.
    """
    def onpress(event):
        if not event.inaxes or event.key != 'shift':
            return
        lon, lat = map(event.xdata, event.ydata, inverse=True)
        dlon = 0.05
        dlat = 0.05
        filter = ('%f < latitude < %f and %f < longitude < %f' %
                   (lat - dlat, lat + dlat, lon - dlon, lon + dlon))
        ev = catalog.filter2(filter)
        if len(ev) == 0:
            print 'No event picked'
            return
        from sito import Stream
        st = Stream()
        print 'Selcet', ev
        for arrival in ev[0].origins[0].arrivals:
            arrival.pick_id.convertIDToQuakeMLURI()
            pick = arrival.pick_id.getReferredObject()
            if not pick:
                print 'FAIL'
                return
            time = pick.time
            seed_id = pick.waveform_id.getSEEDString()
            try:
                st1 = Stream(client.getWaveform(*(seed_id.split('.') + [time - 50, time + 250])))
            except Exception as ex:
                print '%s for %s' % (ex, seed_id)
                continue

            st1.merge()
            print 'load %s %s %.1f' % (seed_id, pick.phase_hint, time - ev[0].origins[0].time)
            st1[0].stats['label'] = '%s %s %.1f' % (seed_id, pick.phase_hint, time - ev[0].origins[0].time)
            st += st1
        st.setHI('filter', '')
        st.filter2(2, None)
        #st.plot(automerge=False, method='fast', type='relative')
        im = st.plot_()
        plt.show()
    if client is None:
        catalog.plot(**kwargs)
    else:
        map, ax = catalog.plot(handle=True, **kwargs)
        plt.connect('button_press_event', onpress)
    plt.show()
