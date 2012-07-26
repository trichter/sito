from obspy.signal.util import utlGeoKm
import numpy as np
from obspy.core.util import gps2DistAzimuth
from sito.util.main import cosd
from sito.util.helper import vectorize_args
import sys
import gc

def migrate(stream, stations, lats, lons, velocity, skip=0, normalize=True, station_splitter='-'):
    """
    migrate reflections of autocorrs and xcorrs to surface
    
    loop over stream
    loop over lons
    loop over lats
    check (dist_st / velocity + skip) * sr < (dist1 + dist2) / velocity * sr < tr.stats.npts  
    """
    ret = np.zeros((len(lons), len(lats)))
    sr = stream[0].stats.sampling_rate
    for l, tr in enumerate(stream):
        st1, st2 = tr.stats.station.split(station_splitter)
        st1, st2 = st1[:-1], st2[:-1]
        dist_st = gps2DistAzimuth(stations[st1].latitude, stations[st1].longitude,
                           stations[st2].latitude, stations[st2].longitude)[0] / 1000.
        i0 = int(round((dist_st / velocity + skip) * sr))
        for x, lon in enumerate(lons):
            if x % 10 == 0:
                sys.stdout.write('Progress[%.2f%%]\r' % (100. * (x + 1) / len(lons) * (l + 1) / len(stream)))
                sys.stdout.flush()
            for y, lat in enumerate(lats):

                dist1 = gps2DistAzimuth(stations[st1].latitude, stations[st1].longitude, lat, lon)[0] / 1000.
                if st1 == st2:
                    dist2 = dist1
                else:
                    dist2 = gps2DistAzimuth(stations[st2].latitude, stations[st2].longitude, lat, lon)[0] / 1000
                time = (dist1 + dist2) / velocity
                i = int(round(time * sr))
                if i0 < i < tr.stats.npts:
                    to_add = tr.data[i]
                    if normalize:
                        to_add /= np.max(np.abs(tr.data[i0:]))
                    ret[x, y] += to_add
            gc.collect(2)
    return ret

def migrate2(stream, stations, lats, lons, velocity, skip=0, normalize=True, station_splitter='-'):
    """
    migrate reflections of autocorrs and xcorrs to surface
    
    loop over stream
    check ellipsis condition
    assumption: flat earth
    """
    ret = np.zeros((len(lons), len(lats)))
    lat0 = np.median(lats)
    lon0 = np.median(lons)
    ys = (lats - lat0) * 111.32 # km
    xs = (lons - lon0) * cosd(lat0) * 111.32  # km
    ys = ys[:, np.newaxis]
    sr = stream[0].stats.sampling_rate
    for l, tr in enumerate(stream):
        sys.stdout.write('Progress[%.2f%%]\r' % (100. * (l + 1) / len(stream)))
        sys.stdout.flush()
        st1, st2 = tr.stats.station.split(station_splitter)
        st1, st2 = st1[:-1], st2[:-1]
        y_st1 = (stations[st1].latitude - lat0) * 111.32
        y_st2 = (stations[st2].latitude - lat0) * 111.32
        x_st1 = (stations[st1].longitude - lon0) * cosd(lat0) * 111.32
        x_st2 = (stations[st2].longitude - lon0) * cosd(lat0) * 111.32
        x0 = (x_st1 + x_st2) / 2
        y0 = (y_st1 + y_st2) / 2
        phi = np.arctan2(y_st2 - y_st1, x_st2 - x_st1)
        xs_rot = (xs - x0) * np.cos(phi) + (ys - y0) * np.sin(phi) + x0
        ys_rot = -(xs - x0) * np.sin(phi) + (ys - y0) * np.cos(phi) + y0
        dist_st = ((x_st1 - x_st2) ** 2 + (y_st1 - y_st2) ** 2) ** 0.5
        i0 = int((dist_st / velocity + skip) * sr)
        if normalize:
            tr_data = tr.data / np.max(np.abs(tr.data[i0:]))
        else:
            tr_data = tr.data
        e = dist_st / 2 # in km
        for i, dat in enumerate(tr_data[i0:]):
            dist1 = (i + i0 - 0.5) / sr * velocity
            dist2 = (i + i0 + 0.5) / sr * velocity
            a1 = dist1 / 2 # in km
            a2 = dist2 / 2 # in km
            b1 = (a1 ** 2 - e ** 2) ** 0.5
            b2 = (a2 ** 2 - e ** 2) ** 0.5
            ind = ((((xs_rot - x0) / a1) ** 2 + ((ys_rot - y0) / b1) ** 2 >= 1) *
                   (((xs_rot - x0) / a2) ** 2 + ((ys_rot - y0) / b2) ** 2 < 1))
            ret[ind.transpose()] += dat
    return ret

@vectorize_args((2, 3))
def utlGeoKm2(orig_lon, orig_lat, lon, lat):
    return utlGeoKm(orig_lon, orig_lat, lon, lat)
def migrate3(stream, stations, lats, lons, velocity, skip=0, normalize=True, station_splitter='-'):
    """
    migrate reflections of autocorrs and xcorrs to surface
    
    loop over stream
    check ellipsis condition
    assumption: flat earth, grid calculation is ok
    """

    ret = np.zeros((len(lons), len(lats)))
    lat0 = np.median(lats)
    lon0 = np.median(lons)
    xs, ys = utlGeoKm2(lon0, lat0, lons, lats)
    ys = ys[:, np.newaxis]
    sr = stream[0].stats.sampling_rate
    for l, tr in enumerate(stream):
        sys.stdout.write('Progress[%.2f%%]\r' % (100. * (l + 1) / len(stream)))
        sys.stdout.flush()
        st1, st2 = tr.stats.station.split(station_splitter)
        st1, st2 = st1[:-1], st2[:-1]
        x_st1, y_st1 = utlGeoKm2(lon0, lat0, stations[st1].longitude, stations[st1].latitude)
        x_st2, y_st2 = utlGeoKm2(lon0, lat0, stations[st2].longitude, stations[st2].latitude)
        x0 = (x_st1 + x_st2) / 2
        y0 = (y_st1 + y_st2) / 2
        phi = np.arctan2(y_st2 - y_st1, x_st2 - x_st1)
        xs_rot = (xs - x0) * np.cos(phi) + (ys - y0) * np.sin(phi) + x0
        ys_rot = -(xs - x0) * np.sin(phi) + (ys - y0) * np.cos(phi) + y0
        dist_st = ((x_st1 - x_st2) ** 2 + (y_st1 - y_st2) ** 2) ** 0.5
        i0 = int((dist_st / velocity + skip) * sr)
        if normalize:
            tr_data = tr.data / np.max(np.abs(tr.data[i0:]))
        else:
            tr_data = tr.data
        e = dist_st / 2 # in km
        for i, dat in enumerate(tr_data[i0:]):
            dist1 = (i + i0 - 0.5) / sr * velocity
            dist2 = (i + i0 + 0.5) / sr * velocity
            a1 = dist1 / 2 # in km
            a2 = dist2 / 2 # in km
            b1 = (a1 ** 2 - e ** 2) ** 0.5
            b2 = (a2 ** 2 - e ** 2) ** 0.5
            ind = ((((xs_rot - x0) / a1) ** 2 + ((ys_rot - y0) / b1) ** 2 >= 1) *
                   (((xs_rot - x0) / a2) ** 2 + ((ys_rot - y0) / b2) ** 2 < 1))
            ret[ind.transpose()] += dat
    return ret


def _mig(lats, lons, perc, ret, stations, st1, st2, i0, tr, velocity, sr, normalize):
    for x, lon in enumerate(lons):
        if x % 10 == 0:
            sys.stdout.write('Progress[%.2f%%]\r' % (100. * (x + 1) / len(lons) * perc))
            sys.stdout.flush()
        for y, lat in enumerate(lats):

            dist1 = gps2DistAzimuth(stations[st1].latitude, stations[st1].longitude, lat, lon)[0] / 1000.
            if st1 == st2:
                dist2 = dist1
            else:
                dist2 = gps2DistAzimuth(stations[st2].latitude, stations[st2].longitude, lat, lon)[0] / 1000
            time = (dist1 + dist2) / velocity
            i = int(round(time * sr))
            if i0 < i < tr.stats.npts:
                to_add = tr.data[i]
                if normalize:
                    to_add /= np.max(np.abs(tr.data[i0:]))
                ret[x, y] += to_add


def migrate_sep(stream, stations, lats, lons, velocity, skip=0, normalize=True, station_splitter='-'):
    ret = np.zeros((len(lons), len(lats)))
    sr = stream[0].stats.sampling_rate
    for l, tr in enumerate(stream):
        st1, st2 = tr.stats.station.split(station_splitter)
        st1, st2 = st1[:-1], st2[:-1]
        dist_st = gps2DistAzimuth(stations[st1].latitude, stations[st1].longitude,
                           stations[st2].latitude, stations[st2].longitude)[0] / 1000.
        i0 = int(round((dist_st / velocity + skip) * sr))
        run_in_separate_process(_mig, lats, lons, (l + 1.) / len(stream), ret, stations, st1, st2, i0, tr, velocity, sr, normalize)
    return ret

import os, cPickle
#http://code.activestate.com/recipes/511474-run-function-in-separate-process/
def run_in_separate_process(func, *args, **kwds):
    pread, pwrite = os.pipe()
    pid = os.fork()
    if pid > 0:
        os.close(pwrite)
        with os.fdopen(pread, 'rb') as f:
            status, result = cPickle.load(f)
        os.waitpid(pid, 0)
        if status == 0:
            return result
        else:
            raise result
    else:
        os.close(pread)
        try:
            result = func(*args, **kwds)
            status = 0
        except Exception, exc:
            result = exc
            status = 1
        with os.fdopen(pwrite, 'wb') as f:
            try:
                cPickle.dump((status, result), f, cPickle.HIGHEST_PROTOCOL)
            except cPickle.PicklingError, exc:
                cPickle.dump((2, exc), f, cPickle.HIGHEST_PROTOCOL)
        os._exit(0)

#an example of use
def treble(x):
    return 3 * x

def main():
    #calling directly
    print treble(4)
    #calling in separate process
    print run_in_separate_process(treble, 4)
