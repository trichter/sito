# by TR

#import logging
from helper import vectorize_args, add_doc
from main import arcsind, r_earth
import numpy as np

try:
    import seis.geo #@UnusedImport
    import seis.ttt
    # functions for use in util namespace
    from seis.geo import xandaz2coord as dist2gps
    from seis.fe import region as feregion#@UnusedImport
    #import seis.model as model
    #iasp91 = model.Iasp91(5., 1000.)
    from seis._rf import mocorr, zmigr #@UnusedImport
except ImportError:
    print ('WARNING: Seispy not found!')
    from _pspier import pspier as pspier  # adapted Fortran code from XY @UnresolvedImport @UnusedImport
    import _psmout, imaging
    def psmout(stream):
        xpos = imaging.getDataWindow(stream, None, None)
        slow = stream.getHI('slowness')
        st = stream[0].stats
        ypos = _psmout.psmout(xpos, slow, #@UndefinedVariable
                              st.ponset - st.starttime, st.endtime -
                              st.starttime, 1. / st.sampling_rate, 0)
        print ypos
        # TO DO: write y data to stream
    import warnings
    warnings.warn('Unable to import seis module. ')
    feregion = None

#phase_dict = {'S':-2, 'P':-1, 'Ps':0, 'Sp':0, 'Ppps':1, 'Pps':1, 'Ppss':2, 'Pss':2, 'Psss':3,
#    'Ppp':4, 'Pppp':4, 'Spp':4, 'Sspp':4, 'Ssss':5}
phase_dict = {'S':-2, 'P':-1, 'Ps':0, 'Sp':0, 'Ppps':1, 'Sssp':1, 'Ppss':2, 'Psss':3, #, 'Ssss':2
    'Pppp':4, 'Sspp':4, 'Sppp':5, 'Ssss':6}

try:
    from seis.model import EarthModel1D as EarthModel1D_Seis
except ImportError:
    pass
    #EarthModel1D_Seis = object

class EarthModel1D(EarthModel1D_Seis):
    def trace3(self, p=6.4, phase='P', max_z=None, till_turn=False):
        """
        Traces one ray from station downwards.        
        """
        if phase in 'pP':
            v = self.vp
        else:
            v = self.vs
        dz = self.z[1] - self.z[0]
        sini = p * v / (r_earth - self.z) * 180 / np.pi
#        from pylab import plot, show
#        plot(sini)
#        show()

        findN = np.nonzero(np.logical_or(sini >= 1, v == 0.))
        if len(findN[0]) == 0:
            N = len(self.z)
        else:
            N = findN[0][0] - 1
        if max_z is not None:
            N = min(N, self.layer(max_z))
        sini = sini[:N]
        v = v[:N]
        r = r_earth - self.z[:N]
        dx = dz * sini / (1 - sini ** 2) ** 0.5
        phi = np.hstack((0, np.cumsum(dx / r)))[:-1]
        dt = dz / (1 - sini ** 2) ** 0.5 / v
        t = np.hstack((0, np.cumsum(dt)))[:-1]
        if not till_turn:
            phi = np.hstack((phi, phi[-1] + np.cumsum(dx[::-1] / r[::-1])))
            t = np.hstack((t, t[-1] + np.cumsum(dt[::-1])))
            r = np.hstack((r, r[::-1]))
        return t, r, phi

    @vectorize_args((1, 2, 3, 4, 5))
    def pspier(self, depth, slat, slon, slowness, azi, phase='P'):
        """"
        Calculate piercing points of events.
    
        :param depth: depth of pp
        :param slat: station latitude
        :param slat: station longitude
        :param slowness: slowness of wave
        :param azi: azimuth of wave
        :param model: use this model (instance of seis.model.EarthModel1D)
        :return: (distance, latitude of pp, longitude of pp"""
        phi = self.trace3(slowness, phase=phase, till_turn=True)[2] #@UnusedVariable
        x = phi[self.layer(depth)] * r_earth
        return (x,) + dist2gps(x, azi, slat, slon)



def Iasp91(dz=5., zmax=1000.):
    from seis._model import iasp91_arr
    z = dz * np.arange(zmax / dz + 0.5)
    z, vp, vs, rh = iasp91_arr(z)
    mod = EarthModel1D(z=z, vp=vp, vs=vs, rh=rh)
    return mod

def TwoLayer(zl=20., vp1=5.8, vp2=8.0, vs1=3.4, vs2=4.5, dz=1, zmax=200):
    if dz is not None:
        z = 1.*np.arange(int(zmax / dz)) * dz
        N = len(z)
        mod = EarthModel1D(z=z, vp=np.empty(N), vs=np.empty(N), n=np.ones(N))
        mod['vp'] = np.where(z < zl, vp1, vp2)
        mod['vs'] = np.where(z < zl, vs1, vs2)
    else:
        mod = EarthModel1D(z=np.array((0., 1.*zl)), vp=np.array((vp1, vp2)),
                           vs=np.array((vs1, vs2)), n=np.array((1, 0)))
    return mod



iasp91 = Iasp91(5., 1000.)

#@vectorize_args((0,1,2,3,4))
@np.vectorize
def pspier(depth, slat, slon, slowness, azi, model=iasp91):
    """"
    Calculate piercing points of events.

    :param depth: depth of pp
    :param slat: station latitude
    :param slat: station longitude
    :param slowness: slowness of wave
    :param azi: azimuth of wave
    :param model: use this model (instance of seis.model.EarthModel1D)
    :return: (distance, latitude of pp, longitude of pp"""
    tp, ts, xp, xs = model.trace(slowness) #@UnusedVariable
    xs = xs[model.layer(depth)]
    return (xs,) + dist2gps(xs, azi, slat, slon)

def time2depth(t, phase='Ps', slowness=6.4, model=iasp91):
    """
    Calculate depth of conversion for a time in receiver function.

    :param t: time after p-onset for 'P+' phases or
        time after s-onset for 'S+' phases or
        time from depth for 'P' or 'S'
    :param phase: phase from phasedict or
        P*', 'S*'', '*'
    :param slowness: slowness (6.4 deg/s)
    :param model: used velocity model (instance of seis.model.EarthModel1D)
    """
    tp, ts, xp, xs = model.trace(slowness) #@UnusedVariable
    #           Ps     Ppps  Pppss Psss     Pppp    Sppp  Ssss   S   P
    t_list = [ts - tp, tp + ts, 2 * ts, 3 * ts - tp, 2 * tp, 3 * tp - ts, 2 * ts, ts, tp]
    if phase in phase_dict:
        return interpolate(t, t_list[phase_dict[phase]], model.z)
    elif isinstance(phase, tuple or list):
        phase_list = phase
    elif phase[0] in 'PS':
        phase_list = [ph for ph in phase_dict.keys() if ph.startswith(phase[0])]
    else:
        phase_list = phase_dict.keys()
    return dict([(ph, interpolate(t, t_list[phase_dict[ph]], model.z)) for ph in phase_list])

def depth2time(z, phase='Ps', slowness=6.4, model=iasp91):
    """
    Calculate time of conversion in receiver function for specified depth.

    :param z: depth of conversion
    :param phase: phase from phasedict or 'P*', 'S*'', '*'
    :param slowness: slowness (6.4 deg/s)
    :param model: used velocity model (instance of seis.model.EarthModel1D)
    """
    tp, ts, xp, xs = model.trace(slowness) #@UnusedVariable
    tp = interpolate(z, model.z, tp)
    ts = interpolate(z, model.z, ts)
    #           Ps     Ppps  Ppss Psss     Pppp    Sppp, Ssss   S   P
    t_list = [ts - tp, tp + ts, 2 * ts, 3 * ts - tp, 2 * tp, 3 * tp - ts, 2 * ts, ts, tp]
    if phase in phase_dict:
        return t_list[phase_dict[phase]]
    elif isinstance(phase, tuple or list):
        phase_list = phase
    elif phase[0] in 'PS':
        phase_list = [ph for ph in phase_dict.keys() if ph.startswith(phase[0])]
    else:
        phase_list = phase_dict.keys()
    return dict([(ph, t_list[phase_dict[ph]]) for ph in phase_list])

@vectorize_args((0,))
def multiples(mm, mmpers, phase_list=('Ps', 'Pppp', 'Ppps')):
    """
    Return multiples for time of converted S-phase after P-onset.

    :param mm: time in mm
    :param mmpers: for conversion or param mm to sec
    :param phase_list: names of multiplies to be returned
    :return: dictonary of depth, times and mms
    """
    time0 = mm / mmpers
    depth = time2depth(time0)
    times = depth2time(depth, phase=phase_list)
    mms = dict([(ph, times[ph] * mmpers) for ph in phase_list])
    return {'depth': depth, 'times': times, 'mms': mms}

@vectorize_args((0,))
def interpolate(t, t_arr, z_arr):
    """
    Interpolate depth with given time.

    :param t: time
    :param t_arr: time array
    :param z_arr: depth array
    :return: depth
    """
    i1 = np.flatnonzero(t_arr <= t)
    if len(i1) > 0:
        i1 = i1[-1]
        if i1 < len(t_arr) - 1:
            i2 = i1 + 1
        else:
            i2 = i1
    else:
        i1 = 0
        i2 = 0
    if i1 == i2:
        z = z_arr[i1]
    else:
        z = z_arr[i1] + (z_arr[i2] - z_arr[i1]) / (t_arr[i2] - t_arr[i1]) * (t - t_arr[i1])
    return z

try:
    @add_doc(seis.ttt.compute)
    def ttt(dist, depth, surface=False):
        """
        Return arrival list.

        Doc of seis.ttt.compute:"""
        arrlist = seis.ttt.compute(dist, depth)
        for arr in arrlist:
            if arr.phase.find('S') != -1 and arr.phase.rfind('S') > arr.phase.rfind('P'):
                velocity = 3.36
            else:
                velocity = 5.8
            arr.inci = arcsind(arr.slow / r_earth * velocity * 180 / np.pi)
        if surface:
            arr = seis.ttt.Arrival(phase='Surf', time=2 * np.pi * r_earth * dist / 360 / 4.4)
            arrlist.append(arr)
        return arrlist
except:
    ttt = None

def correct_ttt(phase='P', slowness=6.4, newmodel=iasp91, oldmodel=iasp91):
    """
    Correct travel times of taup for other (crust) model
    """
    #if oldmodel == newmodel or phase not in 'PS':
    #    return 0.
    tpo, tso, xp, xs = oldmodel.trace(slowness) #@UnusedVariable
    tpn, tsn, xp, xs = newmodel.trace(slowness) #@UnusedVariable
    z = min(oldmodel.z[-1], newmodel.z[-1])
    if phase == 'P':
        return interpolate(z, newmodel.z, tpn) - interpolate(z, oldmodel.z, tpo)
        #tpn[newmodel.layer(z)] - tpo[oldmodel.layer(z)]
    else:
        return interpolate(z, newmodel.z, tsn) - interpolate(z, oldmodel.z, tso)
