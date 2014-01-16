#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR
import matplotlib as mpl
import pylab as plt
import numpy as np
from scipy.integrate import quad
import scipy.interpolate
import scipy.optimize
from math import pi, log10
from scipy.special import erf
from numpy import sqrt, exp, angle, real, absolute
from numpy import linspace
from functools import wraps
import tempfile
import os.path
from IPython import embed
import argparse

parser = argparse.ArgumentParser(description='Set velocity and mean free path '
                                 'for calculation of K3D_rt and Kz_rt.')
parser.add_argument('-c', type=float, help='velocity in m/s')
parser.add_argument('-l', type=float, help='mean free path in m')
parser.add_argument('--show', '-s', action='store_true', help='show plots')
cmd_args = parser.parse_args()
cmd_show = cmd_args.c is None and cmd_args.l is None or cmd_args.show
limit = 5000

if cmd_show:
    plt.ion()
tempfile.tempdir = '/home/richter/cache/'


def T(z, gamma, T0=1):
    return T0 * exp(-(1 + 1j) * gamma * z)

def P_diff(r, t, D):
    return 1 / (4 * pi * D * t) ** 1.5 * exp(-r ** 2 / (4 * D * t))  # factor 1.5 in exponent!!


def D(r, t, l, c):
    #factor is (4*pi*l*c)**-1.5
    if isinstance(t, (float, int)) and isinstance(r, (float, int)):
        return D_(r, t, l, c) if c * t > r else 0
    elif isinstance(t, (float, int)):
        r3 = r[c * t > r]
        t3 = t
        n = len(r) - len(r3)
        return np.hstack((D_(r, t, l, c), np.zeros(n)))
    elif isinstance(r, (float, int)):
        t3 = t[c * t > r]
        r3 = r
        n = len(t) - len(t3)
        return np.hstack((np.zeros(n), D_(r, t, l, c)))
    raise ValueError

def P_rt(r, t, l, c):  # D = c * l / 3
    return A(t, l, c) * C(l, c) * D(r, t, l, c)

def get_K3D_rt(ts, l, c, N=500, cache=True, calc_again=False):
    rfuncs = []
    #rs = np.logspace(-3, c * t * 0.6, 100)
    #rs = linspace(0.001, 500, 100)
    if cache:
        print '--- Use cache for K3D_rt (tail-tail) min value 0.01 ---'
    for t in ts:
        if cache:
            fname = os.path.join(tempfile.gettempdir(), 'K3D_rt_t=%s_l=%s_c=%s_N=%s.npy' % (t, l, c, N))
        if cache and not calc_again and os.path.exists(fname):
            res = np.load(fname)
            rs = res[:, 0]
            res = res[:, 1:]
        else:
            res = np.empty((2 * N, 2))
            rs1 = np.logspace(log10(0.01), log10(c * t * 0.1), N, endpoint=False)
            rs2 = np.linspace(c * t * 0.1, c * t * 0.6, N)
            rs = np.hstack((rs1, rs2))
            #rs = linspace(0.001, c * t * 0.6, N)
            fac0 = 2 * (exp(-c * t / l) * (4. / 3 * pi * l * c) ** -1.5 * t ** 1.5 *
                    (1 + 2.026 * l / c / t) ** -0.5)
            print 't', ts
            for i, r in enumerate(rs):
                t0 = r / c
                te = t / 2
                if te <= t0:
                    res[i, :] = 0.
                else:
                    f = lambda tp: D_(r, tp, l, c) * D_(r, t - tp, l, c)
                    res[i, :] = quad(f, t0, te)
                print 'r', r, res[i, :]
            res = fac0 * res
            if cache:
                res2 = np.hstack((rs[:, np.newaxis], res))
                np.save(fname, res2)
        rfuncs.append(scipy.interpolate.interp1d(rs, res[:, 0], kind='linear', bounds_error=False, fill_value=0.))
    def K3D(r, t):
        assert t in ts
        for i in range(len(ts)):
            if t == ts[i]:
                return rfuncs[i](r)
    K3D.rp = rs
    K3D.ts = ts
    K3D.funcs = rfuncs
    return K3D

def K3D_diff(r, t, D):
    return 1 / (2 * pi * D * r) * exp(-r ** 2 / (D * t))

def Kz_diff(z, t, D):
    return sqrt(pi * t / D) * (1 - erf(z / sqrt(D * t)))

def dvv_diff(a, const):  # a = Dtgamma**2
    return (1.5 * const / sqrt(a) * exp(-0.25j * pi) *
            (1 - exp(0.5j * a) * (1 - erf((0.5 + 0.5j) * sqrt(a)))))

def complex_quad(func, a, b, **kwargs):
    def real_func(x):
        return scipy.real(func(x))
    def imag_func(x):
        return scipy.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return np.array([real_integral[0] + 1j * imag_integral[0], real_integral[1] + 1j * imag_integral[1]])
epsilon = 0.001
def int_cyl(func, c, tp, l=None, N=500, cache=None, calc_again=False, use_rp=None):
    #print 'calc_again', calc_again
    if cache:
        dir_ = tempfile.gettempdir()
        fname = '%s_t=%s_c=%s_l=%s_N=%s.npy' % (cache, tp, c, l, N)
        fname = os.path.join(dir_, fname)
    if cache and not calc_again and os.path.exists(fname):
        res = np.load(fname)
        zs = res[:, 0]
        res = res[:, 1:]
    else:
        #zs = np.hstack((np.zeros(1), np.logspace(log10(0.01), log10(c * tp / 2), N)))
        zs = np.hstack((np.zeros(1), np.logspace(-10, -4, 7),
                        np.logspace(log10(0.001), log10(0.1 * c * tp), N, endpoint=False),
                        linspace(0.1 * c * tp, c * tp / 2, N)))
        res = np.empty((len(zs), 2))
        for i, z in enumerate(zs):
            b = (c * tp / 2) ** 2 - z ** 2
            if b <= 0:
                res[i, :] = 0
                continue
            b = sqrt(b)
            if use_rp is not None:
                rp = use_rp
                x = sqrt((rp[rp >= z]) ** 2 - z ** 2)
                x = x[x < b]
                x = np.hstack((np.zeros(1), x))
                res[i, :] = np.trapz(2 * 2 * pi * x * func(sqrt(z ** 2 + x ** 2)), x)
            else:
                #if b > 1:
                b0 = 0
                if z < epsilon:
                    rho0 = sqrt(epsilon ** 2 - z ** 2)
                res[i, :] = quad(lambda x: 2 * 2 * pi * x * func(sqrt(z ** 2 + x ** 2)), epsilon, b, limit=limit)
                                # quad(lambda x: 2 * 2 * pi * x * func(sqrt(z ** 2 + x ** 2)), 1, b, limit=limit))
                #else:
                #   res[i, :] = quad(lambda x: 2 * 2 * pi * x * func(sqrt(z ** 2 + x ** 2)), 0, b, limit=limit)
            print z, b, res[i, :]
        res[res == np.inf] = 0
        if cache:
            res2 = np.hstack((zs[:, np.newaxis], res))
            np.save(fname, res2)
    sum_ = np.trapz(res[:, 0], zs)
    func = scipy.interpolate.interp1d(zs, res[:, 0], kind='linear', bounds_error=False, fill_value=0.)
    return func, sum_

def int_cyls(func, c, ts, *args, **kwargs):
    rfuncs = []
    sums = []
    for t in ts:
        rfunc, sum_ = int_cyl(lambda x:func(x, t) if func else func, c, t, *args, **kwargs)
        rfuncs.append(rfunc)
        sums.append(sum_)
    def Kz(z, t):
        assert t in ts
        for i in range(len(ts)):
            if t == ts[i]:
                return rfuncs[i](z)
    def norm(t):
        assert t in ts
        for i in range(len(ts)):
            if t == ts[i]:
                return sums[i]
    Kz.ts = ts
    Kz.funcs = rfuncs
    Kz.norm = norm
    return Kz

def Kz_BB(z, t, l, c):
    n1 = np.count_nonzero(c * t / 2 >= z)
    n2 = len(z) - n1
    return np.hstack((np.ones(n1) * Kz_BB_(t, l, c), np.zeros(n2)))

def G(x): return exp(x) * sqrt(1 + 2.026 / x)
def A(t, l, c): return exp(-c * t / l)
def C(l, c): return (4. / 3 * pi * l * c) ** -1.5
def D_(r, t, l, c):
    r2 = 1 - (r / c / t) ** 2
    assert np.all(c * t > r)  # Heaviside(ct-r)
    return r2 ** 0.125 * t ** -1.5 * G(c * t / l * r2 ** 0.75)
def P0(t, l, c): return C(l, c) * sqrt(1 + 2.026 * l / c / t) * t ** -1.5
def Kz_BB_(t, l, c): return exp(-c * t / l) / P0(t, l, c) * 2 / pi / c ** 4 / t ** 3  # equals K3D_BB/c/t*2
@np.vectorize
def K3D_BT(r, t, l, c):
    return A(t, l, c) * C(l, c) / P0(t, l, c) / (2 * pi * r ** 2 * c) * D(r, t - r / c, l, c)
@np.vectorize
def pi4r2_K3D_BT(r, t, l, c):
    return 2 * A(t, l, c) * C(l, c) / P0(t, l, c) / c * D(r, t - r / c, l, c)

ALT = 1  #1e-5 * 10
SIGMA = 0.2
E = 1  #000
def strain_x(z, gamma, k, alT=ALT, sigma=SIGMA):
    gammac = (1 + 1j) * gamma
    fak = (2 * (1 - sigma) - k * z) * exp(-k * z) - k / gammac * exp(-gammac * z)
    return (1 + sigma) / (1 - sigma) * k / gammac * fak * alT
def strain_z(z, gamma, k, alT=ALT, sigma=SIGMA):
    gammac = (1 + 1j) * gamma
    fak = -k / gammac * (2 * sigma - k * z) * exp(-k * z) + exp(-gammac * z)
    return (1 + sigma) / (1 - sigma) * fak * alT
def strain_x2d(x, z, gamma, k, **kwargs):
    return strain_x(z[:, np.newaxis], gamma, k, **kwargs) * np.cos(k * x)
def strain_z2d(x, z, gamma, k, **kwargs):
    return strain_z(z[:, np.newaxis], gamma, k, **kwargs) * np.cos(k * x)
def stress_x(z, gamma, k, alT=ALT, sigma=SIGMA, E=E):
    args = z, gamma, k, alT, sigma
    gammac = (1 + 1j) * gamma
    fak = ((1 - sigma) * strain_x(*args) + sigma * strain_z(*args)) / (1 + sigma)
    return E / (1 - 2 * sigma) * (fak - alT * exp(-gammac * z))

def stress_y(z, gamma, k, alT=ALT, sigma=SIGMA, E=E):
    args = z, gamma, k, alT, sigma, E
    gammac = (1 + 1j) * gamma
    return sigma * (stress_x(*args) + stress_z(*args)) - E * alT * exp(-gammac * z)

def stress_z(z, gamma, k, alT=ALT, sigma=SIGMA, E=E):
    args = z, gamma, k, alT, sigma
    gammac = (1 + 1j) * gamma
    fak = ((1 - sigma) * strain_z(*args) + sigma * strain_x(*args)) / (1 + sigma)
    return E / (1 - 2 * sigma) * (fak - alT * exp(-gammac * z))
def stress_xaz(z, gamma, k, alT=ALT, sigma=SIGMA, E=E):
    gammac = (1 + 1j) * gamma
    return E / (1 - sigma) * alT * (-exp(-gammac * z)  #approximation - exp(-gammac * z)* k ** 2 / gamma ** 2 / (1 - 2 * sigma) -> 0
                                    + (1 - 1j) * k / gamma * exp(-k * z))
def stress_(z, gamma, k, alT=ALT, sigma=SIGMA, E=E):
    gammac = (1 + 1j) * gamma
    return E / (1 - sigma) * alT * (-2 * exp(-gammac * z)  #approximation - exp(-gammac * z)* k ** 2 / gamma ** 2 / (1 - 2 * sigma) -> 0
                                    + (1 + sigma) * (1 - 1j) * k / gamma * exp(-k * z))




def stress_xz(z, gamma, k, alT=ALT, sigma=SIGMA, E=E):
    gammac = (1 + 1j) * gamma
    fak = (1 - k * z) * exp(-k * z) - exp(-gamma * z)
    return 1j * E / (1 - sigma) * k / gammac * fak * alT
def stress_x2d(x, z, gamma, k, **kwargs):
    return stress_x(z[:, np.newaxis], gamma, k, **kwargs) * np.cos(k * x)
def stress_z2d(x, z, gamma, k, **kwargs):
    return stress_z(z[:, np.newaxis], gamma, k, **kwargs) * np.cos(k * x)

NN = 500

# value vs time at distance r
print 'tempdir %s' % tempfile.gettempdir()
print '--- plot P_diff and P_rt ---'
fig = plt.figure()
fig.suptitle('P_diff vs P_rt')
ax1 = fig.add_subplot(111)
ax1.set_xlabel('time (s)')
#c | 100        | 300 | 50 |
#l | 100, 50, 10 | 100 | 50 |
c = cmd_args.c or 1000.
l = cmd_args.l or 500.
r = 100.
t = linspace(0.001, 50, 1000)
ax1.plot(t, r ** 2 * P_diff(r, t, D=c * l / 3.), 'b--')
t2 = t[t > r / c]
ax1.plot(t2, r ** 2 * P_rt(r, t2, l=l, c=c), 'b')
r = 20.
ax1.plot(t, r ** 2 * P_diff(r, t, D=c * l / 3.), 'g--')
t2 = t[t > r / c]
ax1.plot(t2, r ** 2 * P_rt(r, t2, l=l, c=c), 'g')
#if cmd_show:
plt.draw()

t = 17.5
ts = (7.5, 12.5, 17.5)  #, 25., 50.)

print 'check P_rt: 1 =', quad(lambda x: 4 * pi * x ** 2 * P_rt(x, t, l, c), 0, c * t)
print 'exp(-ct/l) =', exp(-c * t / l)

print '--- integrate to get 3D RT kernel for tail-tail ---'
print 'D(t -> r/c) = ', D(r, r / c + 0.001, l, c), D(r, r / c + 0.0001, l, c)
K3D_TT = get_K3D_rt(ts, l, c, N=NN)
print '--- plot 3D kernels ---'
r = linspace(0.01, 20 * c / 4, 1000)


mpl.rcParams.update({'font.size': 7, 'lines.linewidth':1.})
# JGR
fw = 85 / 25.4
fh = fw / 1.61

fig = plt.figure(figsize=(fw, fh))
ax1 = fig.add_axes([0.14, 0.18, 0.82, 0.77])
ax1.set_xlabel(r'distance $r$ (m)')
ax1.set_ylabel(r'$4\pi r^2 K_{3\rm{d}} / t$ ($10^{-4}\rm{m}^{-1}$)')
colors = 'brgck'

for i, t in enumerate(ts):
    norm3D_TT = quad(lambda x: 4 * pi * x ** 2 * K3D_TT(x, t), 0, c * t / 2)
    norm3D_BT = quad(lambda x: 4 * pi * x ** 2 * K3D_BT(x, t, l, c), 0, c * t / 2)
    norm3D_BB = Kz_BB(np.array([0.1]), t, l, c)[0] * c * t / 2
    norm3D_RT = norm3D_TT[0] + norm3D_BT[0] + norm3D_BB
    print 't =', t
    print 'norm K3D_TT', norm3D_TT
    print 'norm K3D_BT', norm3D_BT
    print 'norm K3D_BB', norm3D_BB
    print 'norm K3D_RT', norm3D_RT
    print
    norm = 1 / norm3D_RT
    ax1.plot(r, 1e4 * norm * 4 * pi * r ** 2 * (K3D_TT(r, t) + K3D_BT(r, t, l, c)), color=colors[i], label='t=%.1fs' % t)
    ax1.plot(r, 1e4 * 4 * pi * r ** 2 * K3D_diff(r, t, D=l * c / 3) / t, ':', color=colors[i])
    ax1.plot(r, 1e4 * 4 * pi * r ** 2 * K3D_BT(r, t, l, c) / t, color=colors[i])
    ax1.plot(r, 1e4 * 4 * pi * r ** 2 * K3D_TT(r, t) / t, '--', color=colors[i])
ax1.plot((0, 0), (0, 0), color='gray', label='rt')
ax1.plot((0, 0), (0, 0), ':', color='gray', label='diff')
ax1.legend()
#fig.savefig('/home/richter/Documents/reports/paper/K3d.pdf')
#plt.close(fig)

print '--- calc Kz_TT ---'
Kz_TT = int_cyls(lambda x, t:K3D_TT(x, t), c, ts, l, cache='Kz_TT', use_rp=K3D_TT.rp, N=NN)
print '--- calc Kz_BT ---'
Kz_BT = int_cyls(lambda x, t:K3D_BT(x, t, l, c), c, ts, l, cache='Kz_BT', calc_again=False, N=NN)
Kz = lambda z, t: t * (Kz_TT(z, t) + Kz_BT(z, t)) / (Kz_TT.norm(t) + Kz_BT.norm(t))

print '--- plot Kz kernels ---'
fig = plt.figure(figsize=(fw, fh))
ax1 = fig.add_axes([0.17, 0.18, 0.79, 0.77])
ax1.set_xlabel(r'depth $z$ (m)')
ax1.set_ylabel(r'$K_z / t$ ($10^{-3}\rm{m}^{-1}$)')
z = linspace(0.1, c * 2, 1000)


for i, t in enumerate(ts):
    norm_TT = Kz_TT.norm(t)
    norm_BT = Kz_BT.norm(t)
    norm_BB = Kz_BB(np.ones(1), t, l, c) * c * t / 2
    norm_RT = norm_BB + norm_BT + norm_TT
    print 't =', t
    print 'norm Kz_diff:', quad(lambda x: Kz_diff(x, t, c * l / 3), 0, np.inf)
    #Kz_diff2 = int_cyl(lambda x:K3D_diff(x, 7.5, D=l * c / 3))
    #print 'norm diff', np.trapz(Kz_diff(z, 7.5, D=3333), z)
    print 'norm Kz_TT', norm_TT
    print 'norm Kz_BT', norm_BT
    print 'norm Kz_BB', norm_BB
    print 'norm Kz_RT', norm_RT
    print 'Kz_BB=', Kz_BB(np.array([0.1]), t, l, c)[0] / norm3D_RT * t
    print 'vs', Kz(0.49 * c * t , t)
    print 'Kz_BT*r for r < epsilon:', quad(lambda x: D(x, t - x / c, l, c), 0, epsilon)[0] * 2 * A(t, l, c) * C(l, c) / P0(t, l, c) / c
    print 'vel r < epsilon:', 100 / norm_RT / t * 1.5 * 5000 * 10e-5 * quad(lambda x: D(x, t - x / c, l, c), 0, epsilon)[0] * 2 * A(t, l, c) * C(l, c) / P0(t, l, c) / c
    norm = t / norm_RT
    ax1.plot(z, 1e3 * Kz(z, t) / t, color=colors[i], label='t=%.1fs' % t)
    ax1.plot(z, 1e3 * Kz_diff(z, t, D=l * c / 3) / t, ':', color=colors[i])
    print 'Kz(0, t)/t =', Kz(0, t) / t
ax1.plot((0, 0), (0, 0), color='gray', label='rt')
ax1.plot((0, 0), (0, 0), ':', color='gray', label='diff')
ax1.legend()
fig.savefig('/home/richter/Documents/reports/paper/Kz.pdf')
plt.close(fig)

if False:
    f = lambda x: quad(lambda x2: np.cos(2 * pi / 10000 * x2) / (1 - x2 ** 2 / x ** 2) ** 0.5, 0, x)[0] * 4 / np.pi / 2 / x
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([0.1, 0.2, 0.85, 0.75])
    z = linspace(0, 10000, 100)
    ax.plot(z, [f(x) for x in z], '--', color='gray')
    ax.legend()
    ax.annotate(r'$2 \int_0^\rho \frac{\cos kx}{\sqrt{ 1-x^2/\rho^2 }} \rm{d}x / ( \pi \rho )$', (5000, 0.4))
    ax.set_xlabel(r'$\rho$ (m)')
    fig.savefig('/home/richter/Documents/reports/paper/coskxterm.pdf')
    plt.draw()


print '--- plot temperature vs depth---'

def gamma_func(period, kappa):
    return sqrt(pi / (period * kappa))
def error_opt(kappa, T0_d, T0_a, z_obs_d, abs_obs_d, phase_obs_d, z_obs_a, abs_obs_a, phase_obs_a):
    kappa = abs(kappa)
    T0_d = abs(T0_d)
    T_d = T(z_obs_d, gamma_func(period_d, kappa), T0_d)
    err1 = np.sum((absolute(T_d) - abs_obs_d) ** 2) / abs_obs_d[0] ** 2
    err2 = np.sum(((angle(T_d) - phase_obs_d + pi) % (2 * pi) - pi) ** 2) / (2 * pi) ** 2
    if T0_a is None:
        return sqrt(err1 + err2)  # 0.96039
    T0_a = abs(T0_a)
    T_a = T(z_obs_a, gamma_func(period_a, kappa), T0_a)
    err3 = np.sum((absolute(T_a) - abs_obs_a) ** 2) / abs_obs_a[0] ** 2
    err4 = np.sum(((angle(T_a) - phase_obs_a + pi) % (2 * pi) - pi) ** 2) / abs(phase_obs_a[-1]) ** 2
    return sqrt(err1 + err2 + err3 + err4)



period_a = 365.25 * 24 * 3600
period_d = 24 * 3600

z_obs_d = np.array([0, 0.25, 0.5, 0.75])
phase_obs_d = 2 * pi * np.array([0, -0.36, -0.89, -0.21])
abs_obs_d = np.array([38, 2.5, 0.1, 0.05]) / 2
z_obs_a = np.array([])
abs_obs_a = np.array([])
phase_obs_a = 2 * pi * np.array([])
args = (None, z_obs_d, abs_obs_d, phase_obs_d, z_obs_a, abs_obs_a, phase_obs_a)
to_opt = lambda x: error_opt(x[0], x[1], *args)
res = scipy.optimize.minimize(to_opt, x0=(4e-7, 37.8), method='Nelder-Mead')
kappa, T0_d = absolute(res.x)
T0_a = None
gamma_d = gamma_func(period_d, kappa)
gamma_a = gamma_func(period_a, kappa)
print res
print 'kappa', kappa
print 'gamma day year', gamma_d, gamma_a
#gamma_d = 0.3
#gamma_a = 6.
print 'T0 day year', T0_d, T0_a
T0_a = 5.77  #11.54
lw = 2
mew = 1
ms = 5

z = np.linspace(0.001, 5, 1000)
# JGR
fw = 85 / 25.4
fh = fw / 1.61
mpl.rcParams.update({'font.size': 7, 'lines.linewidth':1.})
fig = plt.figure(figsize=(fw, fh))
ax1 = fig.add_axes([0.12, 0.2, 0.28, 0.75])
ax2 = fig.add_axes([0.45, 0.2, 0.45, 0.75], sharey=ax1)
ax1.invert_yaxis()

Td = lambda x: T(x, gamma=gamma_d, T0=T0_d)
Ta = lambda x: T(x, gamma=gamma_a, T0=T0_a)
ang_d = angle(Td(z))
#ang_d[1:] = ang_d[1:] % (2 * pi) - 2 * pi
ang_d = ang_d % (2 * pi) - 2 * pi
try:
    i = np.where(ang_d[1:] - ang_d[:-1] > 1)[0][1]
except IndexError:
    i = len(ang_d) - 1
ax1.plot(ang_d[:i + 1], z[:i + 1], 'b')
ax1.plot(angle(Ta(z)), z, 'g')
ax2.plot(abs(Td(z)), z, 'b', label='daily')
ax2.plot(abs(Ta(z)), z, 'g', label='annual')

ax1.scatter(phase_obs_d, z_obs_d, marker='x', color='r', s=30, clip_on=False)
ax2.scatter(abs_obs_d, z_obs_d, marker='x', color='r', s=30, clip_on=False)
ax2.scatter((8.8 / 2,), (0.5,), marker='x', color='r', s=30, label='data CHO2', clip_on=False)

ax1.set_ylabel(r'depth (m)')
ax1.set_xlabel('phase')
ax2.set_xlabel(r'amplitude of temp fluct $T_0 \exp (-\gamma t)$ (K)')
ax2.legend(loc='lower right', scatterpoints=1)
ax2.set_xlim([0, 20])
ax1.set_ylim(2, 0)
ax1.set_xlim(-2 * np.pi, 0)
ax1.set_xticks((-2 * np.pi, -np.pi, 0))
ax1.set_xticklabels(('$-2\pi$', '$-\pi$', ''))
[l2.set_visible(False) for l2 in ax2.get_yticklabels()]
fig.savefig('/home/richter/Documents/reports/paper/temperature.pdf')
plt.close(fig)

print '--- plot stress/strain vs depth and stress/strain field ---'
mpl.rcParams.update({'font.size': 12, 'lines.linewidth':1.})
plot_for_year = True
lam = 10000
x_split = 20 if plot_for_year else 3
k = 2 * pi / lam
z1 = np.linspace(0.01, x_split, 100)
z2 = np.linspace(x_split, lam / 2, 100)
fig = plt.figure(figsize=(15, 8))
ax1 = fig.add_axes([0.05, 0.7, 0.07, 0.25])
ax2 = fig.add_axes([0.05, 0.1, 0.07, 0.55], sharex=ax1)
ax3 = fig.add_axes([0.13, 0.7, 0.13, 0.25], sharey=ax1)
ax4 = fig.add_axes([0.13, 0.1, 0.13, 0.55], sharey=ax2)
ax5 = fig.add_axes([0.28, 0.7, 0.07, 0.25], sharey=ax1, sharex=ax1)
ax6 = fig.add_axes([0.28, 0.1, 0.07, 0.55], sharey=ax2, sharex=ax1)
ax7 = fig.add_axes([0.36, 0.7, 0.13, 0.25], sharey=ax1)
ax8 = fig.add_axes([0.36, 0.1, 0.13, 0.55], sharey=ax2)
ax9 = fig.add_axes([0.52, 0.7, 0.2, 0.25], sharey=ax1)
ax10 = fig.add_axes([0.52, 0.1, 0.2, 0.55], sharey=ax2, sharex=ax9)
ax11 = fig.add_axes([0.75, 0.7, 0.2, 0.25], sharey=ax1)
ax12 = fig.add_axes([0.75, 0.1, 0.2, 0.55], sharey=ax2, sharex=ax11)


gamma = gamma_a if plot_for_year else gamma_d
strainx1 = strain_x(z1, gamma, k)
strainz1 = strain_z(z1, gamma, k)
strainxaz1 = strainx1 + strainz1
strainx2 = strain_x(z2, gamma, k)
strainz2 = strain_z(z2, gamma, k)
strainxaz2 = strainx2 + strainz2

ax1.plot(angle(strainxaz1) % (2 * pi) - 2 * pi, z1, 'b', lw=lw)
ax1.plot(angle(strainx1) % (2 * pi) - 2 * pi, z1, 'r--', lw=lw)
ax1.plot(angle(strainz1) % (2 * pi) - 2 * pi, z1, 'c-.', lw=lw)
ax2.plot(angle(strainxaz2) % (2 * pi) - 2 * pi, z2, 'b', lw=lw)
ax2.plot(angle(strainx2) % (2 * pi) - 2 * pi, z2, 'r--', lw=lw)
ax2.plot(angle(strainz2) % (2 * pi) - 2 * pi, z2, 'c-.', lw=lw)
ax3.plot(abs(strainxaz1), z1, 'b', lw=lw)
ax3.plot(abs(strainx1), z1, 'r--', lw=lw)
ax3.plot(abs(strainz1), z1, 'c-.', lw=lw)
ax4.plot(abs(strainxaz2), z2, 'b', lw=lw, label='x+z')
ax4.plot(abs(strainx2), z2, 'r--', lw=lw, label='x')
ax4.plot(abs(strainz2), z2, 'c-.', lw=lw, label='z')

stressx1 = stress_x(z1, gamma, k)
stressy1 = stress_y(z1, gamma, k)
stressz1 = stress_z(z1, gamma, k)
stress_1 = stress_(z1, gamma, k)
#stressxz1 = stress_xz(z1, gamma, k)
stressx2 = stress_x(z2, gamma, k)
stressy2 = stress_y(z2, gamma, k)
stressz2 = stress_z(z2, gamma, k)
stress_2 = stress_(z2, gamma, k)
#stressxz2 = stress_xz(z2, gamma, k)
fak = np.max(absolute(stress_1)) / np.max(absolute(stress_2))
fak = np.round(fak / 100) * 100
fak2 = np.max(absolute(strainxaz1)) / np.max(absolute(strainxaz2))
fak2 = np.round(fak2 / 100) * 100

ax5.plot(angle(-stress_1) % (2 * pi) - 2 * pi, z1, 'b', lw=lw)
ax5.plot(angle(-stressx1) % (2 * pi) - 2 * pi, z1, 'r--', lw=lw)
ax5.plot(angle(-stressz1) % (2 * pi) - 2 * pi, z1, 'c-.', lw=lw)
#ax5.plot(angle(-stressxz1) % (2 * pi) - 2 * pi, z1, 'g:', lw=lw)
ax6.plot(angle(-stress_2) % (2 * pi) - 2 * pi, z2, 'b', lw=lw)
ax6.plot(angle(-stressx2) % (2 * pi) - 2 * pi, z2, 'r--', lw=lw)
ax6.plot(angle(-stressz2) % (2 * pi) - 2 * pi, z2, 'c-.', lw=lw)
#ax6.plot(angle(-stressxz2) % (2 * pi) - 2 * pi, z2, 'g:', lw=lw)
ax7.plot(abs(stress_1), z1, 'b', lw=lw)
ax7.plot(abs(stressx1), z1, 'r--', lw=lw)
ax7.plot(abs(stressz1), z1, 'c-.', lw=lw)
#ax7.plot(abs(stressxz1), z1, 'g:.', lw=lw)
ax8.plot(abs(stress_2), z2, 'b', lw=lw, label='x+z')
ax8.plot(abs(stressx2), z2, 'r--', lw=lw, label='x')
ax8.plot(abs(stressz2), z2, 'c-.', lw=lw, label='z')
#ax8.plot(abs(stressxz2), z2, 'g:.', lw=lw)

z1_ = np.linspace(0.01, x_split, 10)
z2_ = np.linspace(x_split, lam / 2, 10)
x_ = np.linspace(-lam / 2, lam / 2, 9)
strainx1_2d = np.real(strain_x2d(x_, z1_, gamma, k))
strainz1_2d = np.real(strain_z2d(x_, z1_, gamma, k))
strainx2_2d = np.real(strain_x2d(x_, z2_, gamma, k))
strainz2_2d = np.real(strain_z2d(x_, z2_, gamma, k))
ax9.quiver(x_, z1_, strainx1_2d, strainz1_2d, clip_on=False)
ax10.quiver(x_, z2_, strainx2_2d, strainz2_2d, clip_on=False)

stressx1_2d = np.real(stress_x2d(x_, z1_, gamma, k))
stressz1_2d = np.real(stress_z2d(x_, z1_, gamma, k))
stressx2_2d = np.real(stress_x2d(x_, z2_, gamma, k))
stressz2_2d = np.real(stress_z2d(x_, z2_, gamma, k))
ax11.quiver(x_, z1_, stressx1_2d, stressz1_2d, clip_on=False)
ax12.quiver(x_, z2_, stressx2_2d, stressz2_2d, clip_on=False)

ax1.invert_yaxis()
ax2.invert_yaxis()
ax2.set_xlim(-2 * np.pi, 0)
ax2.set_xticks((-2 * np.pi, -np.pi, 0))
ax2.set_xticklabels(('$-2\pi$', '$-\pi$', '0'))
ax2.set_yticks((x_split, 1000, 2000, 3000, 4000, 5000))

ax2.set_ylabel('depth (m)')
ax2.set_xlabel('phase')


#ax4.set_xlabel(r'abs of strain in $\alpha \mathit{T}_0$')
ax4.set_xlabel(r'amplitude of stress in $\alpha T_0 E$')
ax4.legend(loc='lower right')
ax6.set_xlabel('phase')
ax8.set_xlabel(r'abs of stress in $\alpha \mathit{T}_0 \mathit{E}$')
ax8.legend(loc='lower right')
ax9.set_title('strain field')
ax10.set_xlabel('x (m)')
ax11.set_title('stress field')
ax12.set_xlabel('x (m)')

ax2.yaxis.set_label_coords(-0.5, 0.9)


[l2.set_visible(False) for l2 in ax3.get_yticklabels() + ax4.get_yticklabels() +
 ax5.get_yticklabels() + ax6.get_yticklabels() + ax7.get_yticklabels() +
 ax8.get_yticklabels() + ax9.get_yticklabels() + ax10.get_yticklabels() +
 ax11.get_yticklabels() + ax12.get_yticklabels()]
fig.savefig('/home/richter/Documents/reports/paper/strain+stress_field.pdf')
plt.close(fig)

print '--- plot stress/strain vs depth and stress/strain field 2---'
colors = ('c', 'r', 'k', 'b')
lss = ('-', '--', ':', '-.')
mpl.rcParams.update({'font.size': 7, 'lines.linewidth':1.})
# JGR
fw = 2 * 85 / 25.4
fh = fw / 1.61 / 2
fig = plt.figure(figsize=(fw, fh))
ax1 = fig.add_axes([0.07, 0.7, 0.1, 0.25])
ax2 = fig.add_axes([0.07, 0.17, 0.1, 0.45], sharex=ax1)
ax3 = fig.add_axes([0.2, 0.7, 0.27, 0.25], sharey=ax1)
ax4 = fig.add_axes([0.2, 0.17, 0.27, 0.45], sharey=ax2, sharex=ax3)

ax5 = fig.add_axes([0.57, 0.7, 0.1, 0.25], sharey=ax1, sharex=ax1)
ax6 = fig.add_axes([0.57, 0.17, 0.1, 0.45], sharey=ax2, sharex=ax1)
ax7 = fig.add_axes([0.7, 0.7, 0.27, 0.25], sharey=ax1)
ax8 = fig.add_axes([0.7, 0.17, 0.27, 0.45], sharey=ax2, sharex=ax7)


ax1.plot(angle(-stress_1) % (2 * pi) - 2 * pi, z1, color=colors[0], ls=lss[0])
ax1.plot(angle(-stressx1) % (2 * pi) - 2 * pi, z1, color=colors[1], ls=lss[1])
ax1.plot(angle(-stressy1) % (2 * pi) - 2 * pi, z1, color=colors[2], ls=lss[2])
ax1.plot(angle(-stressz1) % (2 * pi) - 2 * pi, z1, color=colors[3], ls=lss[3])
ax2.plot(angle(-stress_2) % (2 * pi) - 2 * pi, z2, color=colors[0], ls=lss[0])
ax2.plot(angle(-stressx2) % (2 * pi) - 2 * pi, z2, color=colors[1], ls=lss[1])
#ax2.plot(angle(-stressy2) % (2 * pi) - 2 * pi, z2, color=colors[2], ls=lss[2])
#ax2.plot(angle(-stressz2) % (2 * pi) - 2 * pi, z2, color=colors[3], ls=lss[3])
ax3.plot(abs(stress_1), z1, color=colors[0], ls=lss[0])
ax3.plot(abs(stressx1), z1, color=colors[1], ls=lss[1])
ax3.plot(abs(stressy1), z1, color=colors[2], ls=lss[2])
ax3.plot(abs(stressz1), z1, color=colors[3], ls=lss[3])
ax4.plot(abs(fak * stress_2), z2, color=colors[0], ls=lss[0], label=r'$\sigma$')
ax4.plot(abs(fak * stressx2), z2, color=colors[1], ls=lss[1], label=r'$\sigma_{xx}$')
ax4.plot(abs(fak * stressy2), z2, color=colors[2], ls=lss[2], label=r'$\sigma_{yy}$')
ax4.plot(abs(fak * stressz2), z2, color=colors[3], ls=lss[3], label=r'$\sigma_{zz}$')
ax4.annotate('xscale %dx\nexaggerated' % fak, (2.41, 600), size=5, va='top', ha='right')

ax5.plot(angle(strainxaz1) % (2 * pi) - 2 * pi, z1, color=colors[0], ls=lss[0])
ax5.plot(angle(strainx1) % (2 * pi) - 2 * pi, z1, color=colors[1], ls=lss[1])
#ax1.plot(angle(-strainy1) % (2 * pi) - 2 * pi, z1, color=colors[2], ls=lss[2])
ax5.plot(angle(strainz1) % (2 * pi) - 2 * pi, z1, color=colors[3], ls=lss[3])
ax6.plot(angle(strainxaz2) % (2 * pi) - 2 * pi, z2, color=colors[0], ls=lss[0])
ax6.plot(angle(strainx2) % (2 * pi) - 2 * pi, z2, color=colors[1], ls=lss[1])
ax6.plot(angle(strainz2) % (2 * pi) - 2 * pi, z2, color=colors[3], ls=lss[3])

ax7.plot(abs(strainxaz1), z1, color=colors[0], ls=lss[0])
ax7.plot(abs(strainx1), z1, color=colors[1], ls=lss[1])
ax7.plot(abs(strainz1), z1, color=colors[3], ls=lss[3])
ax8.plot(abs(fak2 * strainxaz2), z2, color=colors[0], ls=lss[0], label=r'$\varepsilon$')
ax8.plot(abs(fak2 * strainx2), z2, color=colors[1], ls=lss[1], label=r'$\varepsilon_{xx}$')
ax8.plot(abs(fak2 * strainz2), z2, color=colors[3], ls=lss[3], label=r'$\varepsilon_{zz}$')
ax8.annotate('xscale %dx\nexaggerated' % fak2, (1.91, 600), size=5, va='top', ha='right')


ax1.invert_yaxis()
ax2.invert_yaxis()
#ax5.invert_yaxis()
#ax6.invert_yaxis()
ax2.set_yticks((x_split, 1000, 2000, 3000, 4000, 5000))
ax2.set_xlim(-2 * np.pi, 0)
ax2.set_xticks((-2 * np.pi, -np.pi, 0))
ax2.set_xticklabels(('$-2\pi$', '$-\pi$', '0'))

ax2.set_ylabel('depth (m)')
ax2.yaxis.set_label_coords(-0.5, 0.9)
ax2.set_xlabel('phase')
ax4.set_xlabel(r'amplitude of stress in $\alpha T_0 E$')
ax4.legend(loc='lower right', labelspacing=0)  #, fontsize=6, frameon=False)

ax6.set_xlabel('phase')
ax8.set_xlabel(r'amplitude of strain in $\alpha T_0$')
ax8.legend(loc='lower right', labelspacing=0)  #, fontsize=6, frameon=False)

[l2.set_visible(False) for l2 in ax3.get_yticklabels() + ax4.get_yticklabels() +
 ax1.get_xticklabels() + ax3.get_xticklabels() + ax5.get_xticklabels() +
 ax7.get_xticklabels() + ax7.get_yticklabels() + ax8.get_yticklabels()]
fig.savefig('/home/richter/Documents/reports/paper/stress+strain2.pdf')

plt.close()

print '--- plot stress vs depth ---'
colors = ('c', 'r', 'k', 'b')
lss = ('-', '--', ':', '-.')
mpl.rcParams.update({'font.size': 7, 'lines.linewidth':1.})
# JGR
fw = 85 / 25.4
fh = fw / 1.61
fig = plt.figure(figsize=(fw, fh))
ax1 = fig.add_axes([0.15, 0.7, 0.2, 0.25])
ax2 = fig.add_axes([0.15, 0.17, 0.2, 0.45], sharex=ax1)
ax3 = fig.add_axes([0.4, 0.7, 0.55, 0.25], sharey=ax1)
ax4 = fig.add_axes([0.4, 0.17, 0.55, 0.45], sharey=ax2)

ax1.plot(angle(-stress_1) % (2 * pi) - 2 * pi, z1, color=colors[0], ls=lss[0])
ax1.plot(angle(-stressx1) % (2 * pi) - 2 * pi, z1, color=colors[1], ls=lss[1])
ax1.plot(angle(-stressy1) % (2 * pi) - 2 * pi, z1, color=colors[2], ls=lss[2])
ax1.plot(angle(-stressz1) % (2 * pi) - 2 * pi, z1, color=colors[3], ls=lss[3])
ax2.plot(angle(-stress_2) % (2 * pi) - 2 * pi, z2, color=colors[0], ls=lss[0])
ax2.plot(angle(-stressx2) % (2 * pi) - 2 * pi, z2, color=colors[1], ls=lss[1])
#ax2.plot(angle(-stressy2) % (2 * pi) - 2 * pi, z2, color=colors[2], ls=lss[2])
#ax2.plot(angle(-stressz2) % (2 * pi) - 2 * pi, z2, color=colors[3], ls=lss[3])
ax3.plot(abs(stress_1), z1, color=colors[0], ls=lss[0])
ax3.plot(abs(stressx1), z1, color=colors[1], ls=lss[1])
ax3.plot(abs(stressy1), z1, color=colors[2], ls=lss[2])
ax3.plot(abs(stressz1), z1, color=colors[3], ls=lss[3])
ax4.plot(abs(fak * stress_2), z2, color=colors[0], ls=lss[0], label=r'$\sigma$')
ax4.plot(abs(fak * stressx2), z2, color=colors[1], ls=lss[1], label=r'$\sigma_{xx}$')
ax4.plot(abs(fak * stressy2), z2, color=colors[2], ls=lss[2], label=r'$\sigma_{yy}$')
ax4.plot(abs(fak * stressz2), z2, color=colors[3], ls=lss[3], label=r'$\sigma_{zz}$')
ax4.annotate('xscale %dx\nexaggerated' % fak, (2.41, 600), size=5, va='top', ha='right')

ax1.invert_yaxis()
ax2.invert_yaxis()
ax2.set_yticks((x_split, 1000, 2000, 3000, 4000, 5000))
ax2.set_xlim(-2 * np.pi, 0)
ax2.set_xticks((-2 * np.pi, -np.pi, 0))
ax2.set_xticklabels(('$-2\pi$', '$-\pi$', '0'))

ax2.set_ylabel('depth (m)')
ax2.yaxis.set_label_coords(-0.5, 0.9)
ax2.set_xlabel('phase')
ax4.set_xlabel(r'amplitude of stress in $\alpha T_0 E$')
ax4.legend(loc='lower right', labelspacing=0)  #, fontsize=6, frameon=False)

[l2.set_visible(False) for l2 in ax3.get_yticklabels() + ax4.get_yticklabels() +
 ax1.get_xticklabels() + ax3.get_xticklabels()]
fig.savefig('/home/richter/Documents/reports/paper/stress.pdf')
#plt.close(fig)
# sensitivity kernels
#for i, t in enumerate(ts):
#    ax3.plot(Kz_diff(z, t, D=c * l / 3) / t, z, ':', color=colors[i], lw=lw)
#    ax3.plot(Kz(z, t) / t, z, 'b', lw=lw, color=colors[i], label='%.1fs' % t)
    #kz_bb = Kz_BB(z, t, l, c)
    #ax2.plot(kz_bb, z, 'g', lw=5)

plt.draw()
1 / 0


print '--- final part: calculate dvv ---'
mpl.rcParams.update({'font.size': 12, 'lines.linewidth':1.5})
const = 1.5 * 1e-5 * 5000  # b * alpha * d(rho v**2)/dsigma
const2 = const * 0.8  # * (1- nu)
print 'const =', const
#c2 = min(c, 100)
c2 = c
c2 = c  #min(c2, 100)
c3 = min(c2, 100)
split = 0.06

for t in ts:
    #dvv_d1 = complex_quad(lambda z: 1 / t * const * Td(z) * Kz(z, t), 0, c3 * t / 2, limit=limit)
    #dvv_a1 = complex_quad(lambda z: 1 / t * const * Ta(z) * Kz(z, t), 0, c3 * t / 2, limit=limit)
    #dvv_ddiff = complex_quad(lambda z: 1 / t * const * Td(z) * Kz_diff(z, t, D=c * l / 3), 0, c2 * t / 2, limit=limit)
    #dvv_adiff = complex_quad(lambda z: 1 / t * const * Ta(z) * Kz_diff(z, t, D=c * l / 3), 0, c2 * t / 2, limit=limit)
    #phase_err_day = abs((angle(dvv_d1[0] + dvv_d1[1]) - angle(dvv_d1[0])) / 2 / pi * 24)
    #phase_err_year = abs((angle(dvv_a1[0] + dvv_a1[1]) - angle(dvv_a1[0])) / 2 / pi * 365.125)
    to_int = lambda z: 1 / t * const2 * Td(0) * stress_(z, gamma_d, k) * Kz(z, t)
    dvv_d2 = -(0.3 * complex_quad(to_int, 0., split) +
               complex_quad(to_int, split, c * t / 2))

    to_int = lambda z: 1 / t * const2 * Ta(0) * stress_(z, gamma_a, k) * Kz(z, t)
    dvv_a2 = -(0.3 * complex_quad(to_int, 0, split, limit=limit) +
               complex_quad(to_int, split, c * t / 2, limit=limit))
    #dvv_d3 = -complex_quad(lambda z: const2 * Td(0) * stress_(z, gamma_d, k) * exp(-z * 2 * pi / 100.) * 2 * pi / 100 , 0, c2 * t / 2, limit=limit)
    #dvv_a3 = -complex_quad(lambda z: const2 * Ta(0) * stress_(z, gamma_a, k) * exp(-z * 2 * pi / 100.) * 2 * pi / 100, 0, c2 * t / 2, limit=limit)
    #dvv_d4 = -complex_quad(lambda z: const2 * Td(0) * stress_(z, gamma_d, k) / c2 / t * 2, 0, c2 * t / 2, limit=limit)
    #dvv_a4 = -complex_quad(lambda z: const2 * Ta(0) * stress_(z, gamma_a, k) / c2 / t * 2, 0, c2 * t / 2, limit=limit)


    print 't = %.1fs' % t
    #print 'dvv day      (%.3f +- %.2f)%% - phase (%6.2f +- %5.2f)h' % (absolute(dvv_d1)[0] * 100, absolute(dvv_d1)[1] * 100, angle(dvv_d1[0]) / 2 / pi * 24, phase_err_day)
    print 'dvv day2     (%.3f +- %.2f)%% - phase (%6.2f +-     ?)h' % (absolute(dvv_d2)[0] * 100, absolute(dvv_d2)[1] * 100, angle(dvv_d2[0]) / 2 / pi * 24)
    #print 'dvv dayexp   (%.3f +- %.2f)%% - phase (%6.2f +-     ?)h' % (absolute(dvv_d3)[0] * 100, absolute(dvv_d3)[1] * 100, angle(dvv_d3[0]) / 2 / pi * 24)
    #print 'dvv dayconst (%.3f +- %.2f)%% - phase (%6.2f +-     ?)h' % (absolute(dvv_d4)[0] * 100, absolute(dvv_d4)[1] * 100, angle(dvv_d4[0]) / 2 / pi * 24)
    #print 'dvv daydiff  (%.3f +- %.2f)%% - phase (%6.2f +-     ?)h' % (absolute(dvv_ddiff)[0] * 100, absolute(dvv_ddiff)[1] * 100, angle(dvv_ddiff)[0] / 2 / pi * 24)
    #print 'dvv year     (%.3f +- %.2f)%% - phase (%6.2f +- %5.2f)d' % (absolute(dvv_a1)[0] * 100, absolute(dvv_a1)[1] * 100, angle(dvv_a1[0]) / 2 / pi * 365.125, phase_err_year)
    print 'dvv year2    (%.3f +- %.2f)%% - phase (%6.2f +-     ?)d' % (absolute(dvv_a2)[0] * 100, absolute(dvv_a2)[1] * 100, angle(dvv_a2[0]) / 2 / pi * 365.125)
    #print 'dvv yearexp  (%.3f +- %.2f)%% - phase (%6.2f +-     ?)d' % (absolute(dvv_a3)[0] * 100, absolute(dvv_a3)[1] * 100, angle(dvv_a3[0]) / 2 / pi * 365.125)
    #print 'dvv yearconst(%.3f +- %.2f)%% - phase (%6.2f +-     ?)d' % (absolute(dvv_a4)[0] * 100, absolute(dvv_a4)[1] * 100, angle(dvv_a4[0]) / 2 / pi * 365.125)
    #print 'dvv yeardiff (%.3f +- %.2f)%% - phase (%6.2f +-     ?)d' % (absolute(dvv_adiff)[0] * 100, absolute(dvv_adiff)[1] * 100, angle(dvv_adiff)[0] / 2 / pi * 365.125)
print 'KZ', Kz(1, 7.5)

if False:
    print '--- surface wave phase plot ---'

    #lam = linspace(1, 50, 9)
    fs = np.array((1, 3., 4, 5, 6, 7, 8, 9, 10, 11))
    kls = 2 * pi * fs / 100
    gammas_day = 10.37
    gammas_year = 0.542
    z = linspace(0, 100, 100)
    fig = plt.figure()
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)
    for kl in kls:
        ax1.plot(z, exp(-kl * z) * kl, label='k=%.2f' % kl)
    ax1.legend()
    dvv_d = []
    dvv_a = []
    for kl in kls:
        dvv1 = -complex_quad(lambda z: const2 * Td(0) * stress_x(z, gamma_d, k) * exp(-kl * z) * kl, 0, c3 * t / 2, limit=limit)
        dvv2 = -complex_quad(lambda z: const2 * Ta(0) * stress_x(z, gamma_a, k) * exp(-kl * z) * kl, 0, c3 * t / 2, limit=limit)
        dvv_d.append(dvv1)
        dvv_a.append(dvv2)

    dvv_d = np.array(dvv_d)[:, 0]
    dvv_a = np.array(dvv_a)[:, 0]
    ax2.plot(fs, absolute(dvv_d) * 100, 'b')
    ax3.plot(fs, absolute(dvv_a) * 100, 'g')
    ax4.plot(fs, angle(dvv_d) / 2 / pi * 24, 'b')
    ax5.plot(fs, angle(dvv_a) / 2 / pi * 365.125, 'g')

    plt.draw()

if True:
    print '--- compare different scattering parameters ---'
    #c=100, l=10, 20, ..., 150
    cs = (1000.,)
    ls = (50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000.)
    #ls = (500.,)
    fig = plt.figure()
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)
    z = np.linspace(0.01, 10000, 1000)
    t = 12.5
    dvv_d = []
    dvv_a = []
    for c in cs:
        for l in ls:
            #load kernerls
            Kz_TT = int_cyls(None, c, ts, l, cache='Kz_TT', N=NN)
            Kz_BT = int_cyls(None, c, ts, l, cache='Kz_BT', N=NN)
            Kz = lambda z, t: t * (Kz_TT(z, t) + Kz_BT(z, t)) / (Kz_TT.norm(t) + Kz_BT.norm(t))
            ax1.plot(z, Kz(z, 12.5), label='c=%s, l=%s' % (c, l))
            to_int = lambda z: 1 / t * const2 * Td(0) * stress_(z, gamma_d, k) * Kz(z, t)
            dvv_d2 = -(complex_quad(to_int, 0., split) +
                       complex_quad(to_int, split, c * t / 2))

            to_int = lambda z: 1 / t * const2 * Ta(0) * stress_(z, gamma_a, k) * Kz(z, t)
            dvv_a2 = -(complex_quad(to_int, 0, split, limit=limit) +
                       complex_quad(to_int, split, c * t / 2, limit=limit))
            dvv_d.append(dvv_d2)
            dvv_a.append(dvv_a2)

            #print
            #print l
            #dvv_d2 = -complex_quad(lambda z: 1 / t * const2 * Td(0) * stress_(z, gamma_d, k) * Kz(z, t), 0, c2 * t / 2, limit=limit)
            #dvv_a2 = -complex_quad(lambda z: 1 / t * const2 * Ta(0) * stress_(z, gamma_a, k) * Kz(z, t), 0, c2 * t / 2, limit=limit)
            #print 'dvv day2     (%.3f +- %.2f)%% - phase (%6.2f +-     ?)h' % (absolute(dvv_d2)[0] * 100, absolute(dvv_d2)[1] * 100, angle(dvv_d2[0]) / 2 / pi * 24)
            #print 'dvv year2    (%.3f +- %.2f)%% - phase (%6.2f +-     ?)d' % (absolute(dvv_a2)[0] * 100, absolute(dvv_a2)[1] * 100, angle(dvv_a2[0]) / 2 / pi * 365.125)


    dvv_d = np.array(dvv_d)[:, 0]
    dvv_a = np.array(dvv_a)[:, 0]
    ax2.plot(ls, absolute(dvv_d) * 100, 'b')
    ax3.plot(ls, absolute(dvv_a) * 100, 'g')
    ax4.plot(ls, angle(dvv_d) / 2 / pi * 24, 'b')
    ax5.plot(ls, angle(dvv_a) / 2 / pi * 365.125, 'g')

    fig = plt.figure()
    ax1 = fig.add_axes([0.18, 0.18, 0.77, 0.65])
    ax1.scatter(ls, absolute(dvv_d) * 100, label='daily_mod')
    ax1.scatter(ls, absolute(dvv_a) * 100, color='r', label='annual_mod')

    obs_f = [3, 4, 5, 6, 7, 8]
    obs_a = [0.091359370202681825, 0.16250294823409092, 0.18022184335190908, 0.27550038194609094, 0.352675342078, 0.41218299789250001]
    obs_eq = [0.50494604630454543, 0.60439601703281809, 0.66660273773063639, 0.6951867728108182, 0.73710964999509099, 0.80957191190550004]

    ax2 = ax1.twiny()
    ax2.scatter(obs_f, obs_a, marker='d', color='r', label='annual_obs')
    ax2.scatter(obs_f, obs_eq, marker='d', color='c', label='eq_obs')

    ax1.set_xlabel('l (m)')
    ax1.set_ylabel('dv/v (%)')
    ax2.set_xlabel('f (Hz)')
    ax2.set_xlim((9, 2))
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2)
    fig.savefig('/home/richter/Documents/reports/paper/dv_l.pdf')

if cmd_show:
    plt.ioff()
    embed()

