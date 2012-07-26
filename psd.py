from obspy.signal import cosTaper
from obspy.signal.util import nextpow2
import numpy as np
import pylab as plt


#def downsample(data, smoothie):
#    newdata = data[:N-smoothie+1:smoothie]
#    for i in range(1,smoothie):
#        newdata += data[i:N+i-smoothie+1:smoothie]
#    return newdata/smoothie

def oct_downsample(data, df, fac=2., startf=None, endf=None):
    if startf == None:
        startf = df
    if endf == None:
        endf = df * len(data)
    f1 = startf
    j = 0
    newdata = np.empty(len(data) / 10)
    while True:
        ind1 = int(round(f1 / df))
        ind2 = int(round(fac * f1 / df))
        if 1.2 * f1 >= endf:
            break
        elif fac * f1 > endf:
            ind2 = -1
        elif ind1 == ind2:
            ind2 += 1
        newdata[j] = np.mean(data[ind1:ind2])
        f1 *= 2 ** 0.125
        j += 1
    return newdata[:j]
def get_octfreqs(N, df, fac=2., startf=None):
    if startf == None:
        startf = df
    return 2 ** (0.125 * np.arange(N)) * startf * fac ** 0.5



def plot_psd(psd, freq=None, log=False, ax=None):
    if ax == None:
        plt.figure()
        ax1 = plt.subplot(111)
    else:
        ax1 = ax
    if freq == None:
        freq = np.linspace(0, 10, len(psd))
        print freq[1]
        print psd[1]
    if log:
        ax1.loglog(freq, psd)
    else:
        ax1.plot(freq, psd)

    def func(x):
        return 1. / x

    def update(ax):
        invfunc = func
        xmin, xmax = ax.get_xlim()
        if xmin == 0:
            xmin = xmax / 100.
        if xmax == 0:
            xmax = xmin + 100

        fmin = func(xmin)
        fmax = func(xmax)

        roundto = -int(np.round(np.log10(min(abs(fmin), abs(fmax), abs(fmax - fmin) / 20))))
        flabels = np.around(func(np.linspace(xmin, xmax, 10)), roundto)
        ax2.set_xticks(invfunc(flabels))
        ax2.set_xticklabels([str(i) for i in flabels])

    if not log and ax == None:
        ax2 = ax1.twiny()
        ax1.callbacks.connect('xlim_changed', update)
        ax2.set_xlim(ax1.get_xlim())
        update(ax1)
    return ax1

def test():
    from sito import read
    from obspy.signal.freqattributes import mper#, welch
    #from mtspec import mtspec

    ms = read('/home/richter/Data/Parkfield/raw/PKD_1996_296.mseed')
    #ms.plotTrace()

    print ms[0].stats
    # -*- snip -*-
    data = ms[0].data
    data = data - np.mean(data)
    #data -= np.linspace(0,1,len(data))*(data[-1]-data[0])+data[0]

    N = len(data)
    df = 1. / (ms[0].stats.endtime - ms[0].stats.starttime)
    print N // 2 * df


    spec1 = mper(data, cosTaper(N, 0.05), nextpow2(N))[:N // 2]
    #spec2 =  welch(data, cosTaper(N, 0.05), nextpow2(N), len(data)/10, 0.2)[:N//2]
    spec1_d = oct_downsample(spec1, df, fac=1.3)
    freq1 = get_octfreqs(len(spec1_d), df, fac=1.3)

    ax = plot_psd(spec1, log=True)
    ax = plot_psd(spec1_d, freq1, log=True, ax=ax)
    plt.show()

if __name__ == '__main__':
    test()
