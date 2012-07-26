#!/usr/bin/python
# by TR
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors
import os
import glob

CM_DATA = '/home/richter/Data/cm/'

def combine(cmaps, name, splitters=None, get_cdict=False):
    if not splitters:
        N = len(cmaps)
        splitters = np.linspace(0, 1, N + 1)
    cdict = {}
    for i, m in enumerate(cmaps):
        m = plt.get_cmap(m)
        if hasattr(m, '_segmentdata'):
            m = m._segmentdata
        for color in m:
            m[color] = np.array(m[color])
            m[color][:, 0] = (splitters[i] + (m[color][:, 0] - m[color][0, 0]) /
                              (m[color][-1, 0] - m[color][0, 0]) *
                              (splitters[i + 1] - splitters[i]))
            try:
                cdict[color] = np.concatenate((cdict[color], m[color]))
            except KeyError:
                cdict[color] = m[color]
    if get_cdict:
        return cdict
    else:
        return matplotlib.colors.LinearSegmentedColormap(name, cdict)


def show_colormaps(mode='mpl', path=CM_DATA + '*.gpf', cmaps=None):
    plt.rc('text', usetex=False)
    a = np.outer(np.ones(10), np.arange(0, 1, 0.01))
    plt.figure(figsize=(10, 5))
    plt.subplots_adjust(top=0.99, bottom=0.01, left=0.01, right=0.8)
    if mode == 'mpl':
        cmaps = [m for m in plt.cm.datad if not m.endswith("_r")]
    elif mode == 'local':
        cmaps = [createColormapFromGPF(f) for f in glob.glob(path)]
    elif cmaps is None:
        raise ValueError("Mode has to be 'mpl' or 'local' or cmaps=list of cmaps")
    cmaps.sort()
    l = len(cmaps) + 1
    for i, m in enumerate(cmaps):
        plt.subplot(l, 1, i + 1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=plt.get_cmap(m), origin="lower")
        plt.annotate(m.name if hasattr(m, 'name') else m, (1, 0.5), xycoords='axes fraction',
                     fontsize=10, ha='left', va='center')
    return cmaps


def createColormapFromGPF(file_, get_dict=False):
    data = sp.loadtxt(file_)
    cdict = {'red': np.take(data, (0, 1, 1), axis=1),
             'green': np.take(data, (0, 2, 2), axis=1),
             'blue': np.take(data, (0, 3, 3), axis=1)}
    name = os.path.splitext(os.path.basename(file_))[0]
    if get_dict:
        return cdict
    else:
        return matplotlib.colors.LinearSegmentedColormap(name, cdict)

if __name__ == '__main__':
    pass
