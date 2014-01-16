#!/usr/bin/env python

from sito import read
from glob import glob
input = '/home/richter/Results/IPOC/xcorr/FINAL_filter4-6_1bit_auto_hour/xcorr/hout*.QHD'

for fname in glob(input):
    splits = fname.rsplit('/', 2)
    out = '%s/plots/%s.png' % (splits[0], splits[-1].rsplit('.', 1)[0])
    ms = read(fname)
    ms.plotXcorr(0, 20, dateformatter='%y-%m-%d', show=False,
                                     save=out, vmax=0.1)
