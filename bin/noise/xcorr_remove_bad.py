#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR

from sito.data import IPOC
from sito.noisexcorr import  removeBad
import glob
from sito import read
import os.path

method = 'filter0.01-1_water_env2_whitening_1bit_fft'
data = IPOC(xcorr_append='/Tocopilla/' + method, use_local_LVC=True)
data.setXLogger('_' + method)
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter4-6_water_env2_1bit/stretch_t/'
#path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit/stretch2/'
path = '/home/richter/Results/IPOC/xcorr/Tocopilla/filter0.01-1_water_env2_whitening_1bit_fft/xcorr/'

for file in glob.glob(path + '*.QHD'): #@ReservedAssignment
    if not 'filter' in os.path.basename(file):
        print file
        ms = read(file)
        ms.normalize()
        removeBad(ms, 0.8)
        ms.write(os.path.splitext(file)[0], 'Q')
