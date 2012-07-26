#!/usr/bin/env python
import logging
from sito import read, imaging
import os

logging.basicConfig(level=logging.DEBUG)

start = -20
end = 100

# test waterlevel deconvolution against SH
farm = read(os.path.join(os.path.dirname(__file__), 'data', 'TEST_FARM.QHD'))

dec1 = farm.copy()
dec2 = farm.copy()
dec3 = farm.copy()
shdec1 = farm.copy()
shdec2 = farm.copy()
shdec3 = farm.copy()

dec1.receiverf(water=0.005, gauss=10, tshift=20, pad=0, window='tukey', start= -10, end=20, where='ponset', lenslope=5)
dec2.receiverf(water=0.001, gauss=10, tshift=20, pad=0, window='tukey', start= -10, end=50, where='ponset', lenslope=5)
dec3.receiverf(water=0.0001, gauss=10, tshift=20, pad=0, window='tukey', start= -10, end=100, where='ponset', lenslope=5)
imaging.compareRF([dec1, dec2, dec3], start, end, component='Q')

shdec1.receiverSH(start= -10, end=100, spiking=1, cut1= -20, cut2=150)
imaging.compareRF([shdec1, dec1], start, end, component='Q')

shdec2.receiverSH(start= -10, end=50, spiking=1, cut1= -20, cut2=50)
shdec3.receiverSH(start= -10, end=20, spiking=1, cut1= -20, cut2=50)
imaging.compareRF([shdec1, shdec2, shdec3], start, end, component='Q')

from IPython import embed
embed()
