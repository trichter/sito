#!/usr/bin/env python
from sito import imaging, read
import os

start = -5
end = 20

# test moveoutcorrection against SH
dec = read(os.path.join(os.path.dirname(__file__), 'data', 'TEST_DEC.QHD'))

mout = dec.copy()
mout.moveout()
imaging.compareRF([dec, mout], start, end, component='Q')

from IPython import embed
embed()


# test zmigr - zmigr does not work !

#test1
import seis.model
import numpy as np
import seis._rf
model = seis.model.Iasp91()
data = np.sin(np.linspace(0, 100, 1000))
data2 = seis._rf.zmigr(model.z, model.vp, model.vs,
                        data,
                        4.6, 10., 0., 100, 0.5, 0) #slowness, samplingrate, startz, endz, dz, phase 0 ='Ps'
#data2  = seis._rf.zmigr(model.z, model.vp, model.vs, data, 8.0, 100., 0., 100., 0.5, 0)
#util.ipshell()

#test2
#dec = MS.read('TEST_SHDEC.QHD')[34*3:40*3]
#dec.zmigr()
#st = dec[1].stats
#N = st.ntps
#dz = 1/st.sampling_rate
#plt.plot(arange(N)*dz, dec[1].data)

