#!/usr/bin/env python
# TR
import os
modules = 'stations events data trace stream rf imaging' #xcorr
os.chdir('../')
for module in modules.split():
    command = 'python %s.py' % module
    print(command + ' ...\n')
    os.system(command)

