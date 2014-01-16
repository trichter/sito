#!/usr/bin/env python

import os

#cs = (100., 50., 150., 200., 300)
#ls = [float(10 * l) for l in range(1, 31)]
cs = (1000.,)
ls = (50, 100, 150, 200, 250, 300, 350, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)
ls = (400.,)

for c in cs:
    for l in ls:
        cmd = 'model.py -c %s -l %s' % (c, l)
        print cmd
        os.system(cmd)
