#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert SeiscompEventXML version 0.6 to QuakeML.

arguments: files to convert
100 Events are written in files ./1.xml ./2.xml and so on. 
"""
import sys
from sito.events import readSeisComPEventXML0_6
from obspy.core.event import Catalog

def parse_files(fnames):
    """Parses all given files for seiscomp xml"""
    j = 0
    out = Catalog()
    for i, fname in enumerate(fnames):
            print('read ' + fname)
            out += readSeisComPEventXML0_6(fname)
            if (i + 1) % 100 == 0 or i == len(fnames) - 1:
                out_fname = str(j) + '.xml'
                print('write %d events to %s\n' % (len(out), out_fname))
                out.write(out_fname, 'QUAKEML')
                out = Catalog()
                j += 1

if __name__ == '__main__':
    if len(sys.argv) <= 1 or sys.argv[1] in ('-h', '--help'):
        print(__doc__.strip('\n'))
        sys.exit()
    parse_files(sys.argv[1:])

