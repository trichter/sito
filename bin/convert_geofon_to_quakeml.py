#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert SeiscompEventXML version 0.6 to QuakeML.

arguments: files to convert
write to output_geofon.xml
"""
import sys
from sito.events import read_regex
from obspy.core.event import Catalog

def read_files(fnames):
    out = Catalog()
    for fname in fnames:
        out += read_regex(fname)
    out_fname = 'output_geofon.xml'
    out.write(out_fname, 'QUAKEML')

if __name__ == '__main__':
    if len(sys.argv) <= 1 or sys.argv[1] in ('-h', '--help'):
        print(__doc__.strip('\n'))
        sys.exit()
    read_files(sys.argv[1:])

