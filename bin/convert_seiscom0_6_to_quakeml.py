#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from sito.events import readSeisComPEventXML0_6

def parse_files(fnames):
    """Parses all given files for seiscomp xml"""
    out = None
    for i, fname in enumerate(fnames):
            print fname
            if not out:
                out = readSeisComPEventXML0_6(fname)
            else:
                out += readSeisComPEventXML0_6(fname)
            if (i + 1) % 100 == 0:
                out.write(str(i) + '.xml', 'QUAKEML')
                out = None
    out.write(str(i) + '.xml', 'QUAKEML')

if __name__ == '__main__':
    if len(sys.argv) <= 1 or sys.argv[1] in ('-h', '--help'):
        print(__doc__.strip('\n').replace('\n\n', '\n'))
        sys.exit()
    parse_files(sys.argv[1:])

