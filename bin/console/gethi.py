#!/usr/bin/env python
# by TR

import argparse
import logging
from sito import read
from sito.util.main import isint, isfloat
logging.basicConfig()

parser = argparse.ArgumentParser(description='Get header entry of traces in files. Sea sito.Stream.g(s)etHI')
parser.add_argument('header', help='header e.g. dist')
parser.add_argument('files', nargs='*',
                   help='files')
parser.add_argument('-w', '--write', help='write entry into header')
args = parser.parse_args()
header = args.header
if not args.write:
    entries = []
    for f in args.files:
        ms = read(f)
        entries.extend(ms.getHI(header))
    print ' '.join([str(entry) for entry in entries])
else:
    entry = args.write
    entry = int(entry) if isint(entry) else float(entry) if isfloat(entry) else entry
    for f in args.files:
        print 'loading file' + f
        ms = read(f)
        ms.setHI(header, entry)
        print 'write file...'
        ms.write(f, ms[0].stats._format)
