#!/usr/bin/env python
# by TR

import argparse
import logging
from sito import read
import numpy as np
import os.path
logging.basicConfig()

parser = argparse.ArgumentParser(description='Check data in streams for NAN and INF.')
parser.add_argument('files', nargs='*',
                   help='files')
parser.add_argument('-r', '--remove', action='store_true',
                   help='remove traces with NAN and INF data')
parser.add_argument('-n', '--not-inf', action='store_false',
                   help='check not for infinity')
parser.add_argument('-b', '--bigger', type=float,
                   help='check values bigger than this')
parser.add_argument('-l', '--load-format', default=None,
                   help='file format to speed up loading')
parser.add_argument('-f', '--format', default='Q',
                   help='file format to save stream with removed traces')

args = parser.parse_args()

for file in args.files: #@ReservedAssignment
    print 'checking file' + file
    ms = read(file, format=args.load_format)
    found_sth = False
    for tr_no, tr in enumerate(ms):
        if args.not_inf:
            mask = np.isnan(tr.data)
        else:
            mask = np.logical_or(np.isnan(tr.data), np.isinf(tr.data))
        if args.bigger:
            mask = np.logical_or(mask, np.abs(tr.data) >= args.bigger)
        if np.any(mask):
            found_sth = True
            Nnan = np.count_nonzero(np.isnan(tr.data))
            Ninf = np.count_nonzero(np.isinf(tr.data))
            if args.bigger:
                Nbigger = np.count_nonzero(np.abs(tr.data) >= args.bigger)
                if Nbigger > 0:
                    print '%d entries bigger than %f in trace no. %d in file %s\n     %s' % (Nbigger, args.bigger, tr_no, file, tr)
            if Nnan > 0 or Ninf > 0:
                print '%d NANs and %d INFS in trace no. %d in file %s\n     %s' % (Nnan, Ninf, tr_no, file, tr)
            if args.remove:
                ms.remove(tr)
    if found_sth and args.remove:
        print 'write new file...'
        ms.write(os.path.splitext(file)[0], args.format)
print 'checked %d files' % len(args.files)
