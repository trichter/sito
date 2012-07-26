#!/usr/bin/env python
# by TR

import argparse
import logging
from sito import read
import os.path
logging.basicConfig()

parser = argparse.ArgumentParser(description='Trim data. Sea sito.stream.Stream.trim2')
parser.add_argument('files', nargs='*',
                   help='files')
parser.add_argument('-s', '--start',
                   help='start time')
parser.add_argument('-e', '--end',
                   help='end time')
parser.add_argument('-r', '--relative', default='starttime',
                   help='trim relative to this value')
parser.add_argument('-l', '--load-format', default='Q',
                   help='file format to speed up loading')
parser.add_argument('-f', '--format', default='Q',
                   help='file format to save stream with removed traces')

args = parser.parse_args()

start = args.start
end = args.end
if start:
    start = eval(start)
if end:
    end = eval(end)

for file in args.files: #@ReservedAssignment
    print 'loading file' + file
    ms = read(file, format=args.load_format)
    ms.trim2(start, end, args.relative)
    print 'write trimmed file...'
    ms.write(os.path.splitext(file)[0], args.format)
print 'trimmed %d files' % len(args.files)
