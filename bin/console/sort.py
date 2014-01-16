#!/usr/bin/env python
# by TR

import argparse
from sito import read
import glob

parser = argparse.ArgumentParser(description='Sort files by stats field and output filenames.')
parser.add_argument('files', help='files surrounded by ""')
parser.add_argument('field', help='stats field')
parser.add_argument('-r', '--reverse', action='store_false', help='sort in reverse order')

args = parser.parse_args()
ms = read(args.files)
fnames = glob.glob(args.files)
ms.setHI('fname', fnames)
ms.sort(keys=args.field.split(), reverse=args.reverse)
fnames = ms.getHI('fname')

print ' '.join(fnames),
