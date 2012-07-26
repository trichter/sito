#!/usr/bin/env python
# by TR

import argparse
from obspy.core import UTCDateTime as UTC
from sito.util import isint

parser = argparse.ArgumentParser(description=
                                 'Convert Julian day date to normal date or '
                                 'backwards.')
parser.add_argument('date')
parser.add_argument('-u', '--underscores', action='store_true',
                   help='use underscores instead of minus in output')

args = parser.parse_args()
date = args.date
under = args.underscores

utc = UTC(date)
if len(date) < 8 or (len(date) == 8 and not isint(date)):
    if under:
        print utc.strftime('%Y_%m_%d'),
    else:
        print utc.date,
else:
    if under:
        print utc.strftime('%Y_%j'),
    else:
        print utc.strftime('%Y-%j'),
