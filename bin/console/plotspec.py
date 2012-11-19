#!/usr/bin/env python
# by TR

import argparse
from obspy.core import UTCDateTime as UTC
import logging
from obspy.core.utcdatetime import UTCDateTime
logging.basicConfig()

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('file_station',
                   help='file to plot or station to plot')
parser.add_argument('date', nargs='?', default=None, type=UTC,
                   help='if first argument is station: date')

parser.add_argument('-a', '--absolute-scale', type=float, default=0.0005,
                   help='display with different scale, default: 0.0005')
parser.add_argument('-r', '--relative-scale', type=float,
                   help='display with different relative scale - '
                        'overwrites ABSOLUTE_SCALE')
parser.add_argument('-s', '--save',
                   help='save plot to this file instead of showing')

parser.add_argument('-x', '--xcorr-append',
                   help='dont plot raw data and pass this argument to Data object')

parser.add_argument('-c', '--component', default='Z',
                   help='component to plot, default: Z')

parser.add_argument('-d', '--downsample', default=1,
                   help='downsample to this sampling rate, default: 1')

parser.add_argument('-o', '--options',
                   help='dictionary with kwargs passed to plotday')

args = parser.parse_args()

if args.relative_scale is not None:
    args.absolute_scale = None
if args.options is None:
    kwargs = {}
else:
    kwargs = eval('dict(%s)' % args.options)

kwargs.update(dict(absolutescale=args.absolute_scale,
                   scale=args.relative_scale,
                   downsample=args.downsample,
                   save=args.save, show=args.save is None))
print kwargs
if args.date is None:
    from sito import read
    from sito.imaging import plotTrace
    stream = read(args.file_station)
    plotTrace(stream, **kwargs)

else:
    from sito.imaging import plotTrace2
    station = args.file_station
    if station.startswith('PB') or station == 'LVC':
        from sito.data import IPOC
        data = IPOC(xcorr_append=args.xcorr_append)
    elif station == 'PKD':
        from sito.data import Parkfield
        data = Parkfield(xcorr_append=args.xcorr_append)
    else:
        raise argparse.ArgumentError('Not a valid station name')

    day = UTCDateTime(args.date)
    if args.xcorr_append is None:
        stream = data.getRawStream(day, station, component=args.component)
    else:
        stream = data.getStream(day, station, component=args.component)
    if stream[0].stats.is_fft:
        stream.ifft()
    plotTrace(stream, component=args.component, **kwargs)

if args.save is None:
    from matplotlib.pyplot import show
    show()

