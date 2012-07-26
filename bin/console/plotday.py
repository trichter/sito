#!/usr/bin/env python
# by TR

import argparse
from obspy.core import UTCDateTime as UTC
import logging
logging.basicConfig()

parser = argparse.ArgumentParser(description='Day plots.')
parser.add_argument('file_station',
                   help='file to plot or station to plot')
parser.add_argument('date', nargs='?', default=None, type=UTC,
                   help='if first argument is station: date')

parser.add_argument('-a', '--absolute-scale', type=float,
                   help='display with different scale')
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

parser.add_argument('-f', '--frequency', nargs='?', default=False, const=None,
                   help='plot frequency spectrum')

args = parser.parse_args()

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
    if not '.' in args.file_station:
        args.file_station = args.file_station + '.QHD'
    stream = read(args.file_station)
    if args.absolute_scale is None and args.relative_scale is None:
        kwargs['scale'] = 1.
else:
    station = args.file_station
    if station.startswith('PB') or station == 'LVC':
        from sito.data import IPOC
        data = IPOC(xcorr_append=args.xcorr_append)
    elif station == 'PKD':
        from sito.data import Parkfield
        data = Parkfield(xcorr_append=args.xcorr_append)
    else:
        raise argparse.ArgumentError('Not a valid station name')
    day = UTC(args.date)
    if args.xcorr_append is None:
        stream = data.getRawStream(day, station, component=args.component)
    else:
        stream = data.getStream(day, station, component=args.component)
    if args.absolute_scale is None and args.relative_scale is None:
        kwargs['absolutescale'] = 0.0005

tr = stream[0]
if args.frequency is False:
    if tr.stats.is_fft:
        tr.ifft()
    print tr
    tr.plotTrace(component=args.component, **kwargs)
else:
#    if not tr.stats.is_fft:
#        tr.fft(1024)
#    if tr.stats.is_fft:
#        if args.frequency is not None:
#            tr.ffttrim(*eval(args.frequency))
    print tr
    tr.plotPSD()

if args.save is None:
    from matplotlib.pyplot import show
    show()
