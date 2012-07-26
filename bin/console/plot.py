#!/usr/bin/env python
# by TR

import argparse
from sito import read
import logging
from sito.util.main import isnumber

logging.basicConfig()

parser = argparse.ArgumentParser(description='Plot data.')
parser.add_argument('file',
                   help='file')
parser.add_argument('args', nargs='*',
                   help='passed to plotting method')

parser.add_argument('-a', '--absolute-scale', type=float, default=None,
                   help='display with different scale, default: 0.0005')
parser.add_argument('-r', '--relative-scale', type=float, default=1.,
                   help='display with different relative scale - '
                        'overwrites ABSOLUTE_SCALE')
parser.add_argument('-s', '--save',
                   help='save plot to this file instead of showing')

parser.add_argument('-x', '--xcorr', action='store_true',
                   help='use plotXcorr method')

parser.add_argument('-y', '--rf', action='store_true',
                   help='use plotRF method')

parser.add_argument('-c', '--component', default='Z',
                   help='component to plot, default: Z')

parser.add_argument('-d', '--downsample', default=1,
                   help='downsample to this sampling rate, default: 1')

parser.add_argument('-o', '--options',
                   help='dictionary with kwargs passed to plotday')

parser.add_argument('-f', '--frequency', action='store_true',
                   help='plot frequency spectrum')

args = parser.parse_args()

if args.options is None:
    kwargs = {}
else:
    kwargs = eval('dict(%s)' % args.options)

kwargs.update(dict(component=args.component,
                   absolutescale=args.absolute_scale,
                   scale=args.relative_scale,
                   downsample=args.downsample,
                   save=args.save, show=args.save is None))
if args.frequency:
    kwargs.pop('absolutescale')
    kwargs.pop('scale')
    kwargs.pop('downsample')
    # TODO: change that save is working
    kwargs.pop('save')
    kwargs.pop('component')
    kwargs.pop('show')


if not '.' in args.file[-4:]:
    args.file = args.file + '.QHD'
stream = read(args.file)

print stream

if args.frequency:
    meth = 'plotPSD'
elif args.rf:
    meth = 'plotRF'
elif args.xcorr:
    meth = 'plotXcorr'
else:
    meth = 'plot_'

if args.frequency is False and stream[0].stats.is_fft:
    stream.ifft()

def convert_strings_to_floats(tup):
    return [float(entry) if isnumber(entry) else entry for entry in tup]
#    tup = list(tup)
#    for i, entry in enumerate(tup):
#        try:
#            tup[i] = float(entry)
#        except:
#            pass
#    return tup
getattr(stream, meth)(*convert_strings_to_floats(args.args), **kwargs)

if args.save is None:
    from matplotlib.pyplot import show
    show()
