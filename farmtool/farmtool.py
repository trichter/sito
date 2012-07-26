#!/usr/bin/env python
# License: GNU LGPL Version 3
# see http://www.gnu.org/licenses/lgpl.html
# Copyright 2012 Tom Richter

from sito import read
from sito.imaging import plotRFmarks
import matplotlib.pyplot as plt

#import sys
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Verbose',
#color_scheme='LightBG', call_pdb=1)
#from IPython import embed

t1, t2, t3, t4, lw = 70, 80, -5, -3, 2
select = True
num_tr = None
patho1 = '/home/richter/Data/IPOC/receiver/2012_events_mag5.5/'
patho2 = '/home/richter/Results/IPOC/receiver/2012_mag5.5/'
paths = '/home/richter/Results/IPOC/receiver/2012_mag5.5/farm/'
fileo1 = '/home/richter/Data/IPOC/receiver/2012_events_mag5.5/LVC_*.QHD'
fileo2 = '/home/richter/Results/IPOC/receiver/2012_mag5.5/LVC_mout.QHD'
files = '/home/richter/Results/IPOC/receiver/2012_mag5.5/temp/*_mout'

#bandpass = (0.03, 0.1, 2, True)
#bandpass = None #(0.5, None, 2, True) # High Pass

sort = ['azi', 'component']
class FarmTool(object):
    def __init__(self, files=None, select=True, bandpass=None):
        self.select = select
        self.bandpass = bandpass
        self.fig = plt.figure()
        self.fileo1 = self.fileo2 = self.files = None
        if files is not None and len(files) > 1:
            self.fileo1 = files[0]
            self.fileo2 = files[1]
            if len(files) > 2:
                self.fileos = files[2]
        self.open_files()
        self.pressed_m = None

    def connect(self):
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
    def open_files(self):
        if self.fileo1 is None:
            from easygui import filesavebox
            file1 = filesavebox(msg='Choose file with RAW data',
                                        default=patho1,
                                        filetypes=['*.QHD', '*.*'])
        else:
            file1 = self.fileo1
        if self.fileo2 is None:
            file2 = filesavebox(msg='Choose file with RF data',
                                        default=patho2,
                                        filetypes=['*.QHD', '*.*'])
        else:
            file2 = self.fileo2
        if file1 is not None and file2 is not None:
            try:
                self.st1 = read(file1)
                self.st2 = read(file2)
            except ValueError:
                print("Error while reading file!")
                return

            self.st2s = self.st2
            if self.select:
                self.st2s = self.st2.select(expr='not st.mark')
            event_ids = self.st2s.getHI('event.id')
            self.st1s = self.st1
            if  len(self.st1s) != len(self.st2s):
                for tr in self.st1s:
                    if tr.stats.event.id not in event_ids:
                        self.st1s.remove(tr)
            if  len(self.st1s) != len(self.st2s):
                raise ValueError('Streams do not have the same lenth')

            if num_tr is None:
                self.ind1 = 0
                self.ind2 = len(self.st2s)
            else:
                self.ind1 = 0
                self.ind2 = min(3 * num_tr, len(self.st2s))
            if self.bandpass is not None:
                self.st2s.filter2(*self.bandpass)
            #self.st2s.select(component='Q').normalize()
            self.st1s.sort(sort)
            self.st2s.sort(sort)
            self.plot_streams()

    def save_file(self):
        if self.files is None:
            from easygui import filesavebox
            file_ = filesavebox(msg='Please omit file ending. Include * in '
                                    'to write one file for each year!',
                                    default=paths)
        else:
            file_ = self.files
        if file_ is not None:
            if '*' in file_:
                file_ = file_.replace('*', '%s')
                self.st2.writex(file_, 'Q', years=False)
            else:
                self.st2.write(file_, 'Q')
            print('Data saved!')

    def mark_trace(self, i):
        if self.st2s[i].stats.mark:
            self.st2s[i].stats.mark = False
            self.st2s[i + 1].stats.mark = False
            self.st2s[i + 2].stats.mark = False
            self.st2s[i].stats.reason = ''
            self.st2s[i + 1].stats.reason = ''
            self.st2s[i + 2].stats.reason = ''
            self.ax1.plot([t1, t2], [i / 3, i / 3], 'w', linewidth=lw)
            self.ax2.plot([t3, t4], [i / 3, i / 3], 'w', linewidth=lw)

        else:
            self.st2s[i].stats.mark = True
            self.st2s[i + 1].stats.mark = True
            self.st2s[i + 2].stats.mark = True
            self.st2s[i].stats.reason = 'eye'
            self.st2s[i + 1].stats.reason = 'eye'
            self.st2s[i + 2].stats.reason = 'eye'
            self.ax1.plot([t1, t2], [i / 3, i / 3], 'r', linewidth=lw)
            self.ax2.plot([t3, t4], [i / 3, i / 3], 'r', linewidth=lw)
        plt.draw()


    def plot_streams(self):
        try:
            xlim1 = self.ax1.get_xlim()
            xlim2 = self.ax2.get_xlim()
            flag = True
        except AttributeError:
            flag = False
        self.fig.clear()
        self.ax1 = self.fig.add_subplot(121)
        self.ax2 = self.fig.add_subplot(122, sharey=self.ax1)
        self.pl1 = self.st1s[self.ind1:self.ind2].plot_(
                       ax=self.ax1, component='Z', plotinfo=('sum',),
                       plotinfowhere=('top',), plotinfodicts=[dict(pad=0,
                                                                   size=0.4)])
        self.pl2 = self.st2s[self.ind1:self.ind2].plotRF(ax=self.ax2)
        plotRFmarks(self.st2s[self.ind1:self.ind2:3], self.ax1,
                    t1=t1, t2=t2, options='r', lw=2)
        plotRFmarks(self.st2s[self.ind1:self.ind2:3], self.ax2,
                    t1=t3, t2=t4, options='r', lw=2)
        if flag:
            self.ax1.set_xlim(xlim1)
            self.ax2.set_xlim(xlim2)

    def onclick(self, event):
        i = event.ydata
        if abs(i - round(i)) > 0.4:
            return
        i = 3 * int(round(i))
        if event.button in (2, 3):
            self.pressed_m = None
        if event.button == 2:
            print self.st2s[i].stats
        elif event.button == 3:
            self.mark_trace(i)
        elif event.button == 1 and self.pressed_m is not None:
            if self.pressed_m is True:
                self.pressed_m = i
            else:
                i1 = min(i, self.pressed_m)
                i2 = max(i, self.pressed_m)
                self.pressed_m = None
                for j in range(i1, i2 + 1, 3):
                    self.mark_trace(j)



    def onkey(self, event):
        key = event.key
        if key == 'i':
            from IPython import embed
            embed()
        elif key == 'left' and num_tr is not None:
            if self.ind1 == 0:
                print('reached beginning of file')
            else:
                self.ind2 = self.ind1
                self.ind1 = max(0, self.ind1 - 3 * num_tr)
                self.plot_streams()
        elif key == 'right' and num_tr is not None:
            if self.ind2 == len(self.st2s):
                print('reached end of file')
            else:
                self.ind1 = self.ind2
                self.ind2 = min(self.ind2 + 3 * num_tr, len(self.st2s))
                self.plot_streams()
        elif key == 'o':
            self.open_files()
        elif key == 'w':
            self.save_file()
        elif key == 'm':
            self.pressed_m = True


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Farm receiver function data. '
                                     'stats.mark of farmed data will be set to '
                                     'True.')
    parser.add_argument('files', nargs='*',
                       help='File to read. First file: raw data. Second file: '
                       'receiver function data. Third file (optional): '
                       'File to write the results.')
    parser.add_argument('-d', '--display', action='store_true',
                       help='Display not only unmarked data but also marked '
                       'data.')
    parser.add_argument('-b', '--band',
                       help='Band pass for receiver function data '
                       'in the form (min_freq, max_freq, '
                       'corners, zero_phase) eg. (0.03, 0.1, 2, True) or for a '
                       ' high pass (0.5, None, 2, True).')
    args = parser.parse_args()
    if len(args.files) == 0:
        args.files = [fileo1, fileo2, files]
    if args.band is not None:
        args.band = eval(args.band)
    ft = FarmTool(files=args.files, select=not args.display,
                  bandpass=args.band)
    ft.connect()
    plt.show()
