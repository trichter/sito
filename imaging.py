#!/usr/bin/env python
# -*- coding: utf-8 -*-
# by TR


from matplotlib import cbook
from matplotlib.mlab import psd
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from obspy.core import UTCDateTime
from obspy.core.util import deprecated
from sito.util import parameters, ttt, add_doc, calculate
from sito.util.imaging import xcorr_cmap, DLogNorm, getDataWindow, getTimeIntervall
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
log = logging.getLogger(__name__)
cc = mpl.colors.ColorConverter()

def plot_spectrograms(stream, **kwargs):
    fig = plt.figure()
    ax1 = None
    num1 = 2
    num2 = (len(stream) - 1) // num1 + 1
    for i, tr in enumerate(stream):
        if not ax1:
            ax1 = ax = fig.add_subplot(num2, num1, i)
        else:
            ax = fig.add_subplot(num2, num1, i, sharex=ax1, sharey=ax1)
        print tr.stats.station
        tr.spectrogram(axes=ax, title=tr.stats.station, **kwargs)
        ax.set_title(tr.stats.station)

#deltrcs = False
#delindex = []

def getPublicationFigure(axes=None, width=10, ratio=0.618, margin=None, fontsize=10, labelsize=8, backend='eps', usetex=False, distiller='ghostscript'):
    """
    Return Figure instance.
    """
    if not margin:
        margin = [1., 0.1, 1., 0.1] #left, rigth, bottom, top
    return getFigure(axes=axes, width=width, ratio=ratio, margin=margin, fontsize=fontsize, labelsize=labelsize, backend=backend, usetex=usetex, distiller=distiller)

def getFigure(axes=None, width=30, ratio=0.618, margin=None, fontsize=20, labelsize=18, backend='png', usetex=False, distiller='ghostscript'):
    """
    Return Figure instance.

    axes: None -> only one axis
    axes: numpy_array: splitted figure (first row relative widths, second row relative heights)
        sum of the row has to be 1 or smaller. if smaller appropriate space is left between the axes
    width: of whole figure in cm
    ratio: heigth / width
    margin: [left, right, bottom ,top] in cm
    fontsize, labelsize: in cm
    backend: 'ps', 'png'
    """
    #fig_width_pt =   # Get this from LaTeX using \showthe\columnwidth
    #inches_per_pt = 1.0/72.27               # Convert pt to inch
    if not margin:
        margin = [2., 1., 2., 1.] #left, rigth, bottom, top
    fig_width = width / 2.54  # width in inches
    fig_height = width / 2.54 * ratio # height in inches
    margin = np.array(margin)
    margin[:2] = margin[:2] / width  # relative to fig size
    margin[2:] = margin[2:] / width / ratio
    fig_size = [fig_width, fig_height]
    params = {'backend': backend,
        'axes.labelsize': fontsize,
        #'axes.unicode_minus': False, # to save labels as text
        'text.fontsize': fontsize,
        'legend.fontsize': fontsize,
        'xtick.labelsize': labelsize,
        'ytick.labelsize': labelsize,
        'font.size': fontsize,
        'text.usetex': usetex,
        'figure.figsize': fig_size,
        'lines.linewidth': 0.5,
        'lines.markeredgewidth'  : 1.2,
        #'path.simplify' : False,
        #'path.simplify_threshold' : 0.1,
        #'ps.useafm'         : True,    # use of afm fonts, results in small files
        'ps.papersize': 'auto',
        'ps.usedistiller': distiller    # can be: None, ghostscript or xpdf

            # Experimental: may produce smaller files.
            # xpdf intended for production of publication quality files,
            # but requires ghostscript, xpdf and ps2eps
#ps.distiller.res  : 6000      # dpi
#ps.fonttype       : 3         # Output Type 3 (Type3) or Type 42 (TrueType)
}
    plt.rcParams.update(params)
    # Generate data
    fig = plt.figure()
    plot_width = 1 - margin[:2].sum()  # relative to fig size
    plot_height = 1 - margin[2:].sum()
    if not axes:
        # ax = fig.add_axes([margin[0],margin[2],plot_width,plot_height])
        fig.add_axes([margin[0], margin[2], plot_width, plot_height])
    else: #only horzontal split
        if not isinstance(axes[0], list):
            axes = [axes, [1]]
        if len(axes) == 2: # horizontal and vertical split
            Nx = len(axes[0])
            Ny = len(axes[1])
            axes[0] = [i * plot_width for i in axes[0]]
            axes[1] = [i * plot_height for i in axes[1]]
            spacex = spacey = 0
            if Nx > 1:
                spacex = (plot_width - sum(axes[0])) / (Nx - 1)
            if Ny > 1:
                spacey = (plot_height - sum(axes[1])) / (Ny - 1)

            startx = [0] + [sum(axes[0][0:i + 1]) + spacex * (i + 1) for i in range(Nx)]
            starty = [0] + [sum(axes[1][0:i + 1]) + spacey * (i + 1) for i in range(Ny)]
            #ax = []
            for j in range(Ny):
                for i in range(Nx):
                    if j > 0 or i == 0:
                        # ax.append(fig.add_axes([margin[0]+startx[i],margin[2]+starty[j],axes[0][i],axes[1][j]]))
                        fig.add_axes([margin[0] + startx[i], margin[2] + starty[j], axes[0][i], axes[1][j]])
                    else:
                        # ax.append(fig.add_axes([margin[0]+startx[i],margin[2]+starty[j],axes[0][i],axes[1][j]], sharey=ax[0]))
                        fig.add_axes([margin[0] + startx[i], margin[2] + starty[j], axes[0][i], axes[1][j]], sharey=fig.axes[0])
        else:
            return None
    return fig

def get_fig(ax=None, positions=(), adicts=None):
    if ax is None:
        ax = plt.figure().add_subplot(111)
    divider = make_axes_locatable(ax)
    if adicts is None:
        adicts = [dict(pad=0, size=0.8) for i in positions]
    for i, adict in enumerate(adicts):
        if (not adict.has_key('sharex')) and (not adict.has_key('sharey')):
            if positions[i] in ('right', 'left'):
                adicts[i]['sharey'] = ax
            else:
                adicts[i]['sharex'] = ax
    add_axes = []
    for i in range(len(positions)):
        add_axes.append(divider.append_axes(positions[i], **(adicts[i])))
    return add_axes
        #divider.append_axes("top", size=1.2, pad=0.1, sharex=ax)


#    pylab.plot(x,y1,'g:',label='$\sin(x)$')
#    pylab.plot(x,y2,'-b',label='$\cos(x)$')
#   pylab.xlabel('$x$ (radians)')
#    pylab.ylabel('$y$')
#    pylab.legend()
#    pylab.savefig('fig1.eps')


def _insert_zeros(stream, data, min_delta=None):
    N = len(stream)
    npts = np.shape(data)[1]
    starttimes = stream.getHI('starttime')
    deltas = np.array([starttimes[i + 1] - starttimes[i]
                       for i in range(N - 1)])
    if min_delta is None:
        min_delta = np.median(deltas)
    indices = np.nonzero(deltas - min_delta >= 1)
    nums = (np.round(deltas[indices] / min_delta) - 1).astype('int')
    counter = 0
    for i in range(len(nums)):
        index = indices[0][i]
        num = nums[i]
        data = np.vstack((data[:counter + index + 1, :],
                          np.zeros((num, npts)), data[counter + index + 1:]))
        counter += num
    return data

def UTC2year(utc):
    import calendar
    year = utc.year
    return year + utc.julday / (365. + calendar.isleap(year))

def plotRFmarks(stream, ax, t1= -20, t2= -10, options='r', lw=2):
    st = stream.select(component='Z')
    if len(st) == 0:
        st = stream.select(component='L')
    for i, tr in enumerate(st):
        if tr.stats.mark == True:
            ax.plot([t1, t2], [i, i], options, linewidth=lw)

def plotPhases(ms, ax, plotphases='some'):
    """
    Plot phases in given axe or axe list.

    ax: axes instances or list of axes instances
    plotphases: 'all', 'some', 'all3', 'some3'.
    """

    if plotphases == True:
        plotphases = 'some'
    for i, trace in enumerate(ms):
        arrivals = ttt(trace.stats.dist, trace.stats.event.depth, True)
        t0 = arrivals[0].time
        for a in arrivals:
            t = a.time - t0
            if a.phase in ['P', 'Pdiff', 'PcP', 'pP', 'sP', 'PP', 'S', 'Surf'] or 'all' in plotphases:
                if type(ax) != list and t > ax.get_xlim()[0] and t < ax.get_xlim()[1]:
                    ax.plot([t, t], [i - 0.5, i + 0.5], 'r')
                    if i == 0 or (i % 3 == 0 and '3' in plotphases):
                        ax.annotate(a.phase, xy=(t, i - 1), color='r')
                if isinstance(ax, list) and i % 3 == 0:
                    which = i // 3
                    if t > ax[which].get_xlim()[0] and t < ax[which].get_xlim()[1]:
                        ax[which].plot([t, t], ax[which].get_ylim(), 'r')
                        if i == 0 or '3' in plotphases:
                            ax[which].annotate(a.phase, xy=(t, ax[which].get_ylim()[0] * 0.95), color='r')




class Plot(object):
    def __init__(self, stream, start=None, end=None, relative='starttime',
                 rel_label='relative', component='all',
                 filter=None, downsample=None, #@ReservedAssignment
                 xaxis='data', yaxis='num', dateformatter='%y-%m-%d',
                 reverse_x=False, reverse_y=False, minor_x=True, minor_y=True,
                 xlabel=None, ylabel=None,
                 color='kk', topcolor='white', botcolor='white', fast=False, #@UnusedVariable
                 scale=1., absolutescale=None, sumscale=2.,
                 imshow=False, cmap=None, colorbar=True, use_dlognorm=False,
                 alpha=None, #@UnusedVariable
                 vmax=None, vmin=None, cmax=1e-5,
                 plotsum=False, order=None, plotphases=False,
                 figtitle='station component sc:scale', title_xpos=0.5,
                 title_horalign='center', title_in_axis=False, fancy_box=False, box_trans='ax',
                 box_ax=None, box_fs=14,
                 show=True, save=False, #publication=False,#delete=False,
                 fig=None, ax=None, connect_event=True,
                 plotinfo=(), usehardticks='', #plotinfo_width=0.1, #@UnusedVariable
                 plotlabel=None, ax_info=None, #@UnusedVariable
                 plotinfowhere=None, plotinfodicts=None, #@UnusedVariable
                 plot_stack=False, stack_lim=None, plot_psd=False, #@UnusedVariable
                 psd_scale='time', psd_prop=(4096, True, None), #@UnusedVariable
                 annotate=None #@UnusedVariable
                 ):
        """
        Plot stream...
                
        @param stream: stream
        @param start: start time relative to param relative
        @param end: end time relative to param relative
        @param relative: time object, see sito.util.getTimeIntervall
        @param rel_label: time object, labeling relative to this time
        @param component: component or 'all'
        @param xaxis: one of ('data', 'num', 'date', 'sum') or header entries
        @param yaxis: one of ('data', 'num', 'date', 'sum') or header entries
        @param dateformatter: formatter string for dates e.g. '%y-%m-%d'
        @param reverse_x: reverse x-axis?
        @param reverse_y: reverse y-axis?
        @param minor_x: minor ticks on x-axis?
        @param minor_y: minor ticks on y-axis?
        @param xlabel: label for x-axis
        @param ylabel: label for y-axis
        @param color: alternat. color for line plot e.g. 'kk' or ('red','blue')
        @param topcolor: color for filling the upper side of line plot
        @param botcolor: color for filling the lower side of line plot
        @param fast: if True sets params topcolor and botcolor to 'white'
        @param scale: relatvie scale
        @param absolutescale: if set use this absolute scale
        @param sumscale: scale for summation trace relative to normal scale 
        @param imshow: dont plot lines but an image
        @param cmap: colormap for image
        @param colorbar: plot the colorbar for image?
        @param use_dlognorm: imshow in logarithmic scale
        @param vmax: scale for imshow, None or float (if None take maximum)
        @param vmin: scale for imshow, None or float (if None take -vmax)
        @param plotsum: plot the summation trace inside axis?
                        (use plot_stack instead)
        @param order: if set use phaseStack for plotsum
        @param plotphases: plotarriving  phases?
                           only possible for param relative='ponset'
        @param figtitle: title of figure
        @param title_xpos: x-position of title
        @param title_horalign: horizontal alignment of title
        @param title_in_axis: display title in axis?
        @param show: show figure?
        @param save: if set saves figure to this filename
        @param fig: use this figure for plots
        @param ax: use this ax for plot
        @param connect_event: connect scaling event for line plot
        @param plotinfo: plot information in additional axes
                         plotinfo should be 'sum' or a header entry or two
                         header entries seperated by one white space
        @param usehardticks: use in this infos given tick positions
                             (see source code)
        @param ax_info: plot information in these axes (if not given taken from
                        ':param fig:'_ or added to ':param ax:'_)
        @param plotinfodicts: dictionaries passed to get_figure
        @param plot_stack: plot a stack of the traces (add to plotinfo & Co)
        @param stack_lim: y ( or xlimit) of stack as tuple or list (None)
        @param plot_psd: plot a PSD of sum (add to plotinfo & Co)
        @param psd_scale: time or freq
        @param publication: not used at the moment
        """
        ### get and store parameters
        if len(stream) == 0:
            log.warning('You try to plot a stream (id: %s) with 0 traces.' %
                        stream.hash)
            return
        log.info('Plot stream %s: %s' % (stream.hash, parameters()))
        mykwargs = parameters(only=['color', 'minor_x', 'minor_y', 'plotsum',
                                    'plotinfo', 'usehardticks', 'plotlabel',
                                    'plotinfowhere', 'stack_lim',
                                    'psd_scale', 'psd_prop', 'annotate',
                                    'alpha'],
                              format_=None)
        vars(self).update(mykwargs)

        self.colors = topcolor, botcolor
        if fast:
            self.colors = 'white', 'white'
        # determine if data is displayed on x-axis
        self.xaxis_data = xaxis == 'data'
        # determine what is displayed on other axis
        self.otheraxis = oa = self.xaxis_data and yaxis or xaxis
        if oa == 'date':
            self.otheraxis = oa = 'starttime'
        self.plotdate = 'time' in oa or 'ponset' in oa


        #if ax_info is None:
        lenplotinfo = len(plotinfo)
        #else:
        #    lenplotinfo = len(ax_info)
        if self.plotlabel is None:
            self.plotlabel = ('',) * lenplotinfo
        if self.plotinfowhere is None:
            self.plotinfowhere = ('right',) * lenplotinfo
        if plotinfodicts is None:
            plotinfodicts = [dict(pad=0, size=0.8) for i in range(len(self.plotinfo))] #@UnusedVariable
        if plot_stack:
            self.plotinfo += ('sum',)
            self.plotlabel += ('stack',)
            if self.xaxis_data:
                self.plotinfowhere += ('top',)
            else:
                self.plotinfowhere += ('right',)
            plotinfodicts += (dict(pad=0, size=0.8),)
        if plot_psd:
            self.plotinfo += ('psd',)
            self.plotlabel += ('PSD vs. ' + self.psd_scale,)
            if self.xaxis_data:
                self.plotinfowhere += ('top',)
            else:
                self.plotinfowhere += ('right',)
            plotinfodicts += (dict(pad=0.05, size=0.8, sharex=None, sharey=None),)
        if not (len(self.plotinfo) == len(self.plotlabel) ==
                len(self.plotinfowhere) == len(plotinfodicts)) and len(self.plotinfo) > 0:
            raise ValueError('plotinfo, plotlabel, plotinfowhere and '
                             'plotinfodicts must be tuples of same lenght')
        ### prepare stream
        if component in 'QLTZNEZRT':
            stream = stream.select(component=component)
        if filter or downsample or plotsum:
            stream = stream.copy()
        if filter:
            stream.filter2(*filter)
        if downsample:
            stream.downsample2(downsample)
        if plotsum:
            if order is None:
                stream.append(stream.simpleStack(component='all'))
            else:
                stream.append(stream.phaseStack(order=order, component='all'))
            stream[-1].data = stream[-1].data * sumscale
        N = len(stream)
        if plotsum and oa != 'num':
            stream[-1].stats[oa] = stream[-2].stats[oa] + (stream[-2].stats[oa]
                                 - stream[0].stats[oa]) / (N - 1) * 3
        self.stream = stream

        if self.annotate is True:
            self.annotate = self.stream.getHI('id')

        ### get plot data
        self.data = data = getDataWindow(stream, start, end, relative)
        if np.any(np.isnan(data)):
            log.warning('NAN Data in stream to plot!')
            self.data = data = np.nan_to_num(data)
        N2 = data.shape[0]
        start, end = getTimeIntervall(stream, start, end, relative, rel_label)
        self.t = np.linspace(start[0], end[0], len(data[0]))
        self.sampling_rate = stream[0].stats.sampling_rate

        ### get scales
        if absolutescale is not None:
            scale = absolutescale
        else:
            maxima = np.sort((np.max(data, axis=1))[:N2 - plotsum])
            if N2 > 30:
                weights = np.concatenate((np.zeros(10),
                                          np.ones(N2 - 20 - plotsum),
                                          np.zeros(10)))
            else:
                weights = np.ones(N2 - plotsum)
            scale0 = 1 / np.average(maxima, weights=weights)
            scale *= scale0
        if imshow:
            if vmax is None:
                self.vmax = np.max(np.abs(data))
            else:
                self.vmax = vmax
            if vmin is None:
                self.vmin = -self.vmax
            else:
                self.vmin = vmin
            if self.stack_lim is None:
                self.stack_lim = (self.vmin, self.vmax)
            if use_dlognorm:
                self.norm = DLogNorm(vmin=self.vmin, vmax=self.vmax, cmin=cmax, cmax=cmax)
            else:
                self.norm = None
        self.scale = scale

        ### get ax and fig
        fig_N = fig is None
        if fig_N:
            if ax:
                fig = ax.get_figure()
            else:
                fig = plt.figure()
        if ax is None:
            if len(fig.axes) > 0:
                ax = fig.axes[0]
            else:
                ax = fig.add_subplot(111)
        if len(self.plotinfo) > 0 and ax_info is None:
#            if fig_N or len(fig.axes) - 1 < len(plotinfo):
#                if publication:
#                    fig = getPublicationFigure([1. - plotinfo_width * len(plotinfo)] + [plotinfo_width] * len(plotinfo))
#                else:
#                    fig = getFigure([1. - plotinfo_width * len(plotinfo)] + [plotinfo_width] * len(plotinfo))

            ax_info = get_fig(ax, positions=self.plotinfowhere, adicts=plotinfodicts)
#            else:
#                ax_info = fig.axes[1:len(plotinfo) + 1]
        self.ax_info = ax_info
        self.ax = ax
        self.fig = fig
        ### plot
        if ax_info:
            self._plot_info()
        if not imshow:
            self._plot()
        else:
            if cmap == None:
                #import matplotlib.cm
                cmap = xcorr_cmap
            self.cmap = cmap
            self._imshow()
            ax.set_aspect('auto')
            if colorbar:
                if self.norm is None:
                    fig.colorbar(self.image, ax=ax)
                else:
                    fig.colorbar(self.image, ax=ax, ticks=self.norm.ticks(), format='%1.0e')
        if plotphases:
            plotPhases(stream, ax, plotphases)
        ### print title
        if figtitle is not None:
            figtitle = figtitle.replace('station', stream[0].stats.station)
            figtitle = figtitle.replace('component', component)
            figtitle = figtitle.replace('scale', '%f' % (1. / scale))#str(scale))
            try:
                starttime = stream[0].stats.starttime + 0.5
                figtitle = figtitle.replace('time', '%s' % starttime)
                figtitle = figtitle.replace('date', '%s' % starttime.date)
                figtitle = figtitle.replace('year', '%d' % starttime.year)
            except:
                pass

            if title_in_axis:
                ax.text(0.1, 1, figtitle, verticalalignment='top',
                        transform=ax.transAxes)
            if fancy_box:
                if box_ax:
                    box_ax = fig.axes[fancy_box]
                else:
                    box_ax = ax
                if box_trans == 'fig':
                    trans = fig.transFigure
                else:
                    trans = box_ax.transAxes
                from matplotlib import transforms
                offset = transforms.ScaledTranslation(-10. / 72, -10. / 72, fig.dpi_scale_trans)
                trans = trans + offset
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                # place a text box in upper left in axes coords
                box_ax.text(1., 1., figtitle, transform=trans, fontsize=box_fs,
                        verticalalignment='top', ha='right', bbox=props)

            if not title_in_axis and not fancy_box:
                ax.set_title(figtitle, x=title_xpos,
                             horizontalalignment=title_horalign)
                #fig.suptitle(figtitle, x=title_xpos,
                #             horizontalalignment=title_horalign)
                #fig.text(title, 0., 0.95, horizontalalignment = 'left' )

        ### set axis ticks and labels
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        temp = self.xaxis_data and ('x', 'y') or ('y', 'x')
        getattr(ax, 'set_%slim' % temp[0])([self.t[0], self.t[-1]])
        if oa == 'num' and not imshow:
            getattr(ax, 'set_%slim' % temp[1])(-2, N + 2)
        #ax_info[-1].set_ylim((-0.05, 0.05))
        if self.plotdate:
            if self.xaxis_data:
                ax.yaxis_date()
            else:
                ax.xaxis_date()
            from matplotlib.dates import DateFormatter
#            exec("ax.%saxis.set_major_formatter(DateFormatter('%s'))" %
#                 (temp[1], dateformatter))
            getattr(ax, '%saxis' % temp[1]).set_major_formatter(DateFormatter(dateformatter))
            ticks = getattr(ax, 'get_%sticks' % temp[1])()
            if len(ticks) > 20:
                ticks = ticks[::len(ticks) // 20]
                getattr(ax, 'set_%sticks' % temp[1])(ticks)
        if reverse_x:
            ax.set_xlim(ax.get_xlim()[::-1])
        if reverse_y:
            ylim_temp = ax.get_ylim()
            if ylim_temp[1] > ylim_temp[0]:
                ax.set_ylim(ax.get_ylim()[::-1])
        if minor_x:
            ax.xaxis.set_minor_locator(AutoMinorLocator())
        if minor_y:
            ax.yaxis.set_minor_locator(AutoMinorLocator())
        ### connect, save, show
        #self.fig = fig
        if connect_event and not imshow:
            print fig.canvas.mpl_connect('key_press_event', self._onKey_plot)
        elif connect_event:
            fig.canvas.mpl_connect('key_press_event', self._onKey_imshow)
#        if tight_layout:
#            fig.tight_layout()
        if save:
            fig.savefig(save)
        if show:
            fig.show()
        elif save:
            plt.close(fig)
        self.fig = fig

    def _get_offsets(self):
        if self.otheraxis == 'num':
            ret = np.arange(len(self.stream))
        else:
            ret = np.array(self.stream.getHI(self.otheraxis))
            # workaround for starttimes 0.002 seconds before midnight
            if self.otheraxis == 'starttime':
                ret += 0.01
            if self.plotdate:
                ret = np.array([ofs.toordinal() + (ofs - UTCDateTime(ofs.date)) / 24 / 3600 for ofs in ret])
        return ret

    def _imshow(self):
        if self.plotdate:
            data = _insert_zeros(self.stream, self.data)
        else:
            data = self.data

        ax = self.ax
        t = self.t
        dt = (t[1] - t[0]) / 2.
        t = [t[0] - dt, t[-1] + dt]
        if self.otheraxis == 'num':
            offsets = [-0.5, data.shape[0] - 0.5]
        else:
            offsets = self._get_offsets()
        if self.plotdate:
            dt = (offsets[-1] - offsets[0]) / (data.shape[0] - 1)
            if 1.1 > dt > 0.9 and dt != 1.:
                offsets[0] = int(offsets[0])
            dt = (offsets[-1] - offsets[0]) / (data.shape[0] - 1)
            offsets = np.array([offsets[0] - dt / 2 + 0.5, offsets[-1] + dt / 2 - 0.5])
        offsets = np.array([offsets[0], offsets[-1]])
        if self.xaxis_data:
            extent = (t[0], t[-1], offsets[0] - 0.5, offsets[-1] + 0.5)
        else:
            data = np.swapaxes(data, 0, 1)
            extent = (offsets[0], offsets[-1], t[0], t[-1])
        self.image = ax.imshow(data, cmap=self.cmap, interpolation='nearest',
                               origin='lower', extent=extent,
                               vmax=self.vmax, vmin=self.vmin, norm=self.norm)
        #plt.show()
        #ax.axis(extend)

    def _plot(self, refresh=False):
        data = self.data
        color = self.color
        ax = self.ax
        topcolor, botcolor = self.colors
        scale = self.scale
        xaxis_data = self.xaxis_data
        N = len(self.stream)
        if self.alpha is None:
            alpha = np.ones(N)
        elif isinstance(self.alpha, (int, float)):
            alpha = np.ones(N) * self.alpha
        else:
            alpha_data = np.array(self.stream.getHI(self.alpha[0]))
            alpha = ((1.*alpha_data - self.alpha[2]) /
                     (self.alpha[1] - self.alpha[2]) *
                     (self.alpha[3] - self.alpha[4]) + self.alpha[4])
            alpha[alpha > 1] = 1
        t = self.t
        plots = []
        polys = []
        # delete polys
        if refresh:
            for onepoly in self.polys:
                onepoly.remove()
        offsets = self._get_offsets()
        for i in range(N):
            offset = offsets[i]
            data2 = data[i, :] * scale + offset
            # plot lines
            if not refresh:
                if xaxis_data:
                    plots.append(self.ax.plot(t, data2, color[i % 2], alpha=alpha[i])[0])
                else:
                    plots.append(self.ax.plot(data2, t, color[i % 2], alpha=alpha[i])[0])
            else:
                # update lines
                if xaxis_data:
                    self.plots[i].set_ydata(data2)
                else:
                    self.plots[i].set_xdata(data2)
            # plot polys
            if cc.to_rgb(topcolor) != (1, 1, 1):
                if xaxis_data:
                    polys.append(ax.fill_between(t, data2, offset, linewidth=0, where=data[i, :] >= 0, facecolor=topcolor, alpha=alpha[i]))
                else:
                    polys.append(ax.fill_betweenx(t, data2, offset, linewidth=0, where=data[i, :] >= 0, facecolor=topcolor, alpha=alpha[i]))
            if cc.to_rgb(botcolor) != (1, 1, 1):
                if xaxis_data:
                    polys.append(ax.fill_between(t, data2, offset, linewidth=0, where=data[i, :] < 0, facecolor=botcolor, alpha=alpha[i]))
                else:
                    polys.append(ax.fill_betweenx(t, data2, offset, linewidth=0, where=data[i, :] < 0, facecolor=botcolor, alpha=alpha[i]))
            # plot annotations
            if self.annotate is not None:
                if xaxis_data:
                    xy = (t[1], offset)
                else:
                    xy = (offset, t[1])
                ax.annotate(self.annotate[i], xy, backgroundcolor='white')

        self.polys = polys
        if not refresh:
            self.plots = plots

    def _onKey_plot(self, event):
        if event.key == 'up':
            self.scale *= 2.
        elif event.key == 'down':
            self.scale /= 2.
        if self.plotsum:
            if event.key == 'pageup':
                self.data[-1, :] = 2. * self.data[-1, :]
            elif event.key == 'pagedown':
                self.data[-1, :] = 0.5 * self.data[-1, :]
        if event.key in ('up', 'down') or (event.key in ('pageup', 'pagedown')
                                           and self.plotsum):
            print('new scale is ' + str(self.scale))
            self._plot(refresh=True)
            self.fig.show()

    def _onKey_imshow(self, event):
        if event.key == 'up':
            self.vmax /= 2 ** 0.5
            self.vmin /= 2 ** 0.5
        elif event.key == 'down':
            self.vmax *= 2 ** 0.5
            self.vmin *= 2 ** 0.5
#        elif event.key == 's':
#            self.vmin = -self.vmax
#        elif event.key == 'd':
#            self.vmax /= 2 ** 0.5
#        elif event.key == 'c':
#            self.vmax *= 2 ** 0.5
#        elif event.key == 'f':
#            self.vmin /= 2 ** 0.5
#        elif event.key == 'v':
#            self.vmin *= 2 ** 0.5
        if self.plotsum:
            if event.key == 'pageup':
                self.data[-1, :] = 2. * self.data[-1, :]
            elif event.key == 'pagedown':
                self.data[-1, :] = 0.5 * self.data[-1, :]
        if event.key in ('up', 'down') or (event.key in ('pageup', 'pagedown')
                                           and self.plotsum):
            print('new vmax is ' + str(self.vmax))
            self.image.set_clim(vmax=self.vmax, vmin=self.vmin)
            self.stack_lim = (self.vmin, self.vmax)
            self._plot_info()
            self.fig.show()


    def _plot_info(self):
        ax_info = self.ax_info
        stream = self.stream
        usehardticks = self.usehardticks
        ax2 = None
        topcolor, botcolor = self.colors
        polys = []
        for i in range(len(ax_info)):
            #ax_info.append(fig.add_axes([0.95+(i-len(plotinfo))*plotinfo_width, 0.05, plotinfo_width, 0.9]))
            label = self.plotlabel[i]
            info = self.plotinfo[i]
            offsets = self._get_offsets()
            if 'psd' in info:
                tempdata = np.mean(self.data, axis=0)
                # check if autocorrelation
                if tempdata[len(tempdata) // 2] == 1:
                    tempdata = tempdata[len(tempdata) // 2:]
                pxx, freqs = psd(tempdata, NFFT=self.psd_prop[0], Fs=self.sampling_rate,
                scale_by_freq=self.psd_prop[1], pad_to=self.psd_prop[2])
                if self.psd_scale == 'time':
                    pxx = pxx[::-1]
                    freqs = 1. / freqs[::-1]
            if self.plotinfowhere[i] in ('right', 'left'):
                if ' ' in info:
                    alpha = 1.
                    info1, info2 = info.split()
                    ax_info[i].plot(stream.getHI(info1), offsets, '.', mec='b', ms=2, alpha=alpha)
                    ax2 = ax_info[i].twiny()
                    ax2.set_axes_locator(ax_info[i].get_axes_locator())
                    ax2.plot(stream.getHI(info2), offsets, '.r', mec='r', ms=2, alpha=alpha)
                    ax2.xaxis.set_major_locator(MaxNLocator(3))
                    for tl in ax2.get_xticklabels():
#                        from IPython import embed
#                        embed()
#                        if tl.get_text() == '30':
#                            tl.set_ha('left')
                        tl.set_color('r')
                        tl.set_alpha(alpha)
                    for tl in ax_info[i].get_xticklabels():
                        tl.set_color('b')
                        tl.set_alpha(alpha)
                    if '/' in label:
                        label1, label2 = label.split('/')
                        ax_info[i].set_xlabel(label1, color='b', alpha=alpha)
                        ax2.set_xlabel(label2, color='r', alpha=alpha)
                    try:
                        self.twin_axes.append(ax2)
                    except AttributeError:
                        self.twin_axes = [ax2]
                #elif 'azi' in info:
                #    ax_info[i].plot(stream.getHI(info), range(N), lw=1.5)
                    # plot difference between azi and lazi
                    # ax_info[i].plot(abs(np.array(stream.getHI('azi'))-np.array(stream.getHI('lazi'))), range(N))
                elif 'time' in info or 'onset' in info:
                    ax_info[i].plot([UTC2year(j) for j in stream.getHI(info)], offsets, '.', mec='b', ms=2)
                    #ax_info[i].plot_date([j.datetime for j in stream.getHI(info)], range(N))
                elif 'sum' in info or 'mean' in info:
                    temp_mean = np.mean(self.data, axis=0)
                    ax_info[i].plot(temp_mean, self.t, lw=1.5)
                    if self.stack_lim is not None:
                        try:
                            ax_info[i].set_xlim(self.stack_lim)
                        except:
                            mini = np.min(temp_mean)
                            ax_info[i].set_xlim(1.1 * mini, -1.5 * mini)

                elif 'psd' in info:
                    ax_info[i].loglog(pxx, freqs)
                elif self.xaxis_data:
                    ax_info[i].plot(stream.getHI(info), offsets, '.', mec='b', ms=2)
                else:
                    ax_info[i].plot(calculate(self.data, info), self.t, lw=1.5)
                if '/' not in label:
                    ax_info[i].set_xlabel(label)
                if not 'psd' in info:
                    plt.setp(ax_info[i].get_yticklabels(), visible=False)
                ax_info[i].xaxis.grid(color='gray', linestyle='dashed')
                if ' ' in info:
                    ax_info[i].set_xlim(0, 360)
                    ax_info[i].set_xticks([0, 90, 180, 270, 360])
                    ax_info[i].set_xticklabels(['0', '', '180', '', '360'])
                    ax2.set_xlim(30, 90)
                    ax2.set_xticks([30, 45, 60, 75, 90])
                    ax2.set_xticklabels(['30', '', '60', '', '90'])
#                    for tl in ax2.get_xticklabels():
#                        if tl.get_text() == '30':
#                            tl.set_ha('left')
                elif 'azi' in info and 'azi' in usehardticks:
                    ax_info[i].set_xlim(0, 360)
                    ax_info[i].set_xticks([0, 90, 180, 270])
                    ax_info[i].set_xticklabels(['0', '', '180', ''])
                elif info == 'dist' and 'dist' in usehardticks:
                    ax_info[i].set_xlim(0, 100)
                    ax_info[i].set_xticks([0, 25, 50, 75])
                    ax_info[i].set_xticklabels(['0', '', '50', ''])
                elif 'time' in info and 'time' in usehardticks:
                    yeardata = ax_info[i].get_lines()[0].get_xdata()
                    from math import ceil
                    years = range(int(ceil(yeardata[0])), int(yeardata[-1]) + 1)
                    while len(years) > 5:
                        years = years[::2]
                    if len(years) > 0:
                        ax_info[i].set_xticks(years)
                        ax_info[i].set_xticklabels([str(year) for year in years])
                    labels = ax_info[i].get_xticklabels()
                    for label in labels:
                        label.set_rotation(70)
                    ax_info[i].set_xlim((yeardata[0] - 0.3, yeardata[-1] + 0.3))
                #else:
                elif not 'psd' in info:
                    ax_info[i].xaxis.set_major_locator(MaxNLocator(3))
            else:
                ########### start else
                if ' ' in info:
                    info1, info2 = info.split()
                    ax_info[i].plot(offsets, stream.getHI(info1), '.', mec='b', ms=2)
                    ax2 = ax_info[i].twiny()
                    ax2.plot(offsets, stream.getHI(info2), '.', mec='r', ms=2)
                    ax2.yaxis.set_major_locator(MaxNLocator(3))
                    for tl in ax2.get_yticklabels():
                        tl.set_color('r')
                    for tl in ax_info[i].get_yticklabels():
                        tl.set_color('b')
                    if '/' in label:
                        label, label2 = label.split('/')
                        ax2.set_ylabel(label2)
                    try:
                        self.twin_axes.append(ax2)
                    except NameError:
                        self.twin_axes = [ax2]
                elif 'time' in info or 'onset' in info:
                    ax_info[i].plot(offsets, [UTC2year(j) for j in stream.getHI(info)], '.', mec='b', ms=2)
                elif 'sum' in info or 'mean' in info:
                    temp_mean = np.mean(self.data, axis=0)
                    ax_info[i].plot(self.t, temp_mean, 'k', lw=1.5)
                    # plot polys
                    if cc.to_rgb(topcolor) != (1, 1, 1):
                        polys.append(ax_info[i].fill_between(self.t, temp_mean, 0, linewidth=0, where=temp_mean >= 0, facecolor=topcolor))
                    if cc.to_rgb(botcolor) != (1, 1, 1):
                        polys.append(ax_info[i].fill_between(self.t, temp_mean, 0, linewidth=0, where=temp_mean < 0, facecolor=botcolor))
                    if self.stack_lim is not None:
                        try:
                            ax_info[i].set_ylim(self.stack_lim)
                        except:
                            mini = np.min(temp_mean)
                            ax_info[i].set_ylim(1.1 * mini, -1.5 * mini)
                elif 'psd' in info:
                    ax_info[i].loglog(freqs, pxx)
                elif not self.xaxis_data:
                    if 'count' in info:
                        ax_info[i].bar(offsets, stream.getHI(info), width=np.min(np.abs(offsets - np.roll(offsets, 1))), align='center', color='b')
                    else:
                        ax_info[i].plot(offsets, stream.getHI(info), '.', mec='b', ms=2)
                else:
                    ax_info[i].plot(self.t, calculate(self.data, info), lw=1.5)

                ax_info[i].set_ylabel(label)
                ax_info[i].yaxis.grid(color='gray', linestyle='dashed')
                if not 'psd' in info:
                    plt.setp(ax_info[i].get_xticklabels(), visible=False)
                if 'azi' in info and 'azi' in usehardticks:
                    ax_info[i].set_ylim(0, 360)
                    ax_info[i].set_yticks([0, 90, 180, 270])
                    ax_info[i].set_yticklabels(['0', '', '180', ''])
                elif info == 'dist' and 'dist' in usehardticks:
                    ax_info[i].set_ylim(0, 100)
                    ax_info[i].set_yticks([0, 25, 50, 75])
                    ax_info[i].set_yticklabels(['0', '', '50', ''])
                elif 'time' in info and 'time' in usehardticks:
                    yeardata = ax_info[i].get_lines()[0].get_ydata()
                    from math import ceil #@Reimport
                    years = range(int(ceil(yeardata[0])), int(yeardata[-1]) + 1)
                    while len(years) > 5:
                        years = years[::2]
                    if len(years) > 0:
                        ax_info[i].set_yticks(years)
                        ax_info[i].set_yticklabels([str(year) for year in years])
                    labels = ax_info[i].get_yticklabels()
                    for label in labels:
                        label.set_rotation(70)
                elif not 'psd' in info:
                #else:
                    ax_info[i].yaxis.set_major_locator(MaxNLocator(3))

                #### end else
            if not 'psd' in info:
                if self.minor_x:
                    ax_info[i].xaxis.set_minor_locator(AutoMinorLocator())
                if self.minor_y:
                    ax_info[i].yaxis.set_minor_locator(AutoMinorLocator())
            ax2 = None

doc_plots = """
Plot stream with special default kwargs (depends on function).

See source code of this function for default kwargs.
Kwargs are passed to Plot.__init__.
Returns instance of sito.imaging.Plot.

Here is the doc of sito.imaging.Plot.__init__:
"""

@add_doc(doc_plots, Plot.__init__)
def plotRF(stream, *args, **kwargs_in):
    """
    Plot RFs in one axis in time window with filling.
    """
    kwargs = dict(relative='ponset', component='Q',
                  topcolor='black', botcolor='grey',
                  xlabel='time (s)',
                  figtitle='station component',
                  fancy_box=True,
                  plotinfo=('sum', 'azi dist'),
                  plotlabel=('mean', u'azi (°)/dist (°)'),
                  plotinfowhere=('top', 'right'),
                  plotinfodicts=[dict(pad=0, size=0.4), dict(pad=0.1, size=0.8)])
    keys = kwargs_in.keys()
    if 'plotinfo' in keys and 'plotlabel' not in keys:
        kwargs_in['plotlabel'] = None
    if 'plotinfo' in keys and 'plotinfodicts' not in keys:
        kwargs_in['plotinfodicts'] = [dict(pad=0, size=0.4), ]
    kwargs.update(kwargs_in)
    return Plot(stream, *args, **kwargs)

@add_doc(doc_plots, Plot.__init__)
def plotProfile(stream, *args, **kwargs_in):
    """
    Plot RFs as lon profile profile.
    """
    kwargs = dict(relative='ponset', component='Q',
              xaxis='plon', yaxis='data',
              reverse_y=True, scale=0.02,
              figtitle='profile station component sc:scale',
              xlabel=u'longitude (°)', ylabel='time (s)',
              topcolor='red', botcolor='blue',
              plotinfo=('count',), plotinfowhere=('top',),
              plotlabel=('count',)
              )
    kwargs.update(kwargs_in)
    return Plot(stream, *args, **kwargs)

@add_doc(doc_plots, Plot.__init__)
def plotTrace(trace, *args, **kwargs_in):
    """
    Plot trace divided in 24 parts with alternating colors.
    """
    kwargs = dict(num=24, absolutescale=0.0005, topcolor='white', botcolor='white',
                  color='rb', demean=True)
    if kwargs_in.has_key('scale'):
        kwargs.pop('absolutescale')
    kwargs.update(kwargs_in)
    num = kwargs.pop('num')
    demean = kwargs.pop('demean')

    from sito import Stream
    st = Stream()
    start = trace.stats.starttime
    delta = (trace.stats.endtime - start) / num
    for i in range(num):
        st.append(trace.slice(start + delta * i, start + delta * (i + 1)))
    if demean:
        for tr in st:
            tr.data = tr.data - np.mean(tr.data)
    return Plot(st, *args, **kwargs)

@add_doc(plotTrace)
def plotTrace2(data, day, station, component='Z', raw=True, **kwargs):
    """
    Plot first trace of a day file.
    
    @param data: data.Data object
    @param day: UTCDateTime object
    @param station: string
    @param component: string
    @param raw: bool
    """

    day = UTCDateTime(day)
    if raw:
        stream = data.getRawStream(day, station, component=component)
    else:
        stream = data.getStream(day, station, component=component)
    return plotTrace(stream, component=component, **kwargs)


@add_doc(doc_plots, Plot.__init__)
def plotComb(stream, *args, **kwargs_in):
    kwargs = dict(relative='middle',
              yaxis='date',
              figtitle='xcorr station',
              xlabel='lag time (s)', ylabel='xcorr',
              plotinfo=('sum',), plotinfowhere=('top',), plotlabel=('stack',))
    kwargs.update(kwargs_in)
    fig = kwargs.pop('fig', None)
    if fig is None:
        fig = plt.figure()
        ax1 = fig.add_axes((0.1, 0.1, 0.4, 0.8))
        ax2 = fig.add_axes((0.5, 0.1, 0.4, 0.8), sharex=ax1, sharey=ax1)
        ax3 = fig.add_axes((0.92, 0.1, 0.01, 0.8))
    else:
        ax1, ax2 = fig.axes[0], fig.axes[1]
    pl2 = Plot(stream, *args, ax=ax2, imshow=True, colorbar=False, **kwargs)
    pl1 = Plot(stream, *args, ax=ax1, imshow=False, **kwargs)
    for tl in ax2.yaxis.get_ticklabels(): #set_visible(False) #set_yticklabels([], axes=ax2) #visible(False)
        tl.set_visible(False)
    pl2.fig.colorbar(pl2.image, ax3)
    return pl1, pl2

@add_doc(doc_plots, Plot.__init__)
def plotXcorr(stream, *args, **kwargs_in):
    """
    Plot crooss correlations.
    """
    kwargs = dict(relative='middle',
              yaxis='date',
              imshow=True,
              figtitle='xcorr station',
              xlabel='lag time (s)', ylabel='xcorr',
              plot_stack=True)
    kwargs.update(kwargs_in)
    return Plot(stream, *args, **kwargs)

def plotXcorrVsDist(stream, *args, **kwargs_in):
    kwargs = dict(relative='middle',
              yaxis='dist',
              figtitle='xcorr vs dist',
              xlabel='lag time (s)', ylabel='dist (km)')
    kwargs.update(kwargs_in)
    return Plot(stream, *args, **kwargs)



def plotFilterbands(stream, *args, **kwargs_in):
    """
    """
    filterbands = kwargs_in.pop('filterbands', ((0.01, 0.1), (0.1, 0.01)))
    fig = kwargs_in.pop('fig', None)
    show = kwargs_in.pop('show', False)
    save = kwargs_in.pop('save', False)
    kwargs = dict(show=False, save=False, connect_event=False)
    kwargs.update(kwargs_in)

    if filterbands is None:
        filterbands = []
    else:
        filterbands = list(filterbands)
    if len(filterbands) > 0 and not cbook.iterable(filterbands[0]):
            for i in range(len(filterbands) - 1):
                filterbands[i] = (filterbands[i], filterbands[i + 1])
    filterbands.insert(0, (None, None))
    if fig is None:
        fig = plt.figure()

    num_plots = len(filterbands) + 1
    for i, filterband in enumerate(filterbands):
        ax = fig.add_subplot(num_plots, 1, i)
        stream_filter = stream.copy().filter2(*filterband)
        kwargs.update({'ax':ax})
        stream_filter.plot_(*args, **kwargs)
        kwargs.pop('figtitle')

    if save:
        fig.savefig(save)
        if show:
            fig.show()
        elif save:
            plt.close(fig)



@deprecated
def plot2(stream, start=None, end=None, relative='ponset', scale0=1, plotphases=False, show=True):
    """
    Use Plot class instead
    Plot always 3 traces in one axis

    Assuming three components in a row. The order is Z,N,E,Z,N,E, etc.
    Plotting around ponset in time window.
    """

    N = len(stream)
    data = getDataWindow(stream, start, end, relative)
    time = getTimeIntervall(stream, start, end, relative, 'relative')

    tnew = np.linspace(time[0][0], time[1][0], len(data[0]))
    maxima = np.sort(np.max(data, axis=1))
    if N > 30: weights = np.concatenate((np.zeros(10), np.ones(N - 20), np.zeros(10)))
    else: weights = np.ones(N)
    scale2 = 1 / np.average(maxima, weights=weights)
    #scale2 = np.average(maxima)
    scale = scale0 * scale2

    fig = plt.figure()
    fig.suptitle(scale)
    ax_list = []
    Nax = N // 3
    for i in range(N):
        if i % 3 == 0:
            if i == 0:
                ax = fig.add_axes([0.05, 0.05, 0.9, 0.15])
                ax.set_xlabel('time ponset')
            else: ax = fig.add_axes([0.05, 0.20 + 0.75 * (i / 3 - 1) / (Nax - 1), 0.9, 0.75 / (Nax - 1)], sharex=ax_list[0], sharey=ax_list[0])
        ax.plot([tnew[0], tnew[-1]], [0, 0], 'k:')
        ax.plot(tnew, data[i, :] * scale)
        if i % 3 == 2:
            if i != 2:
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
            ax_list.append(ax)
    ax.set_ylim(-2, 2)
    if plotphases: plotPhases(stream, ax_list, plotphases)
    if show: fig.show()
    return fig

@deprecated
def compareRF(streamlist, start=None, end=None, relative='ponset', component='Q', scale0=1, show=True, topcolor='black', botcolor='gray', sumtrcs=True):
    """
    Use Plot class instead
    
    Plot RFs in one axis in time window with filling.
    """

    #mpl.rc('keymap', fullscreen='', home='', back='', forward='', pan='', zoom='',
    #       save='', grid='', yscale='', xscale='', all_axes='')
    global data
    global scale
    global plots
    global polys
    global ax

    fig = plt.figure()
    fig.suptitle('%s_%s ' % (streamlist[0][0].stats.station, component))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    ax.set_xlabel('time (s)')
    plots = []
    colors = ['black', 'blue', 'red', 'green', 'gray', 'orange']


    for i, stream in enumerate(streamlist):
        if component != '':
            stream = stream.select(component=component)
        N = len(stream)

        data = getDataWindow(stream, start, end, relative)
        time = getTimeIntervall(stream, start, end, relative, 'relative')

        tnew = np.linspace(time[0][0], time[1][0], len(data[0]))
        maxima = np.sort(np.max(data, axis=1))
        if N > 30: weights = np.concatenate((np.zeros(10), np.ones(N - 20), np.zeros(10)))
        else: weights = np.ones(N)
        scale2 = 1 / np.average(maxima, weights=weights)
        scale = scale0 * scale2

        for j in range(N):
            # ax.plot([tnew[0], tnew[-1]], [0, 0], 'k:')
            plots.append(ax.plot(tnew, data[j, :] * scale + j, colors[i])[0])
        if sumtrcs:
            global datasum
            global sumplots
            global sumscale
            sumdata = sum(data)
            sumscale = 1 / max(sumdata) * scale0 * 2
            sumplots = []
            sumplots.append(ax.plot(tnew, sumdata * sumscale + N + 3, colors[i])[0])

    def onKey(event):
        global data
        global scale
        global plots
        global polys
        global ax
        if event.key == 'up':
            scale = scale * 2.
        elif event.key == 'down':
            scale = scale / 2.
        if event.key == 'up' or event.key == 'down':
            print('new scale is ' + str(scale))
            for i, oneplot in enumerate(plots):
                oneplot.set_ydata(data[i, :] * scale + i)
            for onepoly in polys:
                onepoly.remove()
        if sumtrcs:
            global datasum
            global sumplots
            global sumscale
            if event.key == 'pageup':
                sumscale *= 2.
            elif event.key == 'pagedown':
                sumscale /= 2.
            if event.key == 'pageup' or event.key == 'pagedown':
                print('new sumscale is ' + str(sumscale))
                sumplots[0].set_ydata(sumdata * sumscale + N + 3)
                if cc.to_rgb(botcolor) != (1, 1, 1): sumplots[-1].remove()
                if cc.to_rgb(topcolor) != (1, 1, 1): sumplots[1].remove()
                sumplots = sumplots[0:1]
                tnew = plots[0].get_xdata()
                if cc.to_rgb(topcolor) != (1, 1, 1): sumplots.append(ax.fill_between(tnew, sumdata * sumscale + N + 3, N + 3, linewidth=0, where=sumdata >= 0, facecolor=topcolor))
                if cc.to_rgb(botcolor) != (1, 1, 1): sumplots.append(ax.fill_between(tnew, sumdata * sumscale + N + 3, N + 3, linewidth=0, where=sumdata < 0, facecolor=botcolor))
        plt.draw()

    fig.canvas.mpl_connect('key_press_event', onKey)
    ax.set_xlim([tnew[0], tnew[-1]])
    ax.set_ylim(-2, N + 1 + 4 * sumtrcs)
    #ax_info = []
    if show: fig.show()
    return fig

@deprecated
def __plot3_old(stream, start=None, end=None, relative='ponset', scale0=1, absolutescale=None, plotphases=False, show=True, #delete=False,
          topcolor='gray', botcolor='white', publication=False, color='kk', zero='left', ax=None, plotdate=False):
    """
    Use Plot class instead
    Plot all traces in one axis around ponset in time window with filling.
    """

    if not stream[0].stats.get(relative):
        relative = 'starttime'
    global data
    global scale
    global plots
    global polys
    global axis

    N = len(stream)
    data = getDataWindow(stream, start, end, relative)
    time = getTimeIntervall(stream, start, end, relative, 'relative')

    if zero == 'left':
        tnew = np.linspace(time[0][0], time[1][0], len(data[0]))
    else:
        tnew = np.linspace(time[0][0], time[1][0], len(data[0])) - (time[1][0] + time[0][0]) / 2

    if absolutescale != None:
        scale = absolutescale
    else:
        maxima = np.sort(np.max(data, axis=1))
        if N > 30: weights = np.concatenate((np.zeros(10), np.ones(N - 20), np.zeros(10)))
        else: weights = np.ones(N)
        scale2 = 1 / np.average(maxima, weights=weights)
        scale = scale0 * scale2
    #scale2 = 1/np.average(maxima)

    if ax == None:
        if publication:
            fig = getPublicationFigure()
        else:
            fig = getFigure()
        axis = ax = fig.axes[0]
        fig.suptitle(scale)#'scale %5.2f' % scale)
        #ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    else:
        fig = None
        ax.text(0.5, 1, scale, verticalalignment='top', transform=ax.transAxes)
    ax.set_xlabel('time (s)')
    plots = []
    polys = []
    for i in range(N):
        # ax.plot([tnew[0], tnew[-1]], [0, 0], 'k:')
        offset = i
        if plotdate:
            offset = stream[i].stats.starttime.toordinal()
        data2 = data[i, :] * scale + offset
        if not plotdate:
            plots.append(ax.plot(tnew, data2, color[i % 2])[0])
        else:
            plots.append(ax.plot_date(tnew, data2, color[i % 2], xdate=False, ydate=True)[0])
        if cc.to_rgb(topcolor) != (1, 1, 1):
            polys.append(ax.fill_between(tnew, data2, offset, linewidth=0, where=data[i, :] >= 0, facecolor=topcolor))
        if cc.to_rgb(botcolor) != (1, 1, 1):
            polys.append(ax.fill_between(tnew, data2, offset, linewidth=0, where=data[i, :] < 0, facecolor=botcolor))

    if plotphases: plotPhases(stream, ax, plotphases)
#    if delete:
#        def mouseClick (event):
#            if event.button == 3: # right mouse button
#                global deltrcs, delindex
#                deltrcs = True
#                index = int(event.ydata / 3)
#                if index not in delindex:
#                    delindex.append(index)
#                print('Select index %s' % delindex)
#        fig.canvas.mpl_connect('button_press_event', mouseClick)
    def onKey(event):
        global data
        global scale
        global plots
        global polys
        global axis
        if event.key == 'up':
            scale = scale * 2.
        elif event.key == 'down':
            scale = scale / 2.
        if event.key == 'up' or event.key == 'down':
            print('new scale is ' + str(scale))
            for i, oneplot in enumerate(plots):
                offset = i
                if plotdate:
                    offset = stream[i].stats.starttime.toordinal()
                oneplot.set_ydata(data[i, :] * scale + offset)
            for onepoly in polys:
                onepoly.remove()
            polys = []
            tnew = plots[0].get_xdata()
            for i in range(len(plots)):
                offset = i
                if plotdate:
                    offset = stream[i].stats.starttime.toordinal()
                if cc.to_rgb(topcolor) != (1, 1, 1):
                    polys.append(axis.fill_between(tnew, data[i, :] * scale + offset, offset, linewidth=0, where=data[i, :] >= 0, facecolor=topcolor))
                if cc.to_rgb(botcolor) != (1, 1, 1):
                    polys.append(axis.fill_between(tnew, data[i, :] * scale + offset, offset, linewidth=0, where=data[i, :] < 0, facecolor=botcolor))
            plt.draw()

    ax.set_xlim([tnew[0], tnew[-1]])
    if not plotdate:
        ax.set_ylim(-2, N + 1)
    else:
        from matplotlib.dates import DateFormatter
        ax.yaxis.set_major_formatter(DateFormatter('%y-%m-%d'))
    if fig != None:
        fig.canvas.mpl_connect('key_press_event', onKey)
        if show: fig.show()
        return fig

@deprecated
def __plotProfile_old(stream, start=None, end=None, relative='ponset', scale0=2,
                    fig=None, show=True, save=False, publication=False, topcolor='red', botcolor='blue', ylabel='time (s)', xlabel=None):
    """
    Use Plot class instead
    Plot RFs as profile.
    """

    # use global variables for event handling
    global data
    global scale
    global plots
    global polys
    global ax
    N = len(stream)
    # get data of given time interval
    data = getDataWindow(stream, start, end, relative)
    # get start and end times of intervall
    time = getTimeIntervall(stream, start, end, relative, 'relative')

    tnew = np.linspace(time[0][0], time[1][0], len(data[0]))
    scale = scale0

    if not fig:
        if publication:
            fig = getPublicationFigure()
        else:
            fig = getFigure()
    title = 'Profile:'
    for tr in stream:
        title += ' ' + tr.stats.station
    fig.suptitle(title)

    ax = fig.axes[0]
    plots = []
    polys = []
    for i in range(N):
        plots.append(ax.plot(data[i, :] * scale + i, tnew, 'k')[0])
    if cc.to_rgb(botcolor) != (1, 1, 1):
        for i in range(N):
            polys.append(ax.fill_betweenx(tnew, data[i, :] * scale + i, i, linewidth=0, where=data[i, :] < 0, facecolor=botcolor))
    if cc.to_rgb(topcolor) != (1, 1, 1):
        for i in range(N):
            polys.append(ax.fill_betweenx(tnew, data[i, :] * scale + i, i, linewidth=0, where=data[i, :] >= 0, facecolor=topcolor))
    def onKey(event):
        global data
        global scale
        global plots
        global polys
        global ax
        if event.key == 'p':
            fig.savefig(save)
        if event.key == 'up':
            scale = scale * 2.
        elif event.key == 'down':
            scale = scale / 2.
        if event.key == 'up' or event.key == 'down':
            print('new scale is ' + str(scale))
            for i, oneplot in enumerate(plots):
                oneplot.set_xdata(data[i, :] * scale + i)
            for onepoly in polys:
                onepoly.remove()
            polys = []
            tnew = plots[0].get_ydata()
            if cc.to_rgb(botcolor) != (1, 1, 1):
                for i in range(N):
                    polys.append(ax.fill_betweenx(tnew, data[i, :] * scale + i, i, linewidth=0, where=data[i, :] < 0, facecolor=botcolor))
            if cc.to_rgb(topcolor) != (1, 1, 1):
                for i in range(N):
                    polys.append(ax.fill_betweenx(tnew, data[i, :] * scale + i, i, linewidth=0, where=data[i, :] >= 0, facecolor=topcolor))
        plt.draw()

    fig.canvas.mpl_connect('key_press_event', onKey)
    ax.set_ylim([tnew[-1], tnew[0]])
    ax.set_ylabel(ylabel)
    ax.set_xlim(-2, N + 1)
    if xlabel:
        ax.set_xlabel(xlabel)
    if save: fig.savefig(save)
    if show: fig.show()
    return fig

class Farm(object):
    """
    Use RFTool instead.
    
    Class for selecting events

    Key:
    left/right: select last/next 3 traces
    home/end: select the 3 traces at the beginning/end of stream
    delete: delete the 3 traces
    backspace: undo (50 entries)
    o: restore original stream
    p: turn on/off of plotting rf
    2: plot2
    3: plot3
    a/y: change waterlevel by factor 3.33
    s/x: change Gauss filter by 0.5Hz
    r: rotate back to ZNE and rotate again to LQT with fresh calculated polarization
    f: apply HP filter
    s: save stream to temp file
    any key: update plot
    h: this help

    use self.ms.write('TEST', 'Q') to save stream
    use self.trcs.filter2(minfreq, maxfreq) to filter the 3 traces
    original stream saved in self.ms_orig
    """
#    3: plot3 (you can select indices with right mouse click there and
#       delete them here with key delete

    def __init__(self, ms, show2=False, show3=False, log=False, winsig=(0., 10.), winnoise=(-20., -10.)):
        from warnings import warn
        warn('Use RfTool instead.', category=DeprecationWarning)
        from matplotlib.widgets import Button, Slider
        #mpl.rcParams['lines.linewidth'] = 2
        #mpl.rcParams['lines.color'] = 'r'
        mpl.rc('keymap', fullscreen='', home='', back='', forward='', pan='', zoom='',
               save='', grid='', yscale='', xscale='', all_axes='')
        self.log = logging.getLogger('Farm')
        if log == True:
            logging.basicConfig(level=logging.DEBUG, format='%(module)s.%(funcName)s - %(levelname)s - %(message)s')
        elif log:
            logging.basicConfig(level=logging.DEBUG, filename=log, format='%(module)s.%(funcName)s - %(levelname)s - %(message)s')
        self.ms = ms
        self.ms_orig = ms.copy()
        self.undo = []
        self.trcs = ms[0:3]
        self.i = 0
        self.plotrf = True
        self.smin = 2.
        self.smax = 3.
        self.wl = 0.01
        self.gauss = 4.
        self.fig2 = None
        self.fig3 = None
        self.winsig = winsig
        self.winnoise = winnoise
        try:
            self.ms[0].stats.signoise
        except KeyError:
            self.ms.signoise(self.winsig, self.winnoise)
        if show2:
            self.fig2 = plot2(ms, show=False)
            self.fig2.show()
        if show3:
            self.fig3 = __plot3_old(ms, show=False, delete=True)
            self.fig3.show()
        fig1 = plt.figure(figsize=(14, 8))
        self.ax1 = fig1.add_subplot(611, autoscale_on=False)
        self.ax2 = fig1.add_subplot(612, sharex=self.ax1, sharey=self.ax1, autoscale_on=False)
        self.ax3 = fig1.add_subplot(613, sharex=self.ax1, sharey=self.ax1, autoscale_on=False)
        self.ax4 = fig1.add_subplot(614, autoscale_on=False)
        self.ax5 = fig1.add_subplot(615, sharex=self.ax4, sharey=self.ax4, autoscale_on=False)
        self.ax6 = fig1.add_subplot(616, sharex=self.ax4, sharey=self.ax4, autoscale_on=False)
        self.plot1, = self.ax1.plot([], [])
        self.plot2, = self.ax2.plot([], [])
        self.plot3, = self.ax3.plot([], [])
        self.plot4, = self.ax4.plot([], [])
        self.plot5, = self.ax5.plot([], [])
        self.plot6, = self.ax6.plot([], [])
        self.ax4.axhline(y=0, xmin=0, xmax=1)
        self.ax5.axhline(y=0, xmin=0, xmax=1)
        self.ax6.axhline(y=0, xmin=0, xmax=1)
        axbut1 = fig1.add_axes([0.01, 0.01, 0.08, 0.03])
        axbut2 = fig1.add_axes([0.11, 0.01, 0.08, 0.03])
        axbut3 = fig1.add_axes([0.21, 0.01, 0.04, 0.03])
        axbut4 = fig1.add_axes([0.26, 0.01, 0.04, 0.03])
        axbut5 = fig1.add_axes([0.31, 0.01, 0.04, 0.03])
        but1 = Button(axbut1, 'auto seis')
        but2 = Button(axbut2, 'auto rf')
        but3 = Button(axbut3, 'LP')
        but4 = Button(axbut4, 'HP')
        but5 = Button(axbut5, 'BP')
        but1.on_clicked(self.but1Push)
        but2.on_clicked(self.but2Push)
        but3.on_clicked(self.but3Push)
        but4.on_clicked(self.but4Push)
        but5.on_clicked(self.but5Push)
        axsl1 = fig1.add_axes([0.45, 0.01, 0.15, 0.03])
        axsl2 = fig1.add_axes([0.75, 0.01, 0.15, 0.03])
        sl1 = Slider(axsl1, 'sec min(LP)', 0.1, 10.0, valinit=2.)
        sl2 = Slider(axsl2, 'max(HP)', 1., 20.0, valinit=3.)
        def slUpdate(value): #@UnusedVariable
            self.smin = sl1.val
            self.smax = sl2.val
        sl1.on_changed(slUpdate)
        sl2.on_changed(slUpdate)
        fig1.canvas.mpl_connect('key_press_event', self.onKey)
        fig1.show()
        self.fig1 = fig1
        print(self.__doc__)
        print(self.info())
        self.update(True, True)
        from IPython import embed
        embed()

    def addUndo(self, remove=False):
        self.undo.append([remove, self.i, self.trcs.copy()])
        if len(self.undo) > 50:
            del self.undo[0]

    def onKey(self, event):
        N = len(self.ms) // 3
        if event.key == 'h':
            print('\n' + self.__doc__)
        elif event.key == 'left':
            self.i = (self.i - 1) % N
        elif event.key == 'right':
            self.i = (self.i + 1) % N
        elif event.key == 'home':
            self.i = 0
        elif event.key == 'end':
            self.i = N - 1
        elif event.key == 'delete':
            self.addUndo(True)
            print('\nDelete 3 traces.')
            self.log.info('Delete %s' % self.info())
            del self.ms[3 * self.i:3 * self.i + 3]
            self.i = self.i % (len(self.ms) // 3)
        elif event.key == 'backspace':
            try:
                un = self.undo.pop()
                if un[0]:
                    self.ms.insert(3 * un[1], un[2])
                    self.log.info('Redo delete of list index %s' % un[1])
                else:
                    self.ms[3 * un[1]:3 * un[1] + 3] = un[2]
                self.i = un[1]
            except IndexError:
                print('\nNo undo possible.')
        elif event.key == 'o':
            print('\nRestore original stream ms_orig.')
            self.ms = self.ms_orig
        elif event.key == '2':
            self.fig2 = plot2(self.ms, show=False)
            self.fig2.show()
        #elif event.key == '3':
        #    self.fig3 = plot3(self.ms, scale=sca, show=False, delete=True)
        #    self.fig3.show()
        elif event.key in 'ayqw':
            if event.key == 'a':
                self.wl = self.wl * 10 ** 0.5
            elif event.key == 'y':
                self.wl = self.wl / 10 ** 0.5
            elif event.key == 's':
                self.gauss = self.gauss + 0.5
            elif event.key == 'x':
                self.gauss = self.gauss - 0.5
            print('\nfrequency deconvolution: water level %s  gauss filter %sHz' % (self.wl, self.gauss))
        elif event.key == 'r':
            print('\nRotate to ZNE and back to LQT.')
            self.addUndo()
            self.trcs.rotateLQT2ZNE()
            self.trcs.rotateZNE2LQT([-2, 10])
        elif event.key == 'f':
            self.addUndo()
            self.trcs.filter2(1 / self.smax, 0)
            self.trcs.signoise(self.winsig, self.winnoise)
        elif event.key == 's':
            self.ms.write('TEMP', 'Q')
        elif event.key == 'p':
            self.plotrf = not self.plotrf
            print('\nPlot receiver functions: %s' % self.plotrf)
            self.plot1.remove()
            self.plot2.remove()
            self.plot3.remove()
            self.ax1.set_visible(False)
            self.ax2.set_visible(False)
            self.ax3.set_visible(False)
            self.ax1, self.ax2, self.ax3
            if self.plotrf:
                self.ax1 = self.fig1.add_subplot(611, autoscale_on=False)
                self.ax2 = self.fig1.add_subplot(612, sharex=self.ax1, sharey=self.ax1, autoscale_on=False)
                self.ax3 = self.fig1.add_subplot(613, sharex=self.ax1, sharey=self.ax1, autoscale_on=False)
                self.ax1.set_visible(True)
                self.ax2.set_visible(True)
                self.ax3.set_visible(True)
                self.ax4.set_visible(True)
                self.ax5.set_visible(True)
                self.ax6.set_visible(True)
            else:
                self.ax1 = self.fig1.add_subplot(311, autoscale_on=False)
                self.ax2 = self.fig1.add_subplot(312, sharex=self.ax1, sharey=self.ax1, autoscale_on=False)
                self.ax3 = self.fig1.add_subplot(313, sharex=self.ax1, sharey=self.ax1, autoscale_on=False)
                self.ax1.set_visible(True)
                self.ax2.set_visible(True)
                self.ax3.set_visible(True)
                self.ax4.set_visible(False)
                self.ax5.set_visible(False)
                self.ax6.set_visible(False)
            self.plot1, = self.ax1.plot([], [])
            self.plot2, = self.ax2.plot([], [])
            self.plot3, = self.ax3.plot([], [])
        self.trcs = self.ms[3 * self.i:3 * self.i + 3]
        print(self.info())
        if event.key in ('left', 'right', 'home', 'end', 'backspace', 'p', 'delete'):
            self.update(True, False)
        else:
            self.update()
##            global deltrcs
##            if deltrcs:
##                global delindex
##                deltrcs = False
##                delindex.sort(reverse=True)
##                print('\nDelete %s traces with index %s' % (3 * len(delindex), delindex))
##                self.log.info('Delete %s traces with index %s' % (3 * len(delindex), delindex))
##                for index in delindex:
##                    del self.ms[3 * index:3 * index + 3]
##            else:
##        elif event.key == 'c':
##            print('\nCopy trcs to undo.')
##            self.addUndo()


    def but1Push(self, event):
        self.update(scaleseis=True)
    def but2Push(self, event):
        self.update(scalerf=True)
    def but3Push(self, event):
        self.addUndo()
        self.trcs.filter2(0, 1 / self.smin)
        self.trcs.signoise(self.winsig, self.winnoise)
        print(self.info())
        self.update()
    def but4Push(self, event):
        self.addUndo()
        self.trcs.filter2(1 / self.smax, 0)
        self.trcs.signoise(self.winsig, self.winnoise)
        print(self.info())
        self.update()
    def but5Push(self, event):
        self.addUndo()
        self.trcs.filter2(1 / self.smax, 1 / self.smin)
        self.trcs.signoise(self.winsig, self.winnoise)
        print(self.info())
        self.update()

    def update(self, scaleseis=False, scalerf=False):
        t1 = self.trcs[0].stats.starttime - self.trcs[0].stats.ponset
        t2 = self.trcs[0].stats.endtime - self.trcs[0].stats.ponset
        dataL = self.trcs[0].data
        dataQ = self.trcs[1].data
        dataT = self.trcs[2].data
        t = np.linspace(t1, t2, len(dataL))
        self.plot1.set_data(t, dataL)
        self.plot2.set_data(t, dataQ)
        self.plot3.set_data(t, dataT)
        if scaleseis:
            maxi = max(max(abs(dataL)), max(abs(dataQ)), max(abs(dataT)))
            self.ax1.set_xlim([t1, t2])
            self.ax1.set_ylim([-maxi, maxi])

        if self.plotrf:
            rf = self.trcs.copy()
            rf.receiverf(self.wl, self.gauss)
            t1 = rf[0].stats.starttime - rf[0].stats.ponset
            t2 = rf[0].stats.endtime - rf[0].stats.ponset
            dataL = rf[0].data
            dataQ = rf[1].data
            dataT = rf[2].data
            t = np.linspace(t1, t2, len(dataL))
            self.plot4.set_data(t, dataL)
            self.plot5.set_data(t, dataQ)
            self.plot6.set_data(t, dataT)
            if scalerf:
                maxi = max(max(abs(dataL)), max(abs(dataQ)), max(abs(dataT)))
                self.ax4.set_xlim((t1, t2))
                self.ax4.set_ylim((-maxi, maxi))
        self.fig1.canvas.draw()

    def info(self):
        ind = 3 * self.i
        tr = self.ms[ind]
        string = '\nlist index tr%s / %s  components %s\n%s\nsignoiseL %.2f\n' + \
                'theo   polar\nazi  %5.1f  %5.1f\n' + \
                'inci %5.1f  %5.1f'
        formating = (self.i, len(self.ms) // 3 - 1,
                       tr.stats.channel[-1] + self.ms[ind + 1].stats.channel[-1] +
                       self.ms[ind + 2].stats.channel[-1],
                       tr.print_(1), tr.stats.signoise,
                       tr.stats.azi, tr.stats.lazi, tr.stats.inci, tr.stats.linci)
        return string % formating


#class FarmRF(object):
#    def __init__(self, stream, start=None, end=None, relative='ponset', component='Q', scale0=2, sumscale0=25,
#           fig=None, show=True, save=False, publication=False, topcolor='black', botcolor='gray',
#           showsum=True, plotinfo=['azi', 'dist'], usehardticks='', plotinfo_width=0.1,
#           plotlabel=['time (s)', u'azi (°)', u'dist (°)'], plotphases=False):
#        """Plot RFs in one axis in time window with filling."""
#
#    #mpl.rc('keymap', fullscreen='', home='', back='', forward='', pan='', zoom='',
#    #       save='', grid='', yscale='', xscale='', all_axes='')
#
#        self.stream_all, self.start, self.end, self.relative = stream, start, end, relative
#        self.botcolor, self.topcolor = botcolor, topcolor
#        self.showsum = showsum
#        if component in 'QLTZNE':
#            self.stream = self.stream_all.select(component=component)
#        else:
#            self.stream = self.stream_all
#        N = len(stream)
#
#        # get data of given time interval
#        self.data = getDataWindow(self.stream, self.start, self.end, self.relative)
#        # get start and end times of intervall
#        self.scale = scale0
#
#        if not fig:
#            if publication:
#                self.fig = getPublicationFigure([1.-plotinfo_width * len(plotinfo)] + [plotinfo_width] * len(plotinfo))
#            else:
#                self.fig = getFigure([1.-plotinfo_width * len(plotinfo)] + [plotinfo_width] * len(plotinfo))
#        self.fig.suptitle('%s %s ' % (self.stream[0].getHI('station'), component))
#        self.ax = self.fig.axes[0]
#        if len(self.fig.axes) > 1:
#            ax_info = self.fig.axes[1:]
#
#        self.ax.set_xlabel(plotlabel[0])
#        self.plots = []
#        self.polys = []
#
#
#
#        time = getTimeIntervall(self.stream, self.start, self.end, self.relative, 'relative')
#        tnew = np.linspace(time[0][0], time[1][0], len(self.data[0]))
#        for i in range(N):
#            # ax.plot([tnew[0], tnew[-1]], [0, 0], 'k:')
#            self.plots.append(self.ax.plot(tnew, self.data[i,:] * self.scale + i, 'k')[0])
#        if cc.to_rgb(self.botcolor) != (1, 1, 1):
#            for i in range(N):
#                self.polys.append(self.ax.fill_between(tnew, self.data[i,:] * self.scale + i, i, linewidth=0, where=self.data[i,:] < 0, facecolor=self.botcolor))
#        if cc.to_rgb(self.topcolor) != (1, 1, 1):
#            for i in range(N):
#                self.polys.append(self.ax.fill_between(tnew, self.data[i,:] * self.scale + i, i, linewidth=0, where=self.data[i,:] >= 0, facecolor=self.topcolor))
#        if self.showsum:
#            self.sumdata = np.mean(data, axis=0)
#            self.sumscale = sumscale0
#            self.sumplots = []
#            self.sumplots.append(self.ax.plot(self.tnew, sumdata * self.sumscale + N + 2, 'k')[0])
#            if cc.to_rgb(self.topcolor) != (1, 1, 1):
#                self.sumplots.append(self.ax.fill_between(tnew, self.sumdata * self.sumscale + N + 2, N + 2, linewidth=0, where=self.sumdata >= 0, facecolor=self.topcolor))
#            if cc.to_rgb(self.botcolor) != (1, 1, 1):
#                self.sumplots.append(self.ax.fill_between(tnew, self.sumdata * self.sumscale + N + 2, N + 2, linewidth=0, where=self.sumdata < 0, facecolor=self.botcolor))
#
#        self.fig.canvas.mpl_connect('key_press_event', self.onKey)
#
#        self.ax.set_xlim([self.tnew[0], self.tnew[-1]])
#        self.ax.set_ylim(-2, N + 1 + 4 * self.showsum)
#        #ax_info = []
#        for i, info in enumerate(plotinfo):
#            #ax_info.append(fig.add_axes([0.95+(i-len(plotinfo))*plotinfo_width, 0.05, plotinfo_width, 0.9]))
#            if 'azi' in info:
#                ax_info[i].plot(self.stream.getHI(info), range(N))
#                # plot difference between azi and lazi
#                # ax_info[i].plot(abs(np.array(stream.getHI('azi'))-np.array(stream.getHI('lazi'))), range(N))
#            elif 'time' in info:
#                ax_info[i].plot([UTC2year(j) for j in self.stream.getHI(info)], range(N))
#                #ax_info[i].plot_date([j.datetime for j in stream.getHI(info)], range(N))
#            else:
#                ax_info[i].plot(self.stream.getHI(info), range(N))
#            ax_info[i].set_xlabel(plotlabel[i + 1])
#            ax_info[i].set_ylim(-2, N + 1 + 5 * self.showsum)
#            #ax_info[i].set_yticklabels(visible=False)
#            plt.setp(ax_info[i].get_yticklabels(), visible=False)
#            ax_info[i].xaxis.grid(color='gray', linestyle='dashed')
#            if 'azi' in info and 'azi' in usehardticks:
#                ax_info[i].set_xlim(0, 360)
#                ax_info[i].set_xticks([0, 90, 180, 270])
#                ax_info[i].set_xticklabels(['0', '', '180', ''])
#            elif info == 'dist' and 'dist' in usehardticks:
#                ax_info[i].set_xlim(0, 100)
#                ax_info[i].set_xticks([0, 25, 50, 75])
#                ax_info[i].set_xticklabels(['0', '', '50', ''])
#            elif 'time' in info and 'time' in usehardticks:
#                #ax_info[i].set_xlim(0.5, 3)
#                ax_info[i].set_xticks([2007, 2008, 2009, 2010])
#                ax_info[i].set_xticklabels(['2007', '2008', '2009', '2010'])
#                labels = ax_info[i].get_xticklabels()
#                for label in labels:
#                    label.set_rotation(70)
#            else:
#                ax_info[i].xaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
#        plt.draw()
#        #if save: fig.savefig(save)
#        #if show: fig.show()
#        #return fig
#
#
#
#    def onKey(self, event):
#        global data
#        global scale
#        global plots
#        global polys
#        global ax
#        N = len(self.stream)
#        #if event.key == 'p':
#        #    self.fig.savefig(save)
#        if event.key == 'up':
#            self.scale *= 2.
#        elif event.key == 'down':
#            self.scale /= 2.
#        if event.key == 'up' or event.key == 'down':
#            print('new scale is ' + str(self.scale))
#            for i, oneplot in enumerate(self.plots):
#                oneplot.set_ydata(self.data[i,:] * self.scale + i)
#            for onepoly in self.polys:
#                onepoly.remove()
#            polys = []
#            tnew = self.plots[0].get_xdata()
#            if cc.to_rgb(self.botcolor) != (1, 1, 1):
#                for i in range(N):
#                    self.polys.append(self.ax.fill_between(tnew, self.data[i,:] * self.scale + i, i, linewidth=0, where=self.data[i,:] < 0, facecolor=self.botcolor))
#            if cc.to_rgb(self.topcolor) != (1, 1, 1):
#                for i in range(N):
#                    self.polys.append(self.ax.fill_between(tnew, self.data[i,:] * self.scale + i, i, linewidth=0, where=self.data[i,:] >= 0, facecolor=self.topcolor))
#        if self.showsum:
#            if event.key == 'pageup':
#                self.sumscale *= 2.
#            elif event.key == 'pagedown':
#                self.sumscale /= 2.
#            if event.key == 'pageup' or event.key == 'pagedown':
#                print('new sumscale is ' + str(self.sumscale))
#                self.sumplots[0].set_ydata(self.sumdata * self.sumscale + N + 2)
#                if cc.to_rgb(self.botcolor) != (1, 1, 1): self.sumplots[-1].remove()
#                if cc.to_rgb(self.topcolor) != (1, 1, 1): self.sumplots[1].remove()
#                self.sumplots = self.sumplots[0:1]
#                tnew = self.plots[0].get_xdata()
#                if cc.to_rgb(self.topcolor) != (1, 1, 1): self.sumplots.append(self.ax.fill_between(tnew, self.sumdata * self.sumscale + N + 2, N + 2, linewidth=0, where=self.sumdata >= 0, facecolor=self.topcolor))
#                if cc.to_rgb(self.botcolor) != (1, 1, 1): self.sumplots.append(self.ax.fill_between(tnew, self.sumdata * self.sumscale + N + 2, N + 2, linewidth=0, where=self.sumdata < 0, facecolor=self.botcolor))
#        plt.draw()


def test():
    try:
        from sito.stream import read
        import os
        msdec = read(os.path.join(os.path.dirname(__file__), 'tests2', 'data', 'TEST_DEC.QHD'))
    except Exception as ex:
        raise IOError('Missing File: ' + str(ex))
    msdec.sort('azi')
    dummy1 = Plot(msdec)
    dummy2 = msdec.plotRF(-10, 50)
    Farm(msdec)

def test2():
    from sito import read
#    ms = read('/home/richter/Results/IPOC/xcorr_test/day/*.QHD')
#    dummy1 = ms.plotComb(-100, 100, relative='middle', yaxis='date')
#    dummy12 = ms.plot_(-100, 100, relative='middle', yaxis='date', imshow=True, plotinfo=('sum',), plotinfowhere=('top',), plotlabel=('sum',), figtitle='station')
#    dummy13 = ms.plotComb(-100, 100, relative='middle', yaxis='date', plotinfo=('sum', 'starttime'), plotinfowhere=('top', 'left'), plotlabel=('sum', 'year'))

    ms2 = read('/home/richter/Results/IPOC/receiver/results/plots/profile_pp_dif0.02.QHD')
    dummy2 = ms2.plotProfile(plotinfo=('count', 'sum'), plotinfowhere=('top', 'right'), plotlabel=('count', 'sum'))
    ms3 = read('/home/richter/Results/IPOC/receiver/results/PB02_2009_mout.QHD')
    dummy3 = ms3.plotRF(-30, 100, plotinfo=('azi dist', 'starttime'), plotlabel=('azi/dist', 'year'), plotsum=True, plotphases=True)
    from IPython import embed
    embed()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    test()
    test2()

#mpl.rc('keymap', fullscreen='', home='', back='', forward='', pan='', zoom='',
#       save='', grid='', yscale='', xscale='', all_axes='')
