#!/usr/bin/env python
# -*- coding: utf-8 -*-

import easygui
import logging
from sito import util, data as mod_data #@UnusedImport
from sito.noisexcorr import FloatingStream
from obspy.signal.psd import PPSD
from obspy.core import UTCDateTime as UTC
from PyQt4 import QtCore
from easygui import exceptionbox, msgbox
from widgets import addConsole, addLogging, addMPLToolBar
import pickle

import sys
from sito.util import TerminalFile

signal = QtCore.SIGNAL
connect = QtCore.QObject.connect
log = logging.getLogger()

def_open = '/home/richter/Results/Parkfield/PPSD/*.pickle'
def_save_data = '/home/richter/Results/Parkfield/PPSD/*.pickle'

def ppsdTool(data, t1, t2, station, component='Z'):
    date = t1.__class__(t1.date)
    enddate = t2.__class__(t2.date)
    ppsd = None
    try:
        daystream = FloatingStream(data, date, station, 1800,
                                      component=component, reserve=20,
                                      use_get_raw=True)
    except:
        log.exception('Can not load FloatingStream use just Data'
                          ' instead')
        daystream = data

    s = TerminalFile(sys.stdin)
    print "Press q to cancel calculation..."

    while date <= enddate:
        if s.getch() == 'q':
            return ppsd
        try:
            stream = daystream.getRawStream(date, station, component)
        except ValueError:
            log.exception('')
            date += 24 * 3600
            continue
        if ppsd == None:
            ppsd = PPSD(stream[0].stats, data.paz)
        ppsd.add(stream, verbose=True)
        date += 24 * 3600
    return ppsd

# TODOni create PPSD plots showing variability

class PPSDTool(object):
    def __init__(self, mainwindow, ui_mainwindow):
        self.win = mainwindow
        self.ui = ui_mainwindow
        self.ui.setupUi(self.win)
        QtCore.pyqtRemoveInputHook() # to be called in __init__ of main-dialog-class
        addMPLToolBar(self.win, self.ui.canv)
        self.canv = self.ui.canv
        self.fig = self.canv.figure

        self.filehandler = self.streamhandler = None
        self.data = None
        self.daystream = None
        self.ppsd = None

        self.timer = QtCore.QTimer()

        connect(self.timer, signal('timeout()'), self.timerEvent)
        connect(self.ui.t1, signal('editingFinished()'), self.set_dates)
        connect(self.ui.t2, signal('editingFinished()'), self.set_dates)
        connect(self.ui.button, signal('clicked()'), self.add_ppsd)
        connect(self.ui.push_open, signal('clicked()'), self.open_ppsd)
        connect(self.ui.push_save, signal('clicked()'), self.save_ppsd)
        connect(self.ui.push_plot, signal('clicked()'), self.plot_ppsd)
        connect(self.ui.actionAdd_Console, signal('triggered()'),
                self.add_console)
        connect(self.ui.actionStart_IPython, signal('triggered()'),
                self.start_ipy)
        connect(self.ui.actionLogging_to_File, signal('triggered()'),
                self.logging_file)
        connect(self.ui.actionLogging_Window, signal('triggered()'),
                self.logging_window)

        self.set_dates()
        self.logging_window()
        log.info('Hallo')


    def add_console(self):
        addConsole(self.win, self.ui, message='self is under namesspace ns',
                    namespace={'ns':self}, multithreaded=False)
    def start_ipy(self):
        from IPython import embed
        embed()

    def open_ppsd(self, file_=None):
        if file_ == None:
            file_ = easygui.fileopenbox(default=def_open,
                                       filetypes=['*.pickle', '*.*'])
        if file_ != None:
            try:
                with open(file_) as f:
                    self.ppsd = pickle.load(f)
            except:
                exceptionbox("Error while reading file!")
            else:
                self.ui.pbar.setValue(0)
                self.ui.button.setText('Add')
                self.ui.project.setEnabled(True)
                self.ui.station.setEnabled(True)
                self.ui.component.setEnabled(True)
                self.data = None
                self.daystream = None
                self.date = self.startdate

    def save_ppsd(self, file_=None):
        if self.ppsd == None:
            msgbox('No PPSD to save.')
            return
        if file_ == None:
            file_ = easygui.filesavebox(default=def_save_data,
                                       filetypes=['*.pickle', '*.*'])
        if file_ != None:
            try:
                self.ppsd.save(file_)
            except:
                exceptionbox("Error while writing file_!")

    def set_dates(self):
        try:
            self.startdate = self.date = UTC(UTC(str(self.ui.t1.text())).date)
            self.enddate = UTC(UTC(str(self.ui.t2.text())).date)
            self.daystream = None
        except:
            log.exception('')

    def add_ppsd(self, dontstop=False):
        if not dontstop and self.timer.isActive():
            self.ui.button.setText('Continue')
            self.timer.stop()
            if self.date > self.enddate:
                self.timer.start(100)
        else:
            self.ui.pbar.setValue(100 * (self.date - self.startdate + 12 * 3600)
                                  / (self.enddate - self.startdate + 24 * 3600))
            self.ui.button.setText('Stop')
            if self.data == None:
                self.data = getattr(mod_data, self.ui.project.text())()
            if self.daystream == None:
                try:
                    self.daystream = FloatingStream(self.data, self.startdate,
                                            str(self.ui.station.text()),
                                            1800, component=
                                                str(self.ui.component.text()),
                                            reserve=20, use_get_raw=True)
                except:
                    log.exception('Can not load FloatingStream use just Data'
                                  ' instead')
                    self.daystream = self.data
            while self.date <= self.enddate:
                try:
                    stream = self.daystream.getRawStream(self.date,
                                                str(self.ui.station.text()),
                                                str(self.ui.component.text()))
                except Exception as ex:
                    log.exception(str(ex))
                    self.date += 24 * 3600
                    stream = None
                else:
                    break
            if stream == None:
                return
            if self.ppsd == None:
                stats = stream[0].stats
                paz = self.data.paz
                self.ppsd = PPSD(stats, paz)
            self.ppsd.add(stream, verbose=True)
            self.date += 24 * 3600
            self.timer.start(100)

    def timerEvent(self):
        if self.date > self.enddate:
            self.timer.stop()
            self.ui.button.setText('Add')
            self.date = self.startdate
            self.ui.pbar.setValue(100)
            return
        self.add_ppsd(dontstop=True)

    def plot_ppsd(self):
        #import numpy as np
        percentiles = [0, 50, 100]
        self.fig.clear()
        self.ppsd.plot(show_coverage=self.ui.check_plot1.isChecked(),
                       show_histogram=self.ui.check_plot2.isChecked(),
                       show_percentiles=self.ui.check_plot3.isChecked(),
                       show_noise_models=self.ui.check_plot4.isChecked(),
                       percentiles=percentiles,
                       fig=self.fig)
        self.canv.draw()

    def logging_file(self):
        file_ = easygui.filesavebox(default='/home/richter/Data/*.log')
        if file_ != None:
            if self.filehandler != None:
                log.removeHandler(self.filehandler)
            self.filehandler = logging.FileHandler(file_, 'a')
            self.filehandler.setLevel(logging.DEBUG)
            #filehandler.setFormatter(logging.Formatter('%(asctime)s
            #        - %(levelname)s - %(name)s.%(funcName)s - %(message)s'))
            log.addHandler(self.filehandler)
    def logging_window(self):
        if self.streamhandler != None:
            log.removeHandler(self.streamhandler)
        self.streamhandler = addLogging(self.win, self.ui, log)

if __name__ == "__main__":
    from PyQt4 import QtGui
    from ppsdtool_gui import Ui_MainWindow
    app = QtGui.QApplication(sys.argv)
    gui = PPSDTool(QtGui.QMainWindow(), Ui_MainWindow())
    gui.win.show()
    sys.exit(app.exec_())
