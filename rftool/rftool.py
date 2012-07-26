#!/usr/bin/env python
# -*- coding: utf-8 -*-

import easygui
import logging
from sito import imaging, read
import numpy as np
from PyQt4 import QtCore
from easygui import exceptionbox, msgbox
from widgets import addConsole, addLogging, addMPLToolBar
from IPython import embed

signal = QtCore.SIGNAL
connect = QtCore.QObject.connect
log = logging.getLogger('rftool')

import sys
from IPython.core import ultratb
sys.excepthook = ultratb.FormattedTB(mode='Verbose',
color_scheme='LightBG', call_pdb=1)

def_open = '/home/richter/Data/Parkfield/receiver/M5.5_events/*.QHD'
def_save_data = '/home/richter/Data/Parkfield/receiver/M5.5_events_handmarked/*'
def_save_rf = '/home/richter/Data/Parkfield/receiver/M5.5_results_theorot/*.QHD'
def_sort = 'ponset,station,component'
class RfTool(object):
    def __init__(self, mainwindow, ui_mainwindow):
        self.win = mainwindow
        self.ui = ui_mainwindow
        self.ui.setupUi(self.win)
        QtCore.pyqtRemoveInputHook() # to be called in __init__ of main-dialog-class
        addMPLToolBar(self.win, self.ui.canv)
        self.canv = self.ui.canv
        self.fig = self.canv.figure
        self.ax1 = self.fig.add_axes([0.03, 0.04, 0.46, 0.91])
        self.ax1_sum = self.fig.add_axes([0.03, 0.95, 0.46, 0.04], sharex=self.ax1)
        self.ax2 = self.fig.add_axes([0.53, 0.04, 0.46, 0.91])
        self.ax2_sum = self.fig.add_axes([0.53, 0.95, 0.46, 0.04])

        self.fill_data = self.fill_rf = 'white-white'
        self.sort_after = def_sort
        self.st_all_old = self.st_all = self.st = self.st2 = None
        self.filehandler = self.streamhandler = None


        connect(self.ui.actionOpen_Data, signal('triggered()'), self.open_data)
        connect(self.ui.actionSave_Data, signal('triggered()'), self.save_data)
        connect(self.ui.actionSave_RF, signal('triggered()'), self.save_rf)
        connect(self.ui.actionFill_Data, signal('triggered()'), self.fill_data_triggered)
        connect(self.ui.actionFill_RF, signal('triggered()'), self.fill_rf_triggered)
        connect(self.ui.actionAdd_Console, signal('triggered()'),
                self.add_console)
        connect(self.ui.actionStart_IPython, signal('triggered()'),
                self.start_ipy)
        connect(self.ui.actionLogging_to_File, signal('triggered()'),
                self.logging_file)
        connect(self.ui.actionLogging_Window, signal('triggered()'),
                self.logging_window)

        connect(self.ui.combo_filter, signal('currentIndexChanged(QString)'),
                self.filter_changed)
        #connect(self.ui.spin_filter1, signal("editingFinished()"),
        #        self.update_raw)
        #connect(self.ui.spin_filter2, signal("editingFinished()"),
        #        self.update_raw)
        connect(self.ui.check_polar, signal('toggled(bool)'),
                self.polar_toggled)
        #connect(self.ui.check_select1, signal('toggled(bool)'),
        #        self.select1_toggled)
        #connect(self.ui.check_select2, signal('toggled(bool)'),
        #        self.select2_toggled)
        connect(self.ui.check_select3, signal('toggled(bool)'),
                self.select3_toggled)
        #connect(self.ui.check_display2, signal('toggled(bool)'),
        #        self.display2_toggled)

        #connect(self.ui.edit_slice, signal('textChanged(QString)'),
        #        self.slice_changed)
        connect(self.ui.push_update1, signal('clicked()'), self.update_raw)
        connect(self.ui.push_update2, signal('clicked()'), self.update_rf)
        connect(self.ui.push_delete, signal('clicked()'), self.delete_selected)
        connect(self.ui.push_select, signal('clicked()'), self.select_stream)
        connect(self.ui.push_sort, signal('clicked()'), self.sort_stream)
        connect(self.ui.push_calc, signal('clicked()'), self.calculate_more)
        self.fig.canvas.mpl_connect('button_press_event', self.button_pressed)

        self.logging_window()
        log.info('Hallo')
        self.open_data()

    def add_console(self):
        addConsole(self.win, self.ui, message='self is under namesspace ns',
                    namespace={'ns':self}, multithreaded=False)
    def start_ipy(self):
        embed()
    def open_data(self, file=None): #@ReservedAssignment
        if file == None:
            file = easygui.fileopenbox(default=def_open, filetypes=['*.QHD', '*.*']) #@ReservedAssignment
        if file != None:
            try:
                self.st_all = read(file)
            except:
                exceptionbox("Error while reading file!")
            if self.st_all[0].stats.event.get('id') == None:
                for tr in self.st_all:
                    tr.stats.event['id'] = str(tr.stats.event.eventno)
            self.handmarked = []
            self.handunmarked = []
            self.st_all.downsample2(self.ui.spin_downsample.value())
            self.st_all.trim2(-50, 300)
            self.st_all.sort(self.sort_after.split(','))
            self.st_all.check()
            self.st_all.setHI('mark', False)
            #self.slice_changed(self.ui.edit_slice.text())
            #self.update_raw()
            self.st_all_old = self.st_all
            self.update_rf()

    def save_data(self, file=None): #@ReservedAssignment
        if self.st == None:
            msgbox('No stream to save.')
            return
        if file == None:
            file = easygui.filesavebox(msg='Please omit file ending.', default=def_save_data) #@ReservedAssignment
        if file != None:
            try:
                if self.ui.check_save.isChecked():
                    self.st.write(file, 'Q')
                else:
                    self.st.select(expr='st.mark==False').write(file, 'Q')
            except:
                exceptionbox("Error while writing file!")
    def save_rf(self, file=None): #@ReservedAssignment
        if self.st3 == None:
            msgbox('No stream to save.')
            return
        if file == None:
            file = easygui.filesavebox(default=def_save_rf) #@ReservedAssignment
        if file != None:
            try:
                if self.ui.check_save.isChecked():
                    self.st3.write(file, 'Q')
                else:
                    self.st3.select(expr='st.mark==False').write(file, 'Q')
            except:
                exceptionbox("Error while saving file!")
    def fill_data_triggered(self):
        self.fill_data = easygui.choicebox(msg='Pick the filling for data.', title=' ', choices=('None', 'black-grey', 'red-blue'))
        if str(self.fill_data) == 'None':
            self.fill_data = 'white-white'
    def fill_rf_triggered(self):
        self.fill_rf = easygui.choicebox(msg='Pick the filling for rf.', title=' ', choices=('None', 'black-grey', 'red-blue'))
        if str(self.fill_rf) == 'None':
            self.fill_rf = 'white-white'
    def logging_file(self):
        file = easygui.filesavebox(default='/home/richter/Data/*.log') #@ReservedAssignment
        if file != None:
            if self.filehandler != None:
                log.removeHandler(self.filehandler)
            self.filehandler = logging.FileHandler(file, 'a')
            self.filehandler.setLevel(logging.DEBUG)
            #filehandler.setFormatter(logging.Formatter('%(asctime)s
            #        - %(levelname)s - %(name)s.%(funcName)s - %(message)s'))
            log.addHandler(self.filehandler)
    def logging_window(self):
        if self.streamhandler != None:
            log.removeHandler(self.streamhandler)
        self.streamhandler = addLogging(self.win, self.ui, log)

    def polar_toggled(self, bool): #@ReservedAssignment
        self.ui.spin_rot1.setEnabled(bool)
        self.ui.spin_rot2.setEnabled(bool)
    def select1_toggled(self, bool): #@ReservedAssignment
        self.ui.spin_select1.setEnabled(bool)
    def select2_toggled(self, bool): #@ReservedAssignment
        self.ui.spin_select2.setEnabled(bool)
        self.ui.spin_select3.setEnabled(bool)
    def select3_toggled(self, bool): #@ReservedAssignment
        if not bool:
            self.handmarked = []
            self.handunmarked = []
    def filter_changed(self, choice):
        if choice in ['LP', 'BP']:
            self.ui.spin_filter1.setEnabled(True)
        else:
            self.ui.spin_filter1.setEnabled(False)
        if choice in ['HP', 'BP']:
            self.ui.spin_filter2.setEnabled(True)
        else:
            self.ui.spin_filter2.setEnabled(False)
    def select_stream(self):
        selection = easygui.enterbox(msg='Select - use st for trace.stats:', default='225 < st.azi < 265.7')
        try:
            self.st_all = self.st_all_old.select(expr=selection)
        except:
            exceptionbox("Error while selecting!")
    def sort_stream(self):
        self.sort_after = easygui.enterbox(msg='Sort after:', default=def_sort)
        if self.sort_after == None:
            self.sort_after = def_sort
        if self.sort_after != def_sort:
            self.st_all.sort(self.sort_after.split(','))

    def calculate_more(self):
        file_load_list = []
        file_save_list = []
        if self.ui.edit_slice.text() != ':' and easygui.ynbox(msg="Maybe you want to change the slice setting to ':'?"):
            self.ui.edit_slice.setText(':')
        while True:
            file_load = easygui.fileopenbox(default='/home/richter/Data/Parkfield/receiver/M5.5_events/*.QHD', filetypes=['*.QHD', '*.*'])
            if file_load == None:
                break
            file_save = easygui.filesavebox(msg='Please omit file ending.', default='/home/richter/Data/Parkfield/receiver/*')
            if file_save == None:
                break
            file_load_list.append(file_load)
            file_save_list.append(file_save)
        for i in range(len(file_load_list)):
            self.open_data(file_load_list[i])
            self.save_rf(file_save_list[i])

    def button_pressed(self, event):
        if not self.ui.check_select3.isChecked() or not event.inaxes:
            return
        i = j = int(round(event.ydata))
        markings = self.st3.getHI('mark')[::3]
        if not self.ui.check_display2.isChecked() :
            if event.inaxes == self.ax2:
                for k, mark in enumerate(markings):
                    if k - sum(markings[:k]) == j and not mark:
                        i = k
                        break
            elif event.inaxes == self.ax1:
                j = i - sum(markings[:i])
        if i >= len(self.st) // 3: return
        if event.button == 1:
            if i in self.handunmarked:
                self.handunmarked.remove(i)
                self.ax1.plot([-20, -10], [i, i], 'r', lw=2)
                if self.ui.check_display2.isChecked():
                    self.ax2.plot([-10, -5], [j, j], 'g', lw=3)
            elif i not in self.handmarked and not self.st3[3 * i].stats.mark:
                self.handmarked.append(i)
                self.ax1.plot([-20, -10], [i, i], 'r', lw=2)
                self.ax2.plot([-10, -5], [j, j], 'g', lw=3)
        elif event.button == 3:
            if i in self.handmarked:
                self.handmarked.remove(i)
                self.ax1.plot([-20, -10], [i, i], 'w', lw=2)
                self.ax2.plot([-10, -5], [j, j], 'w', lw=3)
            elif i not in self.handunmarked and self.st3[3 * i].stats.mark:
                self.handunmarked.append(i)
                self.ax1.plot([-20, -10], [i, i], 'w', lw=2)
                if self.ui.check_display2.isChecked():
                    self.ax2.plot([-10, -5], [j, j], 'w', lw=3)
        # print 'marked %s unmarked %s' %(self.handmarked, self.handunmarked)
        self.canv.draw()


    def update_raw(self):
        if self.st_all == None:
            return
        slice = self.ui.edit_slice.text() #@ReservedAssignment
        if ':' not in slice:
            self.ui.edit_slice.setText(':')
            slice = ':' #@ReservedAssignment
        entry1, entry2 = slice.split(':')
        if entry1 != '' and int(entry1) % 3 != 0:
            entry1 = str(int(entry1) // 3 * 3)
        if entry2 != '' and int(entry2) % 3 != 0:
            entry2 = str(int(entry2) // 3 * 3)
        slice = '%s:%s' % (entry1, entry2) #@ReservedAssignment
        self.ui.edit_slice.setText(slice)
        self.st = eval('self.st_all[%s]' % slice)
        self.st2 = self.st.copy()
        for i in self.handmarked:
            self.st2[3 * i:3 * i + 3].setHI('mark', True)
        # filter
        text = self.ui.combo_filter.currentText()
        value1 = self.ui.spin_filter1.value()
        value2 = self.ui.spin_filter2.value()
        if text == 'LP':
            self.st2.filter2(0, 1. / value1)
        elif text == 'HP':
            self.st2.filter2(1. / value2, 0)
        elif text == 'BP':
            self.st2.filter2(1. / value2, 1. / value1)
        if self.ui.check_display1.isChecked():
            self.draw_raw()
        # rotate
        self.st2.rotateZNE2LQT(self.ui.spin_rot1.value(),
                               self.ui.spin_rot2.value(),
                               usetheo=not self.ui.check_polar.isChecked())
        if not self.ui.check_display1.isChecked():
            self.draw_raw()
        if self.ui.check_select1.isChecked():
            self.st2.afarm(signoise=self.ui.spin_select1.value(), remove=False)
        for i in self.handunmarked:
            self.st2[3 * i:3 * i + 3].setHI('mark', False)
        self.draw_mark('first')
        self.canv.draw()

    def update_rf(self):
        # first seletcion
        self.update_raw()
        self.st3 = self.st2#.copy()
        # deconvolution
        self.st3.receiverf(water=self.ui.spin_dec1.value(),
                           gauss=self.ui.spin_dec2.value(),
                           tshift=self.ui.spin_dec3.value(), pad=0,
                           window='tukey', start=self.ui.spin_win1.value(),
                           end=self.ui.spin_win2.value(), where='ponset',
                           lenslope=self.ui.spin_win3.value())
        # second selection
        if self.ui.check_select2.isChecked():
            self.st3.afarm('rf', signoise=self.ui.spin_select2.value(),
                           signoiseQ=self.ui.spin_select3.value(),
                           maxL=1 / self.ui.spin_select4.value(),
                           sigQ=self.ui.check_select4.isChecked(),
                           broad=self.ui.check_select5.isChecked(),
                           remove=False)
        for i in self.handunmarked:
            self.st3[3 * i:3 * i + 3].setHI('mark', False)
        self.draw_rf()
        self.draw_mark()
        self.canv.draw()
        #self.handmarked = []
        #self.handunmarked = []
    def display2_toggled(self):
        pass
        #self.ax2.shared

    def delete_selected(self):
        markings = self.st3.getHI('mark')
        for i, tr in enumerate(getattr(self, 'st_all[%s]' % self.ui.edit_slice.text())):
            if markings[i]:
                log.info('remove trace %s' % tr)
                self.st_all.remove(tr)

    def draw_raw(self):
        ax = self.ax1
        if ax.lines:
            xlims = list(ax.get_xlim())
            #ylims = list(ax.get_ylim())
        ax.clear()
        if self.ui.check_display1.isChecked():
            comp = 'Z'
        else:
            comp = 'L'
        c1, c2 = self.fill_data.split('-')
        pl = self.st2.plotRF(ax=ax, component=comp, topcolor=c1, botcolor=c2, scale=self.ui.spin_scale1.value(), plotsum=False, show=False, figtitle='')
        self.ax1 = pl.ax
        t = self.ax1.lines[0].get_xdata()
        data = np.mean([i.data for i in self.st2 if i.stats.channel[-1] == comp and (self.ui.check_display3.isChecked() or i.stats.mark == False)], axis=0)
        self.ax1_sum.clear()
        self.ax1_sum.plot(t, data, 'k')
        if c1 != 'white':
            self.ax1_sum.fill_between(t, data, 0, linewidth=0, where=data >= 0, facecolor=c1)
        if c2 != 'white':
            self.ax1_sum.fill_between(t, data, 0, linewidth=0, where=data < 0, facecolor=c2)
        try:
            ax.set_xlim(xlims)
        except UnboundLocalError:
            ax.set_xlim([t[0], t[-1]])
        for tick in self.ax1_sum.get_xticklabels() + self.ax1_sum.get_yticklabels()[:-1]:
            tick.set_visible(False)

    def draw_mark(self, which='both'):
        imaging.plotRFmarks(self.st2, self.ax1)
        if self.ui.check_display2.isChecked() and which in ['both', 'rf']:
            imaging.plotRFmarks(self.st2, self.ax2, -10, -5, 'g', 3)

    def draw_rf(self):
        ax = self.ax2
        if ax.lines:
            xlims = list(ax.get_xlim())
        sharey = None
        if self.ui.check_display2.isChecked():
            sharey = self.ax1
        bounds = ax.get_position().bounds
        self.fig.delaxes(ax)
        ax = self.ax2 = self.fig.add_axes(bounds, sharey=sharey, sharex=self.ax2_sum) # [0.53, 0.04, 0.46, 0.91],

        c1, c2 = self.fill_rf.split('-')
        if not self.ui.check_display2.isChecked():
            st_plot = self.st3.select(expr='not st.mark')
        else:
            st_plot = self.st3
        if len(st_plot) > 0:
            pl = st_plot.plotRF(ax=ax, component='Q', topcolor=c1, botcolor=c2, scale=self.ui.spin_scale3.value(), plotsum=False, show=False, figtitle='')
            self.ax2 = pl.ax
            t = self.ax2.lines[0].get_xdata()
            data = np.mean([i.data for i in self.st2 if i.stats.channel[-1] == 'Q' and (self.ui.check_display3.isChecked() or i.stats.mark == False)], axis=0)
    #        st_plot = self.st3.select(component='Q')
    #        if not self.ui.check_display4.isChecked():
    #            st_plot = st_plot.select(expr='st.mark==False')
    #        data = np.mean(util.imaging.getDataWindow(st_plot, None, None), axis=0)
            self.ax2_sum.clear()
            self.ax2_sum.plot(t, data, 'k')
            if c1 != 'white':
                self.canv.draw()
                self.ax2_sum.fill_between(t, data, 0, linewidth=0, where=data >= 0, facecolor=c1)
            if c2 != 'white':
                self.ax2_sum.fill_between(t, data, 0, linewidth=0, where=data < 0, facecolor=c2)
            try:
                ax.set_xlim(xlims)
            except UnboundLocalError:
                ax.set_xlim([t[0], t[-1]])
            for tick in self.ax2_sum.get_xticklabels() + self.ax2_sum.get_yticklabels()[:-1]:
                tick.set_visible(False)


if __name__ == "__main__":
    import sys
    from PyQt4 import QtGui
    from rftool_gui import Ui_MainWindow
    app = QtGui.QApplication(sys.argv)
    gui = RfTool(QtGui.QMainWindow(), Ui_MainWindow())
    gui.win.show()
    sys.exit(app.exec_())
