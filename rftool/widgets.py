#!/usr/bin/env python

from PyQt4 import QtGui
from PyQt4.QtCore import Qt
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
import logging
import pylab as plt
from spyderlib.widgets.internalshell import InternalShell

class MPLWidget(FigureCanvasQTAgg):
    """
    Class to represent the FigureCanvas widget.
    """
    def __init__(self, parent=None):
        # Standard Matplotlib code to generate the plot
        self.fig = plt.Figure()
        # initialize the canvas where the Figure renders into
        FigureCanvasQTAgg.__init__(self, self.fig)
        self.setParent(parent)

def addMPLToolBar(window, canvas):
    qToolBar = QtGui.QToolBar()
    qToolBar.addWidget(NavigationToolbar2QTAgg(canvas, qToolBar))
    qToolBar.setMovable(False)
    qToolBar.setFloatable(False)
    window.addToolBar(Qt.BottomToolBarArea, qToolBar)

def addConsole(window, ui_window, **kwargs):
   # Create the console widget
    font = QtGui.QFont("Courier new")
    font.setPointSize(10)
    # Note: by default, the internal shell is multithreaded which is safer
    # but not compatible with graphical user interface creation.
    # For example, if you need to plot data with Matplotlib, you will need
    # to pass the option: multithreaded=False
    #Internal Shell(self, parent=None, namespace=None, commands=[], message='', max_line_count=300, font=None, debug=False, exitfunc=None, profile=False, multithreaded=True, light_background=True
    ui_window.console = cons = InternalShell(window, ** kwargs) #debug=True

    # Setup the console widget
    cons.set_font(font)
    cons.set_codecompletion_auto(True)
    cons.set_calltips(True)
    cons.setup_calltips(size=600, font=font)
    cons.setup_completion(size=(300, 180), font=font)
    console_dock = QtGui.QDockWidget("Console", window)
    console_dock.setWidget(cons)

    # Add the console widget to window as a dockwidget
    window.addDockWidget(Qt.BottomDockWidgetArea, console_dock)

class QtStreamHandler(logging.Handler):
    def __init__(self, parent):
        logging.Handler.__init__(self)
        self.textWidget = parent
        #self.formater = logging.Formatter("%(message)s")

    def emit(self, record):
        self.textWidget.appendPlainText(self.format(record))
        self.textWidget.repaint()

def addLogging(window, ui_window, logger):
    ui_window.log_window = QtGui.QPlainTextEdit()
    # Setup the  widget
    streamhandler = QtStreamHandler(ui_window.log_window)
    streamhandler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(name)s.%(funcName)s - %(message)s'))
    #logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(streamhandler)
    log_dock = QtGui.QDockWidget("Logging", window)
    log_dock.setWidget(ui_window.log_window)
    window.addDockWidget(Qt.BottomDockWidgetArea, log_dock)
    return streamhandler
