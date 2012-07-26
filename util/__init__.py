# -*- coding: utf-8 -*-
# TR

from main import (isint, isfloat, isnumber, gps2DistDegree, smooth, filterResp,
                  calculate, yeargen, timegen)
from helper import (NullHandler, checkDir, exha, add_doc,
                    vectorize_args, parameters, pop, runBash, setRootLogger,
                    fillArray, TerminalFile, progress)
from seispy import (dist2gps, iasp91, Iasp91, TwoLayer, mocorr, zmigr, pspier, phase_dict,
                    time2depth, depth2time, multiples, interpolate, ttt,
                    correct_ttt, feregion)
from imaging import getDataWindow, getTimeIntervall, xcorr_cmap, DLogNorm
#from imaging import getWindow, getTimeIntervall, getDataWindow


from _pspier import pspier as pspier2 # adapted Fortran code from XY @UnresolvedImport

from scipy.signal import get_window
from obspy.signal import cosTaper
