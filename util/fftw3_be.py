# by TR

import os
import fftw3
from sito.data import FFTW3_WISDOM_FILE

_plan = None
_plan2 = None

USE_OLD_PLAN = False


class __WisdomHandler(object):
    wisdomfile = FFTW3_WISDOM_FILE

    def __init__(self):
        if os.path.isfile(self.wisdomfile):
            fftw3.import_wisdom_from_file(self.wisdomfile)

    def remove(self):
        if os.path.isfile(self.wisdomfile):
            os.remove(self.wisdomfile)

    def forget(self):
        fftw3.forget_wisdom()

    def setUseOldPlan(self, val):
        """ val 1 == True, 2 = False """
        global USE_OLD_PLAN
        USE_OLD_PLAN = val

    def __del__(self):
        try:
            fftw3.export_wisdom_to_file(self.wisdomfile)
        except TypeError:
            pass

wisdom = __WisdomHandler()


def fft(data, nfft=None, in_place=False, use_old_plan=True, **kwargs_in):
    global _plan
    if USE_OLD_PLAN:
        use_old_plan = (USE_OLD_PLAN == 1)
#    print('USE OLD PLAN %s' % str(use_old_plan))
    kwargs = dict(flags=['measure'])
    kwargs.update(kwargs_in)
    if nfft == None:
        nfft = len(data)
    if _plan is None or (use_old_plan and nfft != len(_plan.inarray)):
        use_old_plan = False
    if not use_old_plan:
        input_ = fftw3.create_aligned_array(nfft)
        if in_place:
            output = None
        else:
            output = fftw3.create_aligned_array(nfft)
        _plan = fftw3.Plan(input_, output, **kwargs)
    _plan.inarray[:len(data)] = data
    _plan.inarray[len(data):] = 0
    _plan.outarray[:] = 0
    _plan()
    if _plan.outarray is None:
        ret = _plan.inarray
    else:
        ret = _plan.outarray
    return ret


def ifft(data, nfft=None, in_place=False, use_old_plan=True, **kwargs_in):
    global _plan2
    if USE_OLD_PLAN:
        use_old_plan = (USE_OLD_PLAN == 1)
    kwargs = dict(direction='backward', flags=['measure'])
    kwargs.update(kwargs_in)
    if nfft == None:
        nfft = len(data)
    if _plan2 is None or(use_old_plan and nfft != len(_plan2.inarray)):
        use_old_plan = False
    if not use_old_plan:
        input_ = fftw3.create_aligned_array(nfft)
        if in_place:
            output = None
        else:
            output = fftw3.create_aligned_array(nfft)
        _plan2 = fftw3.Plan(input_, output, **kwargs)
    _plan2.inarray[:len(data)] = data
    _plan2.inarray[len(data):] = 0
    _plan2()
    if _plan2.outarray is None:
        ret = _plan2.inarray
    else:
        ret = _plan2.outarray
    return ret / nfft
