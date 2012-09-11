# by TR

from select import select
import gc
import inspect
import logging
import numpy as np
import os.path
import shlex
import subprocess
import sys
import tty
try:
    from IPython.display import clear_output #@UnusedImport
    have_ipython = True
except ImportError:
    have_ipython = False
import progressbar #@UnresolvedImport
from functools import wraps

widgets = [progressbar.Bar('>'), ' ', progressbar.ETA(), ' ', progressbar.ReverseBar('<')]
progress = progressbar.ProgressBar(widgets=widgets)


#try:
#    # debugging and ipython utilities for easy import:
#    # import sito; sito.ipshell(); sito.debug()
#    try:
#        import sys
#        from IPython.frontend.terminal.embed import InteractiveShellEmbed
#        import IPython.core.debugger
#        from IPython.core import ultratb
#
#        #from IPython import embed as ipshell
#        ipshell = InteractiveShellEmbed(banner2='Entering IPython.  Press Ctrl-D to exit.',
#                        exit_msg='Leaving Interpreter, back to program.')
#        debug = IPython.core.debugger.Tracer(colors='LightBG')
#        sys.excepthook = ultratb.FormattedTB(mode='Context', #'Verbose',
#                                             color_scheme='LightBG', call_pdb=1)
##        def call_pdb():
##            sys.excepthook = ultratb.FormattedTB(mode='Context', #'Verbose',
##                                             color_scheme='LightBG', call_pdb=1)


#    except:
#        print 'IPython<0.11'
#        import IPython.Shell #@UnresolvedImport
#        ipshell = IPython.Shell.IPShellEmbed(['-pdb'], #'-profile', 'ipshell'], # pylab option does not work proberly
#                                             banner='Entering IPython.  Press Ctrl-D to exit.',
#                                             exit_msg='Leaving Interpreter, back to program.')
#        import IPython.Debugger.Tracer as debug #@UnresolvedImport
#except:
#    ipshell = None
#    debug = None
#    import warnings
#    warnings.warn('Unable to import IPython.Debugger or Ipython.Shell')

log = logging.getLogger(__name__)





def dumpObjects():
    """
    usage:

    class TestClass:
        pass

    testObject1 = TestClass()
    testObject2 = TestClass()

    from wx import Colour
    color = Colour()

    dumpObjects()


    """
    exclude = ["function", "type", "list", "dict", "tuple", #@UnusedVariable
               "wrapper_descriptor",
               "module", "method_descriptor", "member_descriptor",
               "instancemethod",
               "builtin_function_or_method", "frame", "classmethod",
               "classmethod_descriptor", "_Environ", "MemoryError", "_Printer",
               "_Helper", "getset_descriptor", "weakref", "property", "cell",
               "staticmethod", 'set', 'CodecInfo']

    gc.collect()
    oo = gc.get_objects()
    for o in oo:
        if getattr(o, "__class__", None):
            name = o.__class__.__name__
            #print (name)
            #if name not in exclude and 'Error' not in name:
            if name in ('Stream', 'AttribDict', 'Trace', 'Stats') or 'Data' in name:
                print "Object : %s..." % o.__repr__()
                if name == 'Stream':
                    print "__str__: %s..." % o
                print "Class  : %s..." % name
                try:
                    filename = inspect.getabsfile(o.__class__)
                    print "defined: %s\n" % filename
                except TypeError:
                    print "defined: built-in\n"

class NullHandler(logging.Handler):
    """
    NullHandler for use with logging if no logging is wanted.

    Usage:
    import logging, sito.util
    logging.getLogger('').addHandler(sito.util.NullHandler())
    """
    def emit(self, record):
        pass

def checkDir(filename):
    """
    Create directory of filname if it does not exist.
    """
    head = os.path.dirname(filename)
    if head != '' and not os.path.isdir(head):
        os.makedirs(head)

def exha(*posargs):
    """
    Exception handling decorator.

    If an Exception is raised in the decorated function, it get caught
    and logged. (And function is terminated.)
    """
    def wrapper(f):
        if len(posargs) != 0:
            t = tuple(item for item in posargs[0]
                if issubclass(item, Exception)) or (Exception,)
        else:
            t = (Exception,)
        def newfunc(*pargs, **kwargs):
            try:
                f(*pargs, **kwargs)
            except t:
                # Only inform the user that the exception occured
                log.exception('An exception occured:')
        return newfunc
    return wrapper

def add_doc(*docs):
    """
    Decorator to add doc strings of other funtions.
    
    The Docstring of other functions (passed as arguments) or the argument
    itself (if string) are appended to the end of the docstring of the decorated
    function.
    """
    def wrap(func):
        docstrs = list(docs)
        for i, docstr in enumerate(docstrs):
            if not isinstance(docstr, basestring):
                docstrs[i] = docstr.__doc__
                if docstr.__doc__ is None:
                    docstrs[i] = ''
        docstr = ''.join(docstrs)
        if func.__doc__ == None:
            func.__doc__ = docstr
        else:
            func.__doc__ = func.__doc__ + docstr
        return func
    return wrap

def vectorize_args(nums):
    """
    Decorator for vectorization of arguments of a function.

    The positions of the arguments are given in the tuple nums.
    See numpy.vectorize.
    """
    def wrap(func):
        @wraps(func)
        def wrapped(*args, ** kwargs):
            args = list(args)
            for i, arg in enumerate(args):
                if i in nums and type(arg) == list:
                    args[i] = np.array(arg)
            for i, arg in enumerate(args):
                if i in nums and type(arg) == np.ndarray:
                    shape = np.shape(arg)
                    ind = np.transpose(np.ones(shape).nonzero())
                    break
                if i == len(args) - 1:
                    # no need for vectorization as all relevant
                    # arguments are scalars
                    return func(*args, ** kwargs)
            res = np.array([func(
                           * [arg[tuple(j)] if type(arg) == np.ndarray and i in nums else arg for i, arg in enumerate(args)], ** kwargs)
                           for j in ind])
            if np.shape(res) <> shape:
                # func returns more than 1 result, this means the array has to
                # be ordered differently
                res = res.transpose()
            if len(shape) > 1:
                # more than 1D arrays, the shape of the list has to be rearanged
                res = res.reshape((res.shape[0],) + shape)
            return res
        return wrapped
    return wrap


def parameters(only=None, exclude=None, ignore='self', add=None, format_='%s=%s'):
    """
    Return a dictionary or string of the calling functions parameter names and values.
    
    The optional arguments can be used to filter the result:
    :param only: use this to only return parameters from this list of names.
    :param exclude: use this to return every parameter *except* those included
        in this list of names.
    :param ignore: use this inside methods to ignore the calling object's name.
        For convenience, it ignores 'self' by default.
    :return: dictionary of parameters

    http://code.activestate.com/recipes/201195/
    (r2)    
    """
    args, varargs, varkw, defaults = inspect.getargvalues(#@UnusedVariable
                                                          inspect.stack()[1][0])
    if only is None:
        only = args[:]
        if varkw:
            only.extend(defaults[varkw].keys())
            defaults.update(defaults[varkw])
    #ipshell()
    if add:
        only.extend(add.keys())
        defaults.update(add)
    if exclude is None:
        exclude = []
    exclude.append(ignore)
    ret = dict([(attrname, defaults[attrname])
               for attrname in only if attrname not in exclude])
    if format_ is None:
        return ret
    else:
        return ' '.join([format_ % (i, j) for i, j in ret.iteritems()])
# end of http://code.activestate.com/recipes/201195/


def pop(command):
    """
    Execute command, log stout and return the return code of the process.
    """
    sub = subprocess.Popen(shlex.split(command), stdin=None,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           shell=False, close_fds=True)
    stdout = sub.communicate()[0]
    if stdout != '': log.info('stdout: ' + stdout[:-1])
    retcode = sub.poll()
    return retcode

def runBash(command):
    """
    Execute command and return stdout.
    """
    sub = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    return sub.communicate()[0]

def setRootLogger(logfile=None, logdebugfile=None, console=True, fmode='w'):
    """
    Prepare the root logger for use with different handlers.

    :param logfile: logging file for INFO level
    :param logdebugfile: logging file for DEBUG level
    :param console: turn logging to console on/off
    :param fmode: filemode for logging files (w/a)
    """
    rootlog = logging.getLogger()
    rootlog.setLevel(logging.DEBUG)
    for handler in rootlog.handlers[::-1]:
        rootlog.removeHandler(handler)
    if console:
        consolehandler = logging.StreamHandler()
        consolehandler.setLevel(logging.DEBUG)
        consolehandler.setFormatter(logging.Formatter('%(name)s - %(levelname)s - %(message)s'))
        rootlog.addHandler(consolehandler)
    if logfile:
        checkDir(logfile)
        filehandler = logging.FileHandler(logfile, fmode)
        filehandler.setLevel(logging.INFO)
        filehandler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        rootlog.addHandler(filehandler)
    if logdebugfile:
        checkDir(logdebugfile)
        filehandler2 = logging.FileHandler(logdebugfile, fmode)
        filehandler2.setLevel(logging.DEBUG)
        filehandler2.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(name)s.%(funcName)s - %(message)s'))
        rootlog.addHandler(filehandler2)

class NotTTYException(Exception): pass

class TerminalFile:
    """
    Check if key was pressed
    
    Script has to be run in terminal for this.
    
    :usage:
    
    import sys
    from sito.util import TerminalFile
     
    s = TerminalFile(sys.stdin)
    print "Press q to quit..."
    i = 0
    while s.getch() != "q":
        sys.stdout.write("%08d\r" % i)
        i += 1
    print "-- END --"

    http://code.activestate.com/recipes/203830-checking-for-a-keypress-without
    -stop-the-execution/
    """
    def __init__(self, infile):
        if not infile.isatty():
            raise NotTTYException()
        self.file = infile

        #prepare for getch
        self.save_attr = tty.tcgetattr(self.file)
        newattr = self.save_attr[:]
        newattr[3] &= ~tty.ECHO & ~tty.ICANON
        tty.tcsetattr(self.file, tty.TCSANOW, newattr)

    def __del__(self):
        #restoring stdin
        import tty  #this import is required here @Reimport
        tty.tcsetattr(self.file, tty.TCSADRAIN, self.save_attr)

    def getch(self):
        if select([self.file], [], [], 0)[0]:
            c = self.file.read(1)
        else:
            c = ''
        return c

def fillArray(data, mask=None, fill_value=None):
    """
    Fill masked numpy array with value without demasking.
    
    Additonally set fill_value to value.
    If data is not a MaskedArray returns silently data. 
    """
    if mask is not None and mask is not False:
        data = np.ma.MaskedArray(data, mask=mask, copy=False)
    if np.ma.is_masked(data) and fill_value is not None:
        data._data[data.mask] = fill_value
        np.ma.set_fill_value(data, fill_value)
    elif not np.ma.is_masked(data):
        data = np.ma.filled(data)
    return data

#Ipython example
class ProgressBar:
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 40
        self.__update_amount(0)
        if have_ipython:
            self.animate = self.animate_ipython
        else:
            self.animate = self.animate_noipython

    def animate_ipython(self, iter):
        print '\r', self,
        sys.stdout.flush()
        self.update_iteration(iter + 1)

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)

if __name__ == '__main__':
    p = ProgressBar(1000)
    for i in range(1001):
        p.animate(i)
