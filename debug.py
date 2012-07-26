import sys
from IPython.core import ultratb
sys.excepthook = ultratb.FormattedTB(mode='Verbose',
color_scheme='LightBG', call_pdb=1)

from IPython import embed #@UnusedImport
