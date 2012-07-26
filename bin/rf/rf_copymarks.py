#!/usr/bin/env python
# by TR
"""
read header info ('mark') from files in one folder
and write it to other files
"""

import glob
from sito import read
import os.path
glob1 = '/home/richter/Results/IPOC/receiver/2012_mag5.5/*_????_mout.QHD'
glob2 = '/home/richter/Results/IPOC/receiver/2012_mag5.5/*_????_mout.QHD'


merge = False

filelist1 = glob.glob(glob1)
filelist2 = glob.glob(glob2)
N = len(filelist1)
assert len(filelist2) == N
for i in range(N):
    file1 = filelist1[i]
    file2 = filelist2[i]
    st1 = read(file1)
    st2 = read(file2)
    mark1 = st1.getHI('mark')
    reason1 = st1.getHI('reason')
    if merge:
        mark2 = st2.getHI('mark')
        reason2 = st2.getHI('reason')
        mark = [mark1[i] or mark2[i] for i in range(N)]
        reason = reason1.copy()
        for i in range(N):
            if mark2[i]:
                reason[i] = reason2[i]
    else:
        st2.setHI('mark', mark1)
        st2.setHI('reason', reason1)
    files = os.path.splitext(file2)[0]
    st2.write(file2 + 'copymark', 'Q')



