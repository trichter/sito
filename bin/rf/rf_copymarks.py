#!/usr/bin/env python
# by TR
"""
read header info ('mark') from files in one folder
and write it to other files
"""

import glob
from sito import read
import os.path
glob1 = '/home/richter/Results/IPOC/receiver/2012_mag5.5/*_mout.QHD'
file2 = '/home/richter/Results/IPOC/receiver/2012_mag5.5_RT/all/%s_nomout.QHD'
file3 = '/home/richter/Results/IPOC/receiver/2012_mag5.5_RT/%s_nomout'


def v2():
    filelist1 = glob.glob(glob1)
    for file1 in filelist1:
        station = file1.split('/')[-1].split('_')[0]
        print station
        st1 = read(file1)
        st2 = read(file2 % station)
        st1.sort(('event.id', 'channel'))
        st2.sort(('channel',), reverse=True)
        st2.sort(('event.id',))
        print len(st1), len(st2)
        i = 0
        while i < len(st1):
            if st1[i].stats.event.id == st2[i].stats.event.id:
                i += 1
            else:
                st2.pop(i)
        st2 = st2[:i]
        print len(st1), len(st2)
        st2.write(file3 % station, 'Q')

v2()

def v1():
    # not used anymore
    merge = False
    filelist1 = glob.glob(glob1)
    filelist2 = glob.glob(file2)
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



