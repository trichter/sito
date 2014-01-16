#!/usr/bin/env python
# by TR
from obspy.arclink import Client
from obspy.core import UTCDateTime as UTC
import yaml


def save_paz():
    client = Client()
    paz2 = client.getPAZ(network, station, location, channel, t)
    with open(PAZ_FILE % (station, channel), 'w') as f:
        f.write(yaml.dump(dict(paz2)))

def load_paz():
    with open(PAZ_FILE % (station, channel)) as f:
        return yaml.load(f.read())

def save_all_paz():
    client = Client()
    paz = {}
    for st in 'PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PATCX HMBCX PSGCX MNMCX'.split():
        paz[st] = dict(client.getPAZ('CX', st, location, channel, t))
    print paz
    with open(PAZ_FILE % ('ALL', channel), 'w') as f:
        f.write(yaml.dump(paz))



PAZ_FILE = '/home/richter/Data/paz/resp_%s_%s.yaml'
network = 'CX'
station = 'PB10'
channel = 'HHZ'
location = ''
t = UTC('2011-10-01')

save_paz()
#save_all_paz()
#paz = load_paz()
