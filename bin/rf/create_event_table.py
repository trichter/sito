from sito.events import Events
import cPickle

def create_all_events_pickle():
    ev = Events.read('/home/richter/Data/events/2012_03_events_27-93_mag5.5_IPOC.txt')
    with open('./all_events.pickle', 'wb') as f:
        cPickle.dump(ev, f, -1)

def create_all_events_pickle_with_rf():
    with open('./all_events.pickle', 'rb') as f:
        ev2 = cPickle.load(f)
    stations = ('PB01 PB02 PB03 PB04 PB05 PB06 PB07 PB08 PB09 PB10 PB11 PB12 '
                'PB13 PB14 PB15 PB16 HMBCX MNMCX PATCX PSGCX LVC')
    for event in ev2:
        event.rf_usage = '.' * len(stations.split())
    for i_st, station in enumerate(stations.split()):
        print 'Station %s' % station
        ev_st = Events.read('/home/richter/Results/IPOC/receiver/2012_mag5.5/event_maps/events_%s.txt' % station)
        i = 0
        j = 0
        while i < len(ev2) and j < len(ev_st):
            e1 = ev2[i]
            e2 = ev_st[j]
            if e1.id == e2.id:
                e1.rf_usage = e1.rf_usage[:i_st] + 'x' + e1.rf_usage[i_st + 1:]
                j += 1
            i += 1
    with open('./all_events_with_rf_usage.pickle', 'wb') as f:
        cPickle.dump(ev2, f, -1)

def write_all_events():
    with open('./all_events.pickle', 'rb') as f:
        ev = cPickle.load(f)
    ev.write_LATEX('all_events.tex')

def write_all_events_with_rf():
    with open('./all_events_with_rf_usage.pickle', 'rb') as f:
        ev = cPickle.load(f)
    ev.write_LATEX2('all_events_with_rf.tex')

#create_all_events_pickle()
#create_all_events_pickle_with_rf()
#write_all_events()
write_all_events_with_rf()

from IPython import embed
embed()



