#!/usr/bin/env python
# by TR
"""
Convert prepared data for cross correlation to fft
"""
import logging
from obspy.core import UTCDateTime as UTC
from sito.data import IPOC
from sito.util.main import daygen
logging.basicConfig()
log = logging.getLogger('convert_fft')

# import warnings
# warnings.simplefilter("error", np.ComplexWarning)

def main():
    stations = 'PB03 PB04'

    component = 'Z'
    t1 = UTC('2006-01-01')
    t2 = UTC()
#    t1 = UTC('2007-01-01')
#    t2 = UTC('2007-01-03')

    method1 = 'filter0.01-1_water_env2_whitening_1bit'
    method2 = 'filter0.01-1_water_env2_whitening_1bit_fft'

    data1 = IPOC(xcorr_append='/Tocopilla/' + method1, use_local_LVC=True)
    data2 = IPOC(xcorr_append='/Tocopilla/' + method2, use_local_LVC=True)

    for station in stations.split():
        for day in daygen(t1, t2):
            try:
                stream = data1.getStream(day, station, component)
            except:
                log.warning('Could not read stream for day %s station %s' % (day, station))
            else:
                if len(stream) != 1:
                    log.warning('Stream for day %s station %s has wrong length %d' % (day, station, len(stream)))
                elif stream[0].stats.npts / stream[0].stats.sampling_rate < 24 * 3600 * 0.5:
                    log.warning('Stream for day %s station %s has only a coverage of %f  -> discard' % (day, station, 1. * stream[0].stats.npts / stream[0].stats.sampling_rate / 24 / 3600))
                else:
                    stream.fft()
                    stream.write(data2.getDay(station, day), 'Q')

#    day = UTC('2007-01-01')
#    station = 'PB03'
#    stream1 = data1.getStream(day, station, component)
#    stream2 = data2.getStream(day, station, component)
#    stream2.ifft()
#    np.testing.assert_array_almost_equal(stream2[0].data, stream1[0].data)

if __name__ == '__main__':
    main()
