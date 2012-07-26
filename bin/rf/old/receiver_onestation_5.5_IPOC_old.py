from sito import read, util, data as mod_data
import logging
import pylab as plt


ipoc = mod_data.IPOC()
#ipoc.events = '/home/richter/Data/events/events_30-90_mag5.8_IPOC.txt'
ipoc.rf_results = ipoc.results + '/receiver/results5.5/%s_%s'
ipoc.rf_results_dir = ipoc.results + '/receiver/results5.5/'
ipoc.rf_plot_dir = ipoc.rf_results_dir + 'plots/'
data = ipoc

station = 'PB03'

logfile = data.rf_results_dir + 'a_log_' + station + '_5.5_%s.txt'
util.checkDir(logfile)
util.setRootLogger(logfile=logfile % '', logdebugfile=logfile % '_debug')


stream = read('/data/home/richter/Data/IPOC/receiver/M5.5_events/*%s*QHD' % station)
print len(stream)
stream.pspier(60, data.stations)
stream.filter2(0.033, 2.)
stream.trim2(-50, 300)
stream.sort(('event.id', 'station', 'component'))
stream.check()
stream.rotateZNE2LQT(-5, 15, usetheo=True)
stream.afarm(signoise=2.0, remove=False)
stream.receivert()
stream.afarm('rf', signoise=False, signoiseQ=False, maxL=False, sigQ=False, broad=5, remove=False)
stream.writey(data.rf_results % (station, '%s') + '_all_nomout', 'Q')
stream.select(expr='st.mark==False').writey(data.rf_results % (station, '%s') + '_good_nomout', 'Q')
print stream.getReasons()
#
#regions = ((-58.3, -22.0, 800),)
#stream = read('../PB03*_all_nomout.QHD')
#stream.sort('starttime')
#
#for region in regions:
#    str_region = stream.select(around=region)
#    top = 'black'
#    bot = 'gray'
#    fig = str_region.plotRF(-10, 60, show=False, scale=2.2, sumscale=4, topcolor=top, botcolor=bot,
#                            plotinfo=('azi dist', 'starttime'), plotinfo_width=0.1, plotlabel=('time (s)', u'azi/dist', 'year'),
#                            usehardticks='time', figtitle = '%dkm around %.2f, %.2f - %s' % (region[2],region[0], region[1], util.region(region[0], region[1])), x_pos=0.02, horizontalalignment = 'left', showsum=False)
#    #fig.savefig(plotdir + 'rf_PKD_%dkm_around_%.1f_%.1f_Q.png' % (region[2],region[0], region[1]))
#    #plt.close(fig)
#plt.show()
