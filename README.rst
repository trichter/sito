**sito** - Scripts for receiver function analysis and ambient noise cross-correlation - based on ObsPy

| **author**: Tom Eulenfeld
| **license**: MIT license
|

These scripts were written in the scope of my PhD at GFZ Postdam/FU Berlin.
I do not maintain them any longer, but of course you are free to copy&paste code of interest.

The following publicatons are based upon results of these scripts:

Richter, T., C. Sens-Schönfelder, R. Kind, and G. Asch (2014), Comprehensive observation and modeling of earthquake and temperature related seismic velocity changes in northern Chile with passive image interferometry, J. Geophys. Res. Solid Earth, 119, 4747–4765, doi:10.1002/2013JB010695. `pdf1 <http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:823917:3/component/escidoc:828895/823917.pdf>`_

Richter, T. (2014), Temporal Variations of Crustal Properties in Northern Chile Analyzed with Receiver Functions and Passive Image Interferometry,
dissertation, FU Berlin. `pdf2 <http://www.diss.fu-berlin.de/diss/servlets/MCRFileNodeServlet/FUDISS_derivate_000000014929/dissertation_richter.pdf>`_


A lot of alternative Python packages for ambient noise cross-correlation exist.
I refactored the receiver function functionality in a separate repository, see `rf <https://github.com/trichter/rf>`_.

Required:
    - python 2.7
    - numpy
    - scipy
    - matplotlib
    - basemap
    - ipython
    - obspy
    - seispy
    - progressbar

sito consist of the following major modules:

    :data: Class for data handling (You can set up your own data calss here)
    :events: Class for handling event data (better: use new Catalog class in ObsPy)
    :imaging: Customized plotting
    :map: Customized station map plotting
    :noisexcorr: correlation of noise
    :rf: Receiver function calculation
    :station: Class for handling station data
    :stream: Class derived from obspy.trace.Stream with custom methods
    :trace: Class derived from obspy.trace.Trace with custom methods
    :util: Utilities like deconvolution, polarisation, rotation, pspier (mainly imported from other projects)
    :xcorr: Cross correlation functions

minor modules:

    :noise_migration: migrate noise phases in correlations backwards (experimental)
    :debug: ipython debugging functionality
    :seismometer: PAZ of STS-2 and Wood-Anderson