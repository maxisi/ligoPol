'''
Contains all paths used by different modules, making it easy to make changes.

File structure:
files/
    analysis/
        results/
    data/
        ephemerides/
    background/
    templates/
        antennapatterns/
        vectors/
    remote/
        source/
        ephemerides/
    plots/
        
'''
import os

leafList = (
            'analysis/results/',
            'analysis/sigma/',
            'data/ephemerides/',
            'background/',
            'templates/antennapatterns',
            'templates/vectors',
            'plots/'
            )
            
def makestructure():
    print 'Creating file structure:'
    # create files/
    for dir in leafList:
        try:
            os.makedirs('files/' + dir)
            print dir
        except OSError:
            print 'Error: files/%(dir)s already exists' % locals()


# Antenna patterns
psrlist = 'files/templates/psrlist.txt'         # list of pulsars
ap = 'files/templates/antennapatterns/'         # antenna patterns

# Data
rhB = 'files/background/'                       # background heterodynes
originalData = '../../matthew/analyses/S6_all/results' # 'files/remote/source'
importedData = 'files/data/finehet_'         # original data in DataFrame

# Detector
vectors = 'files/templates/vectors/'            # detector and wave vectors

# Source
textfromATNF = 'files/templates/psrcat.txt'     # download from ATNF
psrcat = 'files/templates/pulsar_catalogue'     # pulsar data in DF
psrextra = 'config/psrextra.txt'                 # polarization and inclination angles
psrlist = 'files/templates/psrlist.txt'         # names of PSR names to analyze

# Analysis
sigma = 'files/analysis/sigma/'        # day-long standard deviation
results = 'files/analysis/results/'             # background h_rec and significance
tmat = 'files/analysis/tmat'                    # basis transformation matrix

# Ephemerides
eph = 'files/remote/ephemerides/' #'usr/share/lalpulsar/'
eph_local = 'files/data/ephemerides/'           # imported ephemerides in HDF5 format

# Plots
plots = 'files/plots/'

