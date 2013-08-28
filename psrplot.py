import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import math
import sidereal as sd
from scipy.optimize import fsolve

from sys import exit

# AUXILIARY FUNCTIONS

def lin_fit(lins):
    # assuming input is a series
    
    # get injections
    lins = lins.drop([0])
    
    # put vectors in proper shape
    x = np.reshape(lins.index, (len(lins), 1))
    y = np.reshape(lins, (len(lins), 1))
    
    # fit
    p, _, _, _ = np.linalg.lstsq(x, y)
    
    return np.poly1d([p[0][0], 0])
    
        
def noise_line(d):
    # input can be a dataframe
    # find elements with index 0 (no injection) and take their max over columns and rows
    try:
        max_n = d[0].max()
    except KeyError:
        max_n =d.T[0].max().max()
    # return corresponding constant polynomial
    return np.poly1d([max_n])


def min_det_h(d):
    # assumes input is a series
    # get max noise
    noise = noise_line(d)(1)
    # find max injection below that
    det_injs = d[d<noise]
    try:
        return det_injs.index[np.argmax(det_injs)]
    except ValueError:
        return 0

    
    
def fit_intersect_noise(d):
    # assumes input is a series
    
    noise = noise_line(d)       # noise line
    fit = lin_fit(d)            # fit line
    hdetmin = min_det_h(d)      # min detected h (probably around intersection)
    
    h_intersect = fsolve(lambda x: fit(x) - noise(x), hdetmin)
    
    return h_intersect[0]
 
 
 
# PLOTTING FUNCTIONS
           
def p(hinj=[], hrec=[], s=[], psrname='', detname='', style=sd.default_style, methods=[]):
        
    for method in methods:
        # First Calculate the interquartile range
        #(http://comments.gmane.org/gmane.comp.python.scientific.user/19755)                                                                    
        data = np.sort(hrec)                                                                                                   
        upperQuartile = stats.scoreatpercentile(data,.75)                                                                      
        lowerQuartile = stats.scoreatpercentile(data,.25)                                                                      
        IQR = upperQuartile - lowerQuartile
    
    
        # Get ideal bin size
        #(http://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule)
#         fdsize = 3.49*np.std(data)*len(data)**(-1./3.)
        fdsize = 2 * IQR * len(data)**(-1./3.)
            
        #Get number of bins
        #(http://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram-for-n-where-n-ranges-from-30)
        num_bins = int((np.amax(data) - np.amin(data))/fdsize)

        cumfreqs, lowlim, binsize, _ = stats.cumfreq(data, num_bins)
        pv = [1. - cdf/max(cumfreqs) for cdf in cumfreqs]
        bins = np.linspace(lowlim, num_bins*binsize, num_bins)

        plt.plot(bins, pv, style, color=sd.sd.pltcolor[method], label=method)
        
        plt.yscale('log')

    plt.title(detname + ' PSR ' + psrname)

    plt.xlabel('$h_{rec}$')
    plt.ylabel('1 - CDF (log scale)')

    plt.legend(numpoints=1)
    plt.savefig('plots/p_' + detname + '_' + psrname, bbox_inches='tight')
    
    print 'Plotted and saved in: ',
    print 'plots/p_' + detname + '_' + psrname
    plt.close()
    

# compound plots  
def hinjs(hinj=[], hrec=[], s=[], style=sd.default_style, methods=[]):

    method_plot = plt.plot(hinj, s[methods], style)
        
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=sd.pltcolor[methods[i]], label=methods[i])

    plt.xlabel('$h_{inj}$')
    plt.ylabel('Significance')

    plt.legend(numpoints=1, loc=2)
    

def hinjlins(hinj=[], hrec=[], s=[], style=sd.default_style, methods=[]):
    
    lins = s.applymap(math.sqrt)
    
    method_plot = plt.plot(hinj, lins[methods], style)
    
    slope = []
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=sd.pltcolor[methods[i]], label=methods[i])
        
        # fit
        line = lin_fit(lins[methods[i]])
        slope += [line(1)]
        # plot
        plt.plot(s.index, line(s.index), ls='-',color=sd.pltcolor[methods[i]], linewidth=.05)

    
    # find noise for max slope line
    n_line = noise_line(lins[methods[np.argmax(slope)]])
    plt.plot(s.index, n_line(s.index), '--', color='.95')
    
    linsmax = max([np.amax(s[m].map(lambda x: np.sqrt(x))) for m in methods])

    plt.xlabel('$h_{inj}$')
    plt.ylabel('Linearized significance')
    plt.ylim(ymax=linsmax)
    
    plt.legend(numpoints=1, loc=2)


def hinjrec(hinj=[], hrec=[], s=[], style=sd.default_style, methods=[]):

    method_plot = plt.plot(hinj, hrec[methods], style)
    
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=sd.pltcolor[methods[i]], label=methods[i])

    plt.legend(numpoints=1, loc=2)

    hmax = max([np.amax(hrec[m]) for m in methods])
    
    plt.xlabel('$h_{inj}$')
    plt.ylabel('$h_{rec}$')
    plt.ylim(ymax=hmax)
    

# simple plots
def inj(hinj=[], hrec=[], s=[], style=sd.default_style, methods=[]):
    plt.plot(hinj, style)
    plt.ylabel('$h_{inj}$')
    plt.xlabel('Instantiation')
    

def rec(hinj=[], hrec=[], s=[], style=sd.default_style, methods=[]):

    method_plot = plt.plot(h[methods], style)
    
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=sd.pltcolor[methods[i]], label=methods[i])

    plt.ylabel('$h_{rec}$')
    plt.ylabel('Instantiations')

    plt.legend(numpoints=1)
    

def sig(hinj=[], hrec=[], s=[], style=sd.default_style, methods=[]):
    
    method_plot = plt.plot(s[methods], style)
    
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=sd.pltcolor[methods[i]], label=methods[i])

    plt.ylabel('$Significance$')
    plt.ylabel('Instantiations')

    plt.legend(numpoints=1)
    

def original(detector, psr, location='files/data/'):
    try:
        d = pd.HDFStore(location + 'dataPitkin_' + detector + '.hdf5', 'r')
        a = d[psr].tolist()
        b = [x.real for x in a]
        plt.plot(b[:1000])
    finally:
        d.close()
   
    
def sigma(detector, psrname, location='files/analysis/psrSegmentSigma', compare=False):

    if compare:
        print 'comparing'
        sgS5 = pd.read_table('files/data/source/S5/sigmaS5_' + detector, names=[None, 'matlab'], sep='\s+', header=None, index_col=0)
        sgS5.plot(style='r+')
    else:
        pass
        
#     sg1 = pd.HDFStore(location + '_mat', 'r')
    
    try:
        sg = pd.HDFStore(location, 'r') #location, 'r')
        sg[psrname].plot(style='g+')
    
        plt.legend(numpoints=1)
        plt.title(psrname + detector + 'daily std')
        plt.xlabel('GPS time')
        plt.ylabel('Standard deviation')

        plt.savefig('plots/std_' + detector + '_' + psrname, bbox_inches='tight')
        print 'Plotted and saved in: ',
        print 'plots/' + detector + '_' + psrname
        plt.close()
    
    finally:
        sg.close()
                
        
def p_original(detector, psr, location='files/remote/source/'):
    d = pd.HDFStore(location + 'dataPitkin_' + detector + '.hdf5', 'r')
    a = d[psr].tolist()
    b = [abs(x) for x in a]

    # First Calculate the interquartile range
    #(http://comments.gmane.org/gmane.comp.python.scientific.user/19755)                                                                    
    data = np.sort(d[psr].tolist())                                                                                                   
    upperQuartile = stats.scoreatpercentile(data,.75)                                                                      
    lowerQuartile = stats.scoreatpercentile(data,.25)                                                                      
    IQR = upperQuartile - lowerQuartile
    
    
        # Get ideal bin size
        #(http://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule)
#         fdsize = 3.49*np.std(data)*len(data)**(-1./3.)
    fdsize = 2 * IQR * len(data)**(-1./3.)
        
    #Get number of bins
    #(http://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram-for-n-where-n-ranges-from-30)
    num_bins = int((np.amax(data) - np.amin(data))/fdsize)

    cumfreqs, lowlim, binsize, _ = stats.cumfreq(data, num_bins)
    pv = [1. - cdf/max(cumfreqs) for cdf in cumfreqs]
    bins = np.linspace(lowlim, num_bins*binsize, num_bins)

    plt.plot(bins, pv, style, color=sd.pltcolor[method], label=method)
    
    plt.yscale('log')

    plt.title(detname + ' PSR ' + psrname)

    plt.xlabel('$h_{rec}$')
    plt.ylabel('1 - CDF (log scale)')

    plt.legend(numpoints=1)
    plt.savefig('plots/p_' + detname + '_' + psrname, bbox_inches='tight')
    
    print 'Plotted and saved in: ',
    print 'plots/p_' + detname + '_' + psrname
    plt.close()