import templates as tmp
import matplotlib.pyplot as plt
import numpy as np
import sidereal as sd
import pandas as pd


# '''
# PLOT ANTENNA PATTERNS
# '''
#
# days=1
# t = np.array([x+ 630720013 for x in range(0, int(days*sd.ss), 60)])
# 
# 
# det = 'H1'
# psrlist = ['J0534+2200', 'J0835-4510']
# 
# pols = sd.names # sd.names is a list of all polarization names 
# 
# responseDict = {}
# 
# for psr in psrlist:
#     print psr
#     responseDict[psr] = tmp.Response(psr, det, t, pols)
#     print 'ok'
#     responseDict[psr].get()
# 
# for pol in pols:
#     for psr in psrlist:
#         response = responseDict[psr]
#         getattr(response, pol).plot()
#        
#     plt.gca().legend_=None
#     plt.draw()
#     plt.ylim(ymax=1,ymin=-1)
#     plt.xlim(xmax=t[-1])
#     plt.gca().xaxis.set_major_locator(plt.NullLocator())
# 
#     plt.savefig('files/plots/' + pol +'.pdf', bbox_inches='tight')
#     print 'Plotted and saved in: ',
#     print 'files/plots/response_' + det + '_' + pol
#     plt.close()
    
'''
PLOT FINEHET DATA (real part)
'''
dataloc = 'files/remote/source/S5/' # location of finehet data
psr = 'J0534+2200'

detlist = ['H1', 'L1']

plotloc = 

for det in detlist:
    try:
        d = pd.HDFStore(dataloc + 'finehet_' + psr + '_' + det + '.hdf5', 'r')
        a = d[psr].tolist()
        b = [x.real for x in a]
        plt.plot()
        
    finally:
        d.close()