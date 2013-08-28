import numpy as np
import pandas as pd
from scipy import optimize


# CONSTANTS
ss = 86164.0905             # Seconds in a sidereal day

w = 2*np.pi/ss              # Sidereal angular frequency of Earth

rE = 6378.137e3             # Earth radius (m)


periodLIGO = 1/16384.       # LIGO data sampling period (s), from M. Pitkin

c = 299792458.              # Speed of light (m/s)

###

# CONVERSIONS
def hmsformat(*args):
    if len(args[0]) == 1:
        # Assume hh:mm:ss format
        if type(args) != str:
            argument = args[0][0]
        else:
            argument = args
            
        hms = argument.split(':')
        h, m, s = [float(x) for x in hms]
        
    elif len(args[0]) == 3:
        h, m, s = [x for x in args[0]]
    else:
        print 'ERROR in hmsformat: can\'t take %d arguments' % len(args)    
    return h, m, s
    

def hms_rad(*args):
    # Converts hours, minutes, seconds to radians using the sidereal angular frequency of the Earth
    h, m, s = hmsformat(args)
    sec = s + 60*(m + 60*h)
    return sec*w
    
    
def dms_deg(*args):
    # Converts degrees, minutes, seconds to decimal degrees
    d, m, s = hmsformat(args)
    return d + m/60 + s/(60**2)
    
def masyr_rads(masyr):
    # Converts milliarcseconds/yr to radians/second
    asyr = masyr * 10 ** -3                     # mas/yr to arcseconds/yr 
    radyr = asyr * np.pi / 648000.              # as/yr to rad/yr (Wikipedia)
    rads = radyr / ss                           # rad/yr to rad/s
    return rads
    
def mjd_gps(mjd):
    # Converts MJD time to GPS time (taken from LALBarycenter.c line 749)
    tgps = 86400.*(mjd - 44244.) - 51.184
    return tgps
    
    
###

# PULSAR PARAMETERS
paramNames = [
                '#',
                None,
                'RAS',
                'RAS error',
                'DEC',
                'DEC error',
            ]                

extraParamNames = [None, 'POL', 'POL error', 'INC', 'INC error']

###

# DETECTOR PARAMETERS, all angles in radians. (Source: PRD 58, 063001 p3.)
detectors = pd.DataFrame({
        'LHO': {
                'lat': 0.8107054375513661,
                'lon': -2.084097659806429,
                'x_east': 2.199114857512855,
                'arm_ang': np.pi/2.
                },
    
        'LLO': {
                'lat': 0.5333726194094671, 
                'lon': -1.584235362035253,
                'x_east': 3.4505159311927893,
                'arm_ang': np.pi/2.
                },
    
        'GEO': {
                'lat': 0.9119345341670372,
                'lon': 0.17121679962064373,
                'x_east': 0.37716565135597474,
                'arm_ang': 1.646369083406251
                },
    
        'VIR': {
                'lat': 0.761487152645126,
                'lon': 0.1832595714594046,
                'x_east': 1.2479104151759457,
                'arm_ang': 1.5707963267948966
                },
    
        'TAM': {
                'lat': 0.6227334771115768,
                'lon': 2.4354324382328874,
                'x_east': 3.141592653589793,
                'arm_ang': 1.5707963267948966
                }
        })
        
def detnames(d):
    if d in ['H1', 'H2']:
        det = 'LHO'
    elif d == 'L1':
        det = 'LHO'
    elif d == 'V1':
        det = 'VIR'
    elif d in['LHO', 'LLO', 'VIR']:
        det = d
    else:
        print 'ERROR: %r is an unknown detector name.' % d
        exit()
    return det

###

# BASES
names = ['pl', 'cr', 'xz', 'yz', 'br']              # polarization names

polNames = {
            'pl': 'plus',
            'cr': 'cross',
            'xz': 'vector_x',
            'yz': 'vector_y',
            'br': 'breathing',
            'lo': 'longitudinal'
            }
            
# tuples indicating which vectors need to be multiplied together
# note that detector vectors should be listed second for broadcasting reasons
polComponents = {
                'pl' : [('wx','dx'), ('wx','dy'), ('wy','dx'), ('wy','dy')],
                'cr' : [('wx','dx'), ('wy','dx'), ('wx','dy'), ('wy','dy')],
                'xz' : [('wx','dx'), ('wz','dx'), ('wx','dy'), ('wz','dy')],
                'yz' : [('wy','dx'), ('wz','dx'), ('wy','dy'), ('wz','dy')],
                'br' : [('wx','dx'), ('wx','dy'), ('wy','dx'), ('wy','dy')],
                'lo' : [('wz','dx'), ('wz','dy')]
                }
                
###

# TEMPLATES
tempNames = ['GR', 'G4v', 'AP', 'Sid', 'GRs']       # template names

aps = {
        'GR' : ['pl', 'cr'],
        'G4v': ['xz', 'yz'],
        'GRs' : ['pl', 'cr', 'br'],
        'AP' : ['pl', 'cr', 'xz', 'yz', 'br'],
        'Sid': []
    } 

pcat = {
        'p' : np.pi/2.,
        'm' : - np.pi/2.,
        '0' : 0
        }

def phase(p):
    if isinstance(p, basestring):
        return pcat[p]
    else:
        return p
        
def phase2(pinj):
    if pinj == np.pi/2.:
        return 'p'
    elif pinj == -np.pi/2.:
        return 'm'
    elif pinj == 0:
        return '0'
    else:
        return str(pinj)
        
        
# PLOTS
pltcolor = {
            'GR' : 'g',
            'G4v': 'r',
            'AP' : 'c',
            'Sid': 'b'
            }
            
default_style = '+'


# RESULTS
statkinds = [
            'lin s slope',
            'lin s noise',
            'lin s inter',
            'h rec slope',
            'h rec noise',
            'h rec inter',
            'min inj det'
            ]
            
statlabels = {
            'lin s slope' : 'Lin. sig. fit slope',
            'lin s noise' : 'Lin. sig. noise level',
            'lin s inter' : 'Lin. sig. fit intersect noise level',
            'h rec slope' : '$h_{rec}$ fit slope',
            'h rec noise' : '$h_{rec}$ noise level',
            'h rec inter' : '$h_{rec}$ fit intersect noise level',
            'min inj det' : 'Strength of minimum detection'
            }