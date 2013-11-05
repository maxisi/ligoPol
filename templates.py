import math
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from urllib2 import urlopen
from collections import namedtuple

import sys, os
from sys import argv, exit
import copy
from time import time
from datetime import datetime

import sidereal as sd
import paths

reload(sd)

# Set up logging (from Logging Cookbook, Python online resources)
import logging

# set up logging to file
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='temp/templates_d'+datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+'.log',
                    filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)
    
    
#  vector container to make it easy to move vectors around
Vectors = namedtuple('Vectors', ['dx', 'dy', 'wx', 'wy', 'wz'])


## SOURCE
def currentCatalogue():
    '''
    Returns contents of current pulsar catalogue. Mainly for maintenance reasons.
    '''
    try:
        return pd.read_pickle(paths.psrcat)
    except:
        print 'No pulsar catalogue found.'
        exit()
        

class Source(object):
    '''
    Contains source information for a single pulsar.
    Inputs:
        psr
        loadvec (indicates whether to load src vectors at start, default: True)      [OPT]
    '''
    
    def __init__(self, psr, loadvec=True):
        self.log = logging.getLogger('Source')
        self.log.debug('Initializing source.')
        
        self.psr = psr
        self.npsrs = 1
        self.path = paths.vectors + 'srcVec_' + psr
                     
        # If necessary, rebuild catalogue; otherwise, just load catalogue.
        self.log.debug('Looking for catalogue.')
        try:
            f = open(paths.textfromATNF, 'r')
            pd.read_pickle(paths.psrcat)
        except IOError:
            self.log.warning('No PSR catalogue found.')
            self.build_catalogue()
            f = open(paths.textfromATNF, 'r')
        finally:
            self.log.debug('Checking catalogue.')
            f_text = f.read()
            f.close()
        
        if self.psr not in f_text:
            self.log.debug('PSR not in catalogue.')
            self.build_catalogue()
        
        self.log.debug('Reading catalogue.')
        psrcat = pd.read_pickle(paths.psrcat)
        self.param = psrcat.ix[self.psr]
        
        if loadvec: self.loadVectors()
        

    def build_catalogue(self, extrapsrs=[]):
        '''
        Gets location parameters for pulsars in input list from ATNF online catalogue.
        Creates and pickles corresponding pandas DataFrame 'pulsar_parameters'.
        '''
        self.log.info('Building PSR catalogue.')
        psrs = [self.psr] + extrapsrs
        
        self.log.debug('Creating ATNF url.')
        def atnfurl(psr_list):
            pre = 'http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.47&JName=JName&RaJ=RaJ&DecJ=DecJ&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&'
            post = '&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+errors&no_value=*&nohead=nohead&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=36&table_bottom.y=16'
            names= 'pulsar_names='
    
            for psr in psr_list:
                names+=psr.replace('+', '%2B')
        
                if psr != psr_list[-1]:
                    names+='%0D%0A'
            url = "%(pre)s%(names)s%(post)s" % locals()
    
            if len(url)>2000:
                self.log.warning('WARNING! URL %d characters!' % len(url))
    
            return url

        # Get data
        url = atnfurl(psrs)
        self.log.debug('Parsing webpage.')
        soup = BeautifulSoup(urlopen(url)) # get webpage and parse
        text = str(soup.pre.string)

        self.log.debug('Saving text to file.') 
        f = open(paths.textfromATNF, 'w')
        f.write(text)
        f.close()

        self.log.debug('Create DataFrame from web text.')
        psr = pd.read_table(paths.textfromATNF, sep='\s+', comment='*', names=sd.paramNames, header=None, skiprows=1, index_col=1)
        psrcat=psr.drop("#", axis=1)
    
        # Format
        self.log.debug('Formatting.')
        formatRAS  = lambda x: sd.hms_rad(x)
        formatRASe = lambda x: sd.hms_rad(0., 0., x)
        formatDEC  = lambda y: np.radians(sd.dms_deg(y))
        formatDECe = lambda x: np.radians(sd.dms_deg(0., 0., x))
                
        psrcat['RAS'] = psrcat['RAS'].map(formatRAS)                    # hms -> rad
        psrcat['RAS error'] = psrcat['RAS error'].map(formatRASe)
        psrcat['DEC'] = psrcat['DEC'].map(formatDEC)                    # dms -> rad
        psrcat['DEC error'] = psrcat['DEC error'].map(formatDECe) 
    
    
        self.log.debug('Checking extra parameters.')
        extra = pd.read_table(paths.psrextra, sep=',', header=None, index_col=[0], names=sd.extraParamNames)
        psrcatComplete = pd.merge(psrcat,extra, left_index=True, right_index=True, how='outer')
    
        psrcatComplete['POL error']=psrcatComplete['POL error'].fillna(value=np.pi/4.)
        psrcatComplete['INC error']=psrcatComplete['INC error'].fillna(value=np.pi/4.)

        self.log.debug('Filling missing values.')
        psrcatComplete.fillna(value=0, inplace=True)
        self.log.debug('Saving.')
        psrcatComplete.save(paths.psrcat)
     
     
    def loadVectors(self):
        '''
        Loads source vectors from file.
        '''
        
        self.log.debug('Loading source vectors')
        
        try:
            file = pd.HDFStore(self.path, 'r')
            self.wx = file['wx']
            self.wy = file['wy']
            self.wz = file['wz']
            file.close()
        except IOError:
            self.log.warning('No src vectors found.')
            self.createVectors()
           
            
    def createVectors(self):
        # Return wave vectors for all sources listed.
        
        self.log.debug('Creating src vectors.')
        
        north = np.array([0, 0, 1])
        
        self.log.debug('wz')
        # take source location vector components in celestial coordinates and invert direction multiplying by -1 to get wave vector wz
        wz = [-math.cos(self.param['DEC'])*math.cos(self.param['RAS']), -math.cos(self.param['DEC'])*math.sin(self.param['RAS']), -math.sin(self.param['DEC'])] 
        self.wz = pd.Series(wz, name=self.psr, index=['x', 'y', 'z'])
        
        self.log.debug('wy')
        wy = np.cross(wz, north) 
        wy /= np.sqrt(np.sum(wy ** 2))
        self.wy = pd.Series(wy, name=self.psr, index=['x','y','z'])
        
        self.log.debug('wx')
        wx = np.cross(wy, wz)
        wx /= np.sqrt(np.sum(wx ** 2))
        self.wx = pd.Series(wx, name=self.psr, index=['x','y','z'])
        
        self.log.debug('Saving src vectors.')
        try:
            f = pd.HDFStore(self.path, 'w')
            f['wx'] = self.wx
            f['wy'] = self.wy
            f['wz'] = self.wz
        finally:
            f.close()
         
###         

# DETECTOR
class Detector(object):
    
    def __init__(self, d, t=[]):
        self.log = logging.getLogger('Detector')
        self.log.debug('Initializing detector.')

        self.id = d  
        self.name = sd.detnames(d)
        self.param = sd.detectors[self.name]
        self.t = np.array(t).astype(int)
        self.nentries = len(t)
        self.path = paths.vectors + 'detVec_' + self.name
        
        if self.t!=[]: self.loadVectors()
        

    def fileload(self):
        '''
        Loads detector arm vectors from file.
        '''
        self.log.debug('Loading detector vectors.')
        
        try:
            file = pd.HDFStore(self.path, 'r')
            self.dx = file['dx']
            self.dy = file['dy']
            self.dz = file['dz']
            file.close() 
        except IOError:
            self.log.warning('Detector vectors not found.')
            self.createVectors()


    def loadVectors(self):
        '''
        Loads detector arm vectors if necessary.
        '''
        # Check if vectors are stored in file'
        self.fileload()
        
        self.log.debug('Checking vector health.')
        
        # Check data type'
        try:
            self.dx.columns
        except AttributeError:
            self.log.warning('Wrong data type.')
            self.fileload()
        
#         try:
        if set(self.dx.index)==set(self.t):
            self.log.debug('All times present.')
        else:
            self.log.debug('Incomplete time vector.')
            self.createVectors()

    def createVectors(self):
        '''
        Returns arm vectors in Cartesian sidereal coordinates.
        '''
        
        self.log.debug('Creating detector vectors.')
        northPole = pd.Series(np.array([0, 0, 1]), index=['x', 'y', 'z'])    # Earth center to North pole
        
        self.log.debug('Retrieving detector parameters.')
        lat = self.param.ix['lat']
        lon = self.param.ix['lon']
        x_east = self.param.ix['x_east']
        arm_ang = self.param.ix['arm_ang']
        
        coords = ['x', 'y', 'z']
        
        t = np.array(self.t).astype(int)

        length = self.nentries
        
        self.log.debug('Computing local mean sidereal time.')
        # Angle between detector and Aries (vernal equinox) at time t
        # fiducial GPS time t0=630763213 (12hUT1 1/1/2000, JD245154).
        # See http://aa.usno.navy.mil/faq/docs/GAST.php
        offset = 67310.5484088 * sd.w   # Aries-Greenwich angle at fiducial time (GMST)
        lmst = offset + sd.w*(t-630763213) + lon # (LMST) rows: t, columns: det
        th = pd.Series(lmst, index=t)
        
        self.log.debug('Zenith.')
        z = {
            'x' : np.cos(lat)*np.cos(th),
            'y' : np.cos(lat)*np.sin(th),
            'z' : pd.Series([math.sin(lat)]*len(t), index=t)
            }  # [[x0, ...], [y0, ...], [z0, ...]]
          
        zenith = pd.DataFrame(z) # norm 1 already
        
        z2 = [zenith['x'], zenith['y'], zenith['z']]
        
        self.log.debug('Local East and North.')
        localEast = np.cross(northPole, z2 , axisb=0)
        localNorth = np.cross(z2, localEast, axisa=0)   
        
        self.log.debug('x-arm vector.')
        xArm = np.cos(x_east)*localEast + np.sin(x_east)*localNorth
        xArm /= np.sqrt(np.sum(xArm ** 2., axis=1))[..., None]
        
        self.log.debug('y-arm vector.')
        perp_xz = np.cross(zenith, xArm)
        yArm = xArm*np.cos(arm_ang) + perp_xz*np.sin(arm_ang) # equals perp_xz when angle between arms is 90deg
        
        self.log.debug('Turning vectors into DataFrames.')
        self.dx = pd.DataFrame(xArm, index=t, columns= ['x', 'y', 'z'])
        self.dy = pd.DataFrame(yArm, index=t, columns= ['x', 'y', 'z'])
        self.dz = sd.rE * zenith
        
        self.log.debug('Saving detector vectors.')
        try:
            f = pd.HDFStore(self.path, 'w')
            f['dx'] = self.dx
            f['dy'] = self.dy
            f['dz'] = self.dz
        finally:
            f.close()


###

# ANTENNA PATTERNS
class Polarizations(object):
    '''
    Produces detector response for different polarization given input detector and
    source vectors. Vectors must be in namedtuple form as defined above.
    Effectively computes dyadic product between wave and detector tensors.
    '''
    
    def __init__(self, vectors):
        # assuming vectors is of the Vectors container kind
        self.log = logging.getLogger('Polarizations')
        self.log.debug('Initializing polarizations.')
        self.vec = vectors
    
                    
    def product(self, polKey):
        # check all necessary vector products exist. Otherwise, create them.
        
        self.log.debug('Checking vector products')
        for pair in sd.polComponents[polKey]:
        
            pairName = pair[0] + pair[1]    # 'wxdx' = ('wx','dx')[0] + ('wx','dx')[1]
            
            self.log.debug('Checking if product is already defined.')
            if pairName not in dir(self):
                self.log.debug('Defining: ' + pairName)
                # get vectors
                v0 = getattr(self.vec, pair[0])     # v0 = self.vec.wx
                v1 = getattr(self.vec, pair[1])     # v1 = self.vec.dx
                
                # dot product (mind the order! v1 MUST be detector vector to broadcast)
                setattr(self, pairName, v1.mul(v0, axis='columns').sum(axis=1))
    
    # tensor
    def plus(self):
        self.log.debug('Creating plus polarization.')
        self.product('pl')
        pl = (self.wxdx**2 - self.wxdy**2 - self.wydx**2 + self.wydy**2)/2.
        return pl
        
    def cross(self):
        self.log.debug('Creating cross polarization.')
        self.product('cr')
        cr = self.wxdx*self.wydx - self.wxdy*self.wydy
        return cr
    
    # vector
    def vector_x(self):
        self.log.debug('Creating xz polarization.')
        self.product('xz')
        xz = self.wxdx*self.wzdx - self.wxdy*self.wzdy
        return xz
    
    def vector_y(self):
        self.log.debug('Creating yz polarization.')
        self.product('yz')
        yz = self.wydx*self.wzdx - self.wydy*self.wzdy
        return yz
    
    # scalar
    def breathing(self):
        self.log.debug('Creating breathing polarization')
        self.product('br')
        br = np.sqrt(2)*(self.wxdx**2 - self.wxdy**2 + self.wydx**2 - self.wydy**2)/2.
        return br
        # Added factor of sqrt(2) to distribute power equally among polarizations.
        # Same for longitudinal.
        
    def longitudinal(self):
        self.log.debug('Creating longitudinal polarization')
        self.product('lo')
        lo = (np.sqrt(2)*(self.wzdx**2 - self.wzdy**2))/2.
        return lo
        # Modified:1/2 (based on derivation of dyadic products using tensors shown in 
        # "Gravitational wave polarizations" by Bryant Garcia.
        # The factor of 2 shouldn't be there)


class Response(object):
    '''
    Contains response of 'det' to signals from source 'psr' of kind 'kinds', over 't'.
    
    For default polarization angle, use 'get' method. This will try to recover patterns
    from file and will generate them if not found. It rotates source vectors by
    polarization angles established on file. After 'get' is called, the patterns are
    stored in the class, so there is no need to call again in same session.
    
    To input polarization angle, call 'create' method directly with argument 'psi='.
    This will ignore existing APs (if any) and will fall back on detector vectors to
    compute response (after rotation by 'psi'). Call with 'savefile=False' to prevent
    method from storing results. This routine will compute dyadic products once, so it
    can be called multiple times a session.
    
    To do: Add option to select psi randomly from range?
    '''
    
    def __init__(self, psr, det, t, kinds, loadvectors=False):
        self.log = logging.getLogger('Response')
        self.log.debug('Initializing response')
         
        self.det = sd.detnames(det)
        self.t = np.array(t)
        
        self.psr = psr
        self.path = paths.ap + 'ap' + psr + '_' + self.det
        
        self.log.debug('Determining what template was requested.')
        if kinds in sd.tempNames:
            self.log.debug('Preset template ' + kinds)
            self.kinds = sd.aps[kinds]
            
        elif isinstance(kinds, list) and all([k in sd.tempNames for k in kinds]):
            self.log.debug('List of preset templates.')
            self.kinds = sum([sd.aps[temp] for temp in kinds], [])
            
        elif all([k in sd.names for k in kinds]) and not isinstance(kinds, basestring):
            self.log.debug('List of bases by name.')
            self.kinds = kinds
            
        elif isinstance(kinds, basestring) and kinds in sd.names:
            self.log.debug('Single basis.')
            self.kinds = [kinds]
            
        else:
            # wrong input
            self.log.error('ERROR: %s is not recognized as a valid basis.\nValid bases are:\n\t %r \n\t %r' % (kinds, sd.names, sd.aps))
            exit()
            
        self.hasvectors = False
        self.haspatterns = False
        
        if loadvectors:
            self.log.debug('Loading all vectors.')
            self.src = Source(self.psr)
            self.obs = Detector(self.det, self.t)
            self.hasvectors = True                
        else:
            self.log.debug('Not loading any vectors.')    
                
    def create(self, savefile=True, psi=[]):
    
        self.log.debug('Creating response.')
    
        # Retrieve detector vectors
        if not self.hasvectors:
            self.log.debug('Loading vectors.')
            self.src = Source(self.psr)
            self.obs = Detector(self.det, self.t)
            self.hasvectors = True
        else:
            self.log.debug('Vectors already loaded.')
            
        self.log.debug('Rotating source vectors.')
        if psi==[]:
            self.psi = self.src.param['POL']
        else:
            self.psi = psi

        wxRot = -self.src.wy*np.cos(self.psi) + self.src.wx*np.sin(self.psi)
        wyRot = self.src.wx*np.cos(self.psi) + self.src.wy*np.sin(self.psi) 
        
        self.log.debug('Packing vectors.')
        vecs = Vectors(self.obs.dx, self.obs.dy, wxRot, wyRot, self.src.wz)      
           
        # Get polarizations
        pols = Polarizations(vecs)
        [setattr(self, k, getattr(pols, sd.polNames[k])() ) for k in self.kinds]
            
        # Save if requested
        if savefile:
            self.log.debug('Saving APs.')
            try:
                apF = pd.HDFStore(self.path, 'w')
                for k in self.kinds: apF[k] = getattr(self, k)
            finally:
                apF.close()
                
        self.haspatterns = True
    
    def get(self):
        # Assumes no APs loaded. Otherwise, will re-write.
        self.log.debug('Getting APs.')
        
        self.log.debug('Checking disk.')
        try:
            apFile = pd.HDFStore(self.path, 'r')
            
        except IOError:
            self.log.debug('No patterns on file.')
            self.create()
            
        else:
            self.log.debug('File found.')
            try:
                filePols = [s.strip('/') for s in apFile.keys()]
            
                self.log.debug('Checking file health.')
                allPolsPresent = all([k in filePols for k in self.kinds])
            
                timeCoincides = set(self.t.astype(int)) == set(apFile[filePols[0]].index)
                # if file is empty, raises IndexError
            
                if allPolsPresent and timeCoincides:
                    self.log.debug('Healthy file.')
                    [setattr(self, p, apFile[p]) for p in self.kinds]
                else:
                    self.log.debug('Unhealthy file.')
                    apFile.close()
                    self.create()
                    
            except IndexError:
                self.log.debug('File is empty.')
                apFile.close()
                self.create()
            
            finally:
                apFile.close()
            self.haspatterns = True        
        
    def exportAPmatlab(psr, detname, t):
        p = ap.getAP('LHO', t)
    
    #     psrdict = {
    #                 'pl':pl.T[psr].tolist(),
    #                 'cr':cr.T[psr].tolist(),
    #                 'br':br.T[psr].tolist(),
    #                 'lo':lo.T[psr].tolist(),
    #                 'xz':xz.T[psr].tolist(),
    #                 'yz':yz.T[psr].tolist()
    #                 }
    
        sio.savemat('%(paths.ap)scrab_py' % locals(), p)


###      
        
## SIMULATE            
class Signal(object):
    '''
    Sets absolute properties of a signal: detector, PSR, template and phase difference
    between components. A time vector must also be provided.
    Includes methods to return design matrix and simulated signal given extra inputs
    of polarization and inclination angles.
    '''

    def __init__(self, detector, psr, kind, pdif, t, isloaded=False):
        
        self.log = logging.getLogger('Signal')
    
        # source and detector information
        self.detector = detector
        self.det = sd.detnames(detector)
        self.psr = psr
        
        # time
        self.t = t
        
        # signal info
        self.kind = kind            # kind of signal
        self.basis = sd.aps[kind]   # polarizations composing signal
        self.pdif = sd.phase(pdif)  # phase difference between components
        
        # get antenna patterns
        self.response = Response(psr, detector, t, kind, loadvectors=True)
        
    
    def signalinfo(self, pdif, iota=[], p=0):
        
        self.log.debug('Producing amplitudes and phase differences.')
        
        # determine inclination angle
        if iota==[]:
            self.log.debug('No preset iota, taking default.')
            iota = self.response.src.param['INC']
        else:
            self.log.debug('Taking iota = ' + str(iota))

        # return dataframe with amplitude (h) and phase (phi) for each component
        self.log.debug('Creating DataFrame.')
        info = pd.DataFrame(index=self.basis, columns=['h', 'p'])      

        self.log.debug('Setting amplitudes and phases for ' + self.kind)
        if self.kind == 'GR':
            info['h']['pl'] = (1. + np.cos(iota)**2)/2.
            info['h']['cr'] = np.cos(iota)
            
            info['p']['pl'] = p
            info['p']['cr'] = p + pdif
            
        elif self.kind == 'G4v':
            info['h']['xz'] = np.sin(iota)
            info['h']['yz'] = np.sin(iota) * math.cos(iota)
            
            info['p']['xz'] = p
            info['p']['yz'] = p + pdif
            
        elif self.kind == 'AP':
            for k in sd.names:
                info['h'][k] = 1.
                info['p'][k] = p
            

#         elif self.kind == 'GRs':
#             self.log.debug('GRs')
#             info['h']['pl'] = (1. + np.cos(iota)**2)/2.
#             info['h']['cr'] = np.cos(iota)
#             info['h']['br'] = h_s
# 
#             info['p']['pl'] = p
#             info['p']['cr'] = p + pdif
#             info['p']['br'] = p + pdif_s
            
        else:
            self.log.error('%s is not a recognized template (temp 647)' % self.kind)
            exit()
            
        return info
        

    def design_matrix(self, pol_angle, incl_angle): 
    
        self.log.debug('Building design matrix for ' + self.kind)
        
        if self.kind=='Sid':
            # no need to get antenna patterns
            
            self.log.debug('Building basis set.')
            th = pd.Series(sd.w * np.array(self.t), index=self.t)
            
            basis = [th.map(np.cos), (2.*th).map(np.cos), th.map(np.sin), (2.*th).map(np.sin)]
            
            self.log.debug('Constructing matrix.')
            dm = pd.concat(basis, axis=1, keys=['cos1', 'cos2', 'sin1', 'sin2'])
            dm['cnst']=1

        else:
            self.log.debug('Copying response.')
            response_local = copy.copy(self.response)
            
            self.log.debug('Getting antenna patterns.')
            response_local.create(psi=pol_angle, savefile=False)
            
            self.log.debug('Getting amplitude info.')
            info = self.signalinfo(0, iota=incl_angle)
            
            self.log.debug('Constructing matrix.')
            dmDict = {pol: getattr(response_local, pol) * info['h'][pol] / 2. for pol in self.basis}
            dm = pd.DataFrame(dmDict).dropna(axis=1) # DF cols: comp. names, index: t.
        
        return dm
         
            
    def simulate(self, pol_angle, incl_angle, phase=0):
        '''
        Simulates a signal based given polarization and inclination angles.
        Warning: Does not scale output signal, need to multiply output by h0.
        '''
        
        self.log.debug('Simulating %s signal.' % self.kind)
    
        if self.kind=='Sid':
            self.log.warning('Cannot simulate "Sid". Change Signal.kind')
        else:
            # form design matrix
            dm = self.design_matrix(pol_angle, incl_angle)

            # get phase info
            info = self.signalinfo(self.pdif, iota=1., p=phase)
            
            self.log.debug('Constructing signal.')
            # raise phases to exponent (apply),
            # multiply amplitudes and phases (mul),
            # drop extra columns which come in from phases (dropna),
            # add up columns (sum).
            raisePhi = lambda x: np.exp(1j*x)
            s = dm.mul(info['p'].map(raisePhi)).dropna(axis=1).sum(axis=1)

            return s