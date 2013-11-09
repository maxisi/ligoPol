from sys import exit
import matplotlib.pyplot as plt
import os
import datetime
from time import time
from datetime import datetime
import sys

import numpy as np
import pandas as pd
import random
import math

import templates
import sidereal as sd
import paths
import psrplot

reload(sd)
reload(templates)
reload(psrplot)

# Set up logging (from Logging Cookbook, Python online resources)
import logging

# set up logging to file
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='temp/process'+datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+'.log',
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


def het(vector, f, *arg):
    '''
    Heterodynes vector at frequencies f. Preferred input for vector is series indexed
    over time; f can be a list or an array. Returns DataFrame.
    '''
    hetlog = logging.getLogger('heterodyne')

    hetlog.debug('Ready to heterodyne.')
    
    if len(arg)==0:
        try:
            t = vector.index.tolist()
        except AttributeError:
            hetlog.error('No time vector for heterodyne.')
            exit(0)
            
    elif len(arg)==1:
        t = arg[0]
        
    else:
        hetlog.error('Het needs input time or indexed vector, not %d extra arguments.' % len(arg))
        exit(0)
    
    temp = np.exp(2*np.pi*1j*np.multiply.outer(f, t))
    
    try:
        template = pd.DataFrame(temp, columns=t)
    except ValueError:
        template = pd.Series(temp, index=t)
    
    rh = vector*template
    hetlog.debug('Herodyne succeeded.')
    return rh.T


def listpsrs(detector, directory):
    # Retrieves names of PSRs heterodyned for "detector" in "directory".
    
    pre = 'finehet_'
    pst = '_' + detector

    fnames = os.listdir(directory)

    allpsrs = [f.strip(pre).strip(pst) for f in fnames] 

    return allpsrs


class Data(object):
    '''
    Holds original data and related information.
    '''
    def __init__(self, detector, psr):        
        self.detector = detector
        self.det = sd.detnames(detector)
        self.psr = psr
        
        # data info
        self.datadir = paths.importedData + self.psr + '_' + self.detector + '.hdf5'
        self.seedname = 'finehet_' + self.psr + '_' + self.detector
        
        self.log = logging.getLogger('Data')

    

    def imp(self):
        '''
        Return DF with original data (col: PSR; index: t). Assuming execution on ATLAS.
        '''
        
        self.log.info('Importing seed')
        
        struct = '/data' + self.detector + '/' + self.seedname
        pathOptions = [
                     paths.originalData + struct,
                     paths.originalData + '/' + self.psr + '_' + self.detector + struct
                     ]
        
        try:
            d = pd.HDFStore(self.datadir, 'w')
                
            for p in pathOptions:
                try:
                    dParts = pd.read_table(p, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                except IOError:
                    pass
            
            # check file was found
            try:
                dParts
                self.log.debug('Import success.')
            except NameError:
                self.log.error('Could not find %s data for PSR %s in: %r' % (self.detector, self.psr, p))
                raise IOError

            self.finehet = dParts['Re']+dParts['Im']*1j
            d[self.psr] = self.finehet
            
        finally:
            d.close()


    def get(self):
        '''
        Retrieves original heterodyned data for pulsars in list.
        Imports data from M.Pitkin if necessary.
        '''
        
        self.log.info('Getting data.')
        
        try:
            d = pd.HDFStore(self.datadir, 'r')
        
            try:
                self.finehet = d[self.psr]
            except KeyError:
                self.log.warning('File is empty or corrupted.')
                d.close()
                self.imp()
            else:
                self.log.debug('Data retrieved successfully.')
                d.close()
        
        except IOError:
            self.log.warning('File not found.')
            self.imp()
                

class Background(object):
    '''
    Manages background files for a given detector and source: gets and creates.
    '''
    def __init__(self, detector, psr, freq, filesize=100):
        # data
        self.seed = Data(detector, psr)
        self.seed.get()
        
        # instantiation set info
        self.freq = freq              # frequencies to heterodyne
        self.filesize = filesize      # number of series per file. Adjust!
        
        self.nsets = int(len(freq)/filesize)    # final number of files
        if self.nsets<1: self.nsets = 1         # minimum number is 1
        
        # create frequency sets
        self.fset = {n : freq[n*filesize:min(len(freq),(n+1)*filesize)] for n in range(self.nsets)}
        
        # storing info
        self.dir = paths.rhB + self.seed.det + '/' + psr + '/'
        self.name = 'back_' + psr + '_' + self.seed.det + '_'
        self.path = self.dir + self.name
        
        self.log = logging.getLogger('Background')

    
    def writelog(self):
        now = datetime.datetime.now()
        comments = '# ' + self.seed.detector + '\n# ' + self.seed.psr + '\n# ' + str(now) + '\n'
        fileinfo = 'nsets\tfilesize\n' + str(self.nsets) + '\t' + str(self.filesize)
        
        try:
            f = open(self.dir + 'log.txt', 'w')
            f.write(comments + fileinfo)
        finally:
            f.close()  
            
               
    def create(self):
        '''
        Re heterodynes and saves data at frequencies f. Number of heterodynes is determined by
        f and data can be for more than one pulsar
        '''
        
        self.log.info('Creating background.')
        
        # create background directory
        try:
            os.makedirs(self.dir)
        except OSError:
            pass
            
        # create background files
        for n in range(self.nsets):
            path = self.dir + self.name + str(n)
            
            try:
                rh = pd.HDFStore(path, 'w')
                rh[self.seed.psr] = het(self.seed.finehet, self.fset[n])
            finally:
                rh.close()
        
        self.writelog()
        self.log.info('Background created.')

    def get(self):
        '''
        Checks background required for search exits and creates it if needed.
        Returns filename list.
        '''
        self.log.info('Getting background.')
        # read log
        try:
            self.log.debug('Reading log.')
            readme = pd.read_table(self.dir + 'log.txt', sep='\s+', skiprows=3)
            log_nfiles = readme['nsets'].ix[0]
            log_filesize = readme['filesize'].ix[0]
            log_nfreq = log_nfiles * log_filesize
            
            # get actual number of background files in directory
            files = [name for name in os.listdir(self.dir) if 'back' in name]
            nfiles = len(files)
            
            if nfiles!=log_nfiles or log_nfreq!=len(self.freq) or log_filesize!=self.filesize:
                self.log.warning('Background log inconsistent.')
                self.create()
        except IOError:
            # no log found
            self.log.warning('No background log found.')
            self.create()


class Sigma(object):
    def __init__(self, detector, psr, data, justload=False):
        self.log = logging.getLogger('Sigma')
        self.log.debug('Initializing Sigma.')
        
        self.detector = detector
        self.psr = psr
        
        self.data =  data
        
        self.dir = paths.sigma + '/' + self.detector + '/'
        self.name = 'segsigma_' + self.psr + '_' + self.detector
        self.path = self.dir + self.name
        
        self.justload = justload    # if true, will not compute.
        
        self.get()
        
        
    def create(self):
        '''
        Splits data into day-long segments and returns their standard deviation.
        '''
        self.log.info('Computing segment std.')
        
        data  = self.data
          
        t = data.index
        interval_length= sd.ss
        self.log.debug('Taking std over %f second-long intervals.' % interval_length)

        # Slice up data into day-long bins and get groupby stats (see Ch 9 of Python for Data Analysis).
        bins = np.arange(t[0]-interval_length, t[-1]+interval_length, interval_length)
        slices = pd.cut(t, bins, right=False)
        self.log.debug('Segmented.')

        def getsigma(group):
    #         s = np.std(group)
            g = np.array(group.tolist())
            s = np.std(g)
            return s
            #return group.std(ddof=0) # this is pd unbiased 1/(n-1), should use np.std 1/n?
        
        grouped = data.groupby(slices) # groups by bin
        self.log.debug('Data grouped.')
        
        sigmagroups= grouped.apply(getsigma) # gets std for each bin
        self.log.debug('STD taken.')

        # Create standard deviation time series 
        s = [sigmagroups.ix[slices.labels[t_index]] for t_index in range(0,len(t)) ]
        self.std = pd.Series(s, index=t)
        self.log.debug('Done.')
    
    
    def get(self):
        self.log.info('Retrieving segment standard deviation.' % locals())
        
        try:
            s = pd.HDFStore(self.path)
            try:
                self.std = s[self.psr]
                self.log.debug('File found.')
                
                # check times coincide
                if not self.justload:
                    self.log.debug('Comparing times in data and std.')
                    
                    if not set(self.std.index)==set(self.data.index):
                        self.log.warning('Inconsistent times.')
                        self.create()
                        
                        # save
                        self.log.debug('Saving.')
                        s.close()
                        s = pd.HDFStore(self.path, 'w')
                        s[self.psr] = self.std
                    
            except KeyError:
                self.log.warning('PSR not in file.')
                self.create()
                
                # save
                self.log.debug('Saving.')
                s[self.psr] = self.std
                
        except IOError:
            self.log.warning('Creating std directory.')
            os.makedirs(self.dir)
            self.create()
            # save
            self.log.debug('Saving.')
            s = pd.HDFStore(self.path, 'w')
            s[self.psr] = self.std
            
        finally:
            s.close()
        
        self.log.info('Sigma is ready.')
        
    def plot(self, extra_name=''):
        
        self.std.plot(style='+')
        plt.title('Daily standard deviation for ' + self.detector + ' ' + self.psr + ' data ' + extra_name)
        plt.xlabel('GPS time (s)')
        plt.ylabel('$\sigma$')
        
        # save
        save_dir = paths.plots + '/' + self.detector + '/sigma/'
        save_name = self.name + extra_name + '.png'
        try:
            plt.savefig(save_dir + save_name, bbox_inches='tight')
        except IOError:
            os.makedirs(save_dir)
            plt.savefig(save_dir + save_name, bbox_inches='tight')


class Results(object):
    '''
    Holds search results and contains methods to save them.
    '''
    def __init__(self, detector, psr, methods=[], hinj=[], pdif_s=None, kind=None, pdif=None):
        # system
        self.detector = detector
        self.psr = psr
        
        # search
        self.methods = methods
        
        # injection
        self.hinj = hinj
        self.kind = kind
        self.pdif = pdif
        self.pdif_s = pdif_s
        
        # containers
        self.h = pd.DataFrame(columns = methods, index=range(len(hinj)))
        self.s = pd.DataFrame(columns = methods, index=range(len(hinj)))
        
        self.stats = pd.DataFrame(index=sd.statkinds, columns = methods)
        
        # saving
        self.dir = paths.results + self.psr + '/' + self.detector + '/' 
        self.name = self.psr + '_' + self.detector + '_' + self.kind + '_' + sd.phase2(pdif)
        self.path = self.dir + self.name
        
        self.issaved = False
        self.log = logging.getLogger('Results')
        
    def save(self, extra_name=''):
        
        self.log.info('Saving.')
        
        self.h.index = self.hinj
        self.s.index = self.hinj
        
        self.getstats()
       
        try:
            os.makedirs(self.dir)
        except OSError:
            pass
            
        try:
            f = pd.HDFStore(self.path + extra_name, 'w')
            f['h'] = self.h
            f['s'] = self.s
            f['stats']= self.stats
            self.issaved = True
        except:
            self.log.error("Failed to save (%s)." % self.path + extra_name)
            self.issaved = False
        else:
            f.close()
        
    def load(self):
        self.log.info('Loading.')
        try:
            f = pd.HDFStore(self.path, 'r')
            self.h = f['h']
            self.s = f['s']
        except:
            self.log.error("Failed to load (%s)." % self.path + extra_name)
        finally:
            f.close()            


    def plots(self, pltType, extra_name=''):
    
        self.log.info('Plotting.')
        
        header = self.kind + sd.phase2(self.pdif) + ' injections on ' + self.detector + ' data for ' + self.psr + ' ' + extra_name
          
        getattr(psrplot, pltType)(hinj=self.h.index, hrec=self.h, s=self.s, methods=self.methods)
        
        plt.title(header)
        
        pltdir = paths.plots + self.detector + '/' + self.kind + '/' + pltType + '/'
        pltname = self.detector + '_' + self.kind + '_' + sd.phase2(self.pdif) + '_' + pltType + extra_name
        save_to = pltdir + pltname
        
        try:
            os.makedirs(pltdir)
        except OSError:
            pass
            
        plt.savefig(save_to, bbox_inches='tight')
        plt.close()
        
        self.log.info('Plot saved to:\n %(save_to)s' % locals())
    
        
    def getstats(self, plot=False, store=True):

        self.log.info('Computing statistics.')        
        lins = self.s.applymap(math.sqrt)
        
        for m in self.methods:
            self.stats[m]['min inj det'] = psrplot.min_det_h(lins[m])
            self.stats[m]['lin s slope'] = psrplot.lin_fit(lins[m])(1)
            self.stats[m]['lin s noise'] = psrplot.noise_line(lins[m])(1)
            self.stats[m]['lin s inter'] = psrplot.fit_intersect_noise(lins[m])
            self.stats[m]['h rec noise'] = psrplot.noise_line(self.h[m])(1)
            self.stats[m]['h rec slope'] = psrplot.lin_fit(self.h[m])(1)
            self.stats[m]['h rec inter'] = psrplot.fit_intersect_noise(self.h[m])
                    

class Frequentist(object):
    '''
    Carries out frequentist sensitivity analysis.
    
    Input:
        detector
        psr
        nfreq       (total number of instantiations to sample)
        injkind     (type of injection: 'GR' or 'G4v')
        pdif        (phase difference to form injection template: 'p', 'm' or '0')
        ninj        (number of injections to perform)
        rangeparam  (parameters to swipe over: in ['psi', 'iota', 'phi0'] or 'all')  [OPT]
        frange      (range of frequencies for heterodynes, default [1.0e-7, 1.0e-5]) [OPT]
        hinjrange   (range of injection magnitudes, default [1.0E-27, 1.0E-23])      [OPT]
        filesize    (number of instantiations per background file, default 100)      [OPT]
    '''
    
    def __init__(self, detector, psr, nfreq, injkind, pdif, ninj, rangeparam=[], frange=[1.0e-7, 1.0e-5], hinjrange=[1.0E-27, 1.0E-23], filesize=100):
        # system info
        self.detector = detector
        self.psr = psr
        
        self.log = logging.getLogger('Frequentist')
        self.log.debug('Initializing frequentist sensitivity analysis.')
                
        # data info
        self.freq = np.linspace(frange[0], frange[1], nfreq)

        self.background = Background(detector, psr, self.freq, filesize)
        self.background.get()
        
        self.log.debug('Obtaining time.')
        self.t = self.background.seed.finehet.index
        
        self.log.debug('Obtaining sigma.')
        sigma = Sigma(self.detector, self.psr, self.background.seed.finehet)
        self.sg = sigma.std
        
        self.log.debug('Preparing injection strengths.')
        inj = np.linspace(hinjrange[0], hinjrange[1], ninj)
        injLocations = [int(x) for x in np.linspace(0, nfreq, ninj, endpoint=False)]
        self.hinj = np.zeros(nfreq)
        self.hinj[injLocations] = inj 
               
        self.pdif = pdif
        self.injkind = injkind
        self.injection = templates.Signal(detector, psr, injkind, pdif, self.t)
        
        src = self.injection.response.src

        self.log.debug('Preparing parameter ranges.')
        if 'psi' in rangeparam or rangeparam=='all':
            self.pol_range = [
                            src.param['POL'] - src.param['POL error'],
                            src.param['POL'] + src.param['POL error']
                            ]
        else:
            self.pol_range = [src.param['POL'], src.param['POL']]
        
        if 'iota' in rangeparam or rangeparam=='all':   
            self.inc_range = [
                            src.param['INC'] - src.param['INC error'],
                            src.param['INC'] + src.param['INC error']
                            ]
        else:
            self.inc_range = [src.param['INC'], src.param['INC']]

        if 'phi0' in rangeparam or rangeparam=='all':                     
            self.phi0_range = [0., np.pi/2]
        else:
            self.phi0_range = [0., 0.]
        

    def analyze(self, methods):

        self.log.info('Analyzing %d files.' % self.background.nsets)
    
        self.log.debug('Producing search templates.')
        search = {m: templates.Signal(self.detector, self.psr, m, 0, self.t) for m in methods}

        self.log.debug('Setting up results')
        self.results = Results(self.detector, self.psr, methods=methods, hinj=self.hinj, kind=self.injkind, pdif=self.pdif)
            
        self.log.debug('Looping over files')
        for n in range(self.background.nsets):
            self.log.debug('File %i.' % n)
            try:
                back_file = pd.HDFStore(self.background.path + str(n), 'r')
                data = back_file[self.psr]
            finally:
                back_file.close()
                
            self.log.debug('Looping over instantiations.')
            for inst in data.columns:
                
                inst_number = int(n*self.background.filesize + inst)
                
                self.log.info('%i/%i ' % (inst_number, len(self.hinj)-1))
                
                self.log.debug('Selecting psi, iota, phi0.')
                psi  = random.uniform(self.pol_range[0], self.pol_range[1])
                iota = random.uniform(self.inc_range[0], self.inc_range[1])

                psi_inj  = random.uniform(self.pol_range[0], self.pol_range[1])
                iota_inj = random.uniform(self.inc_range[0], self.inc_range[1])
                phi0 = random.uniform(self.phi0_range[0], self.phi0_range[1])                    

                self.log.debug('Search with POL: %f, INC: %f' % (psi, iota))
                
                self.log.debug('Loop over search methods.')
                # note: important that this follows inst loop to get same psi and iota
                for m in methods:
                    
                    d = data[inst]
                    
                    # inject if necessary
                    h = self.hinj[inst_number]
                    if h != 0:
                        self.log.debug(self.injection.kind + str(self.injection.pdif))
                        self.log.debug('I! %(psi_inj)f %(iota_inj)f %(phi0)f' % locals())
                        d += h * self.injection.simulate(psi_inj, iota_inj, phase=phi0)
                    
                    self.log.debug('Get design matrix.')
                    designMatrix = search[m].design_matrix(psi, iota)
                    
                    A = designMatrix.div(self.sg, axis=0)

                    b = d / self.sg
                    
                    self.log.debug('SVD decomposition.')
                    svd = np.linalg.svd(A, full_matrices=False)
                    
                    U = pd.DataFrame(svd[0], columns=A.columns, index=A.index)
                    W = pd.DataFrame(np.diag(1./svd[1]), index=A.columns, columns=A.columns)
                    V = pd.DataFrame(svd[2], index=A.columns, columns=A.columns)
                    
                    cov = V.T.dot(W**2).dot(V)  # covariance matrix
                    
                    VtW = V.T.dot(W)
                    # need to make U complex before dotting with b
                    Utb = (U + 0j).mul(b, axis=0).sum(axis=0)
                    a = VtW.dot(Utb.T)          # results

                    self.log.debug('Average h0')
                    self.results.h[m][inst_number] = (abs(a).sum()) / len(a)
                    
                    self.log.debug('Significance')
                    self.results.s[m][inst_number] = abs(np.dot(a.conj(), np.linalg.solve(cov, a)))

        ## Save
        self.results.save()
            

## MANY PULSAR ANALYSIS
class ManyPulsars(object):
    '''
    Analyzes sets of multiple pulsars.
    '''
    
    def __init__(self, detector, methods=['GR', 'G4v', 'AP']):
        self.detector = detector
        
        # look for the data in the following directory
        self.dir = paths.originalData + '/data' + detector
        
        self.methods = methods
        
        # get names of all PSRs in directory
        self.allpsrs = listpsrs(self.detector, self.dir)
        
        # book-keeping
        self.hasresults = False
        self.failed = []
        
        self.log = logging.getLogger('Many PSRs')
        self.log.info('Analyzing '+detector+' data with '+str(methods))
    
    def census(self, ratio=[0,1]):
        # Splits list into ratio[1] parts and picks part number ratio[0].
        
        self.log.debug('Performing census.')
        
        names = self.allpsrs
        
        self.log.info('There are %d PSRs on file.' % len(names))
        
        # select subset according to range
        set_choice = ratio[0]
        set_options = ratio[1]
        
        self.log.info('Choosing subset #%d out of %d subsets.' % (set_choice+1, set_options))
       
        setlength = len(names)/int(set_options)
        
        if set_options > len(names): self.log.error('More sets than pulsars!')
        if set_choice < 0 or type(set_choice)!=int: self.log.error('Set number must be a positive integer!')
        
        i0 = set_choice * setlength
        i1 = (set_choice + 1) * setlength
        
        if set_choice < set_options-1: 
            self.psrlist = names[i0: i1]
            
        elif set_choice == set_options-1:       
            self.psrlist = names[i0:]
        else:
            self.psrlist = names
            
    
    def analyze(self, injkind, ratio=[0,1], extra_name=''):
        
        self.log.debug('Beginning MP analysis. Injecting ' + str(injkind))
        
        # get PSR subset
        self.census(ratio)
        
        # setup results
        for m in self.methods:
            setattr(self,'stats'+m, pd.DataFrame(columns=self.psrlist, index=sd.statkinds))
        
        # loop over PSRs
        count = 0
        for psr in self.psrlist:
            count += 1
            self.log.info('Analyzing '+psr+' ('+str(1)+'/'+str(len(self.allpsrs))+')')
            
            try:
                ij = Frequentist(self.detector, psr, 2000, injkind, 'p', 100, rangeparam='all', filesize=200)
            
                ij.analyze(self.methods)
            
                self.log.debug('Recording results.')
            
                for m in self.methods:
                    name = 'stats' + m + '[' + psr + ']'
                    setattr(self,name, ij.results.stats[m])
            except:
                # print error message
                e = sys.exc_info()[0]
                self.log.error("<p>Error: %s</p>" % e)

                self.log.error(psr + ' search failed.')
                self.failed += [psr]
        
        # save stats
        self.save(extra_name=extra_name)
        

    def save(self, extra_name=''):
        self.log.info('Saving results.')
        now = str(datetime.datetime.now())
        path = paths.results + 'manypsr_' + self.detector + '_' + now + '_' + extra_name
        try:
            f = pd.HDFStore(path, 'w')
            
            for m in self.methods:
                f[m] = getattr(self, 'stats' + m )
        except IOError:
            # print error message
            e = sys.exc_info()[0]
            self.log.error("<p>Error: %s</p>" % e)

            self.log.error('Error: cannot save stats, something wrong with directory.\n %s' % path)
        else:
            f.close()
    
            
class MP10gr(ManyPulsars):
    def __init__(self, n, methods=['GR', 'G4v', 'Sid']):
        super(MP10gr, self).__init__('H1', methods=methods)
        self.analyze('GR', [n, 10], extra_name=str(n)+'-9')
    
        
class MP10g4v(ManyPulsars):
    def __init__(self, n, methods=['GR', 'G4v', 'Sid']):
        super(MP10g4v, self).__init__('H1', methods=methods)
        self.analyze('G4v', [n, 10], extra_name=str(n)+'-9')
  
        
class MPstats(object):
    '''
    Loads pulsar search data stored on indicated path and produces histograms of their
    statistics.
    Input:
        detector
        load (indicates whether to load data on start. Default: True)                [OPT]
        pth (path where search data is located. Default is local 'files/analysis/results/a
             tlas/' + detector + '/')                                                [OPT]
    Subfunctions:
        load (loads data)
        hist (histograms)
    '''
    def __init__(self, detector, load=True, pth='files/analysis/results/atlas/' + detector + '/'):
        self.detector = detector
        
        self.path = pth
        self.psrs = os.listdir(self.path)
        self.psrs.remove('.DS_Store')
        
        self.failed = []
        
        self.statdict = {k : k.replace(' ', '_') for k in sd.statkinds}
        
        
    def load(self, injkind, pdif='p', methods=['GR', 'G4v', 'Sid']):
        
        self.injkind = injkind
        self.pdif = pdif
        
        for k in sd.statkinds:
            setattr(self, self.statdict[k], pd.DataFrame(columns=methods, index=self.psrs))
        
        for psr in self.psrs:
            path = self.path + '/' + psr + '/'
            name = psr + '_' + self.detector + '_' + injkind + '_' + pdif
            
            try:
                f = pd.HDFStore(path + name, 'r')
                
                stats = f['stats']

                for k in sd.statkinds:
                    stat_local = getattr(self, self.statdict[k])
                    stat_local.ix[psr] = stats.ix[k]
                
            except IOError:
                print 'No data for ' + psr
                self.failed += [psr]
            else:
                f.close()

    def hist(self, nbins=25, methods=['GR', 'G4v', 'Sid'], kinds=False, log=False, together=False):
        '''
        Plots statistical summary data.
        Input:
            nbins (25)
            methods (['GR', 'G4v', 'Sid'])
            kinds (types of stats to be plotted. Options are (set): 'min inj det',
                  'lin s slope', 'lin s noise', 'lin s inter', 'h rec noise', 'h rec
                  slope', 'h rec inter'. If 'False' (default), takes full set.)
            log (False)
            together (if True, histograms are plotted on single figure. Def: False)
        '''
        
        if not kinds: kinds=self.statdict
        
        for k in kinds:
            # get values to histogram
            stat = getattr(self, self.statdict[k])
            
            stat_label = sd.statlabels[k]
            
            # setup save            
            path = paths.plots + self.detector + '/manypulsars/' + k + '/'

            if log:
                stat_label += ' (log scale)'
                logname = '_log'
            else:
                logname = ''
           
            try:
                os.makedirs(path)
            except:
                pass
            
            plt.figure()
            
            for m in methods:
                
                # histogram
                ax = stat[m].hist(color=sd.pltcolor[m], bins=nbins, label=m, histtype='step')
                
                if log: ax.set_xscale('log')
                
                # format
                plt.title(self.injkind + ' injection on ' + self.detector + ' data ' + k.replace('_', ' ') + ' for ' + str(len(self.psrs)-len(self.failed)) + ' PSRs')
                
                plt.xlabel(stat_label)
                plt.ylabel('Count')
                plt.legend(numpoints=1)
                
                # save
                if not together:
                    plt.savefig(path + k + '_inj' + self.injkind + self.pdif + '_srch' + m + '_' + self.detector + logname, bbox_inches='tight')
                    plt.close()
                    
            if together:
                plt.savefig(path + k + '_inj' + self.injkind + self.pdif + '_srchAll_' + self.detector + logname, bbox_inches='tight')
                plt.close()
                
            print 'Saved in ' + path
                
                
## SPECIAL CASES
class SinglePulsar(object):
    def __init__(self, detector, psr, pd=['p'], methods=['GR', 'G4v']):
        self.detector = detector
        self.psr = psr
        
        self.injection_kinds = ['GR', 'G4v']

        self.search_methods = methods
        
        self.pd = pd

        self.plots = ['hinjrec', 'hinjs', 'hinjlins']


    def scan(self, range='', hinjrange=[1.0E-27, 1.0E-23], extra_name='2'):

        for kind in self.injection_kinds:

            for p in self.pd:
    
                ij = InjSearch(self.detector, self.psr, 2000, kind, p, 100, hinjrange=hinjrange, rangeparam=[range])
        
                ij.analyze(self.search_methods)
        
                for pl in self.plots:
                    ij.results.plots(pl, extra_name=extra_name + '_range'+range)
                    
                ij.results.save(extra_name=extra_name + '_range'+range)
                
                
class Crab(SinglePulsar):
    def __init__(self, paramrange='all', methods=['GR', 'G4v'], extra_name='S6'):
        super(Crab, self).__init__('H1', 'J0534+2200', methods=methods)
        self.scan(range=paramrange, hinjrange=[1.0E-27, 1.0E-24])

      
class Vela(SinglePulsar):
    def __init__(self, paramrange='all', methods=['GR', 'G4v'], extra_name='2'):
        super(Vela, self).__init__('V1', 'J0835-4510', methods=methods)
        self.scan(range=paramrange, hinjrange=[1.0E-27, 1.0E-23])