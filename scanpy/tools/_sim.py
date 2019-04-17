# Author: Alex Wolf (http://falexwolf.de)
"""Simulate Data

Simulate stochastic dynamic systems to model gene expression dynamics and
cause-effect data.

TODO
----
Beta Version. The code will be reorganized soon.
"""

import os
import itertools
import collections
from collections import OrderedDict as odict

import numpy as np
import scipy as sp
from anndata import AnnData

from .. import utils
from .. import readwrite
from .._settings import settings
from .. import logging as logg


def sim(
    model,
    params_file=True,
    tmax=None,
    branching=None,
    nrRealizations=None,
    noiseObs=None,
    noiseDyn=None,
    step=None,
    seed=None,
    writedir=None,
) -> AnnData:
    """Simulate dynamic gene expression data [Wittmann09]_ [Wolf18]_.

    Sample from a stochastic differential equation model built from
    literature-curated boolean gene regulatory networks, as suggested by
    [Wittmann09]_. The Scanpy implementation is due to [Wolf18]_.

    Parameters
    ----------
    model : {'krumsiek11', 'toggleswitch'}
        Model file in 'sim_models' directory.
    params_file : `bool`, (default: `True`)
        Read default params from file.
    tmax : `int`, optional (default: `None`)
        Number of time steps per realization of time series.
    branching : `bool`, optional (default: `None`)
        Only write realizations that contain new branches.
    nrRealizations : int, optional (default: `None`)
        Number of realizations.
    noiseObs : float, optional (default: `None`)
        Observatory/Measurement noise.
    noiseDyn : float, optional (default: `None`)
        Dynamic noise.
    step : int, optional (default: `None`)
        Interval for saving state of system.
    seed : int, optional (default: `None`)
        Seed for generation of random numbers.
    writedir: str, optional (default: `None`)
        Path to directory for writing output files.

    Returns
    -------
    Annotated data matrix.

    Examples
    --------
    See this `use case <https://github.com/theislab/scanpy_usage/tree/master/170430_krumsiek11>`__
    """
    params = locals()
    if params_file:
        model_key = os.path.basename(model).replace('.txt', '')
        from .. import sim_models
        pfile_sim = (os.path.dirname(sim_models.__file__)
                     + '/' + model_key + '_params.txt')
        default_params = readwrite.read_params(pfile_sim)
        params = utils.update_params(default_params, params)
    adata = sample_dynamic_data(params)
    adata.uns['iroot'] = 0
    return adata


def add_args(p):
    """
    Update parser with tool specific arguments.

    This overwrites was is done in utils.uns_args.
    """
    # dictionary for adding arguments
    dadd_args = {
        '--opfile': {
            'default': '',
            'metavar': 'f',
            'type': str,
            'help': 'Specify a parameter file '
                    '(default: "sim/${exkey}_params.txt")'}}
    p = utils.add_args(p, dadd_args)
    return p


def sample_dynamic_data(params):
    """
    Helper function.
    """
    model_key = os.path.basename(params['model']).replace('.txt', '')
    if 'writedir' not in params or params['writedir'] is None:
        params['writedir'] = settings.writedir + model_key + '_sim'
    if not os.path.exists(params['writedir']): os.makedirs(params['writedir'])
    readwrite.write_params(params['writedir'] + '/params.txt', params)
    # init variables
    writedir = params['writedir']
    tmax = params['tmax']
    branching = params['branching']
    noiseObs = params['noiseObs']
    noiseDyn = params['noiseDyn']
    nrRealizations = params['nrRealizations']
    step = params['step']  # step size for saving the figure

    nrSamples = 1  # how many files?
    maxRestarts = 1000
    maxNrSamples = 1

    # simple vector auto regressive process or
    # hill kinetics process simulation
    if 'krumsiek11' not in model_key:
        # create instance, set seed
        grnsim = GRNsim(model=model_key, params=params)
        nrOffEdges_list = np.zeros(nrSamples)
        for sample in range(nrSamples):
            # random topology / for a given edge density
            if 'hill' not in model_key:
                Coupl = np.array(grnsim.Coupl)
                for sampleCoupl in range(10):
                    nrOffEdges = 0
                    for gp in range(grnsim.dim):
                        for g in range(grnsim.dim):
                            # only consider off-diagonal edges
                            if g != gp:
                                Coupl[gp, g] = 0.7 if np.random.rand() < 0.4 else 0
                                nrOffEdges += 1 if Coupl[gp, g] > 0 else 0
                            else:
                                Coupl[gp, g] = 0.7
                    # check that the coupling matrix does not have eigenvalues
                    # greater than 1, which would lead to an exploding var process
                    if max(sp.linalg.eig(Coupl)[0]) < 1:
                        break
                nrOffEdges_list[sample] = nrOffEdges
                grnsim.set_coupl(Coupl)
            # init type
            real = 0
            X0 = np.random.rand(grnsim.dim)
            Xsamples = []
            for restart in range(nrRealizations+maxRestarts):
                # slightly break symmetry in initial conditions
                if 'toggleswitch' in model_key:
                    X0 = (np.array([0.8 for i in range(grnsim.dim)])
                          + 0.01*np.random.randn(grnsim.dim))
                X = grnsim.sim_model(tmax=tmax, X0=X0,
                                     noiseDyn=noiseDyn)
                # check branching
                check = True
                if branching:
                    check, Xsamples = _check_branching(X, Xsamples, restart)
                if check:
                    real += 1
                    grnsim.write_data(X[::step], dir=writedir,
                                      noiseObs=noiseObs,
                                      append=(False if restart==0 else True),
                                      branching=branching,
                                      nrRealizations=nrRealizations)
                # append some zeros
                if 'zeros' in writedir and real == 2:
                    grnsim.write_data(noiseDyn*np.random.randn(500,3), dir=writedir,
                                      noiseObs=noiseObs,
                                      append=(False if restart==0 else True),
                                      branching=branching,
                                      nrRealizations=nrRealizations)
                if real >= nrRealizations:
                    break
        if False:
            logg.info('mean nr of offdiagonal edges',nrOffEdges_list.mean(),
                      'compared to total nr',grnsim.dim*(grnsim.dim-1)/2.)

    # more complex models
    else:
        initType = 'random'

        dim = 11
        step = 5

        grnsim = GRNsim(dim=dim,initType=initType,model=model_key,params=params)
        curr_nrSamples = 0
        Xsamples = []
        for sample in range(maxNrSamples):
            # choose initial conditions such that branchings result
            if initType == 'branch':
                X0mean = grnsim.branch_init_model1(tmax)
                if X0mean is None:
                    grnsim.set_coupl()
                    continue
            real = 0
            for restart in range(nrRealizations+maxRestarts):
                if initType == 'branch':
                    # vary initial conditions around mean
                    X0 = X0mean + (0.05*np.random.rand(dim) - 0.025*np.ones(dim))
                else:
                    # generate random initial conditions within [0.3,0.7]
                    X0 = 0.4*np.random.rand(dim)+0.3
                if model_key in [5,6]:
                    X0 = np.array([0.3,0.3,0,0,0,0])
                if model_key in [7,8,9,10]:
                    X0 = 0.6*np.random.rand(dim)+0.2
                    X0[2:] = np.zeros(4)
                if 'krumsiek11' in model_key:
                    X0 = np.zeros(dim)
                    X0[grnsim.varNames['Gata2']] = 0.8
                    X0[grnsim.varNames['Pu.1']] = 0.8
                    X0[grnsim.varNames['Cebpa']] = 0.8
                    X0 += 0.001*np.random.randn(dim)
                    if False:
                        switch_gene = restart - (nrRealizations - dim)
                        if switch_gene >= dim:
                            break
                        X0[switch_gene] = 0 if X0[switch_gene] > 0.1 else 0.8
                X = grnsim.sim_model(tmax,X0=X0,
                                     noiseDyn=noiseDyn,
                                     restart=restart)
                # check branching
                check = True
                if branching:
                    check, Xsamples = _check_branching(X,Xsamples,restart)
                if check:
                    real += 1
                    grnsim.write_data(X[::step],dir=writedir,noiseObs=noiseObs,
                                     append=(False if restart==0 else True),
                                     branching=branching,
                                     nrRealizations=nrRealizations)
                if real >= nrRealizations:
                    break
    import glob
    # load the last simulation file
    filename = glob.glob(writedir + '/sim*.txt')[-1]
    logg.info('reading simulation results', filename)
    adata = readwrite._read(filename, first_column_names=True,
                            suppress_cache_warning=True)
    adata.uns['tmax_write'] = tmax/step
    return adata


def write_data(X, dir='sim/test', append=False, header='',
               varNames={}, Adj=np.array([]), Coupl=np.array([]),
               boolRules={}, model='', modelType='', invTimeStep=1):
    """ Write simulated data.

        Accounts for saving at the same time an ID
        and a model file.
    """
    # check if output directory exists
    if not os.path.exists(dir):
        os.makedirs(dir)
    # update file with sample ids
    filename = dir+'/id.txt'
    if os.path.exists(filename):
        f = open(filename,'r')
        id = int(f.read()) + (0 if append else 1)
        f.close()
    else:
        id = 0
    f = open(filename,'w')
    id = '{:0>6}'.format(id)
    f.write(str(id))
    f.close()
    # dimension
    dim = X.shape[1]
    # write files with adjacancy and coupling matrices
    if not append:
        if False:
            if Adj.size > 0:
                # due to 'update formulation' of model, there
                # is always a diagonal dependence
                Adj = np.copy(Adj)
                if 'hill' in model:
                    for i in range(Adj.shape[0]):
                        Adj[i,i] = 1
                np.savetxt(dir+'/adj_'+id+'.txt',Adj,
                           header=header,
                           fmt='%d')
            if Coupl.size > 0:
                np.savetxt(dir+'/coupl_'+id+'.txt',Coupl,
                           header=header,
                           fmt='%10.6f')
        # write model file
        if varNames and Coupl.size > 0:
            f = open(dir+'/model_'+id+'.txt','w')
            f.write('# For each "variable = ", there must be a right hand side: \n')
            f.write('# either an empty string or a python-style logical expression \n')
            f.write('# involving variable names, "or", "and", "(", ")". \n')
            f.write('# The order of equations matters! \n')
            f.write('# \n')
            f.write('# modelType = ' + modelType + '\n')
            f.write('# invTimeStep = '+ str(invTimeStep) + '\n')
            f.write('# \n')
            f.write('# boolean update rules: \n')
            for rule in boolRules.items():
                f.write(rule[0] + ' = ' + rule[1] + '\n')
            # write coupling via names
            f.write('# coupling list: \n')
            names = list(varNames.keys())
            for gp in range(dim):
                for g in range(dim):
                    if np.abs(Coupl[gp,g]) > 1e-10:
                        f.write('{:10} '.format(names[gp])
                                + '{:10} '.format(names[g])
                                + '{:10.3} '.format(Coupl[gp,g]) + '\n')
            f.close()
    # write simulated data
    # the binary mode option in the following line is a fix for python 3
    # variable names
    if varNames:
        header += '{:>2} '.format('it')
        for v in varNames.keys():
            header += '{:>7} '.format(v)
    f = open(dir+'/sim_'+id+'.txt','ab' if append else 'wb')
    np.savetxt(f,np.c_[np.arange(0,X.shape[0]),X],header=('' if append else header),
               fmt=['%4.f']+['%7.4f' for i in range(dim)])
    f.close()

class GRNsim:
    """
    Simlulation of stochastic dynamic systems.

    Main application: simulation of gene expression dynamics.

    Also standard models are implemented.
    """

    availModels = collections.OrderedDict([
        ('krumsiek11',
             ('myeloid progenitor network, Krumsiek et al., PLOS One 6, e22649, \n      '
              'equations from Table 1 on page 3, doi:10.1371/journal.pone.0022649 \n')),
        ('var','vector autoregressive process \n'),
        ('hill','process with hill kinetics \n')])

    writeOutputOnce = True

    def __init__(self,dim=3,model='ex0',modelType='var',
                 initType='random',show=False,verbosity=0,
                 Coupl=None,params={}):
        """
            model : either string for predefined model, or directory with
                    a model file and a coupl matrix file
        """
        self.dim = dim if Coupl is None else Coupl.shape[0] # number of nodes / dimension of system
        self.maxnpar = 1 # maximal number of parents
        self.p_indep = 0.4 # fraction of independent genes
        self.model = model
        self.modelType = modelType
        self.initType = initType # string characterizing a specific initial
        self.show = show
        self.verbosity = verbosity
        # checks
        if initType not in ['branch','random']:
            raise RuntimeError('initType must be either: branch, random')
        read = False
        if model not in self.availModels.keys():
            message = 'model not among predefined models \n'
        # read from file
        from .. import sim_models
        model = os.path.dirname(sim_models.__file__) + '/' + model + '.txt'
        if not os.path.exists(model):
            if not os.path.exists(model):
                message = '  cannot read model from file ' + model
                message += '\n as the directory does not exist'
                raise RuntimeError(message)
        self.model = model
        # set the coupling matrix, and with that the adjacency matrix
        self.set_coupl(Coupl=Coupl)
        # seed
        np.random.seed(params['seed'])
        # header
        self.header  = 'model = ' + self.model+ ' \n'
        # params
        self.params = params

    def sim_model(self,tmax,X0,noiseDyn=0,restart=0):
        """ Simulate the model.
        """
        self.noiseDyn = noiseDyn
        #
        X = np.zeros((tmax,self.dim))
        X[0] = X0 + noiseDyn*np.random.randn(self.dim)
        # run simulation
        for t in range(1,tmax):
            if self.modelType == 'hill':
                Xdiff = self.Xdiff_hill(X[t-1])
            elif self.modelType == 'var':
                Xdiff = self.Xdiff_var(X[t-1])
            #
            X[t] = X[t-1] + Xdiff
            # add dynamic noise
            X[t] += noiseDyn*np.random.randn(self.dim)
        return X

    def Xdiff_hill(self,Xt):
        """ Build Xdiff from coefficients of boolean network,
            that is, using self.boolCoeff. The employed functions
            are Hill type activation and deactivation functions.

            See Wittmann et al., BMC Syst. Biol. 3, 98 (2009),
            doi:10.1186/1752-0509-3-98 for more details.
        """
        verbosity = self.verbosity>0 and self.writeOutputOnce
        self.writeOutputOnce = False
        Xdiff = np.zeros(self.dim)
        for ichild,child in enumerate(self.pas.keys()):
            # check whether list of parents is non-empty,
            # otherwise continue
            if self.pas[child]:
                Xdiff_syn = 0 # synthesize term
                if verbosity > 0:
                    Xdiff_syn_str = ''
            else:
                continue
            # loop over all tuples for which the boolean update
            # rule returns true, these are stored in self.boolCoeff
            for ituple,tuple in enumerate(self.boolCoeff[child]):
                Xdiff_syn_tuple = 1
                Xdiff_syn_tuple_str = ''
                for iv,v in enumerate(tuple):
                    iparent = self.varNames[self.pas[child][iv]]
                    x = Xt[iparent]
                    threshold = 0.1/np.abs(self.Coupl[ichild,iparent])
                    Xdiff_syn_tuple *= self.hill_a(x,threshold) if v else self.hill_i(x,threshold)
                    if verbosity > 0:
                        Xdiff_syn_tuple_str += (('a' if v else 'i')
                                                +'('+self.pas[child][iv]+','+'{:.2}'.format(threshold)+')')
                Xdiff_syn += Xdiff_syn_tuple
                if verbosity > 0:
                    Xdiff_syn_str += ('+' if ituple != 0 else '') + Xdiff_syn_tuple_str
            # multiply with degradation term
            Xdiff[ichild] = self.invTimeStep*(Xdiff_syn - Xt[ichild])
            if verbosity > 0:
                Xdiff_str = (child+'_{+1}-' + child + ' = ' + str(self.invTimeStep)
                             + '*('+Xdiff_syn_str+'-'+child+')' )
                settings.m(0,Xdiff_str)
        return Xdiff

    def Xdiff_var(self,Xt,verbosity=0):
         """
         """
         # subtract the current state
         Xdiff = -Xt
         # add the information from the past
         Xdiff += np.dot(self.Coupl,Xt)
         return Xdiff

    def hill_a(self,x,threshold=0.1,power=2):
        """ Activating hill function. """
        x_pow = np.power(x,power)
        threshold_pow = np.power(threshold,power)
        return x_pow / (x_pow + threshold_pow)

    def hill_i(self,x,threshold=0.1,power=2):
        """ Inhibiting hill function.

            Is equivalent to 1-hill_a(self,x,power,threshold).
        """
        x_pow = np.power(x,power)
        threshold_pow = np.power(threshold,power)
        return threshold_pow / (x_pow + threshold_pow)

    def nhill_a(self,x,threshold=0.1,power=2,ichild=2):
        """ Normalized activating hill function. """
        x_pow = np.power(x,power)
        threshold_pow = np.power(threshold,power)
        return x_pow / (x_pow + threshold_pow) * (1 + threshold_pow)

    def nhill_i(self,x,threshold=0.1,power=2):
        """ Normalized inhibiting hill function.

            Is equivalent to 1-nhill_a(self,x,power,threshold).
        """
        x_pow = np.power(x,power)
        threshold_pow = np.power(threshold,power)
        return threshold_pow / (x_pow + threshold_pow) * (1 - x_pow)

    def read_model(self):
        """ Read the model and the couplings from the model file.
        """
        if self.verbosity > 0:
            settings.m(0,'reading model',self.model)
        # read model
        boolRules = []
        for line in open(self.model):
            if line.startswith('#') and 'modelType =' in line:
                keyval = line
                if '|' in line:
                    keyval, type = line.split('|')[:2]
                self.modelType = keyval.split('=')[1].strip()
            if line.startswith('#') and 'invTimeStep =' in line:
                keyval = line
                if '|' in line:
                    keyval, type = line.split('|')[:2]
                self.invTimeStep = float(keyval.split('=')[1].strip())
            if not line.startswith('#'):
                boolRules.append([s.strip() for s in line.split('=')])
            if line.startswith('# coupling list:'):
                break
        self.dim = len(boolRules)
        self.boolRules = collections.OrderedDict(boolRules)
        self.varNames = collections.OrderedDict([(s, i)
            for i, s in enumerate(self.boolRules.keys())])
        names = self.varNames
        # read couplings via names
        self.Coupl = np.zeros((self.dim, self.dim))
        boolContinue = True
        for line in open(self.model):  # open(self.model.replace('/model','/couplList')):
            if line.startswith('# coupling list:'):
                boolContinue = False
            if boolContinue:
                continue
            if not line.startswith('#'):
                gps, gs, val = line.strip().split()
                self.Coupl[int(names[gps]), int(names[gs])] = float(val)
        # adjancecy matrices
        self.Adj_signed = np.sign(self.Coupl)
        self.Adj = np.abs(np.array(self.Adj_signed))
        # build bool coefficients (necessary for odefy type
        # version of the discrete model)
        self.build_boolCoeff()

    def set_coupl(self, Coupl=None):
        """ Construct the coupling matrix (and adjacancy matrix) from predefined models
            or via sampling.
        """
        self.varNames = collections.OrderedDict([(str(i), i) for i in range(self.dim)])
        if (self.model not in self.availModels.keys()
            and Coupl is None):
            self.read_model()
        elif 'var' in self.model:
            # vector auto regressive process
            self.Coupl = Coupl
            self.boolRules = collections.OrderedDict(
                              [(s, '') for s in self.varNames.keys()])
            names = list(self.varNames.keys())
            for gp in range(self.dim):
                pas = []
                for g in range(self.dim):
                    if np.abs(self.Coupl[gp,g] > 1e-10):
                        pas.append(names[g])
                self.boolRules[
                    names[gp]] = ''.join(pas[:1]
                                         + [' or ' + pa for pa in pas[1:]])
                self.Adj_signed = np.sign(Coupl)
        elif self.model in ['6','7','8','9','10']:
            self.Adj_signed = np.zeros((self.dim,self.dim))
            n_sinknodes = 2
#             sinknodes = np.random.choice(np.arange(0,self.dim),
#                                              size=n_sinknodes,replace=False)
            sinknodes = np.array([0,1])
            # assume sinknodes have feeback
            self.Adj_signed[sinknodes,sinknodes] = np.ones(n_sinknodes)
#             # allow negative feedback
#             if self.model == 10:
#                 plus_minus = (np.random.randint(0,2,n_sinknodes) - 0.5)*2
#                 self.Adj_signed[sinknodes,sinknodes] = plus_minus
            leafnodes = np.array(sinknodes)
            availnodes = np.array([i for i in range(self.dim) if i not in sinknodes])
#             settings.m(0,leafnodes,availnodes)
            while len(availnodes) != 0:
                # parent
                parent_idx = np.random.choice(np.arange(0,len(leafnodes)),
                                                    size=1,replace=False)
                parent = leafnodes[parent_idx]
                # children
                children_ids = np.random.choice(np.arange(0,len(availnodes)),
                                                       size=2,replace=False)
                children = availnodes[children_ids]
                settings.m(0,parent,children)
                self.Adj_signed[children,parent] = np.ones(2)
                if self.model == 8:
                    self.Adj_signed[children[0],children[1]] = -1
                if self.model in [9,10]:
                    self.Adj_signed[children[0],children[1]] = -1
                    self.Adj_signed[children[1],children[0]] = -1
                # update leafnodes
                leafnodes = np.delete(leafnodes,parent_idx)
                leafnodes = np.append(leafnodes,children)
                # update availnodes
                availnodes = np.delete(availnodes,children_ids)
#                 settings.m(0,availnodes)
#                 settings.m(0,leafnodes)
#                 settings.m(0,self.Adj)
#                 settings.m(0,'-')
        else:
            self.Adj = np.zeros((self.dim,self.dim))
            for i in range(self.dim):
                indep = np.random.binomial(1,self.p_indep)
                if indep == 0:
                    # this number includes parents (other variables)
                    # and the variable itself, therefore its
                    # self.maxnpar+2 in the following line
                    nr = np.random.randint(1,self.maxnpar+2)
                    j_par = np.random.choice(np.arange(0,self.dim),
                                             size=nr,replace=False)
                    self.Adj[i,j_par] = 1
                else:
                    self.Adj[i,i] = 1
        #
        self.Adj = np.abs(np.array(self.Adj_signed))
        #settings.m(0,self.Adj)

    def set_coupl_old(self):
        """ Using the adjacency matrix, sample a coupling matrix.
        """
        if self.model == 'krumsiek11' or self.model == 'var':
            # we already built the coupling matrix in set_coupl20()
            return
        self.Coupl = np.zeros((self.dim,self.dim))
        for i in range(self.Adj.shape[0]):
            for j,a in enumerate(self.Adj[i]):
                # if there is a 1 in Adj, specify co and antiregulation
                # and strength of regulation
                if a != 0:
                    co_anti = np.random.randint(2)
                    # set a lower bound for the coupling parameters
                    # they ought not to be smaller than 0.1
                    # and not be larger than 0.4
                    self.Coupl[i,j] = 0.0*np.random.rand() + 0.1
                    # set sign for coupling
                    if co_anti == 1:
                        self.Coupl[i,j] *= -1
        # enforce certain requirements on models
        if self.model == 1:
            self.coupl_model1()
        elif self.model == 5:
            self.coupl_model5()
        elif self.model in [6,7]:
            self.coupl_model6()
        elif self.model in [8,9,10]:
            self.coupl_model8()
        # output
        if self.verbosity > 1:
           settings.m(0,self.Coupl)

    def coupl_model1(self):
        """ In model 1, we want enforce the following signs
            on the couplings. Model 2 has the same couplings
            but arbitrary signs.
        """
        self.Coupl[0,0] = np.abs(self.Coupl[0,0])
        self.Coupl[0,1] = -np.abs(self.Coupl[0,1])
        self.Coupl[1,1] = np.abs(self.Coupl[1,1])

    def coupl_model5(self):
        """ Toggle switch.
        """
        self.Coupl = -0.2*self.Adj
        self.Coupl[2,0] *= -1
        self.Coupl[3,0] *= -1
        self.Coupl[4,1] *= -1
        self.Coupl[5,1] *= -1

    def coupl_model6(self):
        """ Variant of toggle switch.
        """
        self.Coupl = 0.5*self.Adj_signed

    def coupl_model8(self):
        """ Variant of toggle switch.
        """
        self.Coupl = 0.5*self.Adj_signed
        # reduce the value of the coupling of the repressing genes
        # otherwise completely unstable solutions are obtained
        for x in np.nditer(self.Coupl,op_flags=['readwrite']):
            if x < -1e-6:
                x[...] = -0.2

    def coupl_model_krumsiek11(self):
        """ Variant of toggle switch.
        """
        self.Coupl = self.Adj_signed

    def sim_model_back_help(self,Xt,Xt1):
        """ Yields zero when solved for X_t
            given X_{t+1}.
        """
        return - Xt1 + Xt + self.Xdiff(Xt)

    def sim_model_backwards(self,tmax,X0):
        """ Simulate the model backwards in time.
        """
        X = np.zeros((tmax,self.dim))
        X[tmax-1] = X0
        for t in range(tmax-2,-1,-1):
            sol = sp.optimize.root(self.sim_model_back_help,
                                   X[t+1],
                                   args=(X[t+1]),method='hybr')
            X[t] = sol.x
        return X

    def branch_init_model1(self,tmax=100):
        # check whether we can define trajectories
        Xfix = np.array([self.Coupl[0,1]/self.Coupl[0,0],1])
        if Xfix[0] > 0.97 or Xfix[0] < 0.03:
            settings.m(0,'... either no fixed point in [0,1]^2! \n' +
                  '    or fixed point is too close to bounds' )
            return None
        #
        XbackUp = grnsim.sim_model_backwards(tmax=tmax/3,X0=Xfix+np.array([0.02,-0.02]))
        XbackDo = grnsim.sim_model_backwards(tmax=tmax/3,X0=Xfix+np.array([-0.02,-0.02]))
        #
        Xup = grnsim.sim_model(tmax=tmax,X0=XbackUp[0])
        Xdo = grnsim.sim_model(tmax=tmax,X0=XbackDo[0])
        # compute mean
        X0mean = 0.5*(Xup[0] + Xdo[0])
        #
        if np.min(X0mean) < 0.025 or np.max(X0mean) > 0.975:
            settings.m(0,'... initial point is too close to bounds' )
            return None
        #
        if self.show and self.verbosity > 1:
            pl.figure()
            pl.plot(XbackUp[:,0],'.b',XbackUp[:,1],'.g')
            pl.plot(XbackDo[:,0],'.b',XbackDo[:,1],'.g')
            pl.plot(Xup[:,0],'b',Xup[:,1],'g')
            pl.plot(Xdo[:,0],'b',Xdo[:,1],'g')
        return X0mean

    def parents_from_boolRule(self,rule):
        """ Determine parents based on boolean updaterule.

            Returns list of parents.
        """
        rule_pa = rule.replace('(','').replace(')','').replace('or','').replace('and','').replace('not','')
        rule_pa = rule_pa.split()
        # if there are no parents, continue
        if not rule_pa:
            return []
        # check whether these are meaningful parents
        pa_old = []
        pa_delete = []
        for pa in rule_pa:
            if pa not in self.varNames.keys():
                settings.m(0,'list of available variables:')
                settings.m(0,list(self.varNames.keys()))
                message = ('processing of rule "' + rule
                             + ' yields an invalid parent: ' + pa
                             + ' | check whether the syntax is correct: \n'
                             + 'only python expressions "(",")","or","and","not" '
                             + 'are allowed, variable names and expressions have to be separated '
                             + 'by white spaces')
                raise ValueError(message)
            if pa in pa_old:
                pa_delete.append(pa)
        for pa in pa_delete:
            rule_pa.remove(pa)
        return rule_pa

    def build_boolCoeff(self):
        ''' Compute coefficients for tuple space.
        '''
        # coefficients for hill functions from boolean update rules
        self.boolCoeff = collections.OrderedDict([(s,[]) for s in self.varNames.keys()])
        # parents
        self.pas = collections.OrderedDict([(s,[]) for s in self.varNames.keys()])
        #
        for key in self.boolRules.keys():
            rule = self.boolRules[key]
            self.pas[key] = self.parents_from_boolRule(rule)
            pasIndices = [self.varNames[pa] for pa in self.pas[key]]
            # check whether there are coupling matrix entries for each parent
            for g in range(self.dim):
                if g in pasIndices:
                    if np.abs(self.Coupl[self.varNames[key],g]) < 1e-10:
                        raise ValueError('specify coupling value for '+str(key)+' <- '+str(g))
                else:
                    if np.abs(self.Coupl[self.varNames[key],g]) > 1e-10:
                        raise ValueError('there should be no coupling value for '+str(key)+' <- '+str(g))
            if self.verbosity > 1:
                settings.m(0,'...'+key)
                settings.m(0,rule)
                settings.m(0,rule_pa)
            # now evaluate coefficients
            for tuple in list(itertools.product([False,True],repeat=len(self.pas[key]))):
                if self.process_rule(rule,self.pas[key],tuple):
                    self.boolCoeff[key].append(tuple)
            #
            if self.verbosity > 1:
                settings.m(0,self.boolCoeff[key])

    def process_rule(self,rule,pa,tuple):
        ''' Process a string that denotes a boolean rule.
        '''
        for i,v in enumerate(tuple):
            rule = rule.replace(pa[i],str(v))
        return eval(rule)

    def write_data(self,X,dir='sim/test',noiseObs=0.0,append=False,
                  branching=False,nrRealizations=1,seed=0):
        header = self.header
        tmax = int(X.shape[0])
        header += 'tmax = ' + str(tmax) + '\n'
        header += 'branching = ' + str(branching) + '\n'
        header += 'nrRealizations = ' + str(nrRealizations) + '\n'
        header += 'noiseObs = ' + str(noiseObs) + '\n'
        header += 'noiseDyn = ' + str(self.noiseDyn) + '\n'
        header += 'seed = ' + str(seed) + '\n'
        # add observational noise
        X += noiseObs*np.random.randn(tmax,self.dim)
        # call helper function
        write_data(X,dir,append,header,
                  varNames=self.varNames,
                  Adj=self.Adj,Coupl=self.Coupl,
                  model=self.model,modelType=self.modelType,
                  boolRules=self.boolRules,invTimeStep=self.invTimeStep)

def _check_branching(X,Xsamples,restart,threshold=0.25):
    """\
    Check whether time series branches.

    Parameters
    ----------
    X (np.array): current time series data.
    Xsamples (np.array): list of previous branching samples.
    restart (int): counts number of restart trials.
    threshold (float, optional): sets threshold for attractor
        identification.

    Returns
    -------
    check : bool
        true if branching realization
    Xsamples
        updated list
    """
    check = True
    if restart == 0:
        Xsamples.append(X)
    else:
        for Xcompare in Xsamples:
            Xtmax_diff = np.absolute(X[-1,:] - Xcompare[-1,:])
            # If the second largest element is smaller than threshold
            # set check to False, i.e. at least two elements
            # need to change in order to have a branching.
            # If we observe all parameters of the system,
            # a new attractor state must involve changes in two
            # variables.
            if np.partition(Xtmax_diff,-2)[-2] < threshold:
                check = False
        if check:
            Xsamples.append(X)
    if not check:
        logg.m('realization {}:'.format(restart), 'no new branch', v=4)
    else:
        logg.m('realization {}:'.format(restart), 'new branch', v=4)
    return check, Xsamples

def check_nocycles(Adj, verbosity=2):
    """\
    Checks that there are no cycles in graph described by adjacancy matrix.

    Parameters
    ----------
    Adj (np.array): adjancancy matrix of dimension (dim, dim)

    Returns
    -------
    True if there is no cycle, False otherwise.
    """
    dim = Adj.shape[0]
    for g in range(dim):
        v = np.zeros(dim)
        v[g] = 1
        for i in range(dim):
            v = Adj.dot(v)
            if v[g] > 1e-10:
                if verbosity > 2:
                    settings.m(0,Adj)
                    settings.m(0,'contains a cycle of length',i+1,
                          'starting from node',g,
                          '-> reject')
                return False
    return True


def sample_coupling_matrix(dim=3,connectivity=0.5):
    """\
    Sample coupling matrix.

    Checks that returned graphs contain no self-cycles.

    Parameters
    ----------
    dim : int
        dimension of coupling matrix.
    connectivity : float
        fraction of connectivity, fully connected means 1.,
        not-connected means 0, in the case of fully connected, one has
        dim*(dim-1)/2 edges in the graph.

    Returns
    -------
    Tuple (Coupl,Adj,Adj_signed) of coupling matrix, adjancancy and
    signed adjacancy matrix.
    """
    max_trial = 10
    check = False
    for trial in range(max_trial):
        # random topology for a given connectivity / edge density
        Coupl = np.zeros((dim,dim))
        n_edges = 0
        for gp in range(dim):
            for g in range(dim):
                if gp != g:
                    # need to have the factor 0.5, otherwise
                    # connectivity=1 would lead to dim*(dim-1) edges
                    if np.random.rand() < 0.5*connectivity:
                        Coupl[gp,g] = 0.7
                        n_edges += 1
        # obtain adjacancy matrix
        Adj_signed = np.zeros((dim,dim),dtype='int_')
        Adj_signed = np.sign(Coupl)
        Adj = np.abs(Adj_signed)
        # check for cycles and whether there is at least one edge
        if check_nocycles(Adj) and n_edges > 0:
            check = True
            break
    if not check:
        raise ValueError('did not find graph without cycles after',
                         max_trial,'trials')
    return Coupl, Adj, Adj_signed, n_edges

class StaticCauseEffect:
    """
    Simulates static data to investigate structure learning.
    """

    availModels = odict([
        ('line', 'y = αx \n'),
        ('noise', 'y = noise \n'),
        ('absline', 'y = |x| \n'),
        ('parabola', 'y = αx² \n'),
        ('sawtooth', 'y = x - |x| \n'),
        ('tanh', 'y = tanh(x) \n'),
        ('combi', 'combinatorial regulation \n'),
    ])

    def __init__(self):

        # define a set of available functions
        self.funcs = {
            'line': lambda x: x,
            'noise': lambda x: 0,
            'absline': lambda x: np.abs(x),
            'parabola': lambda x: x**2,
            'sawtooth': lambda x: 0.5*x - np.floor(0.5*x),
            'tanh': lambda x: np.tanh(2*x),
        }

    def sim_givenAdj(self, Adj: np.array, model='line'):
        """\
        Simulate data given only an adjacancy matrix and a model.

        The model is a bivariate funtional dependence. The adjacancy matrix
        needs to be acyclic.

        Parameters
        ----------
        Adj
            adjacancy matrix of shape (dim,dim).

        Returns
        -------
        Data array of shape (n_samples,dim).
        """
        # nice examples
        examples = [{'func' : 'sawtooth', 'gdist' : 'uniform',
                     'sigma_glob' : 1.8, 'sigma_noise' : 0.1}]

        # nr of samples
        n_samples = 100

        # noise
        sigma_glob = 1.8
        sigma_noise = 0.4

        # coupling function / model
        func = self.funcs[model]

        # glob distribution
        sourcedist = 'uniform'

        # loop over source nodes
        dim = Adj.shape[0]
        X = np.zeros((n_samples,dim))
        # source nodes have no parents themselves
        nrpar = 0
        children = list(range(dim))
        parents = []
        for gp in range(dim):
            if Adj[gp,:].sum() == nrpar:
                if sourcedist == 'gaussian':
                    X[:,gp] = np.random.normal(0,sigma_glob,n_samples)
                if sourcedist == 'uniform':
                    X[:,gp] = np.random.uniform(-sigma_glob,sigma_glob,n_samples)
                parents.append(gp)
                children.remove(gp)

        # all of the following guarantees for 3 dim, that we generate the data
        # in the correct sequence
        # then compute all nodes that have 1 parent, then those with 2 parents
        children_sorted = []
        nrchildren_par = np.zeros(dim)
        nrchildren_par[0] = len(parents)
        for nrpar in range(1,dim):
            # loop over child nodes
            for gp in children:
                if Adj[gp,:].sum() == nrpar:
                    children_sorted.append(gp)
                    nrchildren_par[nrpar] += 1
        # if there is more than a child with a single parent
        # order these children (there are two in three dim)
        # by distance to the source/parent
        if nrchildren_par[1] > 1:
            if Adj[children_sorted[0],parents[0]] == 0:
                help = children_sorted[0]
                children_sorted[0] = children_sorted[1]
                children_sorted[1] = help

        for gp in children_sorted:
            for g in range(dim):
                if Adj[gp,g] > 0:
                    X[:,gp] += 1./Adj[gp,:].sum()*func(X[:,g])
            X[:,gp] += np.random.normal(0,sigma_noise,n_samples)

#         fig = pl.figure()
#         fig.add_subplot(311)
#         pl.plot(X[:,0],X[:,1],'.',mec='white')
#         fig.add_subplot(312)
#         pl.plot(X[:,1],X[:,2],'.',mec='white')
#         fig.add_subplot(313)
#         pl.plot(X[:,2],X[:,0],'.',mec='white')
#         pl.show()

        return X

    def sim_combi(self):
        """ Simulate data to model combi regulation.
        """
        n_samples = 500
        sigma_glob = 1.8

        X = np.zeros((n_samples,3))

        X[:,0] = np.random.uniform(-sigma_glob,sigma_glob,n_samples)
        X[:,1] = np.random.uniform(-sigma_glob,sigma_glob,n_samples)

        func = self.funcs['tanh']

        # XOR type
#         X[:,2] = (func(X[:,0])*sp.stats.norm.pdf(X[:,1],0,0.2)
#                   + func(X[:,1])*sp.stats.norm.pdf(X[:,0],0,0.2))
        # AND type / diagonal
#         X[:,2] = (func(X[:,0]+X[:,1])*sp.stats.norm.pdf(X[:,1]-X[:,0],0,0.2))
        # AND type / horizontal
        X[:,2] = (func(X[:,0])*sp.stats.norm.cdf(X[:,1],1,0.2))

        pl.scatter(X[:,0],X[:,1],c=X[:,2],edgecolor='face')
        pl.show()

        pl.plot(X[:,1],X[:,2],'.')
        pl.show()

        return X

def sample_static_data(model,dir,verbosity=0):
    # fraction of connectivity as compared to fully connected
    # in one direction, which amounts to dim*(dim-1)/2 edges
    connectivity = 0.8
    dim = 3
    n_Coupls = 50
    model = model.replace('static-','')
    np.random.seed(0)

    if model != 'combi':
        n_edges = np.zeros(n_Coupls)
        for icoupl in range(n_Coupls):
            Coupl, Adj, Adj_signed, n_e = sample_coupling_matrix(dim,connectivity)
            if verbosity > 1:
                settings.m(0,icoupl)
                settings.m(0,Adj)
            n_edges[icoupl] = n_e
            # sample data
            X = StaticCauseEffect().sim_givenAdj(Adj,model)
            write_data(X,dir,Adj=Adj)
        settings.m(0,'mean edge number:',n_edges.mean())

    else:
        X = StaticCauseEffect().sim_combi()
        Adj = np.zeros((3,3))
        Adj[2,0] = Adj[2,1] = 0
        write_data(X,dir,Adj=Adj)

if __name__ == '__main__':
    import argparse
#     epilog = ('    1: 2dim, causal direction X_1 -> X_0, constraint signs\n'
#               + '    2: 2dim, causal direction X_1 -> X_0, arbitrary signs\n'
#               + '    3: 2dim, causal direction X_1 <-> X_0, arbitrary signs\n'
#               + '    4: 2dim, mix of model 2 and 3\n'
#               + '    5: 6dim double toggle switch\n'
#               + '    6: two independent evolutions without repression, sync.\n'
#               + '    7: two independent evolutions without repression, random init\n'
#               + '    8: two independent evolutions directed repression, random init\n'
#               + '    9: two independent evolutions mutual repression, random init\n'
#               + '   10: two indep. evol., diff. self-loops possible, mut. repr., rand init\n')
    epilog = ''
    for k,v in StaticCauseEffect.availModels.items():
        epilog += '    static-' + k + ': ' + v
    for k,v in GRNsim.availModels.items():
        epilog += '    ' + k + ': ' + v
    # command line options
    p = argparse.ArgumentParser(
        description=('Simulate stochastic discrete-time dynamical systems,\n'
                     'in particular gene regulatory networks.'),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=('  MODEL: specify one of the following models, or one of \n'
                '    the filenames (without ".txt") in the directory "models" \n'
                + epilog
                ))
    aa = p.add_argument
    aa('--dir',required=True,
        type=str,default='',
        help=('specify directory to store data, '
              + ' must start with "sim/MODEL_...", see possible values for MODEL below '))
    aa('--show',
        action='store_true',
        help='show plots')
    aa('--verbosity',
        type=int,default=0,
        help='specify integer > 0 to get more output [default 0]')
    args = p.parse_args()

    # run checks on output directory
    dir = args.dir
    if not dir.startswith('sim/'):
        raise IOError('prepend "sim/..." to --dir argument,'
                      + '"..." being an arbitrary string')
    else:
        model = dir.split('/')[1].split('_')[0]
        settings.m(0,'...model is: "'+model+'"')
    if os.path.exists(dir) and 'test' not in dir:
        message = ('directory ' + dir +
                   ' already exists, remove it and continue? [y/n, press enter]')
        if str(input(message)) != 'y':
            settings.m(0,'    ...quit program execution')
            quit()
        else:
            settings.m(0,'   ...removing directory and continuing...')
            os.system('rm -r ' + dir)

    settings.m(0,model)
    settings.m(0,dir)

    # sample data
    if 'static' in model:
        sample_static_data(model=model,dir=dir,verbosity=args.verbosity)
    else:
        sample_dynamic_data(model=model,dir=dir)
