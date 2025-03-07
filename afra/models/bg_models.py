import numpy as np
import camb
from afra.tools.ps_estimator import pstimator
from afra.tools.icy_decorator import icy


class bgmodel(object):

    def __init__(self, freqlist, estimator):
        self.freqlist = freqlist
        self.estimator = estimator
        self.params = dict()  # base class holds empty dict
        self.paramdft = dict()
        self.paramrange = dict()
        self.paramlist = list()
        self.blacklist = list()  # fixed parameter list

    @property
    def freqlist(self):
        return self._freqlist

    @property
    def estimator(self):
        return self._estimator

    @property
    def nfreq(self):
        return self._nfreq

    @property
    def params(self):
        return self._params

    @property
    def paramdft(self):
        return self._paramdft

    @property
    def paramrange(self):
        return self._paramrange

    @property
    def paramlist(self):
        return self._paramlist

    @property
    def blacklist(self):
        return self._blacklist

    @freqlist.setter
    def freqlist(self, freqlist):
        assert isinstance(freqlist, (list,tuple))
        self._freqlist = freqlist
        self._nfreq = len(self._freqlist)

    @estimator.setter
    def estimator(self, estimator):
        assert isinstance(estimator, pstimator)
        self._estimator = estimator

    @params.setter
    def params(self, params):
        assert isinstance(params, dict)
        self._params = params

    @paramdft.setter
    def paramdft(self, paramdft):
        assert isinstance(paramdft, dict)
        self._paramdft = paramdft

    @paramrange.setter
    def paramrange(self, paramrange):
        assert isinstance(paramrange, dict)
        self._paramrange = paramrange

    @paramlist.setter
    def paramlist(self, paramlist):
        assert isinstance(paramlist, (list,tuple))
        self._paramlist = paramlist

    @blacklist.setter
    def blacklist(self, blacklist):
        assert isinstance(blacklist, (list,tuple))
        for p in blacklist:
            assert (p in self._paramlist)
        self._blacklist = blacklist

    def reset(self, pdict):
        """(re)set parameters"""
        assert isinstance(pdict, dict)
        for name in pdict.keys():
            if name in self._paramlist:
                self._params.update({name: pdict[name]})


@icy
class ncmbmodel(bgmodel):

    def __init__(self, freqlist, estimator): 
        super(ncmbmodel, self).__init__(freqlist,estimator)
        self.paramlist = self.initlist()
        self.paramrange = self.initrange()
        self.paramdft = self.initdft()

    def initlist(self):
        """parameters are set as
        - bandpower "cmb_x", exponential index of amplitude
        """
        plist = list()
        for t in self._estimator._targets:
            for j in range(len(self._estimator._modes)):
                plist.append('cmb_'+t+'{}'.format(j))
        return plist

    def initrange(self):
        """parameter sampling range,
        in python dict
        {param name : [low limit, high limit]
        """
        prange = dict()
        for i in self._paramlist:
            prange[i] = [0.,1.0e+4]
        return prange

    def initdft(self):
        """register default parameter values
        """
        pdft = dict()
        for key in self._paramrange.keys():
            pdft[key] = 0.5*(self._paramrange[key][0] + self._paramrange[key][1])
        self.reset(pdft)  # update self.params
        return pdft

    def bandpower(self):
        """cross-(frequency)-power-spectrum
        
        Parameters
        ----------
        
        freq_list : float
            list of frequency in GHz
        
        freq_ref : float
            synchrotron template reference frequency
        """
        fiducial_bp = np.zeros((self._estimator._ntarget,self._estimator._nmode),dtype=np.float64)
        for t in range(self._estimator._ntarget):
            for l in range(self._estimator._nmode):
                    fiducial_bp[t,l] = self._params['cmb_'+self._estimator._targets[t]+'{}'.format(l)]
        fiducial_bp = self._estimator.filtrans(fiducial_bp)
        bp_out = np.ones((self._estimator._ntarget,self._estimator._nmode,self._nfreq,self._nfreq),dtype=np.float64)
        for t in range(self._estimator._ntarget):
            for l in range(self._estimator._nmode):
                bp_out[t,l] = fiducial_bp[t,l]
        return bp_out


@icy
class acmbmodel(bgmodel):
    
    def __init__(self, freqlist, estimator):
        super(acmbmodel, self).__init__(freqlist,estimator)
        self.paramlist = self.initlist()
        self.paramrange = self.initrange()
        self.paramdft = self.initdft()
        if not ('EE' in self._estimator.targets):
            self._blacklist.append('AE')
        if not ('BB' in self._estimator.targets):
            self._blacklist.append('r')
            self._blacklist.append('AL')
        # calculate camb template CMB PS with default parameters
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5,ombh2=0.022,omch2=0.122,mnu=0.06,omk=0,tau=0.06)
        pars.InitPower.set_params(As=2e-9,ns=0.965,r=0.)
        pars.set_for_lmax(max(4000,self._estimator._lmax),lens_potential_accuracy=2)
        pars.WantTensors = False
        pars.DoLensing = True
        self._template_sl = camb.get_results(pars)
        self._template_ps = self._template_sl.get_lensed_scalar_cls(max(4000,self._estimator._lmax),CMB_unit='muK',raw_cl=True)
        pars.WantTensors = True
        pars.DoLensing = False
        self._template_sl = camb.get_results(pars)

    @property
    def template_sl(self):
        return self._template_sl

    @property
    def template_ps(self):
        return self._template_ps

    def initlist(self):
        return ['r','AL','AE']

    def initrange(self):
        prange = dict()
        prange['r'] = [0.,1.]
        prange['AL'] = [0.,2.]
        prange['AE'] = [0.,2.]
        return prange

    def initdft(self):
        """register default parameter values
        """
        pdft = {'AE':1.,'AL':1.,'r':0.05}
        self.reset(pdft)
        return pdft

    def bandpower(self):
        enum = {('TT',):[0],('EE',):[1],('BB',):[2],('EE','BB'):[1,2]}
        # fast process for recalculating with new r
        inflation_params = camb.initialpower.InitialPowerLaw()
        inflation_params.set_params(As=2e-9,ns=0.965,r=self._params['r'])
        self._template_sl.power_spectra_from_transfer(inflation_params) #warning OK here, not changing scalars
        fiducial_ps = np.transpose(self._template_sl.get_total_cls(max(4000,self._estimator._lmax),CMB_unit='muK',raw_cl=True)*[1.,self._params['AE'],1.,1.]+self._template_ps*[0.,0.,self._params['AL'],0.])[enum[self._estimator._targets],:self._estimator._lmax+1]
        # from PS to BP
        fiducial_bp = self._estimator.bpconvert(fiducial_ps)
        # impose filtering
        fiducial_bp = self._estimator.filtrans(fiducial_bp)
        # form BP matrix
        bp_out = np.ones((self._estimator._ntarget,self._estimator._nmode,self._nfreq,self._nfreq),dtype=np.float64)
        for t in range(self._estimator._ntarget):
            for l in range(self._estimator._nmode):
                bp_out[t,l] = fiducial_bp[t,l]
        return bp_out
