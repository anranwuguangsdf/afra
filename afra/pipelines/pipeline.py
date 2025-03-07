import numpy as np
import matplotlib.pyplot as plt
from afra.models.fg_models import *
from afra.models.bg_models import *
from afra.tools.ps_estimator import pstimator
from afra.tools.aux import gvec, empcov


class pipe(object):

    def __init__(self, data, noises=None, mask=None, beams=None, targets=('TT',),
                 fiducials=None, fiducial_beams=None,
                 templates=None, template_noises=None, template_beams=None,
                 foreground=None, background=None,
                 likelihood='gauss', solver='dynesty', filt=None):
        """
        Parameters
        ----------

        data : dict
            Measured data maps,
            should be arranged in form {frequency (GHz): (map #, pixel #)}.

        noises : dict
            Simulated noise map samples,
            should be arranged in form: {frequency (GHz): (sample #, map #, pixel #)}.

        fiducials : numpy.array
            Simulated fiducial map samples,
            should be arranged in form {frequency: (sample #, map #, pixel #)}.

        mask : numpy.ndarray
            Single mask map,
            should be arranged in form (1, pixel #).

        beams : dict
            FWHM of gaussian beams for each frequency.

        fiducials : dict
            Fiducial map dict,
            should be arranged in form {frequency: (sample #, map #, pixel #)}.

        fiducial_beams : dict
            Fiducial map fwhm dict,
            should be arranged in form {frequency: fwhm}.

        templates : dict
            Template map dict,
            should be arranged in form {frequency: (map #, pixel #)}.

        template_noises : dict
            Template noise map dict,
            should be arranged in form {frequency: (sample #, map #, pixel #)}.

        template_beams : dict
            Template map fwhm dict,
            should be arranged in form {frequency: fwhm}.

        targets : tuple
            Tuple consists of 'TT', 'TE', 'TB', 'EE', 'EB', 'BB'.
            'TT': TT spectra
            'EE': EE spectra
            'BB': BB spectra
            'TE': TE spectra
            'TB': TB spectra
            'EB': EB spectra

        likelihood : str
            likelihood type, should be either 'gauss' or 'hl'.

        foreground : str
            foreground model name, chosen among "dust", "sync", "syncdust".

        background : str
            background model name, chosen among "acmb", "ncmb".

        filt : dict
            filtering matrix for CMB band-power (from original to filted),
            entry name should contain "targets".

        solver : str
            Bayesian solver name.
        """
        # measurements
        self.data = data
        self.data_bp = None  # data BP
        self.noises = noises
        self.noise_bp = None  # noise BP
        self.mask = mask
        self.beams = beams
        self.targets = targets
        # fiducials
        self.fiducials = fiducials
        self.fiducial_beams = fiducial_beams
        self.fiducial_bp = None  # fiducial BP
        # templates
        self.templates = templates
        self.template_noises = template_noises
        self.template_beams = template_beams
        self.template_bp = None  # template BP
        # choose fore-/back-ground models
        self._foreground_catalog = {'async':asyncmodel,'tsync':tsyncmodel,
                                    'adust':adustmodel,'tdust':tdustmodel,
                                    'asyncadust':asyncadustmodel,'tsynctdust':tsynctdustmodel}
        self._background_catalog = {'acmb':acmbmodel,'ncmb':ncmbmodel}
        self.foreground = foreground  # after catalog
        self.background = background  # after catalog
        self.foreground_obj = None
        self.background_obj = None
        # choose likelihood method
        self.likelihood = likelihood
        # init parameter list
        self.paramlist = list()
        self.paramrange = dict()
        # ps estimator, to be assigned
        self.estimator = None
        # solver name
        self.solver = solver
        # filtering matrix dict
        self.filt = filt
        # covariance matrix
        self.covmat = None

    @property
    def data(self):
        return self._data

    @property
    def data_bp(self):
        return self._data_bp

    @property
    def noises(self):
        return self._noises

    @property
    def noise_bp(self):
        return self._noise_bp

    @property
    def noise_nsamp(self):
        return self._noise_nsamp

    @property
    def noise_flag(self):
        return self._noise_flag

    @property
    def mask(self):
        return self._mask
 
    @property
    def freqlist(self):
        return self._freqlist

    @property
    def nfreq(self):
        return self._nfreq

    @property
    def nside(self):
        return self._nside

    @property
    def npix(self):
        return self._npix

    @property
    def targets(self):
        return self._targets

    @property
    def ntarget(self):
        return self._ntarget

    @property
    def beams(self):
        return self._beams

    @property
    def filt(self):
        return self._filt

    @property
    def fiducials(self):
        return self._fiducials

    @property
    def fiducial_nsamp(self):
        return self._fiducial_nsamp

    @property
    def fiducial_flag(self):
        return self._fiducial_flag

    @property
    def fiducial_beams(self):
        return self._fiducial_beams

    @property
    def fiducial_bp(self):
        return self._fiducial_bp

    @property
    def templates(self):
        return self._templates

    @property
    def template_flag(self):
        return self._template_flag

    @property
    def template_noises(self):
        return self._template_noises

    @property
    def template_nsamp(self):
        return self._template_nsamp

    @property
    def template_beams(self):
        return self._template_beams

    @property
    def template_freqlist(self):
        return self._template_freqlist

    @property
    def template_nfreq(self):
        return self._template_nfreq

    @property
    def template_bp(self):
        return self._template_bp

    @property
    def background(self):
        return self._background

    @property
    def background_catalog(self):
        return self._background_catalog

    @property
    def foreground(self):
        return self._foreground

    @property
    def foreground_catalog(self):
        return self._foreground_catalog

    @property
    def likelihood(self):
        return self._likelihood

    @property
    def paramlist(self):
        return self._paramlist

    @property
    def paramrange(self):
        return self._paramrange

    @property
    def estimator(self):
        return self._estimator

    @property
    def background_obj(self):
        return self._background_obj

    @property
    def foreground_obj(self):
        return self._foreground_obj

    @property
    def covmat(self):
        return self._covmat

    @property
    def solver(self):
        return self._solver

    @data.setter
    def data(self, data):
        assert isinstance(data, dict)
        self._freqlist = sorted(data.keys())
        self._nfreq = len(self._freqlist)
        assert (data[next(iter(data))].ndim == 2)
        assert (data[next(iter(data))].shape[0] == 3)
        self._npix = data[next(iter(data))].shape[1]
        self._nside = int(np.sqrt(self._npix//12))
        self._data = data.copy()

    @noises.setter
    def noises(self, noises):
        if noises is not None:
            assert isinstance(noises, dict)
            assert (noises.keys() == self._data.keys())
            assert (noises[next(iter(noises))].ndim == 3)
            self._noise_nsamp = noises[next(iter(noises))].shape[0]
            assert (noises[next(iter(noises))].shape[1] == 3)
            assert (noises[next(iter(noises))].shape[2] == self._npix)
            self._noises = noises.copy()
            self._noise_flag = True
        else:
            self._noises = None
            self._noise_flag = False

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = np.ones(self._npix,dtype=np.float64)
        else:
            assert isinstance(mask, np.ndarray)
            assert (len(mask) == self._npix)
            self._mask = mask.copy()
        # clean up input maps with mask
        for f in self._freqlist:
            self._data[f][:,self._mask==0.] = 0.
            if self._noise_flag:
                    self._noises[f][:,:,self._mask==0.] = 0.

    @beams.setter
    def beams(self, beams):
        if beams is not None:
            assert isinstance(beams, dict)
            assert (beams.keys() == self._data.keys())
            self._beams = beams.copy()
        else:
            self._beams = dict()
            for f in self._freqlist:
                self._beams[f] = None

    @targets.setter
    def targets(self, targets):
        assert isinstance(targets, tuple)
        self._targets = targets
        self._ntarget = len(self._targets)

    @fiducials.setter
    def fiducials(self, fiducials):
        if fiducials is not None:
            assert isinstance(fiducials, dict)
            assert (fiducials.keys() == self._data.keys())
            assert (fiducials[next(iter(fiducials))].ndim == 3)
            self._fiducial_nsamp = fiducials[next(iter(fiducials))].shape[0]
            assert (fiducials[next(iter(fiducials))].shape[1] == 3)
            assert (fiducials[next(iter(fiducials))].shape[2] == self._npix)
            self._fiducial_flag = True
            self._fiducials = fiducials.copy()
        else:
            self._fiducial_flag = False
            self._fiducials = None

    @fiducial_beams.setter
    def fiducial_beams(self, fiducial_beams):
        if fiducial_beams is not None:
            assert isinstance(fiducial_beams, dict)
            assert (fiducial_beams.keys() == self._fiducials.keys())
            self._fiducial_beams = fiducial_beams.copy()
        elif self._fiducial_flag:
            self._fiducial_beams = dict()
            for name in list(self._fiducials.keys()):
                self._fiducial_beams[name] = None
        else:
            self._fiducial_beams = None

    @templates.setter
    def templates(self, templates):
        if templates is not None:
            assert isinstance(templates, dict)
            self._template_freqlist = sorted(templates.keys())
            self._template_nfreq = len(self._template_freqlist)
            assert (self._template_nfreq==1 or self._template_nfreq==2)
            assert (templates[next(iter(templates))].shape == (3,self._npix))
            self._template_flag = True
            self._templates = templates.copy()
        else:
            self._template_flag = False
            self._template_freqlist = list()
            self._template_nfreq = 0
            self._templates = None

    @template_noises.setter
    def template_noises(self, template_noises):
        if template_noises is not None:
            assert isinstance(template_noises, dict)
            assert (template_noises.keys() == self._templates.keys())
            assert (len(template_noises) == self._template_nfreq)
            self._template_nsamp = len(template_noises[next(iter(template_noises))])
            self._template_noises = template_noises.copy()
        else:
            self._template_noises = None
            self._template_nsamp = 0

    @template_beams.setter
    def template_beams(self, template_beams):
        if template_beams is not None:
            assert isinstance(template_beams, dict)
            assert (template_beams.keys() == self._templates.keys())
            self._template_beams = template_beams.copy()
        else:
            self._template_beams = dict()
            if self._template_flag:
                for name in self._template_freqlist:
                    self._template_beams[name] = None

    @filt.setter
    def filt(self, filt):
        if filt is not None:
            assert isinstance(filt, dict)
            assert (self._targets in filt)
        self._filt = filt

    @likelihood.setter
    def likelihood(self, likelihood):
        assert isinstance(likelihood, str)
        self._likelihood = likelihood

    @background.setter
    def background(self, background):
        if background is None:
            self._background = None
        else:
            assert (background in self._background_catalog)
            self._background = self._background_catalog[background]

    @foreground.setter
    def foreground(self, foreground):
        
        if foreground is None:
            self._foreground = None
        else:
            assert (foreground in self._foreground_catalog)
            self._foreground = self._foreground_catalog[foreground]

    @paramlist.setter
    def paramlist(self, paramlist):
        assert isinstance(paramlist, list)
        self._paramlist = paramlist

    @paramrange.setter
    def paramrange(self, paramrange):
        assert isinstance(paramrange, dict)
        self._paramrange = paramrange

    @data_bp.setter
    def data_bp(self, data_bp):
        if data_bp is not None:
            assert (data_bp.shape == (self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq))
        self._data_bp = data_bp

    @noise_bp.setter
    def noise_bp(self, noise_bp):
        if noise_bp is not None:
            assert (noise_bp.shape == (self._noise_nsamp,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq))
        self._noise_bp = noise_bp

    @fiducial_bp.setter
    def fiducial_bp(self, fiducial_bp):
        if fiducial_bp is not None:
            assert (fiducial_bp.shape == (self._fiducial_nsamp,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq))
        self._fiducial_bp = fiducial_bp

    @template_bp.setter
    def template_bp(self, template_bp):
        if template_bp is not None:
            assert isinstance(template_bp, dict)
        self._template_bp = template_bp

    @estimator.setter
    def estimator(self, estimator):
        if estimator is not None:
            assert isinstance(estimator, pstimator)
        self._estimator = estimator

    @background_obj.setter
    def background_obj(self, background_obj):
        if background_obj is not None:
            assert isinstance(background_obj, bgmodel)
        self._background_obj = background_obj

    @foreground_obj.setter
    def foreground_obj(self, foreground_obj):
        if foreground_obj is not None:
            assert isinstance(foreground_obj, fgmodel)
        self._foreground_obj = foreground_obj

    @covmat.setter
    def covmat(self, covmat):
        if covmat is not None:
            assert (covmat.ndim == 2)
            assert (covmat.shape[0] == self._ntarget*self._estimator.nmode*self._nfreq*(self._nfreq+1)//2)
            assert (covmat.shape[1] == covmat.shape[0])
        self._covmat = covmat

    @solver.setter
    def solver(self, solver):
        assert (solver in ('minuit','dynesty','emcee'))
        self._solver = solver

    def preprocess(self, aposcale, psbin, lmin=None, lmax=None):
        """
        preprocess routine, converts maps into band-powers.

        Parameters
        ----------

        aposcale : float
            Apodization scale in deg.

        psbin : integer
            Number of bins,
            for conducting pseudo-PS estimation.

        lmin/lmax : integer
            Lower/Upper multipole limit.
        """
        assert isinstance(aposcale, float)
        assert isinstance(psbin, int)
        assert (psbin > 0)
        assert (aposcale > 0)
        # STEP I
        # init PS estimator
        self.estimator = pstimator(nside=self._nside,mask=self._mask,aposcale=aposcale,psbin=psbin,lmin=lmin,lmax=lmax,targets=self._targets,filt=self._filt)
        # STEP II
        # template PS estimations (auto corr. only)
        if self._template_flag:
            self.template_bp = dict()
            for i in range(self._template_nfreq):
                # allocate for template
                data_bp = np.zeros((self._ntarget,self._estimator.nmode),dtype=np.float64)
                noise_bp = np.zeros((self._noise_nsamp,self._ntarget,self._estimator.nmode),dtype=np.float64)
                _fi = self._template_freqlist[i]
                # template workspace
                twsp = self._estimator.autoWSP(self._templates[_fi],beams=self._template_beams[_fi])
                # template auto-corr.
                stmp = self._estimator.autoBP(self._templates[_fi],wsp=twsp,beams=self._template_beams[_fi])
                data_bp = np.array(stmp[1:1+self._ntarget])
                # template noise auto corr.
                for s in range(self._template_nsamp):
                    ntmp = self._estimator.autoBP(self._template_noises[_fi][s],wsp=twsp,beams=self._template_beams[_fi])
                    noise_bp[s] = np.array(ntmp[1:1+self._ntarget])
                # mean noise subtraction
                self._template_bp[_fi] = data_bp - np.mean(noise_bp,axis=0)
        # STEP III
        # prepare model, parameter list generated during init models
        if self._background is not None:
            self._background_obj = self._background(self._freqlist,self._estimator)
        if self._foreground is not None:
            self._foreground_obj = self._foreground(self._freqlist,self._estimator,self._template_bp)
        # STEP IV-A
        # data PS estimations (with workspace)
        # allocate
        wsp_dict = dict()  # data nmt workspace
        fwsp_dict = dict()  # fiducial nmt workspace
        self.data_bp = np.zeros((self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
        for i in range(self._nfreq):
            _fi = self._freqlist[i]
            # auto corr.
            wsp_dict[(i,i)] = self._estimator.autoWSP(self._data[_fi],beams=self._beams[_fi])
            stmp = self._estimator.autoBP(self._data[_fi],wsp=wsp_dict[(i,i)],beams=self._beams[_fi])
            self._data_bp[:,:,i,i] = np.array(stmp[1:1+self._ntarget])
            for j in range(i+1,self._nfreq):
                _fj = self._freqlist[j]
                # cross corr.
                wsp_dict[(i,j)] = self._estimator.crosWSP(np.r_[self._data[_fi],self._data[_fj]],beams=[self._beams[_fi],self._beams[_fj]])
                stmp = self._estimator.crosBP(np.r_[self._data[_fi],self._data[_fj]],wsp=wsp_dict[(i,j)],beams=[self._beams[_fi],self._beams[_fj]])
                self._data_bp[:,:,i,j] = np.array(stmp[1:1+self._ntarget])
                self._data_bp[:,:,j,i] = np.array(stmp[1:1+self._ntarget])
        # STEP IV-B
        # fiducial PS estimation 
        if self._fiducial_flag:
            # allocate
            self.fiducial_bp = np.zeros((self._fiducial_nsamp,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
            for i in range(self._nfreq):
                _fi = self._freqlist[i]
                # auto corr.
                fwsp_dict[(i,i)] = self._estimator.autoWSP(self._fiducials[_fi][0],beams=self._fiducial_beams[_fi])
                for s in range(self._fiducial_nsamp):
                    # auto corr.
                    ftmp = self._estimator.autoBP(self._fiducials[_fi][s],wsp=fwsp_dict[(i,i)],beams=self._fiducial_beams[_fi])
                    self._fiducial_bp[s,:,:,i,i] = np.array(ftmp[1:1+self._ntarget])
                for j in range(i+1,self._nfreq):
                    _fj = self._freqlist[j]
                    fwsp_dict[(i,j)] = self._estimator.crosWSP(np.r_[self._fiducials[_fi][0],self._fiducials[_fj][0]],beams=[self._fiducial_beams[_fi],self._fiducial_beams[_fj]])
                    for s in range(self._fiducial_nsamp):
                        # cross corr.
                        ftmp = self._estimator.crosBP(np.r_[self._fiducials[_fi][s],self._fiducials[_fj][s]],wsp=fwsp_dict[(i,j)],beams=[self._fiducial_beams[_fi],self._fiducial_beams[_fj]])
                        self._fiducial_bp[s,:,:,i,j] = np.array(ftmp[1:1+self._ntarget])
                        self._fiducial_bp[s,:,:,j,i] = np.array(ftmp[1:1+self._ntarget])
        # STEP IV-C
        # noise PS estimations
        if self._noise_flag:
            # allocate
            self.noise_bp = np.zeros((self._noise_nsamp,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
            for s in range(self._noise_nsamp):
                for i in range(self._nfreq):
                    _fi = self._freqlist[i]
                    # auto corr.
                    ntmp = self._estimator.autoBP(self._noises[_fi][s],wsp=wsp_dict[(i,i)],beams=self._beams[_fi])
                    self._noise_bp[s,:,:,i,i] = np.array(ntmp[1:1+self._ntarget])
                    for j in range(i+1,self._nfreq):
                        _fj = self._freqlist[j]
                        # cross corr.
                        ntmp = self._estimator.crosBP(np.r_[self._noises[_fi][s],self._noises[_fj][s]],wsp=wsp_dict[(i,j)],beams=[self._beams[_fi],self._beams[_fj]])
                        self._noise_bp[s,:,:,i,j] = np.array(ntmp[1:1+self._ntarget])
                        self._noise_bp[s,:,:,j,i] = np.array(ntmp[1:1+self._ntarget])
        # STEP V
        # fiducial+noise PS covariance matrix
        # fiducial+noise has to be processed in the pixel doamin, in order to yield a proper cov matrix
        if self._fiducial_flag and self._noise_flag:
            ncom = min(self._noise_nsamp,self._fiducial_nsamp)
            # allocate
            nfid = np.zeros((ncom,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
            for s in range(ncom):
                for i in range(self._nfreq):
                    _fi = self._freqlist[i]
                    # auto corr.
                    ntmp = self._estimator.autoBP(self._fiducials[_fi][s]+self._noises[_fi][s],wsp=fwsp_dict[(i,i)],beams=self._fiducial_beams[_fi])
                    nfid[s,:,:,i,i] = np.array(ntmp[1:1+self._ntarget])
                    for j in range(i+1,self._nfreq):
                        _fj = self._freqlist[j]
                        # cross corr.
                        ntmp = self._estimator.crosBP(np.r_[self._fiducials[_fi][s]+self._noises[_fi][s],self._fiducials[_fj][s]+self._noises[_fj][s]],wsp=fwsp_dict[(i,j)],beams=[self._fiducial_beams[_fi],self._fiducial_beams[_fj]])
                        nfid[s,:,:,i,j] = np.array(ntmp[1:1+self._ntarget])
                        nfid[s,:,:,j,i] = np.array(ntmp[1:1+self._ntarget])
            # full cov
            self.covmat = empcov(gvec(nfid))

    def reprocess(self, data):
        """
        replace data maps only
        """
        # read new data dict
        assert isinstance(data, dict)
        assert (data.keys() == self._data.keys())
        assert (data[next(iter(data))].shape == (3,self._npix))
        assert self._data_bp is not None  # check preprocess
        self.data_bp = np.zeros((self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
        for i in range(self._nfreq):
            _fi = self._freqlist[i]
            stmp = self._estimator.autoBP(data[_fi],beams=self._beams[_fi])
            self._data_bp[:,:,i,i] = np.array(stmp[1:1+self._ntarget])
            for j in range(i+1,self._nfreq):
                _fj = self._freqlist[j]
                stmp = self._estimator.crosBP(np.r_[data[_fi],data[_fj]],beams=[self._beams[_fi],self._beams[_fj]])
                self._data_bp[:,:,i,j] = np.array(stmp[1:1+self._ntarget])
                self._data_bp[:,:,j,i] = np.array(stmp[1:1+self._ntarget])

    def plot_data(self):
        cc = ['firebrick', 'midnightblue', 'darkorange'][:self._ntarget]
        fig = plt.figure(figsize=(5*self._nfreq,5*self._nfreq))
        for i in range(self._nfreq):
            for j in range(i,self._nfreq):
                ax = fig.add_subplot(self._nfreq,self._nfreq,1+j*self._nfreq+i)
                ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[j]),fontsize=15)
                ax.tick_params(axis='both',which='major',labelsize='15')
                ax.set_yscale('log')
                for k in range(self._ntarget):
                    ax.scatter(self._estimator.modes,self._data_bp[k,:,i,j],marker='o',s=20,c=cc[k])
        labels = list()
        for t in self._targets:
            labels.append('{} data'.format(t))
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('data_bp.pdf')

    def plot_fiducial(self):
        cc = ['firebrick', 'midnightblue', 'darkorange'][:self._ntarget]
        fig = plt.figure(figsize=(5*self._nfreq,5*self._nfreq))
        for i in range(self._nfreq):
            for j in range(i,self._nfreq):
                ax = fig.add_subplot(self._nfreq,self._nfreq,1+j*self._nfreq+i)
                ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[j]),fontsize=15)
                ax.tick_params(axis='both',which='major',labelsize='15')
                ax.set_yscale('log')
                for k in range(self._ntarget):
                    ax.errorbar(self._estimator.modes,np.mean(self._fiducial_bp[:,k,:,i,j],axis=0),yerr=np.std(self._fiducial_bp[:,k,:,i,j],axis=0),fmt='o',ms=10,mfc='none',mec=cc[k])
        labels = list()
        for t in self._targets:
            labels.append('{} fiducial'.format(t))
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('fiducial_bp.pdf')

    def plot_noise(self):
        cc = ['firebrick', 'midnightblue', 'darkorange'][:self._ntarget]
        fig = plt.figure(figsize=(5*self._nfreq,5*self._nfreq))
        for i in range(self._nfreq):
            ax = fig.add_subplot(self._nfreq,self._nfreq,1+i*(self._nfreq+1))
            ax.tick_params(axis='both',which='major',labelsize='15')
            ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[i]),fontsize=15)
            ax.set_yscale('log')
            for k in range(self._ntarget):
                ax.errorbar(self._estimator.modes,np.mean(self._noise_bp[:,k,:,i,i],axis=0),yerr=np.std(self._noise_bp[:,k,:,i,i],axis=0),fmt='o',ms=10,mfc='none',mec=cc[k])
            for j in range(i+1,self._nfreq):
                ax = fig.add_subplot(self._nfreq,self._nfreq,1+j*self._nfreq+i)
                ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[j]))
                for k in range(self._ntarget):
                    ax.errorbar(self._estimator.modes,np.mean(self._noise_bp[:,k,:,i,j],axis=0),yerr=np.std(self._noise_bp[:,k,:,i,j],axis=0),fmt='o',ms=10,mfc='none',mec=cc[k])
        labels = list()
        for t in self._targets:
            labels.append('{} noise'.format(t))
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('noise_bp.pdf')

    def plot_template(self):
        cc = ['firebrick', 'midnightblue', 'darkorange'][:self._ntarget]
        fig = plt.figure(figsize=(5*self._template_nfreq,5))
        for i in range(self._template_nfreq):
            ax = fig.add_subplot(self._template_nfreq,1,1+i)
            ax.set_title(str(self._template_freqlist[i]),fontsize=15)
            ax.tick_params(axis='both',which='major',labelsize='15')
            ax.set_yscale('log')
            for k in range(self._ntarget):
                ax.scatter(self._estimator.modes,self._template_bp[k],marker='o',s=20,c=cc[k])
        labels = list()
        for t in self._targets:
            labels.append('{} template'.format(t))
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('template_bp.pdf')

    def plot_result(self, best_bp):
        cc = ['firebrick', 'midnightblue', 'darkorange'][:self._ntarget]
        fig = plt.figure(figsize=(5*self._nfreq,5*self._nfreq))
        for i in range(self._nfreq):
            ax = fig.add_subplot(self._nfreq,self._nfreq,1+i*(self._nfreq+1))
            ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[i]))
            ax.set_yscale('log')
            for k in range(self._ntarget):
                ax.plot(self._estimator.modes,best_bp[k,:,i,i],ls='--',lw=2,color='k')
                smean = self._data_bp[k,:,i,i]-np.mean(self._noise_bp[:,k,:,i,i],axis=0)
                sstd = np.std(self._fiducial_bp[:,k,:,i,i],axis=0)+np.std(self._noise_bp[:,k,:,i,i],axis=0)
                ax.errorbar(self._estimator.modes,smean,yerr=sstd,fmt='o',ms=10,mfc='none',mec=cc[k])
            for j in range(i+1,self._nfreq):
                ax = fig.add_subplot(self._nfreq,self._nfreq,1+j*self._nfreq+i)
                ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[j]))
                for k in range(self._ntarget):
                    ax.plot(self._estimator.modes,best_bp[k,:,i,j],color='k')
                    smean = self._data_bp[k,:,i,j]-np.mean(self._noise_bp[:,k,:,i,j],axis=0)
                    sstd = np.std(self._fiducial_bp[:,k,:,i,j],axis=0)+np.std(self._noise_bp[:,k,:,i,j],axis=0)
                    ax.errorbar(self._estimator.modes,smean,yerr=sstd,fmt='o',ms=10,mfc='none',mec=cc[k])
        labels = list()
        for t in self._targets:
            labels.append('{} bestfit'.format(t))
        for t in self._targets:
            labels.append('{} signal'.format(t))
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('result_bp.pdf')

    def plot_residule(self, best_bp):
        cc = ['firebrick', 'midnightblue', 'darkorange'][:self._ntarget]
        fig = plt.figure(figsize=(5*self._nfreq,5*self._nfreq))
        for i in range(self._nfreq):
            ax = fig.add_subplot(self._nfreq,self._nfreq,1+i*(self._nfreq+1))
            ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[i]))
            ax.set_yscale('log')
            for k in range(self._ntarget):
                rmean = self._data_bp[k,:,i,i]-np.mean(self._noise_bp[:,k,:,i,i],axis=0)-best_bp[k,:,i,i]
                sstd = np.std(self._fiducial_bp[:,k,:,i,i],axis=0)+np.std(self._noise_bp[:,k,:,i,i],axis=0)
                ax.errorbar(self._estimator.modes,np.mean(self._fiducial_bp[:,k,:,i,i],axis=0),yerr=sstd,fmt='o',ms=10,mfc='none',mec=cc[k])
                ax.plot(self._estimator.modes,rmean,ls='--',lw=2,color='grey')
            for j in range(i+1,self._nfreq):
                ax = fig.add_subplot(self._nfreq,self._nfreq,1+j*self._nfreq+i)
                ax.set_title(str(self._freqlist[i])+'x'+str(self._freqlist[j]))
                for k in range(self._ntarget):
                    rmean = self._data_bp[k,:,i,j]-np.mean(self._noise_bp[:,k,:,i,j],axis=0)-best_bp[k,:,i,j]
                    sstd = np.std(self._fiducial_bp[:,k,:,i,j],axis=0)+np.std(self._noise_bp[:,k,:,i,j],axis=0)
                    ax.errorbar(self._estimator.modes,np.mean(self._fiducial_bp[:,k,:,i,j],axis=0),yerr=sstd,fmt='o',ms=10,mfc='none',mec=cc[k])
                    ax.plot(self._estimator.modes,rmean,ls='--',lw=2,color='grey')
        labels = list()
        for t in self._targets:
            labels.append('{} residual'.format(t))
        for t in self._targets:
            labels.append('{} fiducial'.format(t))
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('residual_bp.pdf')

# end
