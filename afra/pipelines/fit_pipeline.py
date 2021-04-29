import numpy as np
from afra.pipelines.pipeline import pipe
from afra.methods.fit import *
from afra.tools.icy_decorator import icy


@icy
class fitpipe(pipe):
    
    def __init__(self, data, noises=None, mask=None, beams=None, targets=('TT',),
                 fiducials=None, fiducial_beams=None,
                 templates=None, template_noises=None, template_beams=None,
                 foreground=None, background=None,
                 likelihood='gauss', solver='dynesty', filt=None):
        super(fitpipe, self).__init__(data,noises,mask,beams,targets,fiducials,fiducial_beams,templates,template_noises,template_beams,foreground,background,likelihood,solver,filt)
        # analyse select dict
        self._anadict = {'gauss':self.analyse_gauss, 'hl':self.analyse_hl}
        # Bayesian engine to be assigned
        self.engine = None

    @property
    def engine(self):
        return self._engine

    @engine.setter
    def engine(self, engine):
        if engine is not None:
            assert isinstance(engine, fit)
        self._engine = engine

    def run(self, aposcale=6., psbin=20, lmin=None, lmax=None, kwargs=dict()):
        self.preprocess(aposcale,psbin,lmin,lmax)
        result = self.analyse(kwargs)
        # visualise data and result
        bestpar = None
        bestbp = None
        if self._solver == 'dynesty':
            from dynesty import plotting as dyplot
            fig,ax = dyplot.cornerplot(result,labels=self._paramlist,quantiles=[0.025, 0.5, 0.975],color='midnightblue',title_fmt='.3f',show_titles=1,smooth=0.04)
            plt.savefig('posterior.pdf')
            bestpar = result.samples[np.where(result['logl']==max(result['logl']))][0]
        elif self._solver == 'emcee':
            import corner
            fig = corner.corner(result,labels=self._paramlist);
            plt.savefig('posterior.pdf')
            bestpar = np.median(result,axis=0)
        elif self._solver == 'minuit':
            print ('params: {}\n bestfit: {}\n stderr: {}'.format(self._paramlist,result[0],result[1]))
            bestpar = result[0]
        else:
            raise ValueError('unsupported solver: {}'.format(self._solver))
        for i in range(len(bestpar)):
            if self._foreground_obj is not None:
                self._foreground_obj.reset({self._paramlist[i]: bestpar[i]})
            if self._background_obj is not None:
                self._background_obj.reset({self._paramlist[i]: bestpar[i]})
        if self._foreground_obj is None:
            bestbp = self._background_obj.bandpower()
        elif self._background_obj is None:
            bestbp = self._foreground_obj.bandpower()
        else:
            bestbp = self._foreground_obj.bandpower() + self._background_obj.bandpower()
        self.plot_result(bestbp)
        self.plot_residule(bestbp)
        return result

    def analyse(self, kwargs=dict()):
        return self._anadict[self._likelihood](kwargs)

    def analyse_gauss(self, kwargs=dict()):
        # gauss likelihood
        self.engine = gaussfit(self._data_bp,np.mean(self._fiducial_bp,axis=0),np.mean(self._noise_bp,axis=0),self._covmat,self._background_obj,self._foreground_obj,self._solver)
        if (len(self._paramrange)):
            self._engine.rerange(self._paramrange)
        result = self._engine.run(kwargs)
        self._paramlist = sorted(self._engine.activelist)
        return result

    def analyse_hl(self, kwargs=dict()):
        # noise offset improved HL likelihood
        offset_bp = np.mean(self._noise_bp,axis=0)  # effecive offset (1503.01347, 2010.01139)
        for l in range(offset_bp.shape[1]):
            offset_bp[:,l,:,:] *= np.sqrt(self._estimator.modes[l]+0.5)
        self.engine = hlfit(self._data_bp,np.mean(self._fiducial_bp,axis=0),np.mean(self._noise_bp,axis=0),self._covmat,self._background_obj,self._foreground_obj,self._solver,offset_bp)
        if (len(self._paramrange)):
            self._engine.rerange(self._paramrange)
        result = self._engine.run(kwargs)
        self._paramlist = sorted(self._engine.activelist)
        return result

# end
