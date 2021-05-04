from __future__ import division

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt  
import timeit
import sys
import os

import xqml
from xqml.xqml_utils import progress_bar, getstokes
from xqml.simulation import muKarcmin2var
from xqml.simulation import extrapolpixwin

import pymaster as nmt

class QMLestimator(object):

    def __init__(self, clth, qml_nside, aposcale=None, noiseA=None, noiseB=None, pixwin=False, mask=None,  \
                         lmin=None, lmax=None, lbin=None, targets='EB',muKarcmin=0.1, fwhm=0.0 ):

        """
        Parameters
        ----------
        mask : numpy.ndarray
            A single-vector array of mask map.
        
        lmin : (positive) integer
            Minimal angular mode.
        
        lmax : (positive) integer
            Maximal angular mode.

        lbin : (positive) integer
            Angular mode bin size for NaMaster calculation.

        lcut : (positive) integer
            Ignoring the first and last number of bins from NaMaster results.
        
        targets : string
            Choosing among 'T', 'E', 'B', 'EB', 'TEB' mode.

        filt : dict
            filtering effect in BP, recording extra mixing/rescaling of BP

        clth:numpy.ndarray
            fudicial model
        
        qml_nside : integer
           nisde of  low multipoles maps

        

        pixwin : boolean, optional
            If True, applies pixel window function to spectra. Default: True
        """
    # # Initialise basic parameters ##############################################################                     
        self.lbin       = lbin
        self.qml_nside  = qml_nside
        self.lmin       = lmin
        self.lmax       = lmax
        self.aposcale   = aposcale
        self.mask       = mask
        self.clth       = clth
        self.fwhm       = fwhm
        self._pixwin     = pixwin
    # # Initialise targets-ralated parameters ######################################################            
        self.targets = targets
        self._modedict = {'T':['TT'],'E':['EE'],'B':['BB'],'EB':['EE','BB'],'TEB':['TT','EE','BB']}        
        self._stokes, self._spec, self._istokes, self._ispecs = getstokes(self._modedict[self._targets])
        self._nstokes = len(self._stokes)
        self._nspec   = len(self._spec)
        if self._nstokes == 1:
            self._pol = False
        else:
            self._pol = True
        toGB = 1024. * 1024. * 1024.
        self._emem = 8.*(self._qml_npix*2*self._qml_npix*2) * ( self._nspec*2 ) / toGB

    # # Initialise Modes ################################################################    
        self._modes=np.arange((3+self._lbin)/2.0, self._lmax, self._lbin)
        self._modes=self._modes[np.where(self._modes>self._lmin)]
        self._nmode=len(self._modes)
        print()
        print("nmode = %d"%self._nmode)
        print("modes = ", self._modes)
    # # Initialise Noise Covariance ################################################################    
        self.noiseA = noiseA
        self.noiseB = noiseB        
        if noiseA != None:
            self.qml_noiseA=degrade_maps(self._qml_nside, self.noiseA, self._pol, self._qml_mask)
            NA=get_maps_covariance(self.qml_noiseA, self._istokes, self._qml_mask)
        else:
            pixvar = muKarcmin2var(muKarcmin, self._qml_nside)
            varmap = np.ones((self._nstokes * self._qml_npix)) * pixvar
            NA = np.diag(varmap)
        if noiseB != None:
            self.qml_noiseB=degrade_maps(self._qml_nside, self.noiseB, self._pol, self._qml_mask)
            NB=get_maps_covariance(self.qml_noiseB, self._istokes, self._qml_mask)
        else:
            NB=NA

        self.estimator=xqml.xQML(self._qml_mask,self._ellbins,self.clth,NA=NA,NB=NB,lmax=3*self._qml_nside,spec=self._spec,pixwin=self._pixwin,bell=self._bell,fwhm=0.0)
    

    @property
    def lbin(self):
        return self._lbin

    @property
    def qml_nside(self):
        return self._qml_nside

    @property
    def aposcale(self):
        return self._aposcale

    @property
    def mask(self):
        return self._mask

    @property
    def qml_mask(self):
        return self._qml_mask

    @property
    def lmin(self):
        return self._lmin

    @property
    def lmax(self):
        return self._lmax

    @property
    def fwhm(self):
        return self._fwhm

    @property
    def modes(self):
        return self._modes
        
    @property
    def nmode(self):
        return self._nmode

    @property
    def pol(self):
        return self._pol

    @property
    def targets(self):
        return self._targets

    @property
    def modedict(self):
        return self._modedict

    @lbin.setter
    def lbin(self, lbin):
        if lbin is None:
            self._lbin = 5
        else:
            assert isinstance(lbin, int)
            assert (lbin > 0)
            self._lbin = lbin

    @qml_nside.setter
    def qml_nside(self, qml_nside):
        assert isinstance(qml_nside, int)
        assert (qml_nside > 0)
        self._qml_nside = qml_nside
        self._ellbins = np.arange(2, 3*self._qml_nside + 2, self._lbin)
        self._ellbins[-1] = 3*self._qml_nside+1
        print("ellbins  = ",self._ellbins)
    
    @aposcale.setter
    def aposcale(self, aposcale):
        if aposcale is None:
            self._aposcale = 0.0
        else:
            assert (aposcale > 0)
            self._aposcale = aposcale    

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = np.ones(hp.nsdie2pix(self._qml_nside),dtype=np.bool)
        else:
            assert isinstance(mask, np.ndarray)
            self._mask = mask.copy()
            self._mask = nmt.mask_apodization(self._mask, self._aposcale, apotype='C2')
            self._qml_mask=degrade_mask(self._qml_nside, self._mask)
            self._fsky=np.mean(self._qml_mask)
            self._qml_npix = sum(self._qml_mask)
    

    @lmin.setter
    def lmin(self, lmin):
        if lmin is None:
            self._lmin = 2
        else:
            assert isinstance(lmin, int)
            assert (lmin <= 3*self._qml_nside)
            self._lmin = lmin

    @lmax.setter
    def lmax(self, lmax):
        if lmax is None:
            self._lmax = 2*self._qml_nside
        else:
            assert isinstance(lmax, (int,np.int64))
            if lmax <= 2*self._qml_nside:
                self._lmax = lmax
            else:
                print("lmax has been changed to 2*qml_nside")
                self._lmax = 2*self._qml_nside

    @fwhm.setter
    def fwhm(self, fwhm):
        if fwhm == 0.0:
            self._fwhm = fwhm
            bl=np.ones(3*self._qml_nside+1)
            l=np.arange(self._qml_nside+1,3*self._qml_nside+1)
            bl[self._qml_nside+1:3*self._qml_nside+1]=0.5*(1+np.sin(l*np.pi/2/self._qml_nside))
            self._bell = bl
        else:
            self._fwhm = fwhm
            bl=np.ones(3*self._qml_nside+1)
            l=np.arange(self._qml_nside+1,3*self._qml_nside+1)
            bl[self._qml_nside+1:3*self._qml_nside+1]=0.5*(1+np.sin(l*np.pi/2/self._qml_nside))
            self._bell = bl*hp.gauss_beam(np.deg2rad(self._fwhm), lmax=3*self._qml_nside)

    @targets.setter
    def targets(self, targets):
        assert isinstance(targets, str)
        self._targets = targets

        
        """"没有用到的变量"""
        #self.filt = filt
        #self.lcut = lcut
        #self.aposcale=None
    

    def print_memory_fsky(self):
        print()
        print("patch: fsky=%.2g %% (npix=%d)" % (100*self._fsky,self._qml_npix))
        print("mem=%.2g Gb" % self._emem)

    def get_stokes(spec):
        """
        Get the Stokes parameters number and name(s)

        Parameters
        ----------
        spec : bool
            If True, get Stokes parameters for polar (default: True)

        Returns
        ----------
        stokes : list of string
            Stokes variables names
        spec : int
            Spectra names
        istokes : list
            Indexes of power spectra

        Example
        ----------
        >>> getstokes(['EE','BB'])
        (['Q', 'U'], ['EE', 'BB'], [1, 2])
        >>> getstokes(['TT','EE','BB','TE'])
        (['I', 'Q', 'U'], ['TT', 'EE', 'BB', 'TE'], [0, 1, 2, 3])
        >>> getstokes(['TT', 'EE', 'BB', 'TE', 'EB', 'TB'])
        (['I', 'Q', 'U'], ['TT', 'EE', 'BB', 'TE', 'EB', 'TB'], [0, 1, 2, 3, 4, 5])
        """
        _temp  = "TT" in spec or "TE" in spec or "TB" in spec
        _polar = "EE" in spec or "BB" in spec or "TE" in spec or "TB" in spec or "EB" in spec
        _corr  = "TE" in spec or "TB" in spec or "EB" in spec
        if not _temp and not _polar and not _corr:
            print("invalid spectra list and/or no options")
        
        stokes = []
        if _temp:
            stokes.extend(["I"])
        if _polar:
            stokes.extend(["Q", "U"])
        
        ispecs = [['TT', 'EE', 'BB', 'TE', 'TB', 'EB'].index(s) for s in spec]
        istokes = [['I', 'Q', 'U'].index(s) for s in stokes]

        return stokes, spec, istokes, ispecs  


    def get_spectra(self, mapsA, mapsB=()):
        """
        Return the unbiased spectra

        Parameters
        ----------
        map1 : 1D array
            Pixel map number 1
        map2 : 2D array
            Pixel map number 2

        Returns
        ----------
        cl : array or sequence of arrays
            Returns cl or a list of cl's (TT, EE, BB, TE, EB, TB)
        
        # Should compute auto-spectra if map2 == None ?
        # Define conditions based on the map size
        """
        print()
        print("get spectrum from downgraded maps. ")
        allcl = []
        nsimu=min(mapsA.shape[0],mapsB.shape[0]) if mapsB != () else mapsA.shape[0]
        if mapsB != ():
            for n in np.arange(nsimu):
                mapA=mapsA[n]
                mapB=mapsB[n]
                allcl.append(self.estimator.get_spectra(mapA,mapB))
        else:
            for n in np.arange(nsimu):
                mapA=mapsA[n]
                allcl.append(self.estimator.get_spectra(mapA))

        hcl = np.mean(allcl, 0)

        return hcl
    def autoBP(self, mapsA):
        """
        Return the unbiased spectra

        Parameters
        ----------
        map1 : 1D array
            Pixel map number 1
        map2 : 2D array
            Pixel map number 2

        Returns
        ----------
        cl : array or sequence of arrays
            Returns cl or a list of cl's (TT, EE, BB, TE, EB, TB)
        
        # Should compute auto-spectra if map2 == None ?
        # Define conditions based on the map size
        """
        print()
        print("get spectrum from downgraded maps. ")
        allcl = []
        nsimu = mapsA.shape[0] 
        
        for n in np.arange(nsimu):
            mapA=mapsA[n]
            allcl.append(self.estimator.get_spectra(mapA))

        hcl = np.mean(allcl, 0)

        return hcl

    def crosBP(self, mapsA, mapsB=()):
        """
        Return the unbiased spectra

        Parameters
        ----------
        map1 : 1D array
            Pixel map number 1
        map2 : 2D array
            Pixel map number 2

        Returns
        ----------
        cl : array or sequence of arrays
            Returns cl or a list of cl's (TT, EE, BB, TE, EB, TB)
        
        # Should compute auto-spectra if map2 == None ?
        # Define conditions based on the map size
        """
        print()
        print("get spectrum from downgraded maps. ")
        allcl = []
        nsimu=min(mapsA.shape[0],mapsB.shape[0]) if mapsB != () else mapsA.shape[0]
        if mapsB != ():
            for n in np.arange(nsimu):
                mapA=mapsA[n]
                mapB=mapsB[n]
                allcl.append(self.estimator.get_spectra(mapA,mapB))
        else:
            print("crosBP need another map.")

        hcl = np.mean(allcl, 0)

        return hcl


def degrade_maps(qml_nside, maps_in, pol, qml_mask=None):
    """cut high multipoles of input map and degrade it to low resolution
    Parameters
    ----------
    nsimu：integer
        mumber of samples

    qml_nside : integer
       nisde of  low multipoles maps

    mapsA, mapsB: numpy.ndarray
             
        mapsA noise maps of chanel A
        mapsB noise maps of chanel B（optional）

    Returns
    ----------
    maps:
        Return a degrade maps with Nside = qml_nside
    """
    print()
    print("down grade maps to Nside = %d"%qml_nside)
    print(maps_in.shape)
    if len(maps_in.shape) ==2:
        nsimu=1
        maps_in=np.array([maps_in])
    else:
        nsimu=maps_in.shape[0]
    nside=hp.npix2nside(maps_in.shape[-1])
    Slmax_old=3*nside
    Slmax_new=3*qml_nside
    print("nsimu = %d"%nsimu)
    print("nside = %d"%nside)

    npix= hp.nside2npix(qml_nside)
    if qml_mask == ():
        qml_mask=np.ones(npix,bool)
    else:
        if npix != len(qml_mask):
            print("The nside of qml_mask inconsistent with qml_nside.")

    ALM = hp.Alm
    maps=[]
    for n in np.arange(nsimu):
        progress_bar(n, nsimu)
        n1=n+1
        #map_out = hp.read_map(maps_dir+"%d.fits"%n1,field=(0,1,2),dtype=float,verbose=0)
        #map_in = hp.read_map(maps_dir+"%d.fits"%n1,field=(0,1,2),dtype=np.float64,verbose=0)
        map_in=maps_in[n]
        
        alm_old=hp.map2alm(map_in,lmax=Slmax_old,pol=pol,use_weights=True,iter=3)
        alm_new=np.zeros_like(np.array(alm_old)[:,:ALM.getsize(Slmax_new)])
        for l in np.arange(Slmax_new+1) :
            for m in np.arange(l+1) :
                idx_new = ALM.getidx(Slmax_new, l, m)
                idx_old = ALM.getidx(Slmax_old, l, m)
                if l<=qml_nside:
                    alm_new[:,idx_new]=alm_old[:,idx_old]
                elif l<=Slmax_new:
                    alm_new[:,idx_new]=alm_old[:,idx_old]*0.5*(1+np.sin(l*np.pi/2/qml_nside))
        map_out = hp.alm2map(alm_new, nside= qml_nside, pixwin=False) 
        map_out = map_out*qml_mask


        # hp.write_map(output_dir+"/map_out_%d.fits"%n1,map_out,dtype=np.float64)
        # hp.mollview(map_out[1], cmap=plt.cm.jet, title='Q map output')
        # plt.savefig(output_dir+'/map_out_%d.eps'%n1,bbox_inches='tight',pad_inches=0.1)
        
        maps.append(map_out)


    maps=np.array(maps)
    return maps


def degrade_mask(qml_nside,  mask_in):
    """down grade mask to low resolution
    Parameters
    ----------
    qml_nside : integer
        nisde of  low multipoles maps

    mapsA, mapsB: numpy.ndarray
        mapsA noise maps of chanel A
        mapsB noise maps of chanel B（optional）

    Returns
    ----------
    map:
        Return a degrade map with Nside = qml_nside
    """
    print()
    print("degrade mask to Nside = %d"%qml_nside)
    mask_out = hp.ud_grade(mask_in,qml_nside)
    index=np.where(mask_out<0.999)
    mask_out[index]=0
    index=np.where(mask_out>=0.999)
    mask_out[index]=1
    #hp.write_map("qml_mask.fits",mask_out,dtype=bool)
    mask_out=np.asarray(mask_out,bool)
    return mask_out

def get_maps_covariance (maps, istokes, qml_mask):
    """get maps pixel covariance
    Parameters
    ----------
    nsimu：integer
        mumber of samples
        
    maps：
        degraded maps

    mode:

    targets:

    Returns
    ----------
    map:
        Return a degrade map with Nside = qml_nside
    """
    print()
    print("get pixel covariance of maps.")
    nsimu = maps.shape[0]
    masked_maps=maps[:,istokes][:,:, qml_mask]
    print("nsimu = %d"%nsimu)
    print("maps.shape = ",maps.shape)
    print("masked_maps.shape = ",masked_maps.shape)
    reshaped_masked_maps = masked_maps.reshape(masked_maps.shape[0],masked_maps.shape[1]*masked_maps.shape[2])
    print("reshaped_masked_maps.shape = ", reshaped_masked_maps.shape)
    rank=reshaped_masked_maps.shape[1]
    NoiseVar=np.cov(np.array(reshaped_masked_maps).T)

    for i in np.arange(rank):
         for j in np.arange(rank):
            if j<i-2 :
                NoiseVar[i,j]=0
               
            else:
                NoiseVar[i,j]=NoiseVar[i,j]

    print("NoiseVar.shape = ", NoiseVar.shape)
    
    return NoiseVar
    
def read_maps(nsimu,  maps_dir, pol=True):
    """cut high multipoles of input map and degrade it to low resolution
    Parameters
    ----------
    nsimu：integer
        mumber of samples

    qml_nside : integer
       nisde of  low multipoles maps

    mapsA, mapsB: numpy.ndarray
             
        mapsA noise maps of chanel A
        mapsB noise maps of chanel B（optional）

    Returns
    ----------
    maps:
        Return a degrade maps with Nside = qml_nside
    """
    print()
    print("Read maps from %s"%maps_dir)

    maps=[]
    for n in np.arange(nsimu):
        progress_bar(n, nsimu)
        n1=n+1
        if pol == True:
            map_in = hp.read_map(maps_dir+"%d.fits"%n1,field=(0,1,2),dtype=float,verbose=0)
        else:
            map_in = hp.read_map(maps_dir+"%d.fits"%n1,field=0,dtype=float,verbose=0)
        maps.append(map_in)
    maps=np.array(maps)
    return maps

    
####################### Start test #########################
print("xQML pipeline  test")
############################################################

#################### Make output dir #######################
output_dir = "output"
os.system("rm -rf "+output_dir)
os.system("mkdir "+output_dir)
nsimu = 10
############################################################



##################### Basic parameters #####################
targets='EB'
nside=512
qml_nside=32
mapsA_dir = "/home/jm/data/Data_common/TQUmaps/r0_05_lmax1024/cmb"
lmax = 3* qml_nside
lbin=11
aposcale=2
#not used parameters

#noiseA_dir =  "/home/jm/data/QML/fits2cl_QML/test1/degrade_map2/ilc7_210_16_"
############################################################



####################### Input model #########################
MODELFILE =  "/home/jm/data/Data_common/Input_cls/PCP18_r0.05wl_K2.fits"
clth = np.array(hp.read_cl(MODELFILE))
clth = np.array( list(clth) + list(clth[0:2]*0.))
############################################################



####################### Input mask #########################
mask=hp.read_map( "/home/jm/data/Data_common/Input_mask/AliCPT/I_Noise_95_G_512_6ukarcmin.fits",field=0,dtype=bool,verbose=0)
############################################################



###################### Make bell ###########################
# bl=np.ones(lmax+1)
# l=np.arange(qml_nside+1,3*qml_nside+1)
# bl[qml_nside+1:3*qml_nside+1]=0.5*(1+np.sin(l*np.pi/2/qml_nside))
# bell = bl
############################################################

################### Initialise xqml class ##################
xqml_estimator=QMLestimator(clth=clth, qml_nside=qml_nside,noiseA=None, mask=mask, aposcale=aposcale, \
                            lmax=lmax, lbin=lbin, targets=targets,  fwhm=0.0)
xqml_estimator.print_memory_fsky()
############################################################


####################### Construct Ps #######################
# start = timeit.default_timer()
# print(xqml_estimator.nsimu,xqml_estimator.nside, xqml_estimator.qml_nside, mapsA_dir)
# print(xqml_estimator.pol, xqml_estimator.qml_mask)
mapsA_in=read_maps(nsimu, mapsA_dir)
mapsA=degrade_maps(xqml_estimator.qml_nside, mapsA_in, xqml_estimator.pol, xqml_estimator.qml_mask)
# print()
# print("mapsA.shape = ",mapsA.shape)
#cl=xqml_estimator.get_spectra(mapsA)
cl=xqml_estimator.autoBP(mapsA)
np_cl=np.array(cl)
np.savetxt(output_dir+"/cl.txt" ,np.transpose(np_cl) ,'%10.6e')
# s1 = timeit.default_timer()
# print( "Construct Ps: %d sec" % (s1-start))
############################################################

# def preprocess(self, aposcale, psbin, lmin=None, lmax=None):
#         """
#         preprocess routine, converts maps into band-powers.

#         Parameters
#         ----------

#         aposcale : float
#             Apodization scale.

#         psbin : integer
#             Number of angular modes in each bin,
#             for conducting pseudo-PS estimation.

#         lmin/lmax : integer
#             Lower/Upper multipole limit.
#         """
#         assert isinstance(aposcale, float)
#         assert isinstance(psbin, int)
#         assert (psbin > 0)
#         assert (aposcale > 0)



#         # STEP I
#         # init PS estimator
#         # self.estimator = pstimator(nside=self._nside,mask=self._mask,aposcale=aposcale,psbin=psbin,lmin=lmin,lmax=lmax,targets=self._targets,filt=self._filt)
#         # self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#         #                     lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
        

#         # STEP II
#         # template PS estimations (auto corr. only)
#         if self._template_flag:
#             self.template_bp = dict()
#             for i in range(self._template_nfreq):
                
#                 # allocate for template
#                 data_bp = np.zeros((self._ntarget,self._estimator.nmode),dtype=np.float64)
#                 noise_bp = np.zeros((self._noise_nsamp,self._ntarget,self._estimator.nmode),dtype=np.float64)
#                 _fi = self._template_freqlist[i]
                

#                 # template workspace
#                 #twsp = self._estimator.autoWSP(self._templates[_fi],beams=self._template_beams[_fi])

#                 #init template PS estimator 
#                 self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)

#                 # template auto-corr.
#                 #stmp = self._estimator.autoBP(self._templates[_fi],wsp=twsp,beams=self._template_beams[_fi])
#                 stmp = self._estimator.autoBP(self._templates[_fi])
#                 data_bp = np.array(stmp[1:1+self._ntarget])
                
#                 # template noise auto corr.
#                 # for s in range(self._template_nsamp):
#                 #     ntmp = self._estimator.autoBP(self._template_noises[_fi][s],wsp=twsp,beams=self._template_beams[_fi])
#                 #     noise_bp[s] = np.array(ntmp[1:1+self._ntarget])
                
#                 self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)

#                 ntmp = self._estimator.autoBP(self._template_noises[_fi][s])
#                 noise_bp[s] = np.array(ntmp[1:1+self._ntarget])
#                 # mean noise subtraction
#                 self._template_bp[_fi] = data_bp - np.mean(noise_bp,axis=0)


#         # STEP III
#         # prepare model, parameter list generated during init models
#         if self._background is not None:
#             self._background_obj = self._background(self._freqlist,self._estimator)
#         if self._foreground is not None:
#             self._foreground_obj = self._foreground(self._freqlist,self._estimator,self._template_bp)
        
#         # STEP IV-A
#         # data PS estimations (with workspace)
#         # allocate
#         wsp_dict = dict()  # data nmt workspace
#         fwsp_dict = dict()  # fiducial nmt workspace

#         self.data_bp = np.zeros((self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
#         for i in range(self._nfreq):
#             _fi = self._freqlist[i]
#             # auto corr.
#             #wsp_dict[(i,i)] = self._estimator.autoWSP(self._data[_fi],beams=self._beams[_fi])
#             #stmp = self._estimator.autoBP(self._data[_fi],wsp=wsp_dict[(i,i)],beams=self._beams[_fi])
#             self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#             stmp = self._estimator.autoBP(self._data[_fi])
#             self._data_bp[:,:,i,i] = np.array(stmp[1:1+self._ntarget])
#             for j in range(i+1,self._nfreq):
#                 _fj = self._freqlist[j]
#                 # cross corr.
#                 #wsp_dict[(i,j)] = self._estimator.crosWSP(np.r_[self._data[_fi],self._data[_fj]],beams=[self._beams[_fi],self._beams[_fj]])
#                 #stmp = self._estimator.crosBP(np.r_[self._data[_fi],self._data[_fj]],wsp=wsp_dict[(i,j)],beams=[self._beams[_fi],self._beams[_fj]])
#                 self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#                 stmp = self._estimator.crosBP(np.r_[self._data[_fi],self._data[_fj]])
#                 self._data_bp[:,:,i,j] = np.array(stmp[1:1+self._ntarget])
#                 self._data_bp[:,:,j,i] = np.array(stmp[1:1+self._ntarget])
#         # STEP IV-B
#         # fiducial PS estimation 
#         if self._fiducial_flag:
#             # allocate
#             self.fiducial_bp = np.zeros((self._fiducial_nsamp,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
#             for i in range(self._nfreq):
#                 _fi = self._freqlist[i]
#                 # auto corr.
#                 #fwsp_dict[(i,i)] = self._estimator.autoWSP(self._fiducials[_fi][0],beams=self._fiducial_beams[_fi])
#                 self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#                 for s in range(self._fiducial_nsamp):
#                     # auto corr.
#                     #ftmp = self._estimator.autoBP(self._fiducials[_fi][s],wsp=fwsp_dict[(i,i)],beams=self._fiducial_beams[_fi])
#                     ftmp = self._estimator.autoBP(self._fiducials[_fi][s])
#                     self._fiducial_bp[s,:,:,i,i] = np.array(ftmp[1:1+self._ntarget])
#                 for j in range(i+1,self._nfreq):
#                     _fj = self._freqlist[j]
#                     #fwsp_dict[(i,j)] = self._estimator.crosWSP(np.r_[self._fiducials[_fi][0],self._fiducials[_fj][0]],beams=[self._fiducial_beams[_fi],self._fiducial_beams[_fj]])
#                     self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#                     for s in range(self._fiducial_nsamp):
#                         # cross corr.
#                         #ftmp = self._estimator.crosBP(np.r_[self._fiducials[_fi][s],self._fiducials[_fj][s]],wsp=fwsp_dict[(i,j)],beams=[self._fiducial_beams[_fi],self._fiducial_beams[_fj]])
#                         ftmp = self._estimator.crosBP(np.r_[self._fiducials[_fi][s],self._fiducials[_fj][s]])
#                         self._fiducial_bp[s,:,:,i,j] = np.array(ftmp[1:1+self._ntarget])
#                         self._fiducial_bp[s,:,:,j,i] = np.array(ftmp[1:1+self._ntarget])
#         # STEP IV-C
#         # noise PS estimations
#         if self._noise_flag:
#             # allocate
#             self.noise_bp = np.zeros((self._noise_nsamp,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
#             for s in range(self._noise_nsamp):
#                 for i in range(self._nfreq):
#                     _fi = self._freqlist[i]
#                     # auto corr.
#                     #ntmp = self._estimator.autoBP(self._noises[_fi][s],wsp=wsp_dict[(i,i)],beams=self._beams[_fi])
#                     self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#                     ntmp = self._estimator.autoBP(self._noises[_fi][s])
#                     self._noise_bp[s,:,:,i,i] = np.array(ntmp[1:1+self._ntarget])
#                     for j in range(i+1,self._nfreq):
#                         _fj = self._freqlist[j]
#                         # cross corr.
#                         #ntmp = self._estimator.crosBP(np.r_[self._noises[_fi][s],self._noises[_fj][s]],wsp=wsp_dict[(i,j)],beams=[self._beams[_fi],self._beams[_fj]])
#                         self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#                         ntmp = self._estimator.crosBP(np.r_[self._noises[_fi][s],self._noises[_fj][s]],wsp=wsp_dict[(i,j)])
#                         self._noise_bp[s,:,:,i,j] = np.array(ntmp[1:1+self._ntarget])
#                         self._noise_bp[s,:,:,j,i] = np.array(ntmp[1:1+self._ntarget])
#         # STEP V
#         # fiducial+noise PS covariance matrix
#         # fiducial+noise has to be processed in the pixel doamin, in order to yield a proper cov matrix
#         if self._fiducial_flag and self._noise_flag:
#             ncom = min(self._noise_nsamp,self._fiducial_nsamp)
#             # allocate
#             nfid = np.zeros((ncom,self._ntarget,self._estimator.nmode,self._nfreq,self._nfreq),dtype=np.float64)
#             for s in range(ncom):
#                 for i in range(self._nfreq):
#                     _fi = self._freqlist[i]
#                     # auto corr.
#                     #ntmp = self._estimator.autoBP(self._fiducials[_fi][s]+self._noises[_fi][s],wsp=fwsp_dict[(i,i)],beams=self._fiducial_beams[_fi])
#                     self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#                     ntmp = self._estimator.autoBP(self._fiducials[_fi][s]+self._noises[_fi][s])
#                     nfid[s,:,:,i,i] = np.array(ntmp[1:1+self._ntarget])
#                     for j in range(i+1,self._nfreq):
#                         _fj = self._freqlist[j]
#                         # cross corr.
#                         #ntmp = self._estimator.crosBP(np.r_[self._fiducials[_fi][s]+self._noises[_fi][s],self._fiducials[_fj][s]+self._noises[_fj][s]],wsp=fwsp_dict[(i,j)],beams=[self._fiducial_beams[_fi],self._fiducial_beams[_fj]])
#                         self.estimator = QMLestimator(nsimu=self._noise_nsamp,clth=qml_clth, nside=self._nside, qml_nside=qml_nside,noiseA=None, mask=self._mask, aposcale=aposcale\
#                             lmin=qml_lmin, lmax=qml_lmax, lbin=qml_lbin, targets=self._targets, pixwin=qml_pixwin, bell=qml_bell)
#                         ntmp = self._estimator.crosBP(np.r_[self._fiducials[_fi][s]+self._noises[_fi][s],self._fiducials[_fj][s]+self._noises[_fj][s]])
#                         nfid[s,:,:,i,j] = np.array(ntmp[1:1+self._ntarget])
#                         nfid[s,:,:,j,i] = np.array(ntmp[1:1+self._ntarget])
#             # full cov
#             self.covmat = empcov(gvec(nfid),self._ntarget)

