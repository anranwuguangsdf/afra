import numpy as np
import healpy as hp
import pymaster as nmt
from afra.tools.icy_decorator import icy


@icy
class pstimator(object):

    def __init__(self, nside, mask=None, aposcale=None, psbin=None, lmin=None, lmax=None, lbin=None, lcut=None, targets=('TT',), filt=None):
        """
        Parameters
        ----------
        nside : integer
            HEALPix Nside.
        
        mask : numpy.ndarray
            A single-vector array of mask map.
        
        aposcale : float
            Apodization size in deg.
        
        psbin : (positive) integer
            Number of PS bins.
        
        lmin : (positive) integer
            Minimal angular mode.
        
        lmax : (positive) integer
            Maximal angular mode.

        lbin : (positive) integer
            Angular mode bin size for NaMaster calculation.

        lcut : (positive) integer
            Ignoring the first and last number of bins from NaMaster results.

	targets : tuple
            Tuple consists of 'TT', 'TE', 'TB', 'EE', 'EB', 'BB'.

        filt : dict
            filtering effect in BP, recording extra mixing/rescaling of BP
        """
        self.nside = nside
        self.aposcale = aposcale
        self.lbin = None
        self.lcut = None
        self.lmin = lmin
        self.lmax = lmax
        self.mask = mask
        self.b = None
        self.psbin = psbin
        self.targets = targets
        self.filt = filt
        self._autowdict = {('TT',):self.autoWSP_TT,('EE',):self.autoWSP_EE,('BB',):self.autoWSP_BB,('EE','BB'):self.autoWSP_EnB}
        self._croswdict = {('TT',):self.crosWSP_TT,('EE',):self.crosWSP_EE,('BB',):self.crosWSP_BB,('EE','BB'):self.crosWSP_EnB}
        self._autodict = {('TT',):self.autoBP_TT,('EE',):self.autoBP_EE,('BB',):self.autoBP_BB,('EE','BB'):self.autoBP_EnB}
        self._crosdict = {('TT',):self.crosBP_TT,('EE',):self.crosBP_EE,('BB',):self.crosBP_BB,('EE','BB'):self.crosBP_EnB}

    @property
    def nside(self):
        return self._nside

    @property
    def npix(self):
        return self._npix

    @property
    def mask(self):
        return self._mask

    @property
    def apomask(self):
        return self._apomask

    @property
    def b(self):
        return self._b

    @property
    def aposcale(self):
        return self._aposcale

    @property
    def nmode(self):
        return self._nmode

    @property
    def lbin(self):
        return self._lbin

    @property
    def lcut(self):
        return self._lcut

    @property
    def psbin(self):
        return self._psbin

    @property
    def lmin(self):
        return self._lmin

    @property
    def lmax(self):
        return self._lmax

    @property
    def modes(self):
        return self._modes

    @property
    def targets(self):
        return self._targets

    @property
    def ntarget(self):
        return self._ntarget

    @property
    def filt(self):
        return self._filt

    @nside.setter
    def nside(self, nside):
        assert isinstance(nside, int)
        assert (nside > 0)
        self._nside = nside
        self._npix = 12*self._nside**2

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = np.ones(self._npix,dtype=np.float64)
            self._apomask = np.ones(self._npix,dtype=np.float64)
        else:
            assert isinstance(mask, np.ndarray)
            assert (len(mask) == self._npix)
            self._mask = mask.copy()
            self._apomask = nmt.mask_apodization(self._mask, self._aposcale, apotype='C2')

    @aposcale.setter
    def aposcale(self, aposcale):
        if aposcale is None:
            self._aposcale = 0.0
        else:
            assert (aposcale > 0)
            self._aposcale = aposcale

    @lbin.setter
    def lbin(self, lbin):
        "# of bin width of original bin"
        if lbin is None:
            self._lbin = 5
        else:
            assert isinstance(lbin, int)
            assert (lbin > 0)
            self._lbin = lbin

    @lcut.setter
    def lcut(self, lcut):
        "# of original bins on each side get trimed out of [lmin,lmax]"
        if lcut is None:
            self._lcut = 5
        else:
            assert isinstance(lcut, int)
            assert (lcut > 0)
            self._lcut = lcut

    @lmin.setter
    def lmin(self, lmin):
        "lower limit of multipole"
        if lmin is None:
            self._lmin = 2+25 # consider default lbin*lcut
        else:
            assert isinstance(lmin, int)
            assert (lmin < 3*self._nside)
            self._lmin = lmin

    @lmax.setter
    def lmax(self, lmax):
        "upper limit of multipole"
        if lmax is None:
            self._lmax = self._nside
        else:
            assert isinstance(lmax, (int,np.int64))
            assert (lmax < 3*self._nside)
            self._lmax = lmax
    
    @b.setter
    def b(self, b):
        if b is None:
            """customize NaMaster multipole band object"""
            ell_ini = np.arange((self._lmax-self._lmin+2*self._lbin*self._lcut)//self._lbin)*self._lbin + (self._lmin - self._lbin*self._lcut)
            assert (ell_ini[0] > 1)
            ell_end = ell_ini + self._lbin
            self._b = nmt.NmtBin.from_edges(ell_ini, ell_end, is_Dell=True)
            self.lmax = self._b.lmax  # correct lmax
        else:
            self._b = b

    @psbin.setter
    def psbin(self, psbin):
        if psbin is None:
            self._psbin = 1
            self.modes = self.rebinning(self._b.get_effective_ells())
        else:
            assert isinstance(psbin, int)
            assert (psbin > 0 and psbin <= (self._lmax-self._lmin)//self._lbin)
            self._psbin = psbin
            self.modes = self.rebinning(self._b.get_effective_ells())

    @modes.setter
    def modes(self, modes):
        assert isinstance(modes, (list,tuple,np.ndarray))
        self._modes = modes
        self._nmode = len(modes)

    @targets.setter
    def targets(self, targets):
        assert isinstance(targets, tuple)
        self._targets = targets
        self._ntarget = len(targets)

    @filt.setter
    def filt(self, filt):
        if filt is not None:
            assert isinstance(filt, dict)
            assert (self._targets in filt)
        self._filt = filt

    def rebinning(self, bp):
        bp_trim = bp[self._lcut:-self._lcut]
        bbp = np.empty(self._psbin, dtype=np.float64)
        idx_ini = np.arange(self._psbin)*(len(bp_trim)//self._psbin)
        idx_end = idx_ini + len(bp_trim)//self._psbin
        for i in range(self._psbin):
            bbp[i] = np.mean(bp_trim[idx_ini[i]:idx_end[i]])
        return bbp

    def bpconvert(self, ps):
        """
        "top-hat" window function matrix
        for converting PS into band-powers

        Parameters
        ----------
            input power-spectrum in multipole range (lmin:lmax+1)
        
        Return
        ----------
            band-power converting matrix in shape (# eff-ell)
        """
        raw_conv = self._b.bin_cell(ps)
        if (raw_conv.ndim == 1):
            return self.rebinning(raw_conv)
        else:
            fine_conv = np.empty((raw_conv.shape[0],self._psbin),dtype=np.float64)
            for i in range(fine_conv.shape[0]):
                fine_conv[i] = self.rebinning(raw_conv[i])
            return fine_conv

    def filtrans(self, bp):
        """
        apply filtering effect on band-powers.

        Parameters
        ----------
        
        bp : numpy.ndarray
            band-power in shape (# targets, # modes).

        Returns
        -------
            filtered band-power in shape (# targets, # modes).
        """
        if self._filt is None:
            return bp
        assert (bp.shape == (self._ntarget,self._nmode))
        transmat = self._filt[self._targets]
        assert (transmat.shape[0] == (bp.shape[0]*bp.shape[1]))
        return (transmat.dot(bp.reshape(-1,1))).reshape(self._ntarget,-1)

    def autoWSP(self, maps, beams=None):
        assert isinstance(maps, np.ndarray)
        assert (maps.shape == (3,self._npix))
        _cleaned = maps.copy()
        _cleaned[:,self._mask==0.] = 0.  # !!!
        return self._autowdict[self._targets](_cleaned,beams)

    def autoWSP_TT(self, maps, beams=None):
        # assemble NaMaster fields
        if beams is None:
            f0 = nmt.NmtField(self._apomask, [maps[0]])
        else:
            f0 = nmt.NmtField(self._apomask, [maps[0]], beam=hp.gauss_beam(beams, 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f0, f0, self._b)
        return w

    def autoWSP_EE(self, maps, beams=None):
        # assemble NaMaster fields
        if beams is None:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False)
        else:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False, beam=hp.gauss_beam(beams, 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f2, f2, self._b)
        return w

    def autoWSP_BB(self, maps, beams=None):
        # assemble NaMaster fields
        if beams is None:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True)
        else:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True, beam=hp.gauss_beam(beams, 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f2, f2, self._b)
        return w

    def autoWSP_EnB(self, maps, beams=None):
        # assemble NaMaster fields
        if beams is None:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True)
        else:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True, beam=hp.gauss_beam(beams, 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f2, f2, self._b)
        return w

    def crosWSP(self, maps, beams=[None,None]):
        assert isinstance(maps, np.ndarray)
        assert (maps.shape == (6,self._npix))
        assert (len(beams) == 2)
        _cleaned = maps.copy()
        _cleaned[:,self._mask==0.] = 0.  # !!!
        return self._croswdict[self._targets](_cleaned,beams)

    def crosWSP_TT(self, maps, beams=[None,None]):
        # assemble NaMaster fields
        if beams[0] is None:
            f01 = nmt.NmtField(self._apomask, [maps[0]])
        else:
            f01 = nmt.NmtField(self._apomask, [maps[0]], beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f02 = nmt.NmtField(self._apomask, [maps[3]])
        else:
            f02 = nmt.NmtField(self._apomask, [maps[3]], beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f01, f02, self._b)
        return w

    def crosWSP_EE(self, maps, beams=[None,None]):
        # assemble NaMaster fields
        if beams[0] is None:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False)
        else:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False, beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=False)
        else:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=False, beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f21, f22, self._b)
        return w

    def crosWSP_BB(self, maps, beams=[None,None]):
        # assemble NaMaster fields
        if beams[0] is None:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True)
        else:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True, beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=False, purify_b=True)
        else:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=False, purify_b=True, beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f21, f22, self._b)
        return w

    def crosWSP_EnB(self, maps, beams=[None,None]):
        # assemble NaMaster fields
        if beams[0] is None:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True)
        else:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True, beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=True)
        else:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=True, beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # prepare workspace
        w = nmt.NmtWorkspace()
        w.compute_coupling_matrix(f21, f22, self._b)
        return w

    def autoBP(self, maps, wsp=None, beams=None):
        """
        Auto BP
        
        Parameters
        ----------
        
        maps : numpy.ndarray
            A single-row array of a single map.
        
        wsp : (PS-estimator-defined) workspace
            A template of mask-induced mode coupling matrix.
        
        beams : float
            FWHM of gaussian beams

        Returns
        -------
        
        pseudo-PS results : tuple of numpy.ndarray
            (ell, XX, wsp(if input wsp is None))
        """
        assert isinstance(maps, np.ndarray)
        assert (maps.shape == (3,self._npix))
        _cleaned = maps.copy()
        _cleaned[:,self._mask==0.] = 0.  # !!!
        return self._autodict[self._targets](_cleaned,wsp,beams)

    def crosBP(self, maps, wsp=None, beams=[None,None]):
        """
        Cross BP
        
        Parameters
        ----------
        
        maps : numpy.ndarray
            A single-row array of a single map.
        
        wsp : (PS-estimator-defined) workspace
            A template of mask-induced mode coupling matrix.
        
        beams : float
            FWHM of gaussian beams.
        
        Returns
        -------
        
        pseudo-PS results : tuple of numpy.ndarray
            (ell, XX, wsp(if input wsp is None))
        """
        assert isinstance(maps, np.ndarray)
        assert (maps.shape == (6,self._npix))
        assert (len(beams) == 2)
        _cleaned = maps.copy()
        _cleaned[:,self._mask==0.] = 0.  # !!!
        return self._crosdict[self._targets](_cleaned,wsp,beams)

    def autoBP_TT(self, maps, wsp=None, beams=None):
        dat = maps[0]
        # assemble NaMaster fields
        if beams is None:
            f0 = nmt.NmtField(self._apomask, [dat])
        else:
            f0 = nmt.NmtField(self._apomask, [dat], beam=hp.gauss_beam(beams, 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl00 = nmt.compute_full_master(f0, f0, self._b)
            return (self._modes, self.rebinning(cl00[0]))
        else:
            cl00c = nmt.compute_coupled_cell(f0, f0)
            cl00 = wsp.decouple_cell(cl00c)
            return (self._modes, self.rebinning(cl00[0]))

    def crosBP_TT(self, maps, wsp=None, beams=[None,None]):
        dat1 = maps[0]
        dat2 = maps[3]
        # assemble NaMaster fields
        if beams[0] is None:
            f01 = nmt.NmtField(self._apomask, [dat1])
        else:
            f01 = nmt.NmtField(self._apomask, [dat1], beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f02 = nmt.NmtField(self._apomask, [dat2])
        else:
            f02 = nmt.NmtField(self._apomask, [dat2], beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl00 = nmt.compute_full_master(f01, f02, self._b)
            return (self._modes, self.rebinning(cl00[0]))
        else:
            cl00c = nmt.compute_coupled_cell(f01, f02)
            cl00 = wsp.decouple_cell(cl00c)
            return (self._modes, self.rebinning(cl00[0]))

    def autoBP_EE(self, maps, wsp=None, beams=None):
        # assemble NaMaster fields
        if beams is None:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False)
        else:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False, beam=hp.gauss_beam(beams, 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl22 = nmt.compute_full_master(f2, f2, self._b)
            return (self._modes, self.rebinning(cl22[0]))
        else:
            cl22c = nmt.compute_coupled_cell(f2, f2)
            cl22 = wsp.decouple_cell(cl22c)
            return (self._modes, self.rebinning(cl22[0]))

    def crosBP_EE(self, maps, wsp=None, beams=[None,None]):
        # assemble NaMaster fields
        if beams[0] is None:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False)
        else:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=False, beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=False)
        else:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=False, beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl22 = nmt.compute_full_master(f21, f22, self._b)
            return (self._modes, self.rebinning(cl22[0]))
        else:
            cl22c = nmt.compute_coupled_cell(f21, f22)
            cl22 = wsp.decouple_cell(cl22c)
            return (self._modes, self.rebinning(cl22[0]))

    def autoBP_BB(self, maps, wsp=None, beams=None):
        # assemble NaMaster fields
        if beams is None:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True)
        else:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True, beam=hp.gauss_beam(beams, 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl22 = nmt.compute_full_master(f2, f2, self._b)
            return (self._modes, self.rebinning(cl22[3]))
        else:
            cl22c = nmt.compute_coupled_cell(f2, f2)
            cl22 = wsp.decouple_cell(cl22c)
            return (self._modes, self.rebinning(cl22[3]))

    def crosBP_BB(self, maps, wsp=None, beams=[None,None]):
        # assemble NaMaster fields
        if beams[0] is None:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True)
        else:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=False, purify_b=True, beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=False, purify_b=True)
        else:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=False, purify_b=True, beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl22 = nmt.compute_full_master(f21, f22, self._b)
            return (self._modes, self.rebinning(cl22[3]))
        else:
            cl22c = nmt.compute_coupled_cell(f21, f22)
            cl22 = wsp.decouple_cell(cl22c)
            return (self._modes, self.rebinning(cl22[3]))

    def autoBP_EnB(self, maps, wsp=None, beams=None):
        # assemble NaMaster fields
        if beams is None:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True)
        else:
            f2 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True, beam=hp.gauss_beam(beams, 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl22 = nmt.compute_full_master(f2, f2, self._b)
            return (self._modes, self.rebinning(cl22[0]), self.rebinning(cl22[3]))
        else:
            cl22c = nmt.compute_coupled_cell(f2, f2)
            cl22 = wsp.decouple_cell(cl22c)
            return (self._modes, self.rebinning(cl22[0]), self.rebinning(cl22[3]))

    def crosBP_EnB(self, maps, wsp=None, beams=[None,None]):
        # assemble NaMaster fields
        if beams[0] is None:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True)
        else:
            f21 = nmt.NmtField(self._apomask, [maps[1], maps[2]], purify_e=True, purify_b=True, beam=hp.gauss_beam(beams[0], 3*self._nside-1))
        if beams[1] is None:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=True)
        else:
            f22 = nmt.NmtField(self._apomask, [maps[4], maps[5]], purify_e=True, purify_b=True, beam=hp.gauss_beam(beams[1], 3*self._nside-1))
        # estimate PS
        if wsp is None:
            cl22 = nmt.compute_full_master(f21, f22, self._b)
            return (self._modes, self.rebinning(cl22[0]), self.rebinning(cl22[3]))
        else:
            cl22c = nmt.compute_coupled_cell(f21, f22)
            cl22 = wsp.decouple_cell(cl22c)
            return (self._modes, self.rebinning(cl22[0]), self.rebinning(cl22[3]))
