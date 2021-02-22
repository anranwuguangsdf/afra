# this is a bandpower visualization script
# designed for visual inspection for data pre-process outputs

import numpy as np
import matplotlib.pyplot as plt


def bpvis(targets, modes, freqs, data_bp, fiducial_bp=None, noise_bp=None, best_bp=None):
    """
    Bandpower visualization, with data, fiducials, noise and best-fit,

    targets: str
        T/E/B mode string

    modes: list/tuple/array
        multipole bin position

    freqs: list/tuple/array
        frequency list

    data_bp: array
        data bandpower matrix

    fiducial_bp: array
        fiducial bandpower matrix

    noise_bp: array
        noise bandpower matrix

    best_bp: array
        best-fit bandpower matrix
    """
    
    assert isinstance(targets, str)
    _ntype = len(targets)
    cc = ['red', 'blue', 'orange'][:_ntype]  # colors
    assert isinstance(modes, (list,tuple,np.ndarray))
    _nmode = len(modes)
    assert isinstance(freqs, (list,tuple,np.ndarray))
    _nfreq = len(freqs)
    assert isinstance(data_bp, np.ndarray)
    assert (_ntype == data_bp.shape[0])
    assert (_nmode == data_bp.shape[1])
    assert (_nfreq == data_bp.shape[2])
    assert (_nfreq == data_bp.shape[3])
    
    # visualize fiducial BPs
    if fiducial_bp is not None:
        assert isinstance(fiducial_bp, np.ndarray)
        assert (_ntype == fiducial_bp.shape[1])
        assert (_nmode == fiducial_bp.shape[2])
        assert (_nfreq == fiducial_bp.shape[3])
        assert (_nfreq == fiducial_bp.shape[4])
        fig = plt.figure(figsize=(5*_nfreq,5*_nfreq))
        for i in range(_nfreq):
            for j in range(i,_nfreq):
                ax = fig.add_subplot(_nfreq,_nfreq,1+j*_nfreq+i)
                for k in range(_ntype):
                    ax.errorbar(modes,np.mean(fiducial_bp[:,k,:,i,j],axis=0),yerr=np.std(fiducial_bp[:,k,:,i,j],axis=0),marker='o',mfc='none',mec=cc[k])
                ax.set_title(str(freqs[i])+'x'+str(freqs[j]))
                ax.set_yscale('log')
        labels = list()
        for t in targets:
            labels.append(t+' fiducial')
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('fiducial_bp.pdf')
    
    # visualize noise BPs
    if noise_bp is not None:
        assert isinstance(noise_bp, np.ndarray)
        assert (_ntype == noise_bp.shape[1])
        assert (_nmode == noise_bp.shape[2])
        assert (_nfreq == noise_bp.shape[3])
        assert (_nfreq == noise_bp.shape[4])
        fig = plt.figure(figsize=(5*_nfreq,5*_nfreq))
        for i in range(_nfreq):
            ax = fig.add_subplot(_nfreq,_nfreq,1+i*(_nfreq+1))
            for k in range(_ntype):
                ax.errorbar(modes,np.mean(noise_bp[:,k,:,i,i],axis=0),yerr=np.std(noise_bp[:,k,:,i,i],axis=0),marker='o',mfc='none',mec=cc[k])
            ax.set_title(str(freqs[i])+'x'+str(freqs[i]))
            ax.set_yscale('log')
            for j in range(i+1,_nfreq):
                ax = fig.add_subplot(_nfreq,_nfreq,1+j*_nfreq+i)
                for k in range(_ntype):
                    ax.errorbar(modes,np.mean(noise_bp[:,k,:,i,j],axis=0),yerr=np.std(noise_bp[:,k,:,i,j],axis=0),marker='o',mfc='none',mec=cc[k])
                ax.set_title(str(freqs[i])+'x'+str(freqs[j]))
        labels = list()
        for t in targets:
            labels.append(t+' noise')
        fig.legend([], labels=labels,loc='upper right',fontsize=15)
        plt.savefig('noise_bp.pdf')
    
    # visualize data BPs
    if best_bp is not None:
        assert (best_bp.shape == data_bp.shape)
    fig = plt.figure(figsize=(5*_nfreq,5*_nfreq))
    for i in range(_nfreq):
        ax = fig.add_subplot(_nfreq,_nfreq,1+i*(_nfreq+1))
        for k in range(_ntype):
            ax.scatter(modes,data_bp[k,:,i,i],color=cc[k],marker='.')
            if best_bp is not None:
                ax.plot(modes,best_bp[k,:,i,i],color='k')
            if fiducial_bp is not None and noise_bp is not None:
                smean = data_bp[k,:,i,i]-np.mean(noise_bp[:,k,:,i,i],axis=0)
                sstd = np.std(fiducial_bp[:,k,:,i,i],axis=0)+np.std(noise_bp[:,k,:,i,i],axis=0)
                ax.errorbar(modes,smean,yerr=sstd,marker='o',mfc='none',mec=cc[k])
        ax.set_title(str(freqs[i])+'x'+str(freqs[i]))
        ax.set_yscale('log')
        for j in range(i+1,_nfreq):
            ax = fig.add_subplot(_nfreq,_nfreq,1+j*_nfreq+i)
            for k in range(_ntype):
                ax.scatter(modes,data_bp[k,:,i,j],color=cc[k],marker='.')
                if best_bp is not None:
                    ax.plot(modes,best_bp[k,:,i,j],color='k')
                if fiducial_bp is not None and noise_bp is not None:
                    smean = data_bp[k,:,i,j]-np.mean(noise_bp[:,k,:,i,j],axis=0)
                    sstd = np.std(fiducial_bp[:,k,:,i,j],axis=0)+np.std(noise_bp[:,k,:,i,j],axis=0)
                    ax.errorbar(modes,smean,yerr=sstd,marker='o',mfc='none',mec=cc[k])
            ax.set_title(str(freqs[i])+'x'+str(freqs[j]))
    labels = list()
    if best_bp is not None:
        for t in targets:
            labels.append(t+' bestfit')
    for t in targets:
        labels.append(t+' data')
    if fiducial_bp is not None and noise_bp is not None:
        for t in targets:
            labels.append(t+' signal')
    fig.legend([], labels=labels,loc='upper right',fontsize=15)
    plt.savefig('data_bp.pdf')
