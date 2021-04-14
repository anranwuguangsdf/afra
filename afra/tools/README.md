## AFRA documentation

### tools

The "tools" directory holds three modules: "aux", "ps\_estimator" and "bp\_vis".
Except the "ps\_estimator" which defines the `pstimator` class for band-power estimation, the other two consist of utility functions.

> "aux" module function list:

| function name | description |
|:--------------|:------------|
| **empcov** | empirical covariance matrix calculator |
| **umap** | unitity mapping function (designed for prior) |
| **gvec** | band-power matrix vectorization for Gaussian likelihood |
| **hvec** | band-power matrix vectorization for HL likelihood |

> "bp\_vis" module function list:

| function name | description |
|:--------------|:------------|
| **bp\_vis** | visualize **data**, **CMB fiducial**, **noise** and best-fit CMB band-powers |

> class `pstimator` attribute list:

| attribute name | description |
|:---------------|:------------|
| **nside** | HEALPix Nside of maps |
| **npix** | HEALPix Npix of maps |
| **mask** | mask map |
| **apomask** | apodized mask map  |
| **b** | NaMaster binning instance |
| **aposcale** | NaMaster apodization scale |
| **lmin** | lower limit of multipole range |
| **lmax** | upper limit of multipole range |
| **lbin** | NaMaster multipole binning width |
| **lcut** | number of NaMaster bins discarded at lower and upper boundaries |
| **psbin** | number of rebinning bins |
| **modes** | list of final multipole bins' central position |
| **nmode** | number of final multipole bins |
| **targets** | band-power estimation targets tuple |
| **ntarget** | number of targets |
| **filt** | CMB hamonic-domain filtering matrix |

> class `pstimator` function list:

| function name | description |
|:--------------|:------------|
| **rebinning** | re-bin the NaMaster band-power outputs |
| **bpconvert** | convert power-spectra into band-powers |
| **filtrans** | apply CMB filtering effect on band-powers |
| **autoWSP** | auto-correlation NaMaster workspace |
| **crosWSP** | cross-correlation NaMaster workspace |
| **autoBP** | calculate auto-correlation bandpower |
| **crosBP** | calculate cross-correlation bandpower |


