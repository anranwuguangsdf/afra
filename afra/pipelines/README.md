## AFRA documentation

### pipelines

The "pipelines" directory holds three modules: `pipeline`, `abs_pipeline`, `fit_pipeline`.

The pipeline module defines the base pipeline class: `pipe`, which defines the in/out-puts and two basic functions: `preprocess` and `reprocess`.
The derived pipeline classes inherit base class' methods and create their own specialized methods.

> base class `pipe` **attribute** list:

| attribute name | description |
|:---------------|:------------|
| **data** | multi-frequency measured TQU maps |
| **data\_bp** | estimated band-power matrix for measured maps |
| **noises** | multi-frequency noise TQU maps |
| **noise\_bp** | estimated band-power matrix for noise maps |
| **noise\_flag** | if have noise maps |
| **noise\_nsamp** | noise maps sample size at each frequency |
| **mask** | mask map |
| **freqlist** | list of frequencies |
| **nside** | HEALPix Nside of maps |
| **npix** | HEALPix Npix of maps |
| **targets** | band-power estimation targets tuple |
| **ntarget** | number of targets |
| **beams** | FWHM of measured maps' beams |
| **filt** | CMB hamonic-domain filtering matrix |
| **fiducals** | simulated multi-frequency fiducial CMB maps |
| **fiducial\_nsamp** | fiducial map sample size at each frequency | 
| **fiducial\_flag** | if have fiducial maps |
| **fiducial\_beams** | FWHM of fiducial maps' beams |
| **fiducial\_bp** | estimated band-power matrix for fiducial maps |
| **templates** | foreground template TQU maps |
| **template\_flag** | if have template maps |
| **template\_noises** | foreground template noise TQU maps |
| **template\_nsamp** | template noise maps sample size at each frequency |
| **template\_beams** | FHWM of template maps' beams |
| **background** | background model name |
| **background\_catalog** | background model name list |
| **background\_obj** | background model instance |
| **foreground** | foreground model name |
| **foreground\_catalog** | foreground model name list |
| **foreground\_obj** | foreground model instance |
| **likelihood** | likelihood type name |
| **paramlist** | parameter name list |
| **paramrange** | parameter sampling ranges |
| **estimator** | band-power estimatior instance |
| **covmat** | covariance matrix for CMB fiducial+noise band-power matrix |

> base class `pipe` **function** list:

| function name | description |
|:--------------|:------------|
| **preprocess** | cook maps into band-power matrices |
| **reprocess** | update output with new **data** |

> derived class `abspipe` **attribute** list:

| attribute name | description |
|:---------------|:------------|
| **engine** | Bayesian analysis engine instance |
| **_anadict** | analysis function dictionary |
| **absrslt** | ABS result |
| **absinfo** | extra information |

> derived class `abspipe` **function** list:

| function name | description |
|:--------------|:------------|
| **run_absonly** | run the ABS only |
| **run** | run preprocess, followed by ABS and then Bayesian analysis for CMB parameters |
| **analyse** | carry out ABS |
| **analyse_quiet** | noise-free ABS method |
| **analyse_noisy** | ABS method with noise |
| **postprocess** | post-fitting for CMB parameters |

> derived class `fitpipe` **attribute** list:

| attribute name | description |
|:---------------|:------------|
| **engine** | Bayesian analysis engine instance |
| **_anadict** | analysis function dictionary |

> derived class `fitpipe` **function** list:

| function name | description |
|:--------------|:------------|
| **run** | run preprocess, followed by Bayesian analysis for model parameters |
| **analyse** | carry out Bayesian analysis |
| **analyse_gauss** | Bayesian analysis with Gauss likelihood |
| **analyse_hl** | Bayesian analysis with HL likelihood |
