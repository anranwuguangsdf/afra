## AFRA documentation

### methods

The "methods" directory holds two classes: `abssep` in "abs.py", `fit` in "fit.py".
For Gaussian and HL likelihood, we design two derived classes `gaussfit` and `hlfit` respectively.

> base class `abssep` **attribute** list:

| attribute name | description |
|:---------------|:------------|
| **data** | multi-frequency band-power matrix |
| **noise** | ensemble averaged multi-frequency noise band-power matrix |
| **noise\_flag** | if have **noise** |
| **sigma** | RMS of ensemble noise auto power-spectrum |
| **shift** | CMB band-power shift |
| **threshold** | signal-to-noise ratio threshold for eigen values |

> base class `abssep` **function** list:

| function name | description |
|:--------------|:------------|
| **run** | run ABS analysis on given **data** and **noise** |
| **run\_info** | run ABS analysis but return eigen values and vectors |

> base class `fit` **attribute** list:

| attribute name | description |
|:---------------|:------------|
| **data** | multi-frequency band-power matrix |
| **fiducial** | multi-frequency fiducial CMB band-power matrix |
| **noise** | ensemble averaged multi-frequency noise band-power matrix |
| **covariance** | covariance matrix for **fiducial** |
| **foreground** | CMB foreground model instance |
| **background** | CMB model instance |
| **params** | model parameter dictionary |
| **paramrange** | parameter sampling range dictionary |
| **activelist** | list of free (to be constrained) parameter names |

> base class `fit` **function** list:

| function name | description |
|:--------------|:------------|
| **rerange** | redefine parameters sampling range |
| **run** | run Bayesian analysis |
| **_core_likelihood** | calculate loglikelihood given sample position |
| **prior** | define Bayesian piror mapping |

> derived class `gaussfit` **function** list:

| function name | description |
|:--------------|:------------|
| **loglikeli** | loglikelihood calculator required by Bayesian sampler |

> derived class `hlfit` **function** list:

| function name | description |
|:--------------|:------------|
| **offset** | modified HL offset for **noise** |
| **loglikeli** | loglikelihood calculator required by Bayesian sampler |
