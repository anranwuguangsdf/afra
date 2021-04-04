## AFRA documentation

### models

The "methods" directory holds two classes: `bgmodel` in "bg\_models.py", `fgmodel` in "fg\_models.py".
Under each base class, we define several derived classes for specific fore-/back-ground models.

> base class `bgmodel` **attribute** list:

| attribute name | description |
|:---------------|:------------|
| **freqlist** | list of frequencies |
| **nfreq** | number of frequencies |
| **estimator** | band-power estimator instance |
| **params** | parameter (present) value dictionary |
| **paramdft** | default parameter values' dictionary |
| **paramrange** | parameter sampling range dictionary |
| **paramlist** | list of parameter names |
| **blacklist** | list of banned parameter names |

> base class `bgmodel` **function** list:

| function name | description |
|:--------------|:------------|
| **reset** | update **params** |

> derived class `ncmbmodel` **function** list:

| function name | description |
|:--------------|:------------|
| **initlist** | initialize **paramlist** |
| **initrange** | initialize **paramrange** |
| **initdft** | initialize **paramdft** |
| **bandpower** | prepare multi-frequency CMB band-power matrix |

> derived class `acmbmodel` **attribute** list:
(notice that these attributes are temporarily designed)

| attribute name | description |
|:---------------|:------------|
| **template\_sl** | tensor-free camb solution template |
| **template\_ps** | camb lensed-scalar template |

> derived class `acmbmodel` **function** list:

| function name | description |
|:--------------|:------------|
| **initlist** | initialize **paramlist** |
| **initrange** | initialize **paramrange** |
| **initdft** | initialize **paramdft** |
| **bandpower** | prepare multi-frequency CMB band-power matrix |

> base class `fgmodel` **attribute** list:

| attribute name | description |
|:---------------|:------------|
| **freqlist** | list of frequencies |
| **nfreq** | number of frequencies |
| **estimator** | band-power estimator instance |
| **params** | parameter (present) value dictionary |
| **paramdft** | default parameter values' dictionary |
| **paramrange** | parameter sampling range dictionary |
| **paramlist** | list of parameter names |
| **blacklist** | list of banned parameter names |
| **template\_bp** | template band-powers |
| **template\_flag** | if have **template\_bp** |
| **template\_freqlist** | frequency list for **template\_bp** |
| **template\_nfreq** | number of frequencies in **template\_freqlist** |

> base class `fgmodel` **function** list:

| function name | description |
|:--------------|:------------|
| **reset** | update **params** |
| **initdft** | initialize **paramdft** |
| **i2cmb** | brightness flux to CMB temperature unit converting factor |

> derived classes `async`, `adust`, `asyncadust`, `tsync`, `tdust`, `tsynctdust` function list:

| function name | description |
|:--------------|:------------|
| **initlist** | initialize **paramlist** |
| **initrange** | initialize **paramrange** |
| **initdft** | initialize **paramdft** |
| **bratio** | thermal dust black-body brightness calculator |
| **bandpower** | prepare multi-frequency CMB band-power matrix |

