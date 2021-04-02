### documentation for AFRA

### design

The AFRA package consists of 4 parts: pipeline modules, method modules, model modules and tools modules.

Pipeline modules contains a base class "pipe" and two derived classes "abspipe" and "fitpipe", 
where "abspipe" is designed for ABS analysis and "fitpipe" is for template/model fitting.

The pipelines are in charge of arranging/assembling the workflows of analysis, while the kernels
for either ABS or fitting are implemented in the "method" modules.
We design respectively two method classes, which are "abssep" and "fit". 

From another view,
we can say that the key functions implemented in the pipelines are related to converting pixel-domain maps
into harmonic-domain band-power matrices and arrange them in the right shape, and to how we can use functions/modules
defined elsewhere to help us find the final results. The method modules belong to the helpers group, along with other
key contributors like the band-power estimator and the models.

The CMB back/fore-ground models are defined in the model modules. There are two base classes "bgmodel" and "fgmodel".
"bgmodel" defines the basic rules for constructing a CMB background model band-power matrices, while "fgmodel" says
how we do for the CMB foregrounds. Several specific models are designed under the base classes: "acmbmodel" stands for
analytic CMB, "ncmbmodel" for numeric CMB, "asyncmodel" for analytic synchrotron emission, "adustmodel" for analytic dust emission,
"asyncadustmodel" for analytic synchrotron-dust emission, "tsyncmodel" for templated synchrotron emission, "tdustmodel" for
templated dust emission, "tsynctdustmodel" for templated synchrotron-dust emission.

Other auxiliary functions are designed as "tools" modules. The most important auxiliary class is "pstimator", which is basically
a wrapper of NaMaster band-power estimation routines. We wrapped up the routine into several easy-to-use functions so in the
pipelines we can call corresponding utilities in a simple way. Other useful methods are defined in either "aux" or "bp_vis"
modules. The "icy_decorator" is a class decorator which can prevent typos in calling/assigning attributes.

For more detailed descriptions please check lower level READMEs or the doc-strings in each module.
