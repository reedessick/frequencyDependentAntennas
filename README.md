# frequencyDependentAntennas

This module implements an expression for the frequency-dependent antenna response functions of ground-based gravitational wave detectors.

## To Do:

  - develop a likelihood functional given a PSD
  - plug the likelihood into emcee or an equivalent to sample over localization parameters
    - need to check for diagonal pp plots
    - if this doesn't work, compute the Fisher matrix for localization parameters and use that
  - set up injection sets and recover both with and without the frequency dependence in the antenna patterns
    - quantify the change in localization areas
    - quantify the biases associated with the omission
    - do this for single detectors as well as a few simple networks
        - aLIGO LHO + aLIGO LLO + aVirgo
        - aLIGO LHO + CE (at LLO)
        - ET + CE?
  - consider how directionality and frequency dependence of anteanna patterns could affect how templated searches work
    - phase involved in antenna patterns (which are different for h+ and hx) could mean we cannot simply search for a single template?
    - this could be compounded by a directional dependence, meaning one has to perform an explicit all sky search?

## References

Expressions are taken from LIGO-T060237 and Anderson, et all PhysRevD 63(04) 2003
