# Pathfinder
A Bayesian-inferred simple climate model.


## How-to

Download a release. *Read The Fine Manual.*

Pathfinder has been developed in Python 3.7 and run preferentially through IPython. Currently, packages required to run it are `numpy` (v1.19.2), `scipy` (v1.5.2) and `xarray` (v0.16.0), and it has hard-coded dependencies on `pymc3` (v3.8) and `theano` (v1.0.4) that are in fact used only for calibration. Newer versions of Python or these packages are likely to work, although they were not tested.


## Known issues

* The model requires a high number of substeps to remain stable under high CO2 (because of the ocean carbon cycle). This can be set using the `nt` argument when calling `run_xarray`.

* The temperature-driven mode (`Tdriven`) is extremely sensitive to its forcings: it can be very difficult to make it transition smoothly from historical to projections. This is unavoidable because mathematically it requires the second derivative of `T` and the first derivative of `ERFx` as input.

* Unclear whether the `my_AR1` class from `cls_calib` is actually needed.


## Changelog

##### v1.0.1
Exact same physical equations and numerical values as v1.0.
* Added: best-guess parameters and outputs (in `internal_data/pyMC_calib/`) for single-configuration runs.
* Improved: README and MANUAL files.

##### v1.0
Exact model described by Bossy et al. (subm).

#### v1
First release!


## References

**v1.0 (full) |** Bossy, T., T. Gasser & P. Ciais. "Pathfinder v1.0: a Bayesian-inferred simple carbon-climate model to explore climate change scenarios." *Geoscientific Model Development* (submitted). [doi:10.5194/egusphere-2022-802](https://doi.org/10.5194/egusphere-2022-802).

