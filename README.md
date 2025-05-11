# About tiktaalik

**tiktaalik** is a package for creating matrices to do ultra-fast x-space evolution
of generalized parton distributions (GPDs).

The package tiktaalik is named after the animal tiktaalik,
a transitional genus of fish with several early features of amphibians.
As an iconic example of biological evolution,
tiktaalik felt like a fitting name for a code package that performs
a different kind of evolution.

### Authors

The lead developer of the code is Adam Freese.
The code was developed with the input of
Daniel Adamiak, Ian Cloët, Jianwei Qiu, Nobuo Sato and Marco Zaccheddu.

The package contains several Fortran90 codes that were previously
developed by other authors (e.g., integration methods and special functions).
Attributions to the original authors are present in the files containing these codes.

Helpful suggestions were additionally provided by Pi-Yuen Chuang and Sylvester Joosten.

### QuantOm Collaboration

This package was created for the
[QUAntum chromodynamics Nuclear TOMography (QuantOm) collaboration](https://quantom-collab.github.io/),
which supported the creation of this code.

### References

- A. Freese, D. Adamiak, I. Cloët, W. Melnitchouk, J.-W. Qiu, N. Sato, and M. Zaccheddu.
  [Computer Physics Communications 311 (2025) 109552](https://inspirehep.net/literature/2860861)

# Installation

The package can be installed using pip.
If you're in the main package directory
(the one containing `pyproject.toml`),
you just need to run:
```
pip install .
```
to install.

### Dependencies

- Python
- A Fortran compiler (ideally gfortran)
- The dependencies listed in `pyproject.toml`

matplotlib is also required to run `examples/examples.py`,
but not to install or run tiktaalik itself.

# Usage

See `examples/examples.py` for examples of usage.

# Recent updates

### May 10, 2025

The build system was refactored to use scikit-build instead of numpy.distutils,
since the latter is discontinued as of Python 3.12.
In addition to the package being future-proofed, it is now easier to install,
and the installation is managed by pip.

### February 18, 2025

tiktaalik now supports next-to-leading order (NLO) evolution of GPDs.
NLO evolution is turned off by default,
and can be activated using the kwarg `nlo=True`
when initializing the evolution matrices, fetching the kernel matrices,
or evolving with an Evolver class object.
Please see `examples.py` for examples of how to use NLO evolution.

tiktaalik now permits two kinds of spacing for the x grid:
1. The default linear spacing.
2. A log-linear-log spacing.
The latter requires nx to be divisible by 4.
Half of the x points are put in the ERBL region,
and the other half are split between the two DGLAP regions.
The points in the ERBL region are linearly spaced,
and the points in the DGLAP regions are geometrically spaced,
with greater concentration closer to the DGLAP-ERBL boundaries.

The default linear spacing works quite well and is recommended for xi above 0.2.
The log-linear-log spacing is recommended for xi below 0.1.
Either spacing is acceptable for 0.1 < xi < 0.2, provided nx is large enough.

The numerical methods used by tiktaalik are unstable when xi is too small.
The minimum xi value you should use with tiktaalik is 0.005.
By contrast, there is no maximum; tiktaalik is reliable all the way up to xi=1.

Since tiktaalik defaults to linear spacing and leading order,
any existing scripts that use tiktaalik will still give the same results as before.
