# About tiktaalik

**tiktaalik** is a package for creating matrices to do ultra-fast x-space evolution
of generalized parton distributions (GPDs).

This package was created by Adam Freese for the
[QUAntum chromodynamics Nuclear TOMography (QuantOm) collaboration](https://quantom-collab.github.io/),
which also supported the creation of this code.
The package contains several Fortran90 codes that were previously
developed by other authors (e.g., integration methods and special functions).
Attributions to the original authors are present in the files containing these codes.

The package tiktaalik is named after the animal tiktaalik,
a transitional genus of fish with several early features of amphibians.
As an iconic example of biological evolution,
tiktaalik felt like a fitting name for a code package that performs
a different kind of evolution.

Some of my code was based on ideas suggested by
Daniel Adamiak, Ian Cloët, Jianwei Qiu, Nobuo Sato and Marco Zaccheddu,
to all of whom I am grateful.
The code would be slower and less accurate without their suggestions.

# Reference

If you use tiktaalik in your research, please cite the paper!
The paper is [Kernel methods for evolution of generalized parton distributions](https://inspirehep.net/literature/2860861),
by A. Freese, D. Adamiak, I. Cloët, W. Melnitchouk, J.-W. Qiu, N. Sato, and M. Zaccheddu.
It is published as *Compututer Physics Communications 311 (2025) 109552*,
and the preprint is available on on arXiv at [2412.13450](https://arxiv.org/abs/2412.13450).

# Dependencies

- cmake version 3.12 or greater
- A Fortran compiler (ideally gfortran)
- Python 3.11 or less
- python-numpy

matplotlib is also optionally required to make and show plots in the example script,
but not to actually use the evolution package.

### Regarding Python version

Because of a regression in Python 3.12 (the deprecation of distutils),
tiktaalik cannot be built in Python>=3.12.
Support for Python>=3.12 is delayed until I can figure out how
to use an alternative build system.

In the meantime, if your system has Python>=3.12 and you would like to use tiktaalik,
I recommend using [pyenv](https://github.com/pyenv/pyenv),
which allows you to install multiple Python versions and switch between them.

# Installation

To install the package:
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/installation
make
make install
```
The `-DCMAKE_INSTALL_PREFIX` may be unnecessary for a system-wide installation.
For a system-wide installation, sudo privileges are probably necessary for the
`make install` step.

If the installation is not system-wide, your `LD_LIBRARY_PATH` and
`PYTHONPATH` environment variables will need to be updated for the tiktaalik
libraries to be found by Python.

# Usage

The user interfaces are all part of the tiktaalik Python package.
One runs:
```
import tiktaalik as tk
```
to import the package in Python.

The examples script `examples.py` provided with the package
provides demonstrations of how to use tiktaalik.

# Recent updates

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
