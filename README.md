# About tiktaalik

**tiktaalik** is a package for creating matrices to do ultra-fast x-space evolution
of generalized parton distributions (GPDs).

This package was created by Adam Freese for the
[QUAntum chromodynamics Nuclear TOMography (QuantOm) collaboration](https://www.anl.gov/phy/quantom),
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
Daniel Adamiak, Ian Cloet, Jianwei Qiu, Nobuo Sato and Marco Zaccheddu,
to all of whom I am grateful.
The code would be slower and less accurate without their suggestions.

# Reference

If you use tiktaalik in your research, please cite the paper!

The paper is in preparation.
Reference information will be added here later.

# Dependencies

- cmake version 3.12 or greater
- A Fortran compiler (ideally gfortran)
- Python
- python-numpy

matplotlib is also optionally required to make and show plots in the example script,
but not to actually use the evolution package.

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
