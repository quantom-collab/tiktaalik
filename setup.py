from __future__ import annotations

from skbuild import setup

setup(
    name="tiktaalik",
    version="0.0.3",
    description="a minimal example package (fortran version)",
    author="The scikit-build team",
    license="BSD 3-clause license",
    packages=["f90src", "tiktaalik"],
    cmake_languages=("C", "Fortran"),
    cmake_minimum_required_version="3.18",
)
