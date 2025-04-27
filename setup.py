from __future__ import annotations

from skbuild import setup

setup(
    packages=["f90src", "tiktaalik"],
    cmake_languages=("C", "Fortran"),
    cmake_minimum_required_version="3.18",
)
