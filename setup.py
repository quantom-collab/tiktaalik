from __future__ import annotations
from skbuild import setup
import re
from pathlib import Path
import shutil

# Mitigation of issue mentioned here:
# https://github.com/scikit-build/scikit-build/issues/521
# The old mitigation is no longer working, so take a more extreme approach...
skpath = Path(__file__).resolve().parent / "_skbuild"
if(skpath.exists()):
    shutil.rmtree(skpath)

setup(
    name="tiktaalik",
    version="0.3.0",
    packages=["f90src", "tiktaalik"],
    cmake_languages=("C", "Fortran"),
    cmake_minimum_required_version="3.18",
)
