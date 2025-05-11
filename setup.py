from __future__ import annotations
from skbuild import setup
import re
from pathlib import Path

# Mitigation of issue mentioned here:
# https://github.com/scikit-build/scikit-build/issues/521
for i in (Path(__file__).resolve().parent / "_skbuild").rglob("CMakeCache.txt"):
    i.write_text(re.sub("^//.*$\n^[^#].*build-env.*$", "", i.read_text(), flags=re.M))

setup(
    name="tiktaalik",
    version="0.3.0",
    packages=["f90src", "tiktaalik"],
    cmake_languages=("C", "Fortran"),
    cmake_minimum_required_version="3.18",
)
