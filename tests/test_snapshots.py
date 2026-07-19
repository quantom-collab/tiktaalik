"""Golden master snapshot tests.

Default settings (linear spacing, leading order) reproduce prior results, so a
default run is a sanctioned regression anchor. Fixed inputs produce fixed
outputs; the outputs are recorded once and future runs must match within
tolerance. This detects change, which is what regression testing needs.

The first run records the snapshot and skips. Commit the recorded files so later
runs compare against them.
"""

from pathlib import Path

import numpy as np
import pytest
import tiktaalik as tk
from conftest import evolve_ns

SNAP = Path(__file__).parent / "_snapshots"

NX = 60
XI = np.array([0.1, 0.2, 0.3])
Q2 = np.array([4.0, 8.0, 12.0, 16.0])
T = np.array([0.0, -0.2, -0.4])


def _check(name, arrays, rtol, atol):
    SNAP.mkdir(exist_ok=True)
    path = SNAP / f"{name}.npz"
    if not path.exists():
        np.savez(path, **arrays)
        pytest.skip(f"recorded new snapshot: {name}")
    reference = np.load(path)
    for key, value in arrays.items():
        np.testing.assert_allclose(value, reference[key], rtol=rtol, atol=atol)


def test_default_non_singlet_run_matches_snapshot():
    x = tk.matrices.pixelspace(NX)
    tk.matrices.initialize_kernels(NX, XI)
    tk.matrices.initialize_evolution_matrices(Q2)
    matrix = tk.matrices.matrix_VNS()
    umin = tk.model.Hu(x, XI, T) + tk.model.Hu(-x, XI, T)
    evolved = evolve_ns(matrix, umin)
    _check("default_ns_umin", {"evolved": evolved}, rtol=1e-9, atol=1e-12)


def test_default_singlet_gluon_run_matches_snapshot():
    x = tk.matrices.pixelspace(NX)
    tk.matrices.initialize_kernels(NX, XI)
    tk.matrices.initialize_evolution_matrices(Q2)
    matrix = tk.matrices.matrix_VSG()
    singlet = tk.model.H_singlet(x, XI, T)
    gluon = tk.model.Hg(x, XI, T)
    sg = np.concatenate((singlet, gluon), axis=0)
    evolved = np.einsum("ijkl,jk...->ik...l", matrix, sg)
    _check("default_sg", {"evolved": evolved}, rtol=1e-9, atol=1e-12)
