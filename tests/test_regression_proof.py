"""Proof that the oracle detects a change.

The numeric core is Fortran, so rather than rebuild it we simulate a regression
at the Python level: scale the evolved result by 0.1 percent, the kind of error a
refactor of the core could introduce. The zero distance identity oracle
separates the clean run from the regressed one at its tolerance, which is what
makes it a useful regression test.
"""

import numpy as np
import tiktaalik as tk
from conftest import evolve_ns


def test_zero_distance_identity_catches_a_small_perturbation():
    nx, xi = 40, 0.3
    x = tk.matrices.pixelspace(nx)
    tk.matrices.initialize_kernels(nx, xi)
    tk.matrices.initialize_evolution_matrices(np.array([4.0, 16.0]))
    matrix = np.array(tk.matrices.matrix_VNS(ns_type=-1))
    gpd = tk.model.Hu(x, xi, 0) + tk.model.Hu(-x, xi, 0)

    clean = evolve_ns(matrix, gpd)[..., 0]
    regressed = clean * 1.001

    assert np.allclose(clean, gpd, rtol=0, atol=1e-9)
    assert not np.allclose(regressed, gpd, rtol=0, atol=1e-9)
