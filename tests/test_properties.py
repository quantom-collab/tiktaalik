"""Property and metamorphic tests for GPD evolution.

These give true or false verdicts without knowing the correct numeric output.
The strongest oracles for evolution code are the identity at zero distance,
linearity, and composition, so those come first. Tolerances are set to the
numerical method, not loosened to pass. The observed errors are far tighter than
the tolerances asserted here (identity is exact, composition is order 1e-6,
linearity is at machine precision).
"""

import numpy as np
import pytest
import tiktaalik as tk
from conftest import evolve_ns
from hypothesis import given
from hypothesis import strategies as st

NX = 40
XI = 0.3
Q2I = 4.0
Q2F = 16.0


def _ns_matrix(q2):
    tk.matrices.initialize_kernels(NX, XI)
    tk.matrices.initialize_evolution_matrices(np.asarray(q2, dtype=float))
    return np.array(tk.matrices.matrix_VNS(ns_type=-1))


@pytest.fixture(scope="module")
def linear_case():
    x = tk.matrices.pixelspace(NX)
    matrix = _ns_matrix([Q2I, Q2F])
    u = tk.model.Hu(x, XI, 0) + tk.model.Hu(-x, XI, 0)
    d = tk.model.Hd(x, XI, 0) + tk.model.Hd(-x, XI, 0)
    return matrix, u, d


def test_zero_distance_identity():
    x = tk.matrices.pixelspace(NX)
    matrix = _ns_matrix([Q2I, Q2F])
    gpd = tk.model.Hu(x, XI, 0) + tk.model.Hu(-x, XI, 0)
    evolved = evolve_ns(matrix, gpd)
    # Q2[0] is the initial scale, so its slice is the identity operator.
    np.testing.assert_allclose(evolved[..., 0], gpd, rtol=0, atol=1e-9)


def test_evolver_zero_distance_identity():
    x = tk.matrices.pixelspace(NX)
    darwin = tk.Evolver(nx=NX, xi=XI, Q2i=Q2I, Q2f=Q2F, nQ2=4)
    gpd = tk.model.Hu(x, XI, 0) + tk.model.Hu(-x, XI, 0)
    evolved = darwin.evolveNS(gpd, ns_type=-1)
    np.testing.assert_allclose(evolved[..., 0], gpd, rtol=0, atol=1e-9)


def test_composition_is_transitive():
    x = tk.matrices.pixelspace(NX)
    gpd = tk.model.Hu(x, XI, 0) + tk.model.Hu(-x, XI, 0)

    direct = evolve_ns(_ns_matrix([Q2I, Q2F]), gpd)[..., 1]

    half = evolve_ns(_ns_matrix([Q2I, 8.0]), gpd)[..., 1]
    composed = evolve_ns(_ns_matrix([8.0, Q2F]), half)[..., 1]

    scale = np.max(np.abs(direct))
    np.testing.assert_allclose(composed, direct, rtol=1e-3, atol=1e-6 * scale)


@given(
    a=st.floats(-5.0, 5.0, allow_nan=False, allow_infinity=False),
    b=st.floats(-5.0, 5.0, allow_nan=False, allow_infinity=False),
)
def test_evolution_is_linear(linear_case, a, b):
    matrix, u, d = linear_case
    lhs = evolve_ns(matrix, a * u + b * d)
    rhs = a * evolve_ns(matrix, u) + b * evolve_ns(matrix, d)
    np.testing.assert_allclose(lhs, rhs, rtol=1e-9, atol=1e-9)


def test_evolver_matches_direct_matrix():
    x = tk.matrices.pixelspace(NX)
    gpd = tk.model.Hu(x, XI, 0) + tk.model.Hu(-x, XI, 0)

    q2 = tk.matrices.Q2space(Q2I, Q2F, 4)
    direct = evolve_ns(_ns_matrix(q2), gpd)

    darwin = tk.Evolver(nx=NX, xi=XI, Q2i=Q2I, Q2f=Q2F, nQ2=4)
    via_evolver = darwin.evolveNS(gpd, ns_type=-1)

    scale = np.max(np.abs(direct))
    np.testing.assert_allclose(via_evolver, direct, rtol=1e-6, atol=1e-8 * scale)


@given(xi=st.floats(5e-3, 0.9, allow_nan=False, allow_infinity=False))
def test_outputs_are_finite_in_valid_range(xi):
    # The GK model is stable down to xi=5e-3, so inputs stay at or above it.
    x = tk.matrices.pixelspace(NX)
    tk.matrices.initialize_kernels(NX, xi)
    tk.matrices.initialize_evolution_matrices(np.array([Q2I, Q2F]))
    matrix = tk.matrices.matrix_VNS(ns_type=-1)
    gpd = tk.model.Hu(x, xi, 0) + tk.model.Hu(-x, xi, 0)
    evolved = evolve_ns(matrix, gpd)
    assert np.all(np.isfinite(evolved))


def test_linear_pixelspace_is_symmetric_and_ordered():
    x = tk.matrices.pixelspace(NX)
    assert x.shape == (NX,)
    assert np.all(x > -1) and np.all(x < 1)
    assert np.all(np.diff(x) > 0)
    # Linear midpoints are symmetric about zero. The grid is single precision,
    # so the tolerance is set to float32, not float64.
    np.testing.assert_allclose(x, -x[::-1], rtol=0, atol=1e-6)
