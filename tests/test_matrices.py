import numpy as np
import pytest
from tiktaalik import pars
from tiktaalik.matrices import (
    initialize_kernels,
    initialize_evolution_matrices,
    pixelspace,
    Q2space,
    passes,
    get_nfl,
    matrix_VNS,
    matrix_VSG,
    matrix_ANS,
    matrix_ASG,
    kernel_VQQ,
    kernel_VQG,
    kernel_VGQ,
    kernel_VGG,
)


def test_pixelspace_linear_grid_midpoints():
    nx = 4
    x = pixelspace(nx, xi=0.5, grid_type=1)
    expected = np.array([-0.75, -0.25, 0.25, 0.75])
    assert x.shape == (nx,)
    assert np.allclose(x, expected, atol=1e-12)


def test_pixelspace_linear_grid_independent_of_xi():
    nx = 6
    x1 = pixelspace(nx, xi=0.3, grid_type=1)
    x2 = pixelspace(nx, xi=0.7, grid_type=1)
    assert np.allclose(x1, x2, atol=1e-12)


def test_pixelspace_linear_grid_symmetry_and_bounds():
    nx = 8
    x = pixelspace(nx, xi=0.5, grid_type=1)
    assert x.shape == (nx,)
    assert np.all(x > -1.0)
    assert np.all(x < 1.0)
    assert np.allclose(x, -x[::-1], atol=1e-12)


def test_pixelspace_log_linear_log_counts():
    nx = 8
    xi = 0.5
    x = pixelspace(nx, xi=xi, grid_type=2)
    assert x.shape == (nx,)
    assert np.all(np.diff(x) > 0)
    assert np.all(x > -1.0)
    assert np.all(x < 1.0)
    assert np.sum(x < -xi) == nx // 4
    assert np.sum(x > xi) == nx // 4
    assert np.sum((x >= -xi) & (x <= xi)) == nx // 2


def test_Q2space_length_sorted_endpoints():
    Q2i = 1.0
    Q2f = 100.0
    nQ2 = 10
    Q2 = Q2space(Q2i, Q2f, nQ2)
    assert Q2.shape == (nQ2,)
    assert np.all(np.diff(Q2) > 0)
    assert np.isclose(Q2[0], Q2i, atol=1e-12)
    assert np.isclose(Q2[-1], Q2f, atol=1e-12)


def test_Q2space_injects_thresholds():
    Q2i = 1.0
    Q2f = 100.0
    nQ2 = 10
    Q2 = Q2space(Q2i, Q2f, nQ2)
    if Q2i < pars.mc2 < Q2f:
        assert np.any(np.isclose(Q2, pars.mc2, atol=1e-12))
    if Q2i < pars.mb2 < Q2f:
        assert np.any(np.isclose(Q2, pars.mb2, atol=1e-12))


def test_passes_strictly_between():
    Q2_array = np.array([1.0, 2.0, 3.0])
    assert passes(Q2_array, 2.0) is True
    assert passes(Q2_array, 1.0) is False
    assert passes(Q2_array, 3.0) is False
    assert passes(Q2_array, 0.5) is False
    assert passes(Q2_array, 3.5) is False


def test_get_nfl_thresholds():
    assert get_nfl(pars.mc2 * 0.5) == 3
    assert get_nfl((pars.mc2 + pars.mb2) * 0.5) == 4
    assert get_nfl(pars.mb2 * 2.0) == 5


def test_initialize_kernels_asserts_nx_even_and_minimum():
    with pytest.raises(AssertionError):
        initialize_kernels(5, 0.5, grid_type=1)
    with pytest.raises(AssertionError):
        initialize_kernels(4, 0.5, grid_type=1)


def test_initialize_evolution_matrices_asserts_nQ2_minimum():
    initialize_kernels(6, 0.5, grid_type=1)
    with pytest.raises(AssertionError):
        initialize_evolution_matrices(np.array([1.0]))


def test_matrix_functions_return_finite_arrays():
    nx = 6
    xi = 0.5
    initialize_kernels(nx, xi, grid_type=1)
    Q2 = np.array([pars.mc2, pars.mb2 * 2.0])
    initialize_evolution_matrices(Q2, nlo=False)

    for M in (matrix_VNS(ns_type=1), matrix_VNS(ns_type=-1),
              matrix_VSG(), matrix_ANS(ns_type=1), matrix_ANS(ns_type=-1),
              matrix_ASG()):
        assert M.size > 0
        assert np.all(np.isfinite(M))


def test_initialize_kernels_accepts_scalar_xi():
    initialize_kernels(6, 0.5, grid_type=1)
    initialize_evolution_matrices(np.array([pars.mc2, pars.mb2 * 2.0]))
    M = matrix_VNS()
    assert np.all(np.isfinite(M))
