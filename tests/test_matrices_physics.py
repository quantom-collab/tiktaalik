import numpy as np


from tiktaalik import matrices, pars


def test_pixelspace_shape():
    x = np.asarray(matrices.pixelspace(20)).ravel()
    assert x.shape == (20,)


def test_pixelspace_strictly_increasing_and_bounded():
    x = np.asarray(matrices.pixelspace(20)).ravel()
    assert np.all(np.diff(x) > 0)
    assert x.min() > -1
    assert x.max() < 1


def test_pixelspace_antisymmetric():
    x = np.asarray(matrices.pixelspace(20)).ravel()
    assert np.allclose(x, -x[::-1], atol=1e-06)


def test_pixelspace_xi_independent_for_linear_grid():
    x_default = np.asarray(matrices.pixelspace(20)).ravel()
    x_xi = np.asarray(matrices.pixelspace(20, xi=0.9)).ravel()
    assert np.allclose(x_default, x_xi, atol=1e-06)


def test_Q2space_contract():
    Q2 = np.asarray(matrices.Q2space(4.0, 16.0, 6)).ravel()
    assert Q2.shape == (6,)
    assert np.all(np.diff(Q2) > 0)
    assert np.isclose(Q2[0], 4.0, atol=1e-06)
    assert np.isclose(Q2[-1], 16.0, atol=1e-06)


def test_Q2space_injects_thresholds():
    Q2 = np.asarray(matrices.Q2space(1.0, 100.0, 8)).ravel()
    assert Q2.shape == (8,)
    assert np.any(np.isclose(Q2, pars.mc2, atol=1e-06))
    assert np.any(np.isclose(Q2, pars.mb2, atol=1e-06))


def test_passes_cases():
    assert matrices.passes(np.array([1.0, 10.0]), pars.mc2) is True
    assert matrices.passes(np.array([4.0, 16.0]), pars.mc2) is False


def test_get_nfl_thresholds():
    assert matrices.get_nfl(1.0) == 3
    assert matrices.get_nfl(4.0) == 4
    assert matrices.get_nfl(100.0) == 5


def test_matrix_shapes():
    matrices.initialize_kernels(20, 0.3)
    matrices.initialize_evolution_matrices(np.array([4.0, 8.0, 12.0, 16.0]))
    Mvns = np.asarray(matrices.matrix_VNS())
    Mvsg = np.asarray(matrices.matrix_VSG())
    assert Mvns.shape == (20, 20, 1, 4)
    assert Mvsg.shape == (40, 40, 1, 4)


def test_matrix_identity_at_initial_scale():
    matrices.initialize_kernels(20, 0.3)
    matrices.initialize_evolution_matrices(np.array([4.0, 8.0, 12.0, 16.0]))
    M = np.asarray(matrices.matrix_VNS())
    assert np.allclose(M[:, :, 0, 0], np.eye(20), atol=1e-06)


def test_pixelspace_strictly_increasing_and_bounded_2():
    x = np.asarray(matrices.pixelspace(20)).ravel()
    assert np.all(np.diff(x) > 0)
    assert x.min() > -1.0
    assert x.max() < 1.0


def test_Q2space_contract_2():
    Q2 = np.asarray(matrices.Q2space(4.0, 16.0, 6)).ravel()
    assert Q2.shape == (6,)
    assert np.all(np.diff(Q2) >= 0)
    assert np.isclose(Q2[0], 4.0, atol=1e-06)
    assert np.isclose(Q2[-1], 16.0, atol=1e-06)


def test_pixelspace_xi_independent_grid_type1():
    x0 = np.asarray(matrices.pixelspace(20)).ravel()
    x1 = np.asarray(matrices.pixelspace(20, xi=0.9)).ravel()
    assert np.allclose(x0, x1, atol=1e-06)


def test_passes_cases_2():
    assert matrices.passes(np.array([1.0, 10.0]), pars.mc2)
    assert not matrices.passes(np.array([4.0, 16.0]), pars.mc2)


def test_pixelspace_xi_independent_for_linear_grid_2():
    x_default = np.asarray(matrices.pixelspace(20)).ravel()
    x_other = np.asarray(matrices.pixelspace(20, xi=0.9)).ravel()
    assert np.allclose(x_default, x_other, atol=1e-06)


def test_matrix_identity_at_initial_scale_2():
    matrices.initialize_kernels(20, 0.3)
    matrices.initialize_evolution_matrices(np.array([4.0, 8.0, 12.0, 16.0]))
    M = np.asarray(matrices.matrix_VNS())
    assert M.shape == (20, 20, 1, 4)
    assert np.allclose(M[:, :, 0, 0], np.eye(20), atol=1e-06)
