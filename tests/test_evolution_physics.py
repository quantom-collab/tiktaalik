import numpy as np


import pytest


import tiktaalik as tk


from tiktaalik.evolution import Evolver, make_evo_basis, make_flv_basis, Tpls_matrix, add_pad, swap_AV


def test_zero_distance_identity():
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    x = tk.matrices.pixelspace(nx)
    g = tk.model.Hu(x, xi, 0.0) + tk.model.Hu(-x, xi, 0.0)
    out = ev.evolveNS(g)
    assert np.allclose(out[..., 0], g, rtol=1e-09, atol=1e-11)


def test_linearity_of_evolution():
    rng = np.random.default_rng(1)
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    f = rng.standard_normal((nx, 1))
    g = rng.standard_normal((nx, 1))
    a, b = (0.7, -1.3)
    lhs = ev.evolveNS(a * f + b * g)
    rhs = a * ev.evolveNS(f) + b * ev.evolveNS(g)
    assert np.allclose(lhs, rhs, rtol=1e-09, atol=1e-11)


def test_t_independence():
    rng = np.random.default_rng(2)
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    g = rng.standard_normal((nx, 1, 3))
    out = ev.evolveNS(g)
    for k in range(3):
        assert np.allclose(out[..., k, :], ev.evolveNS(g[..., k]), rtol=1e-09, atol=1e-11)


def test_evolveNS_output_shape():
    rng = np.random.default_rng(3)
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    g = rng.standard_normal((nx, 1))
    out = ev.evolveNS(g)
    assert out.shape == (nx, 1, 4)


def test_flavor_basis_round_trip():
    rng = np.random.default_rng(0)
    nfl, nx, nxi = (4, 12, 2)
    g = rng.standard_normal((nfl + 1, nx, nxi))
    sg, t_grid, qmin = make_evo_basis(g, nfl=nfl)
    back = make_flv_basis(sg, t_grid, qmin, nfl=nfl)
    assert np.allclose(back, g, rtol=1e-09, atol=1e-11)


def test_tpls_matrix_invertible():
    for nfl in (3, 4, 5):
        M = Tpls_matrix(nfl)
        assert M.shape == (nfl, nfl)
        assert np.allclose(np.linalg.inv(M) @ M, np.eye(nfl), rtol=1e-09, atol=1e-11)


def test_add_pad_contract():
    rng = np.random.default_rng(4)
    g = rng.standard_normal((3, 5, 2))
    n = 2
    out = add_pad(g, n)
    assert out.shape[0] == g.shape[0] + n
    assert np.allclose(out[:g.shape[0]], g, rtol=1e-09, atol=1e-11)
    assert np.allclose(out[g.shape[0]:], 0.0, rtol=1e-09, atol=1e-11)


def test_swap_av():
    assert swap_AV('A') == 'V'
    assert swap_AV('V') == 'A'
    with pytest.raises(ValueError):
        swap_AV('x')


def test_composition_transitivity():
    nx, xi = (20, 0.3)
    evF = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    ev1 = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=8.0, grid_type=1)
    ev2 = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=8.0, Q2f=16.0, grid_type=1)
    x = tk.matrices.pixelspace(nx)
    g = tk.model.Hu(x, xi, 0.0) + tk.model.Hu(-x, xi, 0.0)
    direct = evF.evolveNS(g)[..., -1]
    step = ev2.evolveNS(ev1.evolveNS(g)[..., -1])[..., -1]
    assert np.allclose(direct, step, rtol=0.02, atol=1e-06)


def test_finiteness_lo():
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    x = tk.matrices.pixelspace(nx)
    g = tk.model.Hu(x, xi, 0.0) + tk.model.Hu(-x, xi, 0.0)
    out = ev.evolveNS(g, nlo=False)
    assert np.isfinite(out).all()


def test_finiteness_nlo():
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    x = tk.matrices.pixelspace(nx)
    g = tk.model.Hu(x, xi, 0.0) + tk.model.Hu(-x, xi, 0.0)
    out = ev.evolveNS(g, nlo=True)
    assert np.isfinite(out).all()


def test_helicity_zero_distance():
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    x = tk.matrices.pixelspace(nx)
    g = tk.model.Hu(x, xi, 0.0) + tk.model.Hu(-x, xi, 0.0)
    out = ev.evolveNS(g, gpd_type='A')
    assert np.allclose(out[..., 0], g, rtol=1e-09, atol=1e-11)


def test_evolveGPDs_zero_distance():
    rng = np.random.default_rng(5)
    nx, xi, nfl = (20, 0.3, 5)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=1)
    gpd = rng.standard_normal((nfl + 1, nx, 1))
    out = ev.evolveGPDs(gpd, 'V')
    assert np.allclose(out[..., 0], gpd, rtol=1e-09, atol=1e-11)


def test_grid_type2_runs():
    nx, xi = (20, 0.3)
    ev = Evolver(nx=nx, xi=xi, nQ2=4, Q2i=4.0, Q2f=16.0, grid_type=2)
    x = tk.matrices.pixelspace(nx)
    g = tk.model.Hu(x, xi, 0.0) + tk.model.Hu(-x, xi, 0.0)
    out = ev.evolveNS(g)
    assert np.isfinite(out).all()
