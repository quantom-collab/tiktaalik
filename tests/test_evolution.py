import numpy as np
import pytest
from tiktaalik.evolution import Evolver
from tiktaalik import pars
def _evolver(nx=12, xi=(0.1, 0.2), nQ2=6, grid_type=1):
    return Evolver(nx=nx, xi=np.array(xi), nQ2=nQ2,
                   Q2i=pars.mc2 + pars.epsilon, Q2f=pars.mb2 - pars.epsilon,
                   grid_type=grid_type)
def _random_gpd(ev, seed):
    rng = np.random.default_rng(seed)
    return rng.random((ev.nfl + 1, ev.nx, ev.nxi))


def test_output_shape_matches_declared_grid():
    ev = _evolver(nx=12, xi=(0.1, 0.2, 0.3), nQ2=6)
    gpd = _random_gpd(ev, 10)
    out = ev.evolveGPDs(gpd, 'V', nlo=False)
    assert out.shape == (ev.nfl + 1, ev.nx, ev.nxi, ev.nQ2)


def test_zero_distance_identity_V():
    ev = _evolver()
    gpd = _random_gpd(ev, 0)
    out = ev.evolveGPDs(gpd, 'V', nlo=False)
    assert np.allclose(out[..., 0], gpd, rtol=1e-12, atol=1e-12)


def test_zero_distance_identity_A():
    ev = _evolver()
    gpd = _random_gpd(ev, 1)
    out = ev.evolveGPDs(gpd, 'A', nlo=False)
    assert np.allclose(out[..., 0], gpd, rtol=1e-12, atol=1e-12)


def test_zero_distance_identity_nlo_V():
    ev = _evolver()
    gpd = _random_gpd(ev, 11)
    out = ev.evolveGPDs(gpd, 'V', nlo=True)
    assert np.allclose(out[..., 0], gpd, rtol=1e-12, atol=1e-12)


def test_zero_distance_identity_nlo_A():
    ev = _evolver()
    gpd = _random_gpd(ev, 12)
    out = ev.evolveGPDs(gpd, 'A', nlo=True)
    assert np.allclose(out[..., 0], gpd, rtol=1e-12, atol=1e-12)


def test_zero_distance_identity_single_xi():
    ev = _evolver(xi=(0.3,))
    gpd = _random_gpd(ev, 13)
    out = ev.evolveGPDs(gpd, 'V', nlo=False)
    assert np.allclose(out[..., 0], gpd, rtol=1e-12, atol=1e-12)


def test_zero_distance_identity_grid_type_2():
    ev = _evolver(nx=12, xi=(0.1, 0.2), grid_type=2)
    gpd = _random_gpd(ev, 14)
    out = ev.evolveGPDs(gpd, 'V', nlo=False)
    assert np.allclose(out[..., 0], gpd, rtol=1e-12, atol=1e-12)


def test_linearity_V():
    ev = _evolver()
    f = _random_gpd(ev, 2)
    g = _random_gpd(ev, 3)
    a, b = 0.7, -1.3
    combo = ev.evolveGPDs(a * f + b * g, 'V', nlo=False)
    lin = a * ev.evolveGPDs(f, 'V', nlo=False) + b * ev.evolveGPDs(g, 'V', nlo=False)
    assert np.allclose(combo, lin, rtol=1e-10, atol=1e-12)


def test_linearity_A():
    ev = _evolver()
    f = _random_gpd(ev, 4)
    g = _random_gpd(ev, 5)
    a, b = -0.4, 2.1
    combo = ev.evolveGPDs(a * f + b * g, 'A', nlo=False)
    lin = a * ev.evolveGPDs(f, 'A', nlo=False) + b * ev.evolveGPDs(g, 'A', nlo=False)
    assert np.allclose(combo, lin, rtol=1e-10, atol=1e-12)


def test_linearity_nlo_V():
    ev = _evolver()
    f = _random_gpd(ev, 15)
    g = _random_gpd(ev, 16)
    a, b = 1.5, -0.9
    combo = ev.evolveGPDs(a * f + b * g, 'V', nlo=True)
    lin = a * ev.evolveGPDs(f, 'V', nlo=True) + b * ev.evolveGPDs(g, 'V', nlo=True)
    assert np.allclose(combo, lin, rtol=1e-10, atol=1e-12)


def test_composition_transitivity_V():
    Q0 = pars.mc2 + pars.epsilon
    Qf = pars.mb2 - pars.epsilon
    Qmid = np.sqrt(Q0 * Qf)
    xi = np.array([0.1, 0.2])
    evF = Evolver(nx=12, xi=xi, nQ2=6, Q2i=Q0, Q2f=Qf, grid_type=1)
    ev1 = Evolver(nx=12, xi=xi, nQ2=6, Q2i=Q0, Q2f=Qmid, grid_type=1)
    ev2 = Evolver(nx=12, xi=xi, nQ2=6, Q2i=Qmid, Q2f=Qf, grid_type=1)
    gpd = _random_gpd(evF, 6)
    direct = evF.evolveGPDs(gpd, 'V', nlo=False)[..., -1]
    mid = ev1.evolveGPDs(gpd, 'V', nlo=False)[..., -1]
    chain = ev2.evolveGPDs(mid, 'V', nlo=False)[..., -1]
    assert np.allclose(direct, chain, rtol=2e-2, atol=1e-6)


def test_composition_transitivity_A():
    Q0 = pars.mc2 + pars.epsilon
    Qf = pars.mb2 - pars.epsilon
    Qmid = np.sqrt(Q0 * Qf)
    xi = np.array([0.1, 0.2])
    evF = Evolver(nx=12, xi=xi, nQ2=6, Q2i=Q0, Q2f=Qf, grid_type=1)
    ev1 = Evolver(nx=12, xi=xi, nQ2=6, Q2i=Q0, Q2f=Qmid, grid_type=1)
    ev2 = Evolver(nx=12, xi=xi, nQ2=6, Q2i=Qmid, Q2f=Qf, grid_type=1)
    gpd = _random_gpd(evF, 7)
    direct = evF.evolveGPDs(gpd, 'A', nlo=False)[..., -1]
    mid = ev1.evolveGPDs(gpd, 'A', nlo=False)[..., -1]
    chain = ev2.evolveGPDs(mid, 'A', nlo=False)[..., -1]
    assert np.allclose(direct, chain, rtol=2e-2, atol=1e-6)


def test_finiteness_V_nlo():
    ev = _evolver(nx=12, xi=(0.05, 0.25, 0.5), nQ2=6)
    gpd = _random_gpd(ev, 8)
    out = ev.evolveGPDs(gpd, 'V', nlo=True)
    assert np.isfinite(out).all()


def test_finiteness_A_nlo():
    ev = _evolver(nx=12, xi=(0.05, 0.25, 0.5), nQ2=6)
    gpd = _random_gpd(ev, 9)
    out = ev.evolveGPDs(gpd, 'A', nlo=True)
    assert np.isfinite(out).all()
