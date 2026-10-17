import numpy as np
import tiktaalik.model as m
def _val(r):
    # the GK model returns a length-1 array for scalar inputs
    return float(np.asarray(r).ravel()[0])
def _phys(rng):
    # one physical kinematic point: x in [-0.9,0.9], xi in (0,0.5], t in [-0.6,0]
    x = rng.uniform(-0.9, 0.9)
    xi = rng.uniform(0.05, 0.5)
    t = -rng.uniform(0.0, 0.6)
    return x, xi, t


def test_singlet_definition():
    rng = np.random.default_rng(10)
    for _ in range(50):
        x, xi, t = _phys(rng)
        expected = (
            _val(m.Hu(x, xi, t)) - _val(m.Hu(-x, xi, t))
            + _val(m.Hd(x, xi, t)) - _val(m.Hd(-x, xi, t))
            + _val(m.Hs(x, xi, t)) - _val(m.Hs(-x, xi, t))
        )
        assert np.isclose(_val(m.H_singlet(x, xi, t)), expected, rtol=1e-9, atol=1e-11)


def test_singlet_odd_in_x():
    rng = np.random.default_rng(11)
    for _ in range(50):
        x, xi, t = _phys(rng)
        assert np.isclose(
            _val(m.H_singlet(-x, xi, t)), -_val(m.H_singlet(x, xi, t)), rtol=1e-9, atol=1e-11
        )


def test_sea_odd_in_x():
    rng = np.random.default_rng(12)
    for _ in range(50):
        x, xi, t = _phys(rng)
        assert np.isclose(
            _val(m.Hu_sea(-x, xi, t)), -_val(m.Hu_sea(x, xi, t)), rtol=1e-9, atol=1e-11
        )
        assert np.isclose(
            _val(m.Hd_sea(-x, xi, t)), -_val(m.Hd_sea(x, xi, t)), rtol=1e-9, atol=1e-11
        )


def test_gluon_even_in_x():
    rng = np.random.default_rng(13)
    for _ in range(50):
        x, xi, t = _phys(rng)
        assert np.isclose(_val(m.Hg(-x, xi, t)), _val(m.Hg(x, xi, t)), rtol=1e-9, atol=1e-11)


def test_valence_vanishes_erbl():
    rng = np.random.default_rng(14)
    for _ in range(50):
        xi = rng.uniform(0.05, 0.5)
        t = -rng.uniform(0.0, 0.6)
        x = rng.uniform(-0.95, -xi - 0.01)  # strictly below -xi
        assert np.isclose(_val(m.Hu_val(x, xi, t)), 0.0, atol=1e-11)
        assert np.isclose(_val(m.Hd_val(x, xi, t)), 0.0, atol=1e-11)


def test_sea_flavor_equal():
    rng = np.random.default_rng(15)
    for _ in range(50):
        x, xi, t = _phys(rng)
        assert np.isclose(
            _val(m.Hu_sea(x, xi, t)), _val(m.Hd_sea(x, xi, t)), rtol=1e-9, atol=1e-11
        )


def test_sea_from_strange():
    rng = np.random.default_rng(16)
    for _ in range(50):
        x, xi, t = _phys(rng)
        assert np.isclose(
            _val(m.Hu_sea(x, xi, t)), 1.68 * _val(m.Hs(x, xi, t)), rtol=1e-9, atol=1e-11
        )


def test_quark_val_plus_sea():
    rng = np.random.default_rng(17)
    for _ in range(50):
        x, xi, t = _phys(rng)
        assert np.isclose(
            _val(m.Hu(x, xi, t)),
            _val(m.Hu_val(x, xi, t)) + _val(m.Hu_sea(x, xi, t)),
            rtol=1e-9,
            atol=1e-11,
        )
        assert np.isclose(
            _val(m.Hd(x, xi, t)),
            _val(m.Hd_val(x, xi, t)) + _val(m.Hd_sea(x, xi, t)),
            rtol=1e-9,
            atol=1e-11,
        )
