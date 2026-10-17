import numpy as np


from tiktaalik.qcd import alphaQCD


def _arr(r):
    return np.asarray(r).ravel()


def test_alpha_positive():
    Q2 = np.geomspace(4.0, 10000.0, 40)
    a = _arr(alphaQCD(Q2))
    assert np.all(a > 0)


def test_alpha_finite():
    Q2 = np.geomspace(4.0, 10000.0, 40)
    a = _arr(alphaQCD(Q2))
    assert np.all(np.isfinite(a))


def test_alpha_monotone_decreasing():
    Q2 = np.geomspace(4.0, 10000.0, 40)
    a = _arr(alphaQCD(Q2))
    assert np.all(np.diff(a) < 0)


def test_alpha_scalar_equals_array():
    scalar = _arr(alphaQCD(100.0))
    array = _arr(alphaQCD(np.array([100.0])))
    assert scalar.shape == array.shape
    assert np.allclose(scalar, array, atol=1e-12)


def test_alpha_golden_master():
    assert np.isclose(float(_arr(alphaQCD(10.0))[0]), 0.246303788, atol=1e-06)
    assert np.isclose(float(_arr(alphaQCD(100.0))[0]), 0.1778742689, atol=1e-06)
    assert np.isclose(float(_arr(alphaQCD(1000.0))[0]), 0.1405352341, atol=1e-06)


def test_alpha_scalar_equals_array_2():
    scalar = _arr(alphaQCD(100.0))
    array = _arr(alphaQCD(np.array([100.0])))
    assert scalar.shape == array.shape
    assert np.allclose(scalar, array, atol=1e-12, rtol=0.0)
