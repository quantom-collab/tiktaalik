import numpy as np
import pytest
from tiktaalik.model import (
    H_singlet,
    Hu,
    Hd,
    Hs,
    Hg,
    Hu_val,
    Hd_val,
    Hu_sea,
    Hd_sea,
)
X = np.array([0.15, 0.3, 0.5, 0.7, 0.85])
XI = 0.1
T = 10.0
WRAPPERS = [Hu, Hd, Hs, Hg, Hu_val, Hd_val, Hu_sea, Hd_sea]


def test_up_valence_plus_sea_equals_total():
    total = np.ravel(np.asarray(Hu(X, XI, T)))
    parts = np.ravel(np.asarray(Hu_val(X, XI, T))) + np.ravel(np.asarray(Hu_sea(X, XI, T)))
    assert np.allclose(total, parts, rtol=1e-10, atol=1e-12)


def test_down_valence_plus_sea_equals_total():
    total = np.ravel(np.asarray(Hd(X, XI, T)))
    parts = np.ravel(np.asarray(Hd_val(X, XI, T))) + np.ravel(np.asarray(Hd_sea(X, XI, T)))
    assert np.allclose(total, parts, rtol=1e-10, atol=1e-12)


def test_up_valence_plus_sea_equals_total_at_single_x():
    x = 0.3
    total = np.ravel(np.asarray(Hu(x, XI, T)))
    parts = np.ravel(np.asarray(Hu_val(x, XI, T))) + np.ravel(np.asarray(Hu_sea(x, XI, T)))
    assert np.allclose(total, parts, rtol=1e-10, atol=1e-12)


def test_down_valence_plus_sea_equals_total_at_single_x():
    x = 0.3
    total = np.ravel(np.asarray(Hd(x, XI, T)))
    parts = np.ravel(np.asarray(Hd_val(x, XI, T))) + np.ravel(np.asarray(Hd_sea(x, XI, T)))
    assert np.allclose(total, parts, rtol=1e-10, atol=1e-12)


def test_up_valence_plus_sea_equals_total_at_negative_x():
    total = np.ravel(np.asarray(Hu(-X, XI, T)))
    parts = np.ravel(np.asarray(Hu_val(-X, XI, T))) + np.ravel(np.asarray(Hu_sea(-X, XI, T)))
    assert np.allclose(total, parts, rtol=1e-10, atol=1e-12)


def test_down_valence_plus_sea_equals_total_at_negative_x():
    total = np.ravel(np.asarray(Hd(-X, XI, T)))
    parts = np.ravel(np.asarray(Hd_val(-X, XI, T))) + np.ravel(np.asarray(Hd_sea(-X, XI, T)))
    assert np.allclose(total, parts, rtol=1e-10, atol=1e-12)


def test_singlet_equals_component_combination():
    expected = (
        np.ravel(np.asarray(Hu(X, XI, T))) - np.ravel(np.asarray(Hu(-X, XI, T)))
        + np.ravel(np.asarray(Hd(X, XI, T))) - np.ravel(np.asarray(Hd(-X, XI, T)))
        + np.ravel(np.asarray(Hs(X, XI, T))) - np.ravel(np.asarray(Hs(-X, XI, T)))
    )
    result = np.ravel(np.asarray(H_singlet(X, XI, T)))
    assert np.allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_singlet_equals_component_combination_at_single_x():
    x = 0.3
    expected = (
        np.ravel(np.asarray(Hu(x, XI, T))) - np.ravel(np.asarray(Hu(-x, XI, T)))
        + np.ravel(np.asarray(Hd(x, XI, T))) - np.ravel(np.asarray(Hd(-x, XI, T)))
        + np.ravel(np.asarray(Hs(x, XI, T))) - np.ravel(np.asarray(Hs(-x, XI, T)))
    )
    result = np.ravel(np.asarray(H_singlet(x, XI, T)))
    assert np.allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_singlet_equals_component_combination_at_other_xi():
    xi = 0.25
    expected = (
        np.ravel(np.asarray(Hu(X, xi, T))) - np.ravel(np.asarray(Hu(-X, xi, T)))
        + np.ravel(np.asarray(Hd(X, xi, T))) - np.ravel(np.asarray(Hd(-X, xi, T)))
        + np.ravel(np.asarray(Hs(X, xi, T))) - np.ravel(np.asarray(Hs(-X, xi, T)))
    )
    result = np.ravel(np.asarray(H_singlet(X, xi, T)))
    assert np.allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_singlet_equals_component_combination_at_other_t():
    t = 5.0
    expected = (
        np.ravel(np.asarray(Hu(X, XI, t))) - np.ravel(np.asarray(Hu(-X, XI, t)))
        + np.ravel(np.asarray(Hd(X, XI, t))) - np.ravel(np.asarray(Hd(-X, XI, t)))
        + np.ravel(np.asarray(Hs(X, XI, t))) - np.ravel(np.asarray(Hs(-X, XI, t)))
    )
    result = np.ravel(np.asarray(H_singlet(X, XI, t)))
    assert np.allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_singlet_output_length_matches_input_length():
    result = np.ravel(np.asarray(H_singlet(X, XI, T)))
    assert result.shape == X.shape


def test_singlet_is_odd_in_x_vectorized():
    pos = np.ravel(np.asarray(H_singlet(X, XI, T)))
    neg = np.ravel(np.asarray(H_singlet(-X, XI, T)))
    assert np.allclose(neg, -pos, rtol=1e-10, atol=1e-12)
