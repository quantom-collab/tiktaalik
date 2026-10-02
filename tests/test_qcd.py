import numpy as np
import pytest
from hypothesis import given, settings, seed, strategies as st, HealthCheck
from tiktaalik.qcd import alphaQCD
Q2_LOW = 1.0
Q2_HIGH = 1.0e4


def test_scalar_input_is_wrapped_into_array():
    result = alphaQCD(10.0)
    assert isinstance(result, np.ndarray)
    assert result.shape == (1,)


def test_scalar_input_matches_single_element_array():
    scalar_result = alphaQCD(10.0)
    array_result = alphaQCD(np.array([10.0]))
    assert scalar_result.shape == array_result.shape
    assert np.allclose(scalar_result, array_result, rtol=1e-10, atol=0.0)


def test_output_length_matches_input_length():
    Q2 = np.array([2.0, 5.0, 10.0, 100.0])
    result = alphaQCD(Q2)
    assert result.shape == Q2.shape


def test_batch_call_agrees_with_elementwise_calls():
    Q2 = np.array([1.0, 2.5, 10.0, 50.0, 1000.0])
    batch = alphaQCD(Q2)
    for i in range(Q2.size):
        single = alphaQCD(float(Q2[i]))
        assert np.allclose(batch[i], single[0], rtol=1e-10, atol=0.0)


def test_repeated_values_in_batch_are_consistent():
    Q2 = np.array([7.3, 7.3, 7.3])
    result = alphaQCD(Q2)
    # Same input scale must give the same coupling, element by element.
    assert np.allclose(result, result[0], rtol=1e-12, atol=0.0)


def test_coupling_is_positive_and_finite_over_valid_range():
    Q2 = np.logspace(0, 4, 41)
    result = alphaQCD(Q2)
    assert np.all(np.isfinite(result))
    assert np.all(result > 0.0)


def test_asymptotic_freedom_coupling_decreases_with_scale():
    # QCD coupling decreases as the scale increases; check on a coarse grid
    # so the strict decrease is not masked by roundoff.
    Q2 = np.array([1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 1000.0, 10000.0])
    result = alphaQCD(Q2)
    assert np.all(np.diff(result) < 0.0)


def test_coupling_is_of_perturbative_size():
    # In the valid range the coupling should stay O(1), not blow up.
    Q2 = np.logspace(0, 4, 21)
    result = alphaQCD(Q2)
    assert np.all(result < 10.0)
