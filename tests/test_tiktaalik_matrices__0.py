import numpy as np
import pytest
from hypothesis import given, strategies as st, settings, HealthCheck, seed

from tiktaalik import pars
from tiktaalik.matrices import Q2space


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@seed(0)
@given(
    Q2i=st.floats(min_value=0.1, max_value=5.0, allow_nan=False, allow_infinity=False),
    Q2f=st.floats(min_value=5.1, max_value=100.0, allow_nan=False, allow_infinity=False),
    nQ2=st.integers(min_value=3, max_value=20),
)
def test_Q2space_matches_geomspace_when_no_thresholds(Q2i, Q2f, nQ2):
    """If no quark‑mass thresholds lie between Q2i and Q2f the result must be
    identical to a plain geomspace."""
    thresholds = [m2 for m2 in (pars.mc2, pars.mb2) if Q2i < m2 < Q2f]
    assume = len(thresholds) == 0
    if not assume:
        pytest.skip("generated range contains a threshold")
    q2 = Q2space(Q2i, Q2f, nQ2)
    expected = np.geomspace(Q2i, Q2f, nQ2)
    np.testing.assert_allclose(q2, expected, rtol=0, atol=1e-12)


def test_Q2space_includes_known_thresholds_and_is_sorted():
    """When thresholds lie inside the interval they must appear in the output,
    the array must be sorted and have the requested length."""
    # Choose a range that definitely contains both charm and bottom thresholds.
    Q2i = pars.mc2 - 1.0
    Q2f = pars.mb2 + 1.0
    # nQ2 must be at least the number of thresholds + 1.
    nQ2 = 7

    q2 = Q2space(Q2i, Q2f, nQ2)

    # Length must equal the requested number of points.
    assert q2.shape == (nQ2,)

    # The array must be strictly increasing.
    assert np.all(np.diff(q2) > 0)

    # All thresholds that lie inside the interval must be present.
    thresholds = [m2 for m2 in (pars.mc2, pars.mb2) if Q2i < m2 < Q2f]
    for thr in thresholds:
        assert np.any(np.isclose(q2, thr, atol=0, rtol=0))

    # The result must be a sorted version of the geomspace with thresholds appended.
    base = np.geomspace(Q2i, Q2f, nQ2 - len(thresholds))
    combined = np.sort(np.append(base, thresholds))
    np.testing.assert_allclose(q2, combined, rtol=0, atol=1e-12)


def test_Q2space_excludes_thresholds_outside_interval():
    """Thresholds outside the (Q2i, Q2f) interval must not appear in the output."""
    # Pick a range that lies entirely below both thresholds.
    Q2i = 0.1
    Q2f = pars.mc2 - 0.5
    nQ2 = 5

    q2 = Q2space(Q2i, Q2f, nQ2)

    thresholds = [m2 for m2 in (pars.mc2, pars.mb2) if Q2i < m2 < Q2f]
    # No thresholds should be inside the interval.
    assert thresholds == []

    # Ensure the output matches a pure geomspace.
    expected = np.geomspace(Q2i, Q2f, nQ2)
    np.testing.assert_allclose(q2, expected, rtol=0, atol=1e-12)
