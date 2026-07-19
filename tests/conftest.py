import numpy as np
from hypothesis import settings

# Deterministic Hypothesis: a fixed derangement seed and no example database, so
# repeated runs give the same result. max_examples is kept modest because each
# example that reinitializes kernels does real numerical work.
settings.register_profile(
    "tiktaalik", derandomize=True, database=None, deadline=None, max_examples=15
)
settings.load_profile("tiktaalik")


def evolve_ns(matrix, gpd):
    """Apply a non-singlet evolution matrix to a GPD, following examples.py."""
    return np.einsum("ijkl,jk...->ik...l", matrix, gpd)
