"""Live demonstration that a property test catches an injected regression.

Run directly: python tests/prove_regression.py

The numeric core is Fortran, so the regression is injected at the Python level by
scaling the evolution matrix by 0.1 percent. The zero distance identity oracle
holds on the clean code and fails on the regressed code.
"""

import numpy as np
import tiktaalik as tk

nx, xi = 40, 0.3
x = tk.matrices.pixelspace(nx)
tk.matrices.initialize_kernels(nx, xi)
tk.matrices.initialize_evolution_matrices(np.array([4.0, 16.0]))
gpd = tk.model.Hu(x, xi, 0) + tk.model.Hu(-x, xi, 0)


def identity_holds(matrix):
    evolved = np.einsum("ijkl,jk...->ik...l", matrix, gpd)
    return np.allclose(evolved[..., 0], gpd, rtol=0, atol=1e-9)


clean = np.array(tk.matrices.matrix_VNS(ns_type=-1))
regressed = clean * 1.001

print("clean code:     zero-distance identity holds:", identity_holds(clean))
print("regressed code: zero-distance identity holds:", identity_holds(regressed))

assert identity_holds(clean)
assert not identity_holds(regressed)
print("PASS: the zero-distance identity oracle catches the injected regression")
