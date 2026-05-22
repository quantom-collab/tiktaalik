"""
polynomiality.py

part of the tiktaalik package for GPD evolution
by Adam Freese

Methods to obtain coefficients of polynomiality relations.
"""

import numpy as np
from scipy.optimize import curve_fit

from . import matrices

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_polynomial(H, s, force_even=True):
    # Currently requires H(x,xi), without t or Q2 dependence
    # TODO: generalize
    Ms = matrices.mellin(s)
    raw_poly = np.einsum('zx,xz->z', Ms, H)
    xi = matrices.get_xi_array()
    if(force_even):
        s_eff = 2*(s//2)
        n = np.linspace(0, s_eff, s_eff//2+1)
        c, cov = fit_even_polynomial(xi, raw_poly, s_eff)
        fit_poly = even_polynomial(xi, s_eff, c)
    else:
        n = np.linspace(0, s, s+1)
        c, cov = fit_polynomial(xi, raw_poly, s)
        fit_poly = polynomial(xi, s, c)
    print("n: ", n)
    print("c: ", c)
    data = {
            'c' : c,
            'cov' : cov,
            'raw polynomial' : raw_poly,
            'fit polynomial' : fit_poly,
            'n' : n
            }
    return data


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def polynomial(x, N, c):
    ''' Evaluates a polynomial in x of maximum degree N.
    Input:
    - x ... float or np.array of floats
    - N ... integer
    - c ... length N+1 array of floats
    Output:
    - float or np.array of floats (same as x)
    '''
    f = 0*x
    for i in range(N+1):
        f += c[i] * x**i
    return f

def even_polynomial(x, N, c):
    ''' Evaluates an polynomial in x of maximum degree N.
    Input:
    - x ... float or np.array of floats
    - N ... even integer
    - c ... length N//2+1 array of floats
    Output:
    - float or np.array of floats (same as x)
    '''
    assert(_is_even(N))
    f = 0*x
    for i in range(N//2+1):
        f += c[i] * x**(2*i)
    return f

def fit_polynomial(x, y, N):
    ''' Fits a collection of data to a degree-N polynomial.
    Input:
    - x ... np.array of floats
    - y ... np.array of floats
    - N ... integer
    Output:
    - c ..... length N+1 array of polynomial coefficients
    - cov ... (N+1)x(N+1) covariance matrix
    '''
    f = lambda x, *p0: polynomial(x, N, p0)
    p0 = np.zeros(N+1)
    c, cov = curve_fit(f, x, y, p0=p0)
    return c, cov

def fit_even_polynomial(x, y, N):
    ''' Fits a collection of data to a degree-N even polynomial.
    Input:
    - x ... np.array of floats
    - y ... np.array of floats
    - N ... even integer
    Output:
    - c ..... length N//2+1 array of polynomial coefficients
    - cov ... (N//2+1)x(N//2+1) covariance matrix
    '''
    assert(_is_even(N))
    f = lambda x, *p0: even_polynomial(x, N, p0)
    p0 = np.zeros(N//2+1)
    c, cov = curve_fit(f, x, y, p0=p0)
    return c, cov

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def _is_even(N):
    return (N//2)*2 == N
