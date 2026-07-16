"""
continuum.py

part of the tiktaalik package for GPD evolution
by Adam Freese

Continuum quantities to use in benchmarks of the finite element method.
"""

import numpy as np
from .f90wrap.testing import dummy as f90src

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Continuum CFFs

def cff_q(xi, Q2, nlo=False):
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_cff_q(xi,Q2,nlo)
    return V

def cff_g(xi, Q2, nlo=False):
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_cff_g(xi,Q2,nlo)
    return V

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Continuum shifts (V-type)

def shift_cNS(x, xi, Q2, nlo=False, nstype=1):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cns(x,xi,Q2,nlo,nstype)
    return V

def shift_cQQ(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cqq(x,xi,Q2,nlo)
    return V

def shift_cQG(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cqg(x,xi,Q2,nlo)
    return V

def shift_cGQ(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cgq(x,xi,Q2,nlo)
    return V

def shift_cGG(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cgg(x,xi,Q2,nlo)
    return V

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Continuum shifts (A-type)

def shift_cNS_tilde(x, xi, Q2, nlo=False, nstype=1):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cns_tilde(x,xi,Q2,nlo,nstype)
    return V

def shift_cQQ_tilde(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cqq_tilde(x,xi,Q2,nlo)
    return V

def shift_cQG_tilde(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cqg_tilde(x,xi,Q2,nlo)
    return V

def shift_cGQ_tilde(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cgq_tilde(x,xi,Q2,nlo)
    return V

def shift_cGG_tilde(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cgg_tilde(x,xi,Q2,nlo)
    return V
