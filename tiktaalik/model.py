"""
model.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This implements a model due to Goloshkokov and Kroll.
The specific reference consulted is was:
  P. Kroll, H. Moutarde, F. Sabatie
  European Physical Journal C (2013) 73:2278
  arxiv:1210.6975
  Kroll:2012sm
see /f90src/model/gk.f90 for more details.
"""

import numpy as np
from .f90wrap.model import dummy as f90src

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model GPDs

def H_singlet(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    f = (
              Hu(x,xi,t) - Hu(-x,xi,t)
            + Hd(x,xi,t) - Hd(-x,xi,t)
            + Hs(x,xi,t) - Hs(-x,xi,t)
            )
    return f

def Hu(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hu_wrap(x,xi,t)
    return H

def Hd(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hd_wrap(x,xi,t)
    return H

def Hs(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hs_wrap(x,xi,t)
    return H

def Hg(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hg_wrap(x,xi,t)
    return H

def Hu_val(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hu_val_wrap(x,xi,t)
    return H

def Hd_val(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hd_val_wrap(x,xi,t)
    return H

def Hu_sea(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hu_sea_wrap(x,xi,t)
    return H

def Hd_sea(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    H = f90src.hd_sea_wrap(x,xi,t)
    return H
