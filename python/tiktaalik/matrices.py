"""
matrices.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This file contains methods to construct evolution matrices.
"""

import numpy as np
from matrices_wrap import dummy as f90src

from . import pars

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All-important initialization routines.

def initialize_kernels(nx, xi, grid_type=1):
    ''' Initializes all the kernels to be used to make evolution matrices.
    This **MUST** be called before any other methods:
    - Before methods to get kernels.
    - Before methods to get evolution matrices.
    - Before methods to initialize evolution matrices.
    INPUT:
    - nx .......... integer, number of x points.
                    Must be even, and at least 6.
                    Spacing depends on grid_type
    - xi .......... numpy.array, with the xi points to be used; if a scalar
                    is passed, it will be turned into a one-component array.
    - grid_type ... integer, specifying the grid type.
                    1 : nx linearly spaced points from -1+1/nx to 1-1/nx
                    2 : nx/4 log-spaced points in each of the DGLAP regions,
                        and nx/2 linearly-spaced points in the ERBL region
    OUTPUT:
    - None
    NOTES:
    - If the user wishes to change either nx or xi, just call this method again.
      Calling this method deallocates any existing evolution matrices.
    '''
    assert(nx==(nx//2)*2)
    assert(nx>=6)
    if(np.isscalar(xi)):
        xi = np.array([xi])
    nxi = xi.shape[0]
    f90src.make_kernels_wrap(nx, xi, grid_type)
    return

def initialize_evolution_matrices(Q2, nlo=False):
    ''' Initializes the evolution matrices.
    This **MUST** be called before any method to get the evolution matrices.
    Moreover, initialize_kernels (above) **MUST** be called first.
    nx and nxi should be consistent with those used to initialize the kernels.
    INPUT:
    - Q2 .... numpy.array, with *at least* two values
         first value is the Q2 value at which evolutuon starts
    - nlo ... boolean, False by default; whether to include NLO corrections.
    OUTPUT:
    - None
    NOTES:
    - If the user wishes to change the Q2 array, this must be called again.
    - If the user wishes to change the nx or xi, initialize_kernels
      must be called to do this.
    - If the user wishes to toggle NLO, this must be called again.
    - If initialize_kernels has been called, any initialized evolution matrices
      have been erased and this routine must be called again.
    '''
    nQ2 = Q2.shape[0]
    assert(nQ2 >= 2)
    f90src.make_matrices_wrap(Q2, nlo)
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spaces used

def pixelspace(nx, xi=0.5, grid_type=1):
    ''' The x space used by tiktaalik. '''
    # TODO EXPLAIN
    #x_end = np.linspace(-1, 1, nx+1)
    #x_min = (x_end[0]   +x_end[1] )/2
    #x_max = (x_end[nx-1]+x_end[nx])/2
    #return np.linspace(x_min,x_max,nx)
    x = f90src.pixelspace_wrap(nx, xi, grid_type)
    return x

def Q2space(Q2i, Q2f, nQ2):
    ''' Creates something that's mostly a geomspace,
    but injects any quark mass thresholds if they're present.
    '''
    nThresh = int(Q2i < pars.mc2 and Q2f > pars.mc2) + int(Q2i < pars.mb2 and Q2f > pars.mb2)
    Q2 = np.geomspace(Q2i, Q2f, nQ2-nThresh)
    if(Q2i < pars.mc2 and Q2f > pars.mc2):
        Q2 = np.append(Q2, pars.mc2)
    if(Q2i < pars.mb2 and Q2f > pars.mb2):
        Q2 = np.append(Q2, pars.mb2)
    Q2 = np.sort(Q2)
    return Q2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q2 related methods

def passes(Q2_array, m2):
    # TODO DOCSTRING
    if(Q2_array[0] < m2-pars.epsilon and Q2_array[-1] > m2+pars.epsilon):
        return True
    return False

def get_nfl(Q2):
    # TODO DOCSTRING
    if(Q2 < pars.mc2):
        return 3
    elif(Q2 < pars.mb2):
        return 4
    return 5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to create evolution matrices

def matrix_VNS():
    # TODO DOCSTRING
    ''' Evolution matrix for *non-singlet* Q->Q, helicity-independent (V-type).
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_vns_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_VSG():
    # TODO DOCSTRING
    ''' Evolution matrix for singlet sector, helicity-independent (V-type).
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_vsg_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_ANS():
    # TODO DOCSTRING
    ''' Evolution matrix for *non-singlet* Q->Q, helicity-dependent (A-type).
    The pom parameter must be 1 or -1, and distinguishes between
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar.
    The pom parameter only matters for NLO.
    Note that nxi = nx//2. This also requires nx to be even.
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_ans_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_ASG():
    # TODO DOCSTRING
    ''' Evolution matrix for singlet sector, helicity-dependent (A-type).
    The column matrix that this evolves should consist of nx quark singlet
    GPD values, followed by nx gluon GPD values.
    Note that nxi = nx//2. This also requires nx to be even.
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_asg_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to create kernels

def kernel_VQQ(Q2=pars.mc2, nfl=4, nlo=False, nstype=1):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqq_wrap(Q2, nx, nxi, nfl, nlo, nstype)
    return K

def kernel_VQG(Q2=pars.mc2, nfl=4, nlo=False):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGQ(Q2=pars.mc2, nfl=4, nlo=False):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vgq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGG(Q2=pars.mc2, nfl=4, nlo=False):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vgg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AQQ(Q2=pars.mc2, nfl=4, nlo=False, nstype=1):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqq_wrap(Q2, nx, nxi, nfl, nlo, nstype)
    return K

def kernel_AQG(Q2=pars.mc2, nfl=4, nlo=False):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGQ(Q2=pars.mc2, nfl=4, nlo=False):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGG(Q2=pars.mc2, nfl=4, nlo=False):
    # TODO DOCSTRING
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TESTING AREA (TODO REMOVE BEFORE RELEASE)

def dilog(x):
    if(np.isscalar(x)):
        x = np.array([x])
    return f90src.dilog_wrap(x)
