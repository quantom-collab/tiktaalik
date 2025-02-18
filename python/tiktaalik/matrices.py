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
    ''' The x space used by tiktaalik.
    1. grid_type=1, linear grid (default)
       The x values are the central values of nx intervals evenly dividing
       the domain [-1,1]. For instance, if nx=4, then [-1,1] is divided into the
       intervals [-1,-0.5], [-0.5,0], [0,0.5] and [0.5,1]. The midpoints of these
       are -0.75, -0.25, 0.25 and 0.75. These four midpoints are used as x values.
       Independent of xi.
    2. grid_type=2, log-linear-log
       x is broken down into the two DGLAP regions (x < -xi and x > xi) and the
       ERBL region (-xi < x < xi). nx/4 points are placed in each of the DGLAP
       regions, and are geometrically spaced, more closely around the x=-xi or
       x=xi endpoint. The other nx/2 points are linearly spaced in the ERBL region.
       The placement of the points in each region is done using midpoints
       (geometric rather than arithmetic in the DGLAP region) of interals,
       similarly to the linear grid.
    '''
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
    ''' Routine to check if m2 is located between the minimum and maximum
    values in Q2_array. Assumes Q2_array is monotonically ordered.
    Used to check whether a mass threshold is passed during evolution.
    '''
    if(Q2_array[0] < m2-pars.epsilon and Q2_array[-1] > m2+pars.epsilon):
        return True
    return False

def get_nfl(Q2):
    ''' Gets the effective number of flavors. Returns either 3, 4 or 5,
    under the assumption that Q2 below the strange threshold or above the
    truth threshold is never passed.
    This routine  is used by the Evolver class.
    '''
    if(Q2 < pars.mc2):
        return 3
    elif(Q2 < pars.mb2):
        return 4
    return 5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to create evolution matrices

def matrix_VNS(ns_type=1):
    ''' Evolution matrix for *non-singlet* Q->Q, helicity-independent (V-type).
    The ns_type parameter must be 1 or -1, and distinguishes between
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar.
    The ns_type parameter only matters for NLO.
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_vns_wrap(nx, nxi, nQ2, ns_type)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_VSG():
    ''' Evolution matrix for singlet sector, helicity-independent (V-type).
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_vsg_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_ANS(ns_type=1):
    ''' Evolution matrix for *non-singlet* Q->Q, helicity-dependent (A-type).
    The ns_type parameter must be 1 or -1, and distinguishes between
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar.
    The ns_type parameter only matters for NLO.
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_ans_wrap(nx, nxi, nQ2, ns_type)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_ASG():
    ''' Evolution matrix for singlet sector, helicity-dependent (A-type).
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

def kernel_VQQ(Q2=pars.mc2, nfl=4, nlo=False, ns_type=1):
    ''' Evolution kernel for Q->Q, helicity-independent (V-type).
    The ns_type parameter can be -1, 0, or 1. The +1 and -1 values are for
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar, respectively.
    The 0 value is for the singlet Q->Q kernel.
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO too.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqq_wrap(Q2, nx, nxi, nfl, nlo, ns_type)
    return K

def kernel_VQG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->Q, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGQ(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for Q->G, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vgq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->G, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vgg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AQQ(Q2=pars.mc2, nfl=4, nlo=False, ns_type=1):
    ''' Evolution kernel for Q->Q, helicity-dependent (A-type).
    The ns_type parameter can be -1, 0, or 1. The +1 and -1 values are for
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar, respectively.
    The 0 value is for the singlet Q->Q kernel.
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO too.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqq_wrap(Q2, nx, nxi, nfl, nlo, ns_type)
    return K

def kernel_AQG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->Q, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGQ(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for Q->G, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->G, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agg_wrap(Q2, nx, nxi, nfl, nlo)
    return K
