"""
matrices.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This file contains methods to construct evolution matrices.
"""

import numpy as np
from .f90wrap.matrices import dummy as f90src

from . import pars

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A global dictionary containing information about what's been initialized.
# It's initialized with some sane default values.
# This should be treated like a hidden variable, and only modified by methods
# in this file.

_matrix_dict = {
        'nQ2' : 2,
        'Q2'  : np.array([pars.mc2, pars.mb2]),
        'nxi' : 5,
        'xi'  : np.linspace(0.1, 0.5, 5),
        'nx'  : 80,
        'grid_type' : 1,
        'lagrange_order' : 6,
        'kernel_init' : False,
        'evomat_init' : False,
        'wilson_init' : False,
        'nlo_evo' : False
        }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The autmatically-run initialization routine, called when tiktaalik starts

def first_initialization():
    ''' Method is run when tiktaalik is imported. This ensures some reasonable
    x, xi and Q2 arrays exist in memory.
    '''
    f90src.initialize_x_xi_wrap(_matrix_dict['nx'],
                                _matrix_dict['xi'],
                                _matrix_dict['grid_type'],
                                _matrix_dict['lagrange_order']
                                )
    f90src.initialize_q2_wrap(_matrix_dict['Q2'])
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for the user to retrieve the x, xi and Q2 grids in use

def get_x_grid():
    ''' Returns the x grid in use.
    The returned grid is two-dimensional, with dimensions (nx,nxi).
    The reason for this is that, depending on the grid type,
    the x spacing may be xi-dependent.
    Thus, this is effectively an array of nxi x arrays.
    '''
    nx  = _matrix_dict['nx']
    nxi = _matrix_dict['nxi']
    x = f90src.get_x_wrap(nx, nxi)
    return x

def get_xi_array():
    nxi = _matrix_dict['nxi']
    xi = f90src.get_xi_wrap(nxi)
    return xi

def get_Q2_array():
    nQ2 = _matrix_dict['nQ2']
    Q2 = f90src.get_q2_wrap(nQ2)
    return Q2

def get_nx():
    nx = _matrix_dict['nx']
    return nx

def get_nxi():
    nxi = _matrix_dict['nxi']
    return nxi

def get_nQ2():
    nQ2 = _matrix_dict['nQ2']
    return nQ2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for the user to set x, xi and Q2 grids

def set_x_xi_grids(nx, xi, grid_type, lagrange_order=5):
    ''' Set the x and xi grids used internally to construct kernel and evolution
    matrices.
    ------
    Input:
        - nx (integer, >= 6)
            number of x points
        - xi (float or np.array)
            xi points
        - grid_type (integer)
            specify how the x grid should be constructed
            grid_type=1 for linearly-spaced grid
            grid_type=2 for hybrid log-linear grid with x=+/-xi on the grid
        - lagrange_order (integer, optional)
            order of piecewise lagrange interpolation
    '''
    assert(nx>=6)
    # Cache the data sent by the user
    _matrix_dict['nx']  = nx
    if(np.isscalar(xi)):
        xi = np.array([xi])
    _matrix_dict['nxi'] = xi.shape[0]
    _matrix_dict['xi']  = xi
    _matrix_dict['grid_type'] = grid_type
    # Set the x and xi grids in the Fortran code
    f90src.initialize_x_xi_wrap(nx, xi, grid_type, lagrange_order)
    # Mark all matrices as uninitialized, since their shape may now be invalid
    _matrix_dict['kernel_init'] = False
    _matrix_dict['evomat_init'] = False
    _matrix_dict['wilson_init'] = False
    return

def set_Q2_grid(Q2):
    ''' Set the Q2 grid used internally to construct evolution matrices.
    ------
    Input:
        - Q2 (np.array, size >= 2)
            specific Q2 points
    '''
    assert(Q2.shape[0] >= 2)
    # Cache the data sent by the user
    _matrix_dict['nQ2'] = Q2.shape[0]
    _matrix_dict['Q2']  = Q2
    # Set the Q2 grid in the Fortran code
    f90src.initialize_q2_wrap(Q2)
    # Mark matrices as uninitialized that depend on the Q2 shape/range
    _matrix_dict['evomat_init'] = False
    _matrix_dict['wilson_init'] = False
    return

def do_lo_evolution():
    ''' Tell tiktaalik to do evolution at leading order (LO). '''
    _matrix_dict['nlo_evo'] = False
    _matrix_dict['evomat_init'] = False
    return

def do_nlo_evolution():
    ''' Tell tiktaalik to do evolution at next-to-leading order (NLO). '''
    _matrix_dict['nlo_evo'] = True
    _matrix_dict['evomat_init'] = False
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Method to get an interpixel

def interpixel(i_pixel, x, xi):
    ''' Given arbitrary x (potentially off of the grid),
    return the ith interpixel of the current x grid.
    i_pixel uses 0-indexing, so 0 <= i_pixel < n_pixels.
    '''
    n_pixels = f90src.get_nx_wrap()
    assert(i_pixel >= 0)
    assert(i_pixel < n_pixels)
    grid_type = _matrix_dict['grid_type']
    return f90src.interpixel_wrap(n_pixels, i_pixel, x, xi, grid_type)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to obtain evolution kernel matrices

def kernel_VQQ(Q2=pars.mc2, nfl=4, nlo=False, ns_type=1):
    ''' Evolution kernel for Q->Q, helicity-independent (V-type).
    The ns_type parameter can be -1, 0, or 1. The +1 and -1 values are for
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar, respectively.
    The 0 value is for the singlet Q->Q kernel.
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO too.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    _check_kernel_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqq_wrap(Q2, nx, nxi, nfl, nlo, ns_type)
    return K

def kernel_VQG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->Q, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    _check_kernel_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGQ(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for Q->G, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    _check_kernel_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vgq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->G, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    _check_kernel_initialization()
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
    _check_kernel_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqq_wrap(Q2, nx, nxi, nfl, nlo, ns_type)
    return K

def kernel_AQG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->Q, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    _check_kernel_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGQ(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for Q->G, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    _check_kernel_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->G, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    _check_kernel_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to obtain evolution matrices

def matrix_VNS(ns_type=1):
    ''' Evolution matrix for *non-singlet* Q->Q, helicity-independent (V-type).
    The ns_type parameter must be 1 or -1, and distinguishes between
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar.
    The ns_type parameter only matters for NLO.
    '''
    _check_kernel_initialization()
    _check_evomat_initialization()
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
    _check_kernel_initialization()
    _check_evomat_initialization()
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
    _check_kernel_initialization()
    _check_evomat_initialization()
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
    _check_kernel_initialization()
    _check_evomat_initialization()
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_asg_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to obtain Wilson coefficient matrices

def dvcs_Cq(nlo=False):
    # TODO: docstring
    # Get the matrix dimensions
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    # Make sure the Wilson coefficient matrices are initialized
    if(_matrix_dict['wilson_init']==False):
        f90src.make_wilson_wrap(nQ2)
        _matrix_dict['wilson_init'] = True
    # Get the coefficients
    C = f90src.dvcs_cq_wrap(nx, nxi, nQ2, nlo)
    return C

def dvcs_Cg(nlo=False):
    # TODO: docstring
    # Get the matrix dimensions
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    # Make sure the Wilson coefficient matrices are initialized
    if(_matrix_dict['wilson_init']==False):
        f90src.make_wilson_wrap(nQ2)
        _matrix_dict['wilson_init'] = True
    # Get the coefficients
    C = f90src.dvcs_cg_wrap(nx, nxi, nQ2, nlo)
    return C

def dvcs_Ctilq(nlo=False):
    # TODO: docstring
    # Get the matrix dimensions
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    # Make sure the Wilson coefficient matrices are initialized
    if(_matrix_dict['wilson_init']==False):
        f90src.make_wilson_wrap(nQ2)
        _matrix_dict['wilson_init'] = True
    # Get the coefficients
    C = f90src.dvcs_ctilq_wrap(nx, nxi, nQ2, nlo)
    return C

def dvcs_Ctilg(nlo=False):
    # TODO: docstring
    # Get the matrix dimensions
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    # Make sure the Wilson coefficient matrices are initialized
    if(_matrix_dict['wilson_init']==False):
        f90src.make_wilson_wrap(nQ2)
        _matrix_dict['wilson_init'] = True
    # Get the coefficients
    C = f90src.dvcs_ctilg_wrap(nx, nxi, nQ2, nlo)
    return C

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mellin moments

def mellin(s):
    # TODO: docstring
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    M = f90src.mellin_wrap(nx, nxi, s)
    return M

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q2 related methods

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
# to deprecate

#def pixelspace(nx, xi=0.5, grid_type=1):
#    # TODO: replace by returning cached x space
#    ''' The x space used by tiktaalik.
#    1. grid_type=1, linear grid (default)
#       The x values are the central values of nx intervals evenly dividing
#       the domain [-1,1]. For instance, if nx=4, then [-1,1] is divided into the
#       intervals [-1,-0.5], [-0.5,0], [0,0.5] and [0.5,1]. The midpoints of these
#       are -0.75, -0.25, 0.25 and 0.75. These four midpoints are used as x values.
#       Independent of xi.
#    2. grid_type=2, log-linear-log
#       x is broken down into the two DGLAP regions (x < -xi and x > xi) and the
#       ERBL region (-xi < x < xi). nx/4 points are placed in each of the DGLAP
#       regions, and are geometrically spaced, more closely around the x=-xi or
#       x=xi endpoint. The other nx/2 points are linearly spaced in the ERBL region.
#       The placement of the points in each region is done using midpoints
#       (geometric rather than arithmetic in the DGLAP region) of interals,
#       similarly to the linear grid.
#    '''
#    x = f90src.pixelspace_wrap(nx, xi, grid_type)
#    return x

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Internal helper routines

def _check_kernel_initialization():
    ''' Check if kernels are initialized, and initialize them if they're not.
    Also flag them as initialized after.
    '''
    if(_matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        _matrix_dict['kernel_init'] = True
    return

def _check_evomat_initialization():
    ''' Check if evolution matrices are initialized, and initialize them if
    they're not. Also flag them as initialized after.
    '''
    if(_matrix_dict['evomat_init']==False):
        Q2 = _matrix_dict['Q2']
        nlo = _matrix_dict['nlo_evo']
        f90src.make_matrices_wrap(Q2, nlo)
        _matrix_dict['evomat_init'] = True
    return
