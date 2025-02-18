"""
evolution.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This file contains a class with a simpler evolution interface.
"""

import numpy as np

from . import matrices
from . import pars
from .matrices import pixelspace, Q2space, passes, get_nfl

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a class to cache matrices and properties into

class Evolver:
    ''' A convenience class for evolving GPDs. Should be initialized as:
        Darwin = Evolver( [optional args] )
    The optional arguments are:
    - nx .......... integer, number of x points; must be even
    - xi .......... float or numpy.array; the xi values to be used
    - nQ2 ......... number of Q2 points; must be at least 4
    - Q2i ......... the initial Q2
    - Q2f ......... the final Q2
    - grid_type ... 1 (linear, default) or 2 (log-linear-log), x spacing type
    Notes:
    - The x space to be used is matrices.pixelspace(nx).
      ...
      1. By default, with grid_type=1, the x values
      are the central values of nx intervals evenly dividing the domain [-1,1].
      For instance, if nx=4, then [-1,1] is divided into the four intervals
      [-1,-0.5], [-0.5,0], [0,0.5] and [0.5,1]. The midpoints of these are
      -0.75, -0.25, 0.25 and 0.75. These four midpoints are the x values given
      by pixelspace(4).
      ...
      2. If grid_type=2, then x values are spaced linearly in the ERBL region,
      but logarithmically in each of the DGLAP regions, with greater
      concentration near the x=+/-xi crossoverpoints.
      nx must be divisible by 4 if grid_type=2 is used.
    - The Q2 space used is matrices.Q2space(nQ2). This is a geomspace,
      except that the charm mass and/or bottom mass may be injected if they
      are between Q2i and Q2f.
    - The minimum nx is 6 because the fifth-order Lagrange interpolation
      used to build the finite elements won't work with smaller nx.
    - It's a good idea to use at least nQ2>=4, just in case the charm and
      bottom masses both need to be injected between Q2i and Q2f.
    '''
    # Matrices up to mc2 (if needed)
    Mevo_VNS_LO_mc2 = None
    Mevo_ANS_LO_mc2 = None
    Mevo_VSG_LO_mc2 = None
    Mevo_ASG_LO_mc2 = None
    Mevo_VNS_NLO_mc2 = None # plus-type NS
    Mevo_ANS_NLO_mc2 = None # minus-type NS
    Mevo_VSG_NLO_mc2 = None
    Mevo_ASG_NLO_mc2 = None
    # Matrices up to mb2 (if needed)
    Mevo_VNS_LO_mb2 = None
    Mevo_ANS_LO_mb2 = None
    Mevo_VSG_LO_mb2 = None
    Mevo_ASG_LO_mb2 = None
    Mevo_VNS_NLO_mb2 = None # plus-type NS
    Mevo_ANS_NLO_mb2 = None # minus-type NS
    Mevo_VSG_NLO_mb2 = None
    Mevo_ASG_NLO_mb2 = None
    # Matrices up to final Q2
    Mevo_VNS_LO_Q2f = None
    Mevo_ANS_LO_Q2f = None
    Mevo_VSG_LO_Q2f = None
    Mevo_ASG_LO_Q2f = None
    Mevo_VNS_NLO_Q2f = None # plus-type NS
    Mevo_ANS_NLO_Q2f = None # minus-type NS
    Mevo_VSG_NLO_Q2f = None
    Mevo_ASG_NLO_Q2f = None
    # Flags
    pass_charm  = False
    pass_bottom = False

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __init__(self,
            nx=100,
            xi=0.5,
            nQ2=17,
            Q2i=pars.mc2+pars.epsilon,
            Q2f=pars.mb2-pars.epsilon,
            grid_type=1
            ):
        # First, we need to ensure that there's an even number of x points,
        # and at least 6 of them
        try:
            assert(2*(nx//2)==nx)
            assert(nx>=6)
        except:
            raise ValueError("nx must be even, not nx={:d}".format(nx))
        # Set up x and xi grids
        if(np.isscalar(xi)):
            xi = np.array([xi])
        self.nx  = nx
        self.nxi = xi.shape[0]
        self.x   = pixelspace(self.nx)
        self.xi  = xi
        # Initialize kernels
        matrices.initialize_kernels(nx, xi, grid_type=grid_type)
        # Set up Q2 grid
        self.Q2i = Q2i
        self.Q2f = Q2f
        self.nQ2 = nQ2
        self.Q2 = Q2space(Q2i, Q2f, nQ2)
        # We use (up to) five flaovrs internally
        self.nfl = 5
        # Set up evolution matrices
        self.makeEvolutionMatrices()
        return

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def evolveGPDs(self, gpd_grid, gpd_type='V', nlo=False):
        ''' Main interface to evolve GPDs, given a gpd_grid in terms of the
        physical flavor basis.
        REQUIRED INPUT:
            - gpd_grid ... np.array with dim(nfl,self.nx,self.nxi,...)
              The x values should be spaced as pixelspace(self.nx).
              Even if there is only one xi value, there must be a xi dimension.
              If the gpd_grid is missing a dimension, I will add an extra
              dimension for the single xi value.
        OPTIONAL INPUT:
            - gpd_type ... either 'V' (default) or 'A',
              for helicity-independent (V) or helicity-dependent (A) evolution.
            - nlo ........ either True or False (default)
        OUTPUT:
            - a GPD grid with dim(nfl=5,self.nx,self.nxi,...,self.nQ2)
              if the highest Q2 is below the charm or bottom mass, the charm
              or bottom will just be zero
        NOTES:
            - The x, xi and Q2 grids are set upon class initialization.
        '''
        if(len(gpd_grid.shape) < 3):
            gpd_grid = gpd_grid[...,np.newaxis]
        # Ensure reasonable inputs
        assert(self.nx  == gpd_grid.shape[1])
        assert(self.nxi == gpd_grid.shape[2])
        assert( (gpd_type=='V') or (gpd_type=='A') )
        # Start a new GPD with an extra dimension for Q2
        new_gpd = np.copy(gpd_grid)[...,np.newaxis]
        # Pad with zeroes if missing any flavors
        nfl = new_gpd.shape[0] - 1
        if(nfl < self.nfl):
            needed_fl = self.nfl - nfl
            new_gpd = add_pad(new_gpd, needed_fl)
        # Evolve through thresholds if needed
        if(self.pass_charm):
            if(gpd_type=='V' and nlo==False):
                M_SG  = self.Mevo_VSG_LO_mc2
                M_NSp = self.Mevo_VNS_LO_mc2
                M_NSm = self.Mevo_VNS_LO_mc2
            elif(gpd_type=='V' and nlo==True):
                M_SG  = self.Mevo_VSG_NLO_mc2
                M_NSp = self.Mevo_VNS_NLO_mc2
                M_NSm = self.Mevo_ANS_NLO_mc2
            elif(gpd_type=='A' and nlo==False):
                M_SG  = self.Mevo_ASG_LO_mc2
                M_NSp = self.Mevo_ANS_LO_mc2
                M_NSm = self.Mevo_ANS_LO_mc2
            elif(gpd_type=='A' and nlo==True):
                M_SG  = self.Mevo_ASG_NLO_mc2
                M_NSp = self.Mevo_ANS_NLO_mc2
                M_NSm = self.Mevo_VNS_NLO_mc2
            added_gpd = evolve_fixed_flavor(new_gpd[...,-1], M_SG, M_NSp, M_NSm, nfl=3)
            new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        if(self.pass_bottom):
            if(gpd_type=='V' and nlo==False):
                M_SG  = self.Mevo_VSG_LO_mb2
                M_NSp = self.Mevo_VNS_LO_mb2
                M_NSm = self.Mevo_VNS_LO_mb2
            elif(gpd_type=='V' and nlo==True):
                M_SG  = self.Mevo_VSG_NLO_mb2
                M_NSp = self.Mevo_VNS_NLO_mb2
                M_NSm = self.Mevo_ANS_NLO_mb2
            elif(gpd_type=='A' and nlo==False):
                M_SG  = self.Mevo_ASG_LO_mb2
                M_NSp = self.Mevo_ANS_LO_mb2
                M_NSm = self.Mevo_ANS_LO_mb2
            elif(gpd_type=='A' and nlo==True):
                M_SG  = self.Mevo_ASG_NLO_mb2
                M_NSp = self.Mevo_ANS_NLO_mb2
                M_NSm = self.Mevo_VNS_NLO_mb2
            added_gpd = evolve_fixed_flavor(new_gpd[...,-1], M_SG, M_NSp, M_NSm, nfl=4)
            new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        # Finish the evolution
        if(gpd_type=='V' and nlo==False):
            M_SG  = self.Mevo_VSG_LO_Q2f
            M_NSp = self.Mevo_VNS_LO_Q2f
            M_NSm = self.Mevo_VNS_LO_Q2f
        elif(gpd_type=='V' and nlo==True):
            M_SG  = self.Mevo_VSG_NLO_Q2f
            M_NSp = self.Mevo_VNS_NLO_Q2f
            M_NSm = self.Mevo_ANS_NLO_Q2f
        elif(gpd_type=='A' and nlo==False):
            M_SG  = self.Mevo_ASG_LO_Q2f
            M_NSp = self.Mevo_ANS_LO_Q2f
            M_NSm = self.Mevo_ANS_LO_Q2f
        elif(gpd_type=='A' and nlo==True):
            M_SG  = self.Mevo_ASG_NLO_Q2f
            M_NSp = self.Mevo_ANS_NLO_Q2f
            M_NSm = self.Mevo_VNS_NLO_Q2f
        nfl = get_nfl(self.Q2f)
        added_gpd = evolve_fixed_flavor(new_gpd[...,-1], M_SG, M_NSp, M_NSm, nfl=nfl)
        new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        return new_gpd

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def evolveSinglet(self, singlet_grid, gluon_grid, gpd_type='V', nlo=False):
        ''' A method that evolves only a singlet quark grid and a gluon grid.
        REQUIRED INPUT:
            - singlet_grid ... np.array with dim(self.nx,self.nxi,...)
              The x values should be spaced as pixelspace(self.nx).
              Even if there is only one xi value, there must be a xi dimension.
              If either GPD grid is missing a dimension, I will add an extra
              dimension for the single xi value.
            - gluon_grid ..... np.array with dim(self.nx,self.nxi,...)
              see description for singlet_grid
        OPTIONAL INPUT:
            - gpd_type ..... either 'V' (default) or 'A', for
              helicity-independent (V) or helicity-dependent (A) evolution.
            - nlo .......... either True or False (default)
        OUTPUT:
            - two GPD grids with dim(self.nx,self.nxi,...,self.nQ2),
              the first singlet and the second gluon.
        NOTES:
            - The x, xi and Q2 grids are set upon class initialization.
        '''
        if(len(singlet_grid.shape) < 2):
            singlet_grid = singlet_grid[...,np.newaxis]
        if(len(gluon_grid.shape) < 2):
            gluon_grid = gluon_grid[...,np.newaxis]
        # Ensure reasonable inputs
        assert((self.nx ==singlet_grid.shape[0]) and (self.nx ==gluon_grid.shape[0]))
        assert((self.nxi==singlet_grid.shape[1]) and (self.nxi==gluon_grid.shape[1]))
        assert( (gpd_type=='V') or (gpd_type=='A') )
        # Concatenated grid
        sg_grid = np.concatenate((singlet_grid, gluon_grid), axis=0)
        # Start a new GPD with an extra dimension for Q2
        new_gpd = np.copy(sg_grid)[...,np.newaxis]
        # Evolve through thresholds if needed
        if(self.pass_charm):
            if(gpd_type=='V' and nlo==False):
                M_SG = self.Mevo_VSG_LO_mc2
            elif(gpd_type=='V' and nlo==True):
                M_SG = self.Mevo_VSG_NLO_mc2
            elif(gpd_type=='A' and nlo==False):
                M_SG = self.Mevo_ASG_LO_mc2
            elif(gpd_type=='A' and nlo==True):
                M_SG = self.Mevo_ASG_NLO_mc2
            added_gpd = np.einsum('ijkl,jk...->ik...l', M_SG, new_gpd[...,-1])
            new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        if(self.pass_bottom):
            if(gpd_type=='V' and nlo==False):
                M_SG = self.Mevo_VSG_LO_mb2
            elif(gpd_type=='V' and nlo==True):
                M_SG = self.Mevo_VSG_NLO_mb2
            elif(gpd_type=='A' and nlo==False):
                M_SG = self.Mevo_ASG_LO_mb2
            elif(gpd_type=='A' and nlo==True):
                M_SG = self.Mevo_ASG_NLO_mb2
            added_gpd = np.einsum('ijkl,jk...->ik...l', M_SG, new_gpd[...,-1])
            new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        # Finish the evolution
        if(gpd_type=='V' and nlo==False):
            M_SG = self.Mevo_VSG_LO_Q2f
        elif(gpd_type=='V' and nlo==False):
            M_SG = self.Mevo_VSG_NLO_Q2f
        elif(gpd_type=='A' and nlo==False):
            M_SG = self.Mevo_ASG_LO_Q2f
        elif(gpd_type=='A' and nlo==True):
            M_SG = self.Mevo_ASG_NLO_Q2f
        added_gpd = np.einsum('ijkl,jk...->ik...l', M_SG, new_gpd[...,-1])
        new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        new_singlet = new_gpd[0:self.nx,...]
        new_gluon = new_gpd[self.nx:2*self.nx,...]
        return new_singlet, new_gluon

    def evolveNS(self, ns_grid, gpd_type='V', nlo=False, ns_type=1):
        ''' A method that evolves only a non-singlet quark grid.
        REQUIRED INPUT:
            - ns_grid ... np.array with dim(self.nx,self.nxi,...)
              The x values should be spaced as pixelspace(self.nx).
              Even if there is only one xi value, there must be a xi dimension.
              If the ns_grid is missing a dimension, I will add an extra
              dimension for the single xi value.
        OPTIONAL INPUT:
            - gpd_type ..... either 'V' (default) or 'A', for
              helicity-independent (V) or helicity-dependent (A) evolution.
            - nlo .......... either True or False (default)
        OUTPUT:
            - a GPD grid with dim(self.nx,self.nxi,...,self.nQ2)
        NOTES:
            - The x, xi and Q2 grids are set upon class initialization.
        '''
        if(len(ns_grid.shape) < 2):
            ns_grid = ns_grid[...,np.newaxis]
        # Ensure reasonable inputs
        assert(self.nx ==ns_grid.shape[0])
        assert(self.nxi==ns_grid.shape[1])
        assert( (gpd_type=='V') or (gpd_type=='A') )
        assert( (ns_type==-1) or (ns_type==1) )
        # Take advantage of V+=A-, V-=A+ to simplify code a bit
        if(ns_type==-1):
            gpd_type = swap_AV(gpd_type)
        # Start a new GPD with an extra dimension for Q2
        new_gpd = np.copy(ns_grid)[...,np.newaxis]
        # Evolve through thresholds if needed
        if(self.pass_charm):
            if(gpd_type=='V' and nlo==False):
                M_NS = self.Mevo_VNS_LO_mc2
            elif(gpd_type=='V' and nlo==True):
                M_NS = self.Mevo_VNS_NLO_mc2
            elif(gpd_type=='A' and nlo==False):
                M_NS = self.Mevo_ANS_LO_mc2
            elif(gpd_type=='A' and nlo==True):
                M_NS = self.Mevo_ANS_NLO_mc2
            added_gpd = np.einsum('ijkl,jk...->ik...l', M_NS, new_gpd[...,-1])
            new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        if(self.pass_bottom):
            if(gpd_type=='V' and nlo==False):
                M_NS = self.Mevo_VNS_LO_mb2
            elif(gpd_type=='V' and nlo==True):
                M_NS = self.Mevo_VNS_NLO_mb2
            elif(gpd_type=='A' and nlo==False):
                M_NS = self.Mevo_ANS_LO_mb2
            elif(gpd_type=='A' and nlo==True):
                M_NS = self.Mevo_ANS_NLO_mb2
            added_gpd = np.einsum('ijkl,jk...->ik...l', M_NS, new_gpd[...,-1])
            new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        # Finish the evolution
        if(gpd_type=='V' and nlo==False):
            M_NS = self.Mevo_VNS_LO_Q2f
        elif(gpd_type=='V' and nlo==True):
            M_NS = self.Mevo_VNS_NLO_Q2f
        elif(gpd_type=='A' and nlo==False):
            M_NS = self.Mevo_ANS_LO_Q2f
        elif(gpd_type=='A' and nlo==True):
            M_NS = self.Mevo_ANS_NLO_Q2f
        added_gpd = np.einsum('ijkl,jk...->ik...l', M_NS, new_gpd[...,-1])
        new_gpd = np.concatenate((new_gpd, added_gpd[...,1:]), axis=-1)
        return new_gpd

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def makeEvolutionMatrices(self):
        ''' Make and save evolution matrices to be used by the class.
        If any thresholds are passed, we'll need to build separate sets of evolution matrices
        between the passed thresholds.
        '''
        self.pass_charm = False
        self.pass_bottom = False
        if(passes(self.Q2,pars.mc2)):
            self.pass_charm = True
            subQ2 = self.Q2[self.Q2 <= pars.mc2]
            self.makeEvolutionMatrices_mc2(subQ2)
        if(passes(self.Q2,pars.mb2)):
            self.pass_bottom = True
            subQ2 = self.Q2[(pars.mc2 <= self.Q2) & (self.Q2 <= pars.mb2)]
            self.makeEvolutionMatrices_mb2(subQ2)
        #
        if(self.pass_bottom):
            subQ2 = self.Q2[pars.mb2 <= self.Q2]
            self.makeEvolutionMatrices_Q2f(subQ2)
        elif(self.pass_charm):
            subQ2 = self.Q2[pars.mc2 <= self.Q2]
            self.makeEvolutionMatrices_Q2f(subQ2)
        else:
            self.makeEvolutionMatrices_Q2f(self.Q2)
        self.initialized = True
        return

    def makeEvolutionMatrices_mc2(self,Q2):
        ''' Make evolution matrices that go only up to mc2,
        if the charm threshold is passed.  '''
        # LO first
        matrices.initialize_evolution_matrices(Q2, nlo=False)
        self.Mevo_VNS_LO_mc2 = matrices.matrix_VNS()
        self.Mevo_ANS_LO_mc2 = matrices.matrix_ANS()
        self.Mevo_VSG_LO_mc2 = matrices.matrix_VSG()
        self.Mevo_ASG_LO_mc2 = matrices.matrix_ASG()
        # NLO next
        matrices.initialize_evolution_matrices(Q2, nlo=True)
        self.Mevo_VNS_NLO_mc2 = matrices.matrix_VNS()
        self.Mevo_ANS_NLO_mc2 = matrices.matrix_ANS()
        self.Mevo_VSG_NLO_mc2 = matrices.matrix_VSG()
        self.Mevo_ASG_NLO_mc2 = matrices.matrix_ASG()
        return

    def makeEvolutionMatrices_mb2(self,Q2):
        ''' Make evolution matrices that go only up to mc2,
        if the botton threshold is passed.  '''
        # LO first
        matrices.initialize_evolution_matrices(Q2, nlo=False)
        self.Mevo_VNS_LO_mb2 = matrices.matrix_VNS()
        self.Mevo_ANS_LO_mb2 = matrices.matrix_ANS()
        self.Mevo_VSG_LO_mb2 = matrices.matrix_VSG()
        self.Mevo_ASG_LO_mb2 = matrices.matrix_ASG()
        # NLO first
        matrices.initialize_evolution_matrices(Q2, nlo=True)
        self.Mevo_VNS_NLO_mb2 = matrices.matrix_VNS()
        self.Mevo_ANS_NLO_mb2 = matrices.matrix_ANS()
        self.Mevo_VSG_NLO_mb2 = matrices.matrix_VSG()
        self.Mevo_ASG_NLO_mb2 = matrices.matrix_ASG()
        return

    def makeEvolutionMatrices_Q2f(self,Q2):
        ''' Make evolution matrices that go up to the target Q2. '''
        # LO first
        matrices.initialize_evolution_matrices(Q2, nlo=False)
        self.Mevo_VNS_LO_Q2f = matrices.matrix_VNS()
        self.Mevo_ANS_LO_Q2f = matrices.matrix_ANS()
        self.Mevo_VSG_LO_Q2f = matrices.matrix_VSG()
        self.Mevo_ASG_LO_Q2f = matrices.matrix_ASG()
        # NLO next
        matrices.initialize_evolution_matrices(Q2, nlo=True)
        self.Mevo_VNS_NLO_Q2f = matrices.matrix_VNS()
        self.Mevo_ANS_NLO_Q2f = matrices.matrix_ANS()
        self.Mevo_VSG_NLO_Q2f = matrices.matrix_VSG()
        self.Mevo_ASG_NLO_Q2f = matrices.matrix_ASG()
        return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helpful routines used by the class


def make_evo_basis(gpd_grid, nfl):
    ''' At a fixed nfl, converts a GPD grid with dimensions (nfl+1,nx,nxi,...),
    with the flavor dimension arranged as (glue,up,down,...) into three grids:
        1. sg_grid  .... dim(2*nx,nxi,...),
                         with the nx Singlet values followed by nx gluon values
        2. t_grid ...... dim(nfl-2,nx,nxi,...),
                         with the first dimension arranged as (T3,T8,...)
        3. qmin_grid ... dim(nfl,nx,nxi),
                         with the first dimension arranged as (umin,dmin,...)
    '''
    gpd_work = np.copy(gpd_grid)
    if(gpd_work.shape[0] < nfl+1):
        needed_fl = nfl+1 - gpd_work.shape[0]
        xshape = tuple(gpd_work.shape[1:])
        shape = (needed_fl,) + xshape
        zero_grid = np.zeros((shape))
        gpd_work = np.concatenate((gpd_work, zero_grid), axis=0)
    # Set up qmin and qpls grids
    qmin_grid = gpd_work[1:nfl+1,...] + np.flip(gpd_work[1:nfl+1,...], axis=1)
    qpls_grid = gpd_work[1:nfl+1,...] - np.flip(gpd_work[1:nfl+1,...], axis=1)
    # Create matrix to maximally decouple qpls thingies
    M = Tpls_matrix(nfl)
    st_grid = np.einsum('ij,j...->i...',M,qpls_grid)
    t_grid = st_grid[1:,...]
    s_grid = st_grid[0,...]
    g_grid = gpd_work[0,...]
    sg_grid = np.concatenate((s_grid, g_grid), axis=0)
    return sg_grid, t_grid, qmin_grid


def make_flv_basis(sg_grid, t_grid, qmin_grid, nfl):
    ''' At fixed nfl, converts three evolution basis GPD grids back to a
    flavor basis GPD grid.
    INPUT:
        1. sg_grid  .... dim(2*nx,nxi,...),
                         with the nx Singlet values followed by nx gluon values
        2. t_grid ...... dim(nfl-2,nx,nxi,...),
                         with the first dimension arranged as (T3,T8,...)
        3. qmin_grid ... dim(nfl,nx,nxi),
                         with the first dimension arranged as (umin,dmin,...)
        4. nfl ......... integer, number of flavors
    OUTPUT:
        - gpd_grid, a GPD  grid with dimensions (nfl+1,nx,nxi,...),
          with the flavor dimension arranged as (glue,up,down,...).
    '''
    nx = t_grid.shape[1]
    # First, build the qpls grid
    M = Tpls_matrix(nfl)
    N = np.linalg.inv(M)
    s_grid = sg_grid[np.newaxis,0:nx,...]
    st_grid = np.concatenate((s_grid, t_grid), axis=0)
    qpls_grid = np.einsum('ij,j...->i...',N,st_grid)
    # Next, the quarks
    q_grid = (qpls_grid + qmin_grid) / 2
    # Finally, the gluons
    g_grid = sg_grid[np.newaxis,nx:2*nx,...]
    gpd_grid = np.concatenate((g_grid, q_grid), axis=0)
    return gpd_grid


def evolve_fixed_flavor(gpd_grid, M_SG, M_NSp, M_NSm, nfl=4):
    ''' A routine that uses given singlet and non-singlet evolution matrices
    to evolve a GPD grid in the physical flavor basis. Deals with decoupling
    to the evolution basis, evolving, and recoupling, assuming a fixed
    (user-supplied) number of flavors.
    '''
    # If gpd_grid has more flavor dimensions than nfl,
    # we need to make sure to pad with zeros when we're done.
    extra_nfl = gpd_grid.shape[0] - 1 - nfl
    # Make evolution basis
    sg_ini, t_ini, qmin_ini = make_evo_basis(gpd_grid, nfl=nfl)
    # Evolve the evolution basis stuff
    sg_fin   = np.einsum('ijkq,jk...->ik...q',   M_SG,  sg_ini)
    t_fin    = np.einsum('ijkq,fjk...->fik...q', M_NSp, t_ini)
    qmin_fin = np.einsum('ijkq,fjk...->fik...q', M_NSm, qmin_ini)
    # back to flavor
    gpd_evo = make_flv_basis(sg_fin, t_fin, qmin_fin, nfl=nfl)
    # pad with zeros if needed
    gpd_evo = add_pad(gpd_evo, extra_nfl)
    return gpd_evo


def Tpls_matrix(nfl):
    ''' A matrix for converting the q-plus type distributions into the singlet
    and T-type distributions. Assumes the q-plus type distributions are arranged
    in a column matrix as (u_pls, d_pls, ...). The matrix will convert this to
    a column matrix with (Sigma, T3, T8, ...).
    To convert back, one only need to invert this matrix.
    '''
    M = np.zeros((nfl,nfl))
    M[0,:] = 1
    for ifl in range(1,nfl):
        M[ifl,0:ifl] = 1
        M[ifl,ifl] = -ifl
    return M


def add_pad(grid, nrows):
    ''' Routine to add nrows rows of zeroes to the 0th dimension of a grid.
    Used to pad GPD grids with zeroes for missing flavors.
    '''
    xshape = tuple(grid.shape[1:])
    shape = (nrows,) + xshape
    zero_grid = np.zeros((shape))
    new_grid = np.concatenate((grid, zero_grid), axis=0)
    return new_grid


def swap_AV(gpd_type):
    if(gpd_type=='A'):
        return 'V'
    if(gpd_type=='V'):
        return 'A'
    raise ValueError("gpd_type must be 'V' or 'A', not '{}'.".format(gpd_type))
