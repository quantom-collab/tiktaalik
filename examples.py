"""
examples.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This file contains examples for how to use the tiktaalik package.
"""

# First, import tiktaalik
import tiktaalik as tk

# We also require numpy
import numpy as np

################################################################################
print("="*80)
print("First example section...")

# Suppose we have a model GPD we want to evolve.
# tiktaalik comes with the GK model as an example to use.
# We need to choose some xi values and a number of x values.
# Let's start with one xi value, x=0.5.
xi = 0.5

# Currently, tiktaalik requires the x grid to be spaced in a peculiar way,
# which we call a "pixelspace".
# Basically, the x values are the central values of nx intervals evenly dividing
# the domain [-1,1]. For instance, if nx=4, then [-1,1] is divided into the four
# intervals [-1,-0.5], [-0.5,0], [0,0.5] and [0.5,1]. The midpoints of these are
# -0.75, -0.25, 0.25 and 0.75. These four midpoints are the pixelspace values.

# We can make a pixelspace using the tiktaalik package.
# Let's make one with 20 points to start.
# Note that tiktaalik *requires* the number of x points to be even.
nx = 20
x = tk.matrices.pixelspace(nx)

# Let's create a model GPD we want to evolve
# We'll use the tk-supplied GK model code.
# However, as long as we make sure to use a pixelspace for x,
# we can use any GPD function.
# Let's create a non-singlet combination so it evolves by itself.
# The Hu-minus distribution will do.
umin0 = tk.model.Hu(x, xi, t=0) + tk.model.Hu(-x, xi, t=0)

# Note that the GPD we made has three dimensions, the last one for t.

# Now, we want to make an evolution matrix.
# Before we can do this, we need to call the routines to initialize
# the kernel and evolution matrices that tiktaalik's Fortran code holds in memory.
# (These are held in memory to avoid needing to frequent rebuildings.
# I'll explain this in a little bit.)
# The kernels need to be initialized with info about the number of x points
# and the specific collection of xi values we're using.
tk.matrices.initialize_kernels(nx=nx, xi=xi)

# Next, we need to initialize the evolution matrices.
# We need to come up with some Q2 array.
# The GK model is designed to work at Q2 = 4 GeV**2.
# Let's say we want to evolve to 5, 6 and 7 GeV**2.
# Then we need an array with all four values:
Q2 = np.array([4, 5, 6, 7]) # units of GeV**2
tk.matrices.initialize_evolution_matrices(Q2)

# Now we're allowed to get the evolution matrices!
# Let's get the non-singlet, helicity-independent one.
M_NS = tk.matrices.matrix_VNS()

# The matrix will have dimensions (nx, nx, nxi, nQ2).
# Let's check!
print("Shape of M_NS: ", M_NS.shape)
# Note that the GPD we're using has a shape (nx,nxi,nt).
print("Shape of umin0: ", umin0.shape)
# To get an evolved GPD, we use np.einsum
umin_evo = np.einsum('ijkl,jk...->ik...l', M_NS, umin0)
# The ellipsis makes sure any extra dimensions (t in this case)
# are placed in-between the (x,xi) and (Q2) dimensions.
print("Shape of umin_evo: ", umin_evo.shape)

# We can reuse the same matrix to evolve different non-singlet distributions.
dmin0 = tk.model.Hd(x, xi, t=0) + tk.model.Hd(-x, xi, t=0)
dmin_evo = np.einsum('ijkl,jk...->ik...l', M_NS, dmin0)

# We can also do singlet evolution, but for this we need to concatenate
# the singlet quark and gluon arrays.
sing0 = tk.model.H_singlet(x, xi, t=0)
glue0 = tk.model.Hg(x, xi, t=0)
sg0 = np.concatenate((sing0, glue0), axis=0)
# We grab the singlet-gluon matrix, which has dimensions (2*nx, 2*nx, nxi, nQ2).
M_SG = tk.matrices.matrix_VSG()
print("Shape of M_SG: ", M_SG.shape)
print("Shape of sg0: ", sg0.shape)
# The dimensions match, so we can evolve as usual.
sg_evo = np.einsum('ijkl,jk...->ik...l', M_SG, sg0)
# And the singlet and gluon GPDs can be grabbed by slicing the array
sing_evo = sg_evo[0:nx,...]
glue_evo = sg_evo[nx:2*nx,...]

################################################################################
print("="*80)
print("Second example section...")

# nx=20 is a bit small, and leads to inaccurate evolution.
# Let's do better.
nx = 60
x = tk.matrices.pixelspace(nx)
# We could also evolve at a few different xi values simultaneously.
xi = np.array([0.1, 0.2, 0.3])
# We need to reinitialize the kernels
tk.matrices.initialize_kernels(nx, xi)

# Since we reinitialized the kernels, the evolution matrices have been freed.
# (After all, their dimensions did depend on nx and xi.)
# We need to initialzie them again.
# Let's stick with the same Q2 array as before.
Q2 = np.array([4, 5, 6, 7]) # units of GeV**2
tk.matrices.initialize_evolution_matrices(Q2)

# We should get the matrices again.
M_NS = tk.matrices.matrix_VNS()
M_SG = tk.matrices.matrix_VSG()
print("Shape of M_NS: ", M_NS.shape)
print("Shape of M_SG: ", M_SG.shape)
# We can see from the output that these matrices are bigger (nx=60, nxi=3, nQ2=4)

# Let's also get some GK model GPDs to evolve.
# The GK model also depends on t, so let's get five t values
t = np.array([0, -0.1, -0.2, -0.3, -0.4]) # units of GeV**2
umin0 = tk.model.Hu(x, xi, t) + tk.model.Hu(-x, xi, t)
dmin0 = tk.model.Hd(x, xi, t) + tk.model.Hd(-x, xi, t)
sing0 = tk.model.H_singlet(x, xi, t)
glue0 = tk.model.Hg(x, xi, t)

# Evolution works exactly as before.
# Since GPD evolution is independent of t, the evolution of all five
# t values is done simultaneously.
# (We still have to concatenate and slice to do the singlet-gluon evolution.)
umin_evo = np.einsum('ijkl,jk...->ik...l', M_NS, umin0)
dmin_evo = np.einsum('ijkl,jk...->ik...l', M_NS, dmin0)
sg0 = np.concatenate((sing0, glue0), axis=0)
sg_evo = np.einsum('ijkl,jk...->ik...l', M_SG, sg0)
sing_evo = sg_evo[0:nx,...]
glue_evo = sg_evo[nx:2*nx,...]

# We can see the (nx,nxi,nt) shape is preserved---we just added a Q2 dimension.
print("Shape of umin_evo: ", umin_evo.shape)

# Now, for the reason tiktaalik stores things in memory.
# The most expensive part of the code is computing the kernels.
# If we change our mind about the Q2 grid we want...
Q2 = np.array([4, 5, 6, 7, 8, 9]) # units of GeV**2
# ... then we only need to reinitialize the evolution matrices.
tk.matrices.initialize_evolution_matrices(Q2)
# The same kernels we already computed will be used. This saves time!

################################################################################
print("="*80)
print("Third example section...")

# tiktaalik also offers an Evolver class to make things easy.
# It basically takes care of building matrices, and also does conversion
# between flavor and evolution bases, as well as dealing with flavor thresholds.
nx = 60
x = tk.matrices.pixelspace(nx)
xi = np.array([0.1, 0.2, 0.3])
Q2i = 4 # GeV**2
Q2f = 25 # GeV**2 ... note that we pass the bottom mass threshold
nQ2 = 8
# The Q2 grid used internally by the Evolver class is a custom Q2space:
Q2 = tk.matrices.Q2space(Q2i, Q2f, nQ2)
# ... which is almost a geomspace, except that the charm and/or bottom mass
# is injected if either threshold is passed between Q2i and Q2f.
# To be sure, the total number of points will be nQ2.
# Just in case both mc2 and mb2 are passed, it's a good idea to ensure
# nQ2 is at least 4.

Darwin = tk.Evolver(nx=nx, xi=xi, Q2i=Q2i, Q2f=Q2f, nQ2=nQ2)

# At the model scale, the GK model has gluons, up, and down quarks.
t = np.array([0, -0.1, -0.2, -0.3, -0.4]) # units of GeV**2
g0 = tk.model.Hg(x,xi,t)
u0 = tk.model.Hu(x,xi,t)
d0 = tk.model.Hd(x,xi,t)
s0 = tk.model.Hs(x,xi,t)
# The tk.Evolver can evolve a grid of flavor-basis GPDs,
# with the first axis being a flavor number:
# 0 for glue, 1 for up, 2 for down, etc
# So let's make a grid:
gpd_grid_ini = np.concatenate(
        (g0[np.newaxis,...], u0[np.newaxis,...],
            d0[np.newaxis,...], s0[np.newaxis,...]), axis=0
        )
# The shape is (4, nx, nxi, nt), with 4 for g+u+d+s
print("Shape of gpd_grid_ini: ", gpd_grid_ini.shape)

# The evolved grid is simply:
gpd_evo = Darwin.evolveGPDs(gpd_grid_ini)
print("Shape of gpd_evo: ", gpd_evo.shape)

# We can also do:
gpd_evo = Darwin.evolveGPDs(gpd_grid_ini, gpd_type='V')
# gpd_type can be V for helicity-independent (gamma-plus correlator)
# or A for helicity-dependent (gamma-plus * gamma5 correlator).
# If gpd_type isn't specified, I assume it's helicity-independent.

# The Evolver class also contains easy routines to do singlet or
# non-singlet evolution only.
# This allows us to do evolution without manually calling kernel and matrix
# initialization routine, without manually concatenating and slicing, etc.
# For singlet, for instance:
Sigma0 = u0 + d0 + s0 - np.flip(u0+d0+s0, axis=0)
Sigma_evo, g_evo = Darwin.evolveSinglet(Sigma0, g0)
print("Shape of Sigma0: ", Sigma0.shape)
print("Shape of Sigma_evo: ", Sigma_evo.shape)

# For another exaple, let's do easy evolution of non-singlet
umin0 = u0 + np.flip(u0, axis=0)
umin_evo = Darwin.evolveNS(umin0)
print("Shape of umin0: ", umin0.shape)
print("Shape of umin_evo: ", umin_evo.shape)

################################################################################

# Finally, as a bonus, let's make a plot showing off the evolution we just did.
# This requires matplotlib.
# If we can't import matplotlib, we'll just skip this section.
l_mpl = True
try:
    import matplotlib as mpl
    import matplotlib.pyplot as py
    mpl.rc('font',size=22,weight='normal')
    mpl.rc('text',usetex=True)
    mpl.rc('text.latex')
except:
    print("Couldn't load matplotlib, not making plots.")
    l_mpl = False

if(l_mpl):
    nrows, ncols = 2, 3
    fig = py.figure(figsize=(ncols*8,nrows*6))
    ax1 = py.subplot(nrows,ncols,1)
    ax2 = py.subplot(nrows,ncols,2)
    ax3 = py.subplot(nrows,ncols,3)
    ax4 = py.subplot(nrows,ncols,4)
    ax5 = py.subplot(nrows,ncols,5)
    ax6 = py.subplot(nrows,ncols,6)
    #
    ax1.plot(x, u0[:,2,0],           'x', label=r'$Q^2=4$~GeV$^2$')
    ax1.plot(x, gpd_evo[1,:,2,0,-1], '+', label=r'$Q^2=25$~GeV$^2$')
    ax1.set_title(r'Up quark, $\xi=0.3, t=0$')
    ax1.set_xlabel(r'$x$')
    #
    ax2.plot(x, np.zeros(x.shape),   'x', label=r'$Q^2=4$~GeV$^2$')
    ax2.plot(x, gpd_evo[4,:,2,0,-1], '+', label=r'$Q^2=25$~GeV$^2$')
    ax2.set_title(r'Charm quark, $\xi=0.3, t=0$')
    ax2.set_xlabel(r'$x$')
    #
    ax3.plot(x, g0[:,2,0],           'x', label=r'$Q^2=4$~GeV$^2$')
    ax3.plot(x, gpd_evo[0,:,2,0,-1], '+', label=r'$Q^2=25$~GeV$^2$')
    ax3.set_title(r'Gluon, $\xi=0.3, t=0$')
    ax3.set_xlabel(r'$x$')
    #
    for ix in [31, 37, 42]:
        ax4.plot(Q2, gpd_evo[1,ix,2,0,:], 'x', label=r'$x = {:.2f}$'.format(x[ix]))
    ax4.set_xlabel(r'$Q^2$ (GeV$^2$)')
    ax4.set_title(r'Up quark, $\xi=0.3, t=0$')
    #
    for ix in [31, 37, 42]:
        ax5.plot(Q2, gpd_evo[4,ix,2,0,:], 'x', label=r'$x = {:.2f}$'.format(x[ix]))
    ax5.set_xlabel(r'$Q^2$ (GeV$^2$)')
    ax5.set_title(r'Charm quark, $\xi=0.3, t=0$')
    #
    for ix in [31, 37, 42]:
        ax6.plot(Q2, gpd_evo[0,ix,2,0,:], 'x', label=r'$x = {:.2f}$'.format(x[ix]))
    ax6.set_xlabel(r'$Q^2$ (GeV$^2$)')
    ax6.set_title(r'Gluon, $\xi=0.3, t=0$')
    #
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        _ = ax.legend(prop = { 'size' : 17 })
    fig.tight_layout()
    py.show()
    del fig
