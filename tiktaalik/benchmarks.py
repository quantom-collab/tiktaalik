"""
benchmarks.py

part of the tiktaalik package for GPD evolution
by Adam Freese

Methods to benchmark the accuracy of the finite element method.
"""

import numpy as np
from scipy.optimize import curve_fit

import matplotlib as mpl
import matplotlib.pyplot as plt

from . import continuum, matrices, model, pars, polynomiality

mpl.rc('font',size=30,family='cmr10',weight='normal')
mpl.rc('text',usetex=True)
mpl.rc('text.latex', preamble=r"\usepackage{bm,amsmath,amssymb,amsfonts,mathrsfs}")
plt.rcParams["axes.formatter.use_mathtext"] = True

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Evolution shift benchmark

def shift_benchmark(
        key='NS',
        xi=0.1,
        nx=81,
        pol=False,
        nlo=False,
        ns_type=1,
        grid_type=2
        ):
    ''' Benchmark of the shift induced by a kernel on the right-hand side of the
    GPD evolution equation.
    ------
    Input:
        - key (string)
            'NS', 'qq', 'qg', 'gq', or 'gg'
            What kind of kernel is being benchmarked
        - xi (float)
            skewness value
            For xi<=0.1, grid_type=2 is recommended
            For xi>=0.1, grid_type=1 is recommended
        - nx (int)
            number of x points
            if grid_type=2, then using one more than a multiple of four will
            give x=xi and x=-xi on the grid
        - pol (bool)
            True or helicity-dependent (A-type), False for helicity-independent (V-type)
        - nlo (bool)
            True or NLO evolution, False for LO evolution
        - ns_type (int)
            +1 or -1
            relevant only for NLO evolution; whether it's a q+qbar or
            q-qbar non-singlet mixture that's being evolved.
        - grid_type (int)
            1 : linear x spacing
            2 : log-linear hybrid spacing
    '''
    # Set the grids
    matrices.set_x_xi_grids(nx, xi, grid_type)
    matrices.set_Q2_grid(np.array([4,5]))
    # Retrieve x grid
    x_pixel = matrices.get_x_grid()
    # Retrieve kernel matrix (interpixel method)
    K  = _get_kernel(key=key, pol=pol, ns_type=ns_type, nlo=nlo)
    # Interpixel shift
    H0 = _get_gpd(x_pixel, xi=xi, pol=pol, key=key)
    dH_pixel = np.einsum('ij,j->i', K[:,:,0], H0)
    # Continuum shift ("ground truth")
    x_truth = _make_continuum_x(xi)
    dH_truth = _get_continuum_shift(x_truth, pol=pol, key=key, xi=xi, nlo=nlo, ns_type=ns_type)
    # Error
    dH_truth_2 = _get_continuum_shift(x_pixel, pol=pol, key=key, xi=xi, nlo=nlo, ns_type=ns_type)
    error = 100*abs(dH_pixel - dH_truth_2) / abs(dH_truth_2)
    # Set up the plot
    nrows, ncols = 2, 1
    fig, (ax1, ax2) = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1]},
            figsize=(8,8),
            layout = 'constrained'
            )
    ax1.plot(x_truth, dH_truth, '-', label=r'Truth',     color='tab:orange')
    ax1.plot(x_pixel, dH_pixel, '+', label=r'tiktaalik', color='tab:blue')
    # Error
    ax2.plot(x_pixel, error, '+', label=r'tiktaalik', color='tab:blue')
    # Plot labels etc
    ax2.set_ylim((1e-6,1e3))
    ax2.plot(x_truth, x_truth*0+1, linewidth=1, color='tab:gray')
    for ax in [ax1, ax2]:
        _plot_xi_lines(ax, xi)
        ax.set_xlim((-1,1))
        if(grid_type==2):
            ax.set_xscale('symlog', linthresh=xi)
    if(grid_type==2):
        ax2.set_xticks(ax.get_xticks()[::2])
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$x$')
    ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel(r'$\int \mathrm{d} y \, K(x,y,\xi) H(y)$')
    ax2.set_ylabel(r'percent error')
    legend = ax1.legend(prop = { 'size' : 26 })
    legend.get_frame().set_facecolor('#f8f8f8')
    bbox = dict(facecolor='#f8f8f8', alpha=0.76, edgecolor='gray', boxstyle='round,pad=0.2')
    ax1.annotate(
            r'\textbf{'+key+r'}', xy=(0.03,0.05), xycoords='axes fraction',
            bbox=bbox
            )
    fig.patch.set_alpha(0)
    return fig

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wilson coefficient benchmark

def wilson_benchmark(
        key='q',
        nx = 101,
        t = 0,
        grid_type = 2,
        nlo = False
        ):
    # Set the grids
    xi = np.geomspace(1e-4, 1, 100)
    matrices.set_Q2_grid(np.array([4,5]))
    matrices.set_x_xi_grids(nx, xi, grid_type)
    # Retrieve x grid
    x = matrices.get_x_grid()
    # Retrieve coefficient matrix (interpixel method)
    C = _get_dvcs_coefficient(key=key, nlo=nlo)
    # Interpiixel CFF
    gpd = _get_gpd_dvcs(x, xi, t, key=key)
    cff_pixel = np.einsum('ij,ji...->i...', C, gpd)[:,0]
    # Continuum CFF ("ground truth")
    cff_truth = _get_continuum_cff(xi, key=key, Q2=4, nlo=nlo)
    # Set up plot
    nrows, ncols = 2, 1
    fig, (ax1, ax2) = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1]},
            figsize=(8,8),
            layout = 'constrained'
            )
    ax1.plot(xi, xi*np.real(cff_truth), '-',  label=r'Truth     (real)', color='tab:blue')
    ax1.plot(xi, xi*np.real(cff_pixel), '+',  label=r'tiktaalik (real)', color='tab:orange')
    ax1.plot(xi, xi*np.imag(cff_truth), '--', label=r'Truth     (imag)', color='tab:green')
    ax1.plot(xi, xi*np.imag(cff_pixel), 'x',  label=r'tiktaalik (imag)', color='tab:purple')
    # Error
    ImErr = 100*abs(np.imag(cff_truth-cff_pixel) / np.imag(cff_truth))
    ReErr = 100*abs(np.real(cff_truth-cff_pixel) / np.real(cff_truth))
    ax2.plot(xi, ReErr, '+', label=r'Error (real)', color='tab:orange')
    ax2.plot(xi, ImErr, 'x', label=r'Error (imag)', color='tab:purple')
    # Plot labels etc
    for ax in [ax1, ax2]:
        ax.set_xscale('log')
    ax2.set_xlabel(r'$\xi$')
    ax1.set_ylabel(r'$\mathcal{H}_{'+key+r'}(\xi)$')
    ax2.set_ylabel(r'percent error')
    legend = ax1.legend(prop = { 'size' : 26 })
    legend.get_frame().set_facecolor('#f8f8f8')
    ax1.get_xaxis().set_visible(False)
    if(nlo):
        lo_or_nlo = 'NLO'
    else:
        lo_or_nlo = 'LO'
    bbox = dict(facecolor='#f8f8f8', alpha=0.76, edgecolor='gray', boxstyle='round,pad=0.2')
    ax1.annotate(
            r'\textbf{'+lo_or_nlo+r'}', xy=(0.03,0.05), xycoords='axes fraction',
            bbox=bbox
            )
    fig.patch.set_alpha(0)
    return fig

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Polynomiality benchmark(s)

# TODO: clean up, make more systematic
# - Need singlet/gluon tests too!
# - Coefficients as functions of t? (prolly in polynomiality module)
# - Probably separate plot(s) for the coefficients
# - Maybe nx dependence of mean residual to study how it's a discretization artifact


##def poly(x, N, c):
##    # TODO: clean up...
##    ff = [ c[i]*x**i for i in range(N+1) ]
##    # What if I force polynomial to be even?
##    for i in range(N+1):
##        if( 2*(i//2)!=i ):
##            ff[i] = 0
##    return sum(ff)
##
##def fit_polynomial(xx, yy, Nmax):
##    # TODO: clean up...
##    f = lambda x, *p0: poly(x, Nmax, p0)
##    print(Nmax)
##    c, cov = curve_fit(f, xx, yy, p0 = [0. for _ in range(Nmax+1)])
##    return c

def polynomiality_test_ns(grid_type=2, ns=1, force_even=True):
    # x/xi grids
    nx = 81
    nxi = 80
    if(grid_type==1):
        xi = np.linspace(0.2, 0.8, nxi)
    else:
        #xi = np.geomspace(1e-4, 1-1e-4, nxi)
        xi = np.linspace(1e-2, 1-1e-2, nxi)
    matrices.set_x_xi_grids(nx, xi, grid_type)
    x = matrices.get_x_grid()
    # Q2 grid
    Q2 = np.geomspace(4, 17, 11)
    matrices.set_Q2_grid(Q2)
    # Initial and evolved GPDs
    H0 = (model.Hu(x, xi, 0)[:,:,0] + model.Hu(-x, xi, 0)[:,:,0])/2
    matrices.do_lo_evolution()
    Mev1 = matrices.matrix_VNS(ns_type=-1)
    matrices.do_nlo_evolution()
    Mev2 = matrices.matrix_VNS(ns_type=-1)
    H1 = np.einsum('xyzq,yz->xzq', Mev1, H0)[:,:,-1]
    H2 = np.einsum('xyzq,yz->xzq', Mev2, H0)[:,:,-1]
    # Stuff from polynomiality module
    data0_sets = []
    data1_sets = []
    data2_sets = []
    for i in range(ns):
        s = 1 + 2*i
        data0_sets += [ polynomiality.get_polynomial(H0, s, force_even=force_even) ]
        data1_sets += [ polynomiality.get_polynomial(H1, s, force_even=force_even) ]
        data2_sets += [ polynomiality.get_polynomial(H2, s, force_even=force_even) ]
    # Set up plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nrows, ncols = 3, ns
    fig, axes = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1,1]},
            figsize=(8*ncols,13),
            layout = 'constrained'
            )
    if(ns==1):
        axes = axes[:,np.newaxis]
    # Loop over s values
    for i in range(ns):
        # Stuff from polynomiality fits
        P0_raw = data0_sets[i]['raw polynomial']
        P0_fit = data0_sets[i]['fit polynomial']
        P1_raw = data1_sets[i]['raw polynomial']
        P1_fit = data1_sets[i]['fit polynomial']
        P2_raw = data2_sets[i]['raw polynomial']
        P2_fit = data2_sets[i]['fit polynomial']
        P0_res = P0_raw - P0_fit
        P1_res = P1_raw - P1_fit
        P2_res = P2_raw - P2_fit
        c0 = data0_sets[i]['c']
        c1 = data1_sets[i]['c']
        c2 = data2_sets[i]['c']
        n = data0_sets[i]['n']
        # TODO: errors
        # Plot raw moments
        axes[0,i].plot(xi, P0_raw, 'o',  label=r'Initial', color='tab:blue')
        axes[0,i].plot(xi, P1_raw, 'x',  label=r'LO',      color='tab:orange')
        axes[0,i].plot(xi, P2_raw, '+',  label=r'NLO',     color='tab:green')
        # Plot polynomial fits
        axes[0,i].plot(xi, P0_fit, '-',                    color='tab:blue')
        axes[0,i].plot(xi, P1_fit, '--',                   color='tab:orange')
        axes[0,i].plot(xi, P2_fit, '-.',                   color='tab:green')
        # Plot residuals
        axes[1,i].plot(xi, P0_res, 'o',  label=r'Initial', color='tab:blue')
        axes[1,i].plot(xi, P1_res, 'x',  label=r'LO',      color='tab:orange')
        axes[1,i].plot(xi, P2_res, '+',  label=r'NLO',     color='tab:green')
        # Plot the polynomial coefficients
        axes[2,i].plot(n, c0, 'o',  label=r'Initial', color='tab:blue')
        axes[2,i].plot(n, c1, 'x',  label=r'LO',      color='tab:orange')
        axes[2,i].plot(n, c2, '+',  label=r'NLO',     color='tab:green')
    # Finish up plot decorations
    axes[0,0].set_ylabel(r'$\mathcal{M}(s,\xi)$')
    axes[1,0].set_ylabel(r'Residuals')
    axes[2,0].set_ylabel(r'$A_{s,n}$')
    for i in range(ns):
        axes[0,i].get_xaxis().set_visible(False)
        axes[1,i].set_xlabel(r'$\xi$')
        axes[2,i].set_xlabel(r'$n$')
        #if(grid_type==2):
        #    axes[0,i].set_xscale('log')
        #    axes[1,i].set_xscale('log')
    bbox = dict(facecolor='#f8f8f8', alpha=0.76, edgecolor='gray', boxstyle='round,pad=0.2')
    for i in range(ns):
        axes[0,i].annotate(
                r'$\mathbf{s='+'{:d}'.format(2*i+1)+r'}$', xy=(0.03,0.05), xycoords='axes fraction',
                bbox=bbox
                )
    legend = axes[0,0].legend(prop = { 'size' : 26 })
    legend.get_frame().set_facecolor('#f8f8f8')
    fig.patch.set_alpha(0)
    return fig




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Auxiliary functions

def _get_kernel(key='NS', pol=False, nfl=4, ns_type=1, nlo=False):
    if(key=='NS' and pol==False):
        K = matrices.kernel_VQQ(nfl=nfl, nlo=nlo, ns_type=ns_type)
    elif(key=='qq' and pol==False):
        K = matrices.kernel_VQQ(nfl=nfl, nlo=nlo, ns_type=0)
    elif(key=='qg' and pol==False):
        K = matrices.kernel_VQG(nfl=nfl, nlo=nlo)
    elif(key=='gq' and pol==False):
        K = matrices.kernel_VGQ(nfl=nfl, nlo=nlo)
    elif(key=='gg' and pol==False):
        K = matrices.kernel_VGG(nfl=nfl, nlo=nlo)
    elif(key=='NS' and pol==True):
        K = matrices.kernel_AQQ(nfl=nfl, nlo=nlo, ns_type=ns_type)
    elif(key=='qq' and pol==True):
        K = matrices.kernel_AQQ(nfl=nfl, nlo=nlo, ns_type=0)
    elif(key=='qg' and pol==True):
        K = matrices.kernel_AQG(nfl=nfl, nlo=nlo)
    elif(key=='gq' and pol==True):
        K = matrices.kernel_AGQ(nfl=nfl, nlo=nlo)
    elif(key=='gg' and pol==True):
        K = matrices.kernel_AGG(nfl=nfl, nlo=nlo)
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return K

def _get_dvcs_coefficient(key='q', nlo=False):
    if(key=='q'):
        C = matrices.dvcs_Cq(nlo=nlo)[:,:,0]
    elif(key=='g'):
        C = matrices.dvcs_Cg(nlo=nlo)[:,:,0]
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return C

def _get_gpd(x, xi=0.1, pol=False, key='NS'):
    H = np.zeros(x.shape)
    if(key=='NS' and pol==False):
        H = _ns_gpd(x, xi)
    elif(key[1]=='g' and pol==False):
        H = _gluon_gpd(x, xi)
    elif(key[1]=='q' and pol==False):
        H = _singlet_gpd(x, xi)
    elif(key=='NS' and pol==True):
        H = _ns_gpd_tilde(x, xi)
    elif(key[1]=='g' and pol==True):
        H = _gluon_gpd_tilde(x, xi)
    elif(key[1]=='q' and pol==True):
        H = _singlet_gpd_tilde(x, xi)
    else:
        raise ValueError("Key "+key+" unrecognized.")
    H[np.isnan(H)] = 0
    H[np.isinf(H)] = 0
    return H

def _get_gpd_dvcs(x, xi, t, key='q'):
    H = np.zeros(x.shape)
    if(key=='g'):
        H = _gluon_gpd_with_t(x, xi, t)
    elif(key=='q'):
        H = _dvcs_quark_combo(x, xi, t)
    else:
        raise ValueError("Key "+key+" unrecognized.")
    H[np.isnan(H)] = 0
    H[np.isinf(H)] = 0
    return H

def _get_continuum_shift(x, key='NS', pol=False, xi=0.1, Q2=pars.mc2, nlo=False, ns_type=1):
    if(key=='NS' and pol==False):
        continuum_shift = continuum.shift_cNS(x, xi, Q2, nlo, ns_type)[:,0]
    elif(key=='qq' and pol==False):
        continuum_shift = continuum.shift_cQQ(x, xi, Q2, nlo)[:,0]
    elif(key=='qg' and pol==False):
        continuum_shift = continuum.shift_cQG(x, xi, Q2, nlo)[:,0]
    elif(key=='gq' and pol==False):
        continuum_shift = continuum.shift_cGQ(x, xi, Q2, nlo)[:,0]
    elif(key=='gg' and pol==False):
        continuum_shift = continuum.shift_cGG(x, xi, Q2, nlo)[:,0]
    elif(key=='NS' and pol==True):
        continuum_shift = continuum.shift_cNS_tilde(x, xi, Q2, nlo, ns_type)[:,0]
    elif(key=='qq' and pol==True):
        continuum_shift = continuum.shift_cQQ_tilde(x, xi, Q2, nlo)[:,0]
    elif(key=='qg' and pol==True):
        continuum_shift = continuum.shift_cQG_tilde(x, xi, Q2, nlo)[:,0]
    elif(key=='gq' and pol==True):
        continuum_shift = continuum.shift_cGQ_tilde(x, xi, Q2, nlo)[:,0]
    elif(key=='gg' and pol==True):
        continuum_shift = continuum.shift_cGG_tilde(x, xi, Q2, nlo)[:,0]
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return continuum_shift

def _get_continuum_cff(xi, key='q', Q2=4, nlo=False):
    if(key=='q'):
        continuum_cff = continuum.cff_q(xi, Q2, nlo=nlo)
    elif(key=='g'):
        continuum_cff = continuum.cff_g(xi, Q2, nlo=nlo)
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return continuum_cff

def _ns_gpd(x, xi):
    T3 = model.Hu(x,xi,0) - model.Hd(x,xi,0)
    T3 = T3 - np.flip(T3, axis=0)
    return T3[:,0,0]

def _singlet_gpd(x, xi):
    return model.H_singlet(x,xi,0)[:,0,0]

def _gluon_gpd(x, xi):
    return model.Hg(x,xi,0)[:,0,0]

def _gluon_gpd_with_t(x, xi,t):
    return model.Hg(x,xi,t)

def _ns_gpd_tilde(x, xi):
    T3 = model.Hu_tilde(x,xi,0) - model.Hd_tilde(x,xi,0)
    T3 = T3 + np.flip(T3, axis=0)
    return T3[:,0,0]

def _singlet_gpd_tilde(x, xi):
    Sigma = model.Hu_tilde(x,xi,0) + model.Hd_tilde(x,xi,0)
    Sigma = Sigma - np.flip(Sigma, axis=0)
    return Sigma[:,0,0]

def _gluon_gpd_tilde(x, xi):
    # Sicne Hg_tilde=0 in GK model, need something non-zero
    # with the right symmetries as a stand-in for shift tests...
    G = model.Hu_tilde(x,xi,0) + model.Hd_tilde(x,xi,0)
    G = G + np.flip(G, axis=0)
    return G[:,0,0]

def _dvcs_quark_combo(x, xi, t):
    H = (
            4/9*(model.Hu(x,xi,t) - model.Hu(-x,xi,t))
            +
            1/9*(model.Hd(x,xi,t) - model.Hd(-x,xi,t))
            +
            1/9*(model.Hs(x,xi,t) - model.Hs(-x,xi,t))
            )
    return H

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot utilities

def _make_continuum_x(xi):
    N = 200
    x1 = np.geomspace( -1, -xi, N)
    x2 = np.linspace( -xi,  xi, N)
    x3 = np.geomspace( xi,   1, N)
    x = np.concatenate((x1[:N-1], x2[0:N], x3[1:]))
    return x

def _plot_xi_lines(ax, xi):
    ymin, ymax = ax.get_ylim()
    ax.vlines( xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.vlines(-xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.set_ylim((ymin,ymax))
    return
