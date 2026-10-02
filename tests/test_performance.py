import numpy as np
import pytest
from tiktaalik.evolution import Evolver
from tiktaalik.matrices import pixelspace, initialize_kernels
from tiktaalik.qcd import alphaQCD
from tiktaalik import pars

XI = np.array([0.1, 0.2])


def _build_evolver(nx, nQ2=17):
    return Evolver(nx=nx, xi=XI, nQ2=nQ2,
                   Q2i=pars.mc2 + pars.epsilon, Q2f=pars.mb2 - pars.epsilon, grid_type=1)


@pytest.mark.benchmark(group="setup")
def test_perf_initialize_kernels(benchmark):
    # Dominant setup cost (~seconds); few rounds because each call is heavy.
    benchmark.pedantic(initialize_kernels, args=(100, XI, 1), rounds=3, iterations=1)


@pytest.mark.benchmark(group="setup")
def test_perf_evolver_construction(benchmark):
    # Full Evolver build = kernels + evolution matrices.
    benchmark.pedantic(_build_evolver, args=(100,), kwargs={"nQ2": 17}, rounds=3, iterations=1)


@pytest.mark.benchmark(group="evolve")
def test_perf_evolve_lo(benchmark):
    ev = _build_evolver(100)
    gpd = np.random.default_rng(0).random((ev.nfl + 1, ev.nx, ev.nxi))
    benchmark(ev.evolveGPDs, gpd, 'V', nlo=False)


@pytest.mark.benchmark(group="evolve")
def test_perf_evolve_nlo(benchmark):
    ev = _build_evolver(100)
    gpd = np.random.default_rng(1).random((ev.nfl + 1, ev.nx, ev.nxi))
    benchmark(ev.evolveGPDs, gpd, 'V', nlo=True)


@pytest.mark.benchmark(group="util")
def test_perf_pixelspace(benchmark):
    benchmark(pixelspace, 1000, 0.5, 1)


@pytest.mark.benchmark(group="util")
def test_perf_alphaQCD(benchmark):
    q2 = np.logspace(0, 4, 1000)
    benchmark(alphaQCD, q2)
