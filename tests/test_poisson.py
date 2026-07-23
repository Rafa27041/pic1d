"""Unit tests for the Poisson solvers, via manufactured solutions."""
import numpy as np
import pytest

from pic1d import poisson


def test_periodic_sinusoidal():
    """rho = sin(kx) => phi = sin(kx)/k^2. Max error < 1e-3 at n=256."""
    n, L = 256, 2 * np.pi
    dx = L / n
    x = np.arange(n + 1) * dx
    rho = np.sin(x)
    phi = poisson.solve_periodic(rho, dx)
    assert np.max(np.abs(phi - np.sin(x))) < 1e-3


def test_dirichlet_manufactured():
    """phi_exact = sin(pi x / L) with rho = (pi/L)^2 sin(pi x / L),
    phi(0) = phi(L) = 0. Second-order convergence: error ~ dx^2
    (verify by running at two resolutions)."""
    assert False
