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
    L, N = 1.0, 129
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    phi_exact = np.sin(np.pi * x / L)
    rho = (np.pi / L)**2 * np.sin(np.pi * x / L)
    phi = poisson.solve_dirichlet(rho, dx, 0.0, 0.0)
    assert np.max(np.abs(phi - phi_exact)) < 1e-4
