"""1D Poisson solvers: d^2(phi)/dx^2 = -rho (normalized units).

Two solvers are needed:
- Periodic domain (validation runs): FFT-based solve.
- Dirichlet walls (sheath runs): tridiagonal (Thomas) solve with
  phi(0) = phi_left, phi(L) = phi_right.
"""
import numpy as np


def solve_periodic(rho: np.ndarray, dx: float) -> np.ndarray:
  rho_inner = rho[:-1]
  n = rho_inner.size
  k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
  rho_k = np.fft.fft(rho_inner)
  phi_k = np.zeros_like(rho_k)
  phi_k[1:] = rho_k[1:] / k[1:]**2
  phi = np.real(np.fft.ifft(phi_k))
  return np.append(phi, phi[0])


def solve_dirichlet(rho: np.ndarray, dx: float,
                    phi_left: float = 0.0, phi_right: float = 0.0) -> np.ndarray:
    """Tridiagonal Poisson solve with fixed-potential walls.

    Discretization: (phi[i-1] - 2 phi[i] + phi[i+1]) / dx^2 = -rho[i]
    for interior nodes; boundary values pinned.

    Use scipy.linalg.solve_banded or write the Thomas algorithm
    directly (it is ~10 lines and worth doing once by hand).

    TODO (Week 1). Unit test: manufactured solution, e.g.
    phi_exact = sin(pi x / L) with corresponding rho.
    """
    raise NotImplementedError
