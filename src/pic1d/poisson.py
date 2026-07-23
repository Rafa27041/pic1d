"""1D Poisson solvers: d^2(phi)/dx^2 = -rho (normalized units).

Two solvers are needed:
- Periodic domain (validation runs): FFT-based solve.
- Dirichlet walls (sheath runs): tridiagonal (Thomas) solve with
  phi(0) = phi_left, phi(L) = phi_right.
"""
import numpy as np


def solve_periodic(rho: np.ndarray, dx: float) -> np.ndarray:
    """FFT Poisson solve on a periodic domain.

    Algorithm:
    1. rho_k = fft(rho[:-1])  (drop duplicated end node)
    2. phi_k = rho_k / k^2 for k != 0; phi_k[0] = 0
       (zero-mean potential; the k=0 mode is gauge freedom, and the
       run must be globally neutral or the solve is ill-posed)
    3. phi = real(ifft(phi_k)); append phi[0] to restore end node.

    Use k = 2*pi*fftfreq(n, d=dx). TODO (Week 1).
    """
    raise NotImplementedError


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
