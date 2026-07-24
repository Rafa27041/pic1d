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
    n = rho.size - 2
    a = np.ones(n)
    b = -2.0 * np.ones(n)
    c = np.ones(n)
    d = -rho[1:-1] * dx**2
    d[0]  -= phi_left
    d[-1] -= phi_right

    # forward sweep  (your Thomas code here)
    # back substitution  (your Thomas code here)
    # result: x, length n

    phi = np.empty(rho.size)
    phi[0], phi[-1] = phi_left, phi_right
    phi[1:-1] = x
    return phi
    
