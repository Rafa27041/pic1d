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

def thomas(a, b, c, d):
    """Solve a tridiagonal system. a=sub, b=diag, c=super, d=rhs."""
    n = len(d)
    cp = np.zeros(n)
    dp = np.zeros(n)

    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]

    for i in range(1, n):                        # forward sweep
        m = b[i] - a[i] * cp[i-1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i-1]) / m

    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):                 # back substitution
        x[i] = dp[i] - cp[i] * x[i+1]

    return x
  
def solve_dirichlet(rho: np.ndarray, dx: float,
                    phi_left: float = 0.0, phi_right: float = 0.0) -> np.ndarray:
    n = rho.size - 2
    a = np.ones(n)
    b = -2.0 * np.ones(n)
    c = np.ones(n)
    d = -rho[1:-1] * dx**2
    d[0]  -= phi_left
    d[-1] -= phi_right

    x = thomas(a, b, c, d)

    phi = np.empty(rho.size)
    phi[0], phi[-1] = phi_left, phi_right
    phi[1:-1] = x
    return phi
    
