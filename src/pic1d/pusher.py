"""Leapfrog particle push and field gather (grid -> particles).

The leapfrog scheme staggers v and x by half a timestep:
    v^(n+1/2) = v^(n-1/2) + (q/m) E(x^n) dt
    x^(n+1)   = x^n + v^(n+1/2) dt

It is time-reversible and symplectic: energy drift should stay below
~1% over thousands of steps if implemented correctly. That plot is
your correctness certificate.

Remember the half-step setup: before the first step, pull velocities
back by dt/2 (v -= 0.5 * (q/m) * E * dt) so v and x are staggered.
"""
import numpy as np


def gather_field(grid, x: np.ndarray) -> np.ndarray:
    """Interpolate grid.efield to particle positions.

    MUST use the same linear (CIC) weighting as deposition, otherwise
    the scheme produces a self-force on isolated particles.

    TODO (Week 1).
    """
    raise NotImplementedError


def push(species, grid, dt: float) -> None:
    """One leapfrog step: update species.v then species.x in place.

    TODO (Week 1):
    1. E_p = gather_field(grid, species.x)
    2. species.v += (species.charge / species.mass) * E_p * dt
    3. species.x += species.v * dt
    (Boundary handling lives in boundary.py, applied after the push.)
    """
    raise NotImplementedError
