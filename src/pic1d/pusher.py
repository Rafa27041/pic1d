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


def gather_field(grid, x):
    j = np.floor(x / grid.dx).astype(int)
    f = x / grid.dx - j
    return grid.efield[j] * (1.0 - f) + grid.efield[j + 1] * f


def push(species, grid, dt):
    E_p = gather_field(grid, species.x)
    species.v += (species.charge / species.mass) * E_p * dt
    species.x += species.v * dt
