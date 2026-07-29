"""Unit tests for the leapfrog pusher.

The key physics test: a single particle in a fixed sinusoidal field
must conserve energy to high accuracy over many oscillations.
"""
import numpy as np
import pytest

from pic1d.grid import Grid
from pic1d.particles import Species
from pic1d import pusher


def test_single_particle_energy_conservation():
    L, N = 10.0, 1001
    grid = Grid(length=L, n_nodes=N)
    c = L / 2
    grid.efield = (grid.x - c)          # E = +(x-c) gives restoring accel for q/m=-1

    sp = Species("e", charge=-1.0, mass=1.0,
                 x=np.array([c + 0.5]), v=np.array([0.0]))   # small amplitude
    dt, nsteps = 0.005, 5000

    a = (sp.charge / sp.mass) * gather_field(grid, sp.x)
    sp.v -= 0.5 * a * dt                # leapfrog half-step pullback

    energies = []
    for _ in range(nsteps):
        pusher.push(sp, grid, dt)
        ke = 0.5 * sp.mass * sp.v[0]**2
        pe = 0.5 * (sp.x[0] - c)**2
        energies.append(ke + pe)

    energies = np.array(energies)
    drift = (energies.max() - energies.min()) / energies.mean()
    assert drift < 0.02


def test_gather_matches_deposition_weighting():
    grid = Grid(length=10.0, n_nodes=101)   # dx = 0.1
    grid.efield[30] = 5.0                    # node 30 at x=3.0

    on_node = gather_field(grid, np.array([3.0]))
    assert np.isclose(on_node[0], 5.0)       # particle on node feels full value

    midpoint = gather_field(grid, np.array([3.05]))   # halfway to node 31 (=0)
    assert np.isclose(midpoint[0], 2.5)      # half of 5.0
