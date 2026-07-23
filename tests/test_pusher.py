"""Unit tests for the leapfrog pusher.

The key physics test: a single particle in a fixed sinusoidal field
must conserve energy to high accuracy over many oscillations.
"""
import numpy as np
import pytest

from pic1d.grid import Grid
from pic1d.particles import Species
from pic1d import pusher


@pytest.mark.skip(reason="TODO Week 1: enable once pusher is implemented")
def test_single_particle_energy_conservation():
    """Particle oscillating in E(x) = -sin(2 pi x / L): total energy
    (kinetic + potential) drift < 1e-3 over 5000 steps with dt = 0.05."""
    grid = Grid(length=10.0, n_nodes=101)
    grid.efield = -np.sin(2 * np.pi * grid.x / grid.length)
    p = Species("e", charge=-1.0, mass=1.0,
                x=np.array([2.5]), v=np.array([0.0]))
    # TODO: half-step setup, loop pusher.push, track energy, assert drift.
    assert False


@pytest.mark.skip(reason="TODO Week 1")
def test_gather_matches_deposition_weighting():
    """A particle exactly on a node must feel exactly that node's field;
    a particle at a cell midpoint must feel the average of its two nodes."""
    assert False
