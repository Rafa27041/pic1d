"""Unit tests for cloud-in-cell charge deposition."""
import numpy as np
import pytest

from pic1d.grid import Grid
from pic1d.particles import Species
from pic1d import deposition


def test_total_charge_conserved():
    grid = Grid(length=10.0, n_nodes=101)
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 10.0, size=5000)
    sp = Species("e", charge=-1.0, mass=1.0, x=x, v=np.zeros(5000))
    deposition.deposit_charge(grid, [sp], periodic=True)
    total = np.sum(grid.rho[:-1]) * grid.dx        # drop duplicated end node
    assert np.isclose(total, -5000.0)


def test_particle_on_node():
    grid = Grid(length=10.0, n_nodes=101)          # dx = 0.1, node 30 at x=3.0
    sp = Species("e", charge=-1.0, mass=1.0,
                 x=np.array([3.0]), v=np.array([0.0]))
    deposition.deposit_charge(grid, [sp], periodic=False)
    assert np.isclose(grid.rho[30] * grid.dx, -1.0)   # all charge on node 30
    assert np.isclose(grid.rho[31], 0.0)


def test_particle_at_midpoint():
    grid = Grid(length=10.0, n_nodes=101)          # midpoint of a cell at x=3.05
    sp = Species("e", charge=-1.0, mass=1.0,
                 x=np.array([3.05]), v=np.array([0.0]))
    deposition.deposit_charge(grid, [sp], periodic=False)
    assert np.isclose(grid.rho[30], grid.rho[31])     # 50/50 split
