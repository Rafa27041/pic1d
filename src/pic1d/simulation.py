"""Main PIC loop orchestration.

The canonical loop (per step):
    1. deposit_charge(grid, species)          particles -> grid
    2. solve Poisson                          rho -> phi
    3. grid.compute_efield_from_phi()         phi -> E
    4. push(species, grid, dt)                grid -> particles
    5. apply boundary conditions
    6. diagnostics.record(...)

Stability constraints (check these in __init__ and warn loudly):
    - dt < ~0.2 / w_pe  (resolve plasma oscillations; w_pe = 1 in
      normalized units, so dt <= 0.2)
    - dx < ~0.5 lambda_D (avoid finite-grid instability / self-heating)
    - v_max * dt < dx    (particles must not skip cells)
"""
import numpy as np

from .grid import Grid
from . import deposition, poisson, pusher, boundary


class Simulation:
    def __init__(self, grid: Grid, species_list, dt: float,
                 bc: str = "periodic", diagnostics=None):
        self.grid = grid
        self.species = species_list
        self.dt = dt
        self.bc = bc
        self.diagnostics = diagnostics or []
        self.time = 0.0
        self.step_count = 0
        # TODO (Week 2): validate stability constraints here.

    def initialize(self):
        """Initial field solve + half-step velocity pullback for leapfrog.

        TODO (Week 2).
        """
        raise NotImplementedError

    def step(self):
        """Advance one timestep following the canonical loop. TODO (Week 2)."""
        raise NotImplementedError

    def run(self, n_steps: int):
        for _ in range(n_steps):
            self.step()
