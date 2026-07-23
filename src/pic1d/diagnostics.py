"""Diagnostics: energy history, phase space snapshots, wall fluxes."""
import numpy as np


class EnergyHistory:
    """Accumulate kinetic, field, and total energy each step.

    Field energy (normalized): 0.5 * sum(E^2) * dx.
    Total energy drift < ~1% over the run validates the leapfrog
    implementation. TODO (Week 2).
    """

    def __init__(self):
        self.time, self.kinetic, self.field, self.total = [], [], [], []

    def record(self, t, species_list, grid):
        raise NotImplementedError


class PhaseSpaceRecorder:
    """Store (x, v) snapshots every n steps for animations.

    Two-stream phase-space vortices are the money plot; save frames
    as arrays and render GIFs with matplotlib.animation. TODO (Week 2).
    """

    def __init__(self, every_n_steps: int = 50):
        self.every = every_n_steps
        self.frames = []

    def maybe_record(self, step, species_list):
        raise NotImplementedError


class WallFluxCounter:
    """Accumulate particle count and kinetic energy absorbed at each wall.

    Time-averaged over the steady-state window, this gives:
    - particle flux Gamma = N_absorbed * weight / (A * t_window)
    - energy flux  q_wall = E_absorbed * weight / (A * t_window)
    In normalized units A = 1. q_wall is the input to the Week 4
    transient wall conduction estimate. TODO (Weeks 3-4).
    """

    def __init__(self):
        self.left = {"count": 0.0, "energy": 0.0}
        self.right = {"count": 0.0, "energy": 0.0}

    def add(self, side: str, count: float, energy: float):
        raise NotImplementedError
