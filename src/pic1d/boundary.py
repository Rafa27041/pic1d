"""Particle boundary conditions.

Week 2 needs only periodic wrapping. Week 3 adds absorbing walls and
steady injection for the sheath problem.
"""
import numpy as np


def apply_periodic(species, length: float) -> None:
    """Wrap positions into [0, length): x = x mod L. TODO (Week 2)."""
    raise NotImplementedError


def apply_absorbing(species, length: float, flux_counter=None) -> None:
    """Remove particles with x < 0 or x > length.

    If flux_counter is provided (diagnostics.WallFluxCounter), record
    the count and kinetic energy of absorbed particles per wall. This
    is the raw data for the Week 4 wall heat flux estimate.

    TODO (Week 3).
    """
    raise NotImplementedError


def inject_maxwellian(species, side: str, n_inject: int, v_thermal: float,
                      length: float, dt: float, rng=None) -> None:
    """Inject particles from a wall to sustain the bulk plasma.

    Correct flux injection samples velocities from a *flux* (drifting
    half-Maxwellian) distribution, v ~ v * exp(-v^2 / 2 v_th^2), not a
    plain half-Maxwellian; and positions are placed a random fraction
    of v*dt inside the domain so injection is uniform in time.
    Getting this wrong is the #1 source of sheath run bugs; budget
    debugging time here.

    TODO (Week 3).
    """
    raise NotImplementedError
