"""Particle species container.

Positions in Debye lengths, velocities in electron thermal speeds.
Charge and mass are per-macroparticle in normalized units, e.g.
electrons: q = -1, m = 1; ions: q = +1, m = mass_ratio (use a reduced
ratio like 100 or 400 to keep sheath runs affordable; note this in the
report as standard practice).
"""
import numpy as np


class Species:
    def __init__(self, name: str, charge: float, mass: float,
                 x: np.ndarray, v: np.ndarray, weight: float = 1.0):
        self.name = name
        self.charge = charge
        self.mass = mass
        self.x = np.asarray(x, dtype=float)
        self.v = np.asarray(v, dtype=float)
        self.weight = weight  # macroparticle weight

    @property
    def n_particles(self) -> int:
        return self.x.size

    @classmethod
    def maxwellian(cls, name, charge, mass, n_particles, length,
                   v_thermal, drift=0.0, rng=None, weight=1.0):
        if rng is None:
            rng = np.random.default_rng()
        x = rng.uniform(0.0, length, size=n_particles)
        v = rng.normal(drift, v_thermal, size=n_particles)
        return cls(name, charge, mass, x, v, weight)

    def kinetic_energy(self) -> float:
        """Total kinetic energy: 0.5 * m * w * sum(v^2)."""
        return 0.5 * self.mass * self.weight * np.sum(self.v ** 2)

    def remove(self, mask: np.ndarray) -> dict:
        """Remove particles where mask is True (absorbed at walls).

        Returns a dict of diagnostics for the removed particles
        (count, summed kinetic energy) so flux counters can accumulate
        them. TODO (Week 3).
        """
        raise NotImplementedError
