"""Charge deposition: particles -> grid.

Cloud-in-cell (CIC / linear) weighting: each particle deposits charge
onto its two nearest nodes, weighted by proximity. This is first-order
weighting; it must exactly conserve total charge (see unit test).
"""
import numpy as np


def deposit_charge(grid, species_list, periodic: bool = True) -> None:
    """Accumulate charge density on grid.rho from all species.

    Algorithm (vectorized, no Python loop over particles):
    1. grid.rho[:] = 0
    2. For each species:
       - j = floor(x / dx)               (left node index)
       - f = x / dx - j                  (fractional distance)
       - np.add.at(rho, j,   q * w * (1 - f) / dx)
       - np.add.at(rho, j+1, q * w * f / dx)
    3. Periodic: fold node n-1 onto node 0 (they are the same point),
       or use modular indexing.
    4. Dirichlet/sheath: clamp end contributions; end nodes represent
       half-cells, so divide their density by 2 (volume correction).

    np.add.at is required (not fancy-index +=) because multiple
    particles hit the same node.

    TODO (Week 1).
    """
    raise NotImplementedError
