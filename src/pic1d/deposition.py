"""Charge deposition: particles -> grid.

Cloud-in-cell (CIC / linear) weighting: each particle deposits charge
onto its two nearest nodes, weighted by proximity. This is first-order
weighting; it must exactly conserve total charge (see unit test).
"""
import numpy as np


def deposit_charge(grid, species_list, periodic=True):
    grid.rho[:] = 0.0                       # clear from last step

    for sp in species_list:
        j = np.floor(sp.x / grid.dx).astype(int)   # left node index
        f = sp.x / grid.dx - j                      # fractional distance

        contrib = sp.charge * sp.weight / grid.dx

        np.add.at(grid.rho, j,     contrib * (1.0 - f))   # to left node
        np.add.at(grid.rho, j + 1, contrib * f)           # to right node

    if periodic:
        grid.rho[0]  += grid.rho[-1]        # last node == first node
        grid.rho[-1]  = grid.rho[0]         # keep them equal
