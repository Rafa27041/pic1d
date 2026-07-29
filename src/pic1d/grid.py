"""Spatial grid for the 1D PIC simulation.

All lengths are in units of the Debye length. The grid stores the
cell-centered/node quantities: charge density rho, potential phi,
and electric field E, all defined on nodes for simplicity.
"""
import numpy as np


class Grid:
    """Uniform 1D grid on [0, length] with n_nodes nodes.

    Attributes
    ----------
    x : ndarray, node positions
    dx : float, node spacing
    rho : ndarray, charge density on nodes
    phi : ndarray, electrostatic potential on nodes
    efield : ndarray, electric field on nodes
    """

    def __init__(self, length: float, n_nodes: int):
        self.length = length
        self.n_nodes = n_nodes
        self.dx = length / (n_nodes - 1)
        self.x = np.linspace(0.0, length, n_nodes)
        self.rho = np.zeros(n_nodes)
        self.phi = np.zeros(n_nodes)
        self.efield = np.zeros(n_nodes)

    def compute_efield_from_phi(self, periodic=True):
        self.efield[1:-1] = -(self.phi[2:] - self.phi[:-2]) / (2 * self.dx)
        if periodic:
            self.efield[0]  = -(self.phi[1] - self.phi[-2]) / (2 * self.dx)
            self.efield[-1] = self.efield[0]
        else:
            self.efield[0]  = -(self.phi[1] - self.phi[0]) / self.dx
            self.efield[-1] = -(self.phi[-1] - self.phi[-2]) / self.dx
