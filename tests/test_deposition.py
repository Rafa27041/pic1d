"""Unit tests for cloud-in-cell charge deposition."""
import numpy as np
import pytest

from pic1d.grid import Grid
from pic1d.particles import Species
from pic1d import deposition


@pytest.mark.skip(reason="TODO Week 1")
def test_total_charge_conserved():
    """sum(rho) * dx must equal total particle charge exactly
    (to machine precision), for random particle positions."""
    assert False


@pytest.mark.skip(reason="TODO Week 1")
def test_particle_on_node():
    """A particle sitting exactly on node j deposits all charge on j."""
    assert False


@pytest.mark.skip(reason="TODO Week 1")
def test_particle_at_midpoint():
    """A particle at a cell midpoint splits charge 50/50 between nodes."""
    assert False
