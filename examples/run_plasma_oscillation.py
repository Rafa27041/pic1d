"""Week 2 validation #1: cold plasma oscillations.

Setup: electrons with a small sinusoidal density perturbation on a
uniform neutralizing ion background, periodic domain.
Expected: electric field energy oscillates at 2*w_pe (energy oscillates
at twice the field frequency); measured w_pe should match 1.0 in
normalized units to within a few percent.

Output: report/figures/plasma_oscillation_frequency.png
    (field energy vs time + FFT peak vs theory)

Suggested parameters:
    L = 2*pi * 4        (a few wavelengths)
    n_nodes = 257
    n_particles = 50_000
    dt = 0.05
    perturbation amplitude = 0.01
TODO (Week 2).
"""
