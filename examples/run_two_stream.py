"""Week 2 validation #2: two-stream instability.

Setup: two counter-streaming electron beams (drift = +/- 3 v_th, cold),
uniform ion background, periodic domain, seeded with a tiny sinusoidal
perturbation at the fastest-growing wavenumber.

Measure: exponential growth rate of field energy in the linear phase;
compare against the cold two-stream dispersion relation. Then watch
saturation and phase-space vortex formation.

Outputs:
    report/figures/two_stream_growth_rate.png
    report/figures/two_stream_phase_space.gif   (the money plot)
    report/figures/energy_conservation.png      (total drift < 1%)
TODO (Week 2).
"""
