# pic1d: A 1D Electrostatic Particle-in-Cell Plasma Simulation

A 1D electrostatic PIC code written from scratch in Python (NumPy-vectorized),
validated against classic kinetic plasma problems and applied to self-consistent
plasma sheath formation and wall heat flux estimation, a problem directly
relevant to Hall thruster channel wall erosion and thermal loading.

**Author:** Rafael Cedeno

## Physics roadmap

1. **Cold plasma oscillations**: verify oscillation at the plasma frequency.
2. **Two-stream instability**: verify linear growth rate against theory,
   phase-space vortex formation.
3. **Plasma sheath**: absorbing walls, kinetic electrons and ions, sheath
   potential drop vs. floating potential theory, Bohm criterion check.
4. **Wall heat flux**: particle and energy fluxes to the wall, transient
   wall conduction estimate. Ties plasma results to thermal analysis.

## Normalization

All quantities are dimensionless from day one:

| Quantity  | Normalized by                      |
|-----------|------------------------------------|
| Length    | Debye length (lambda_D)            |
| Time      | Inverse plasma frequency (1/w_pe)  |
| Velocity  | Electron thermal speed (v_th)      |
| Potential | k_B T_e / e                        |
| Density   | Reference density n_0              |

With these choices, the normalized Poisson equation is
`d^2(phi)/dx^2 = -(n_i - n_e)` and the electron equation of motion is
`dv/dt = E` (sign conventions handled in `pusher.py`).

## Repo structure

```
src/pic1d/
    grid.py         # Grid class: domain, spacing, field arrays
    particles.py    # Species class: positions, velocities, charge, mass
    deposition.py   # Cloud-in-cell (linear) charge weighting
    poisson.py      # 1D Poisson solvers (periodic FFT + Dirichlet tridiagonal)
    pusher.py       # Leapfrog particle push + field gather
    boundary.py     # Periodic wrap, absorbing walls, particle injection
    simulation.py   # Main PIC loop orchestration
    diagnostics.py  # Energy history, phase space, flux counters
tests/              # Unit tests for each component (run: pytest)
examples/           # Driver scripts for each physics milestone
report/             # Technical report + generated figures
```

## Quick start

```bash
pip install -r requirements.txt
pytest                                  # unit tests
python examples/run_plasma_oscillation.py
python examples/run_two_stream.py
python examples/run_sheath.py
```

## Development order (matches the 4-week plan)

- [ ] Week 1: implement `pusher.py`, `poisson.py`, `deposition.py`; all unit tests pass
- [ ] Week 2: assemble `simulation.py` with periodic BCs; plasma oscillation + two-stream validation
- [ ] Week 3: absorbing walls + injection in `boundary.py`; sheath formation
- [ ] Week 4: wall heat flux diagnostics; technical report; one-command figure reproduction

## Key references

- Birdsall & Langdon, *Plasma Physics via Computer Simulation*
- L. Brieda, *Plasma Simulations by Example* (practical Python walkthrough)
- Hockney & Eastwood, *Computer Simulation Using Particles*
