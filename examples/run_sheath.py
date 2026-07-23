"""Week 3: self-consistent plasma sheath + Week 4 wall heat flux.

Setup: bounded domain with absorbing walls at phi = 0 (grounded) or
floating; kinetic electrons AND kinetic ions with reduced mass ratio
m_i/m_e = 100 (state this in the report; standard practice). Steady
injection from a central source region or from the boundaries keeps
the bulk plasma sustained.

Validate against theory:
    1. Sheath potential drop vs analytical floating potential
       (depends on mass ratio; derive for m_i/m_e = 100).
    2. Density profiles: quasi-neutral presheath, non-neutral sheath
       (a few Debye lengths thick).
    3. Bohm criterion: ion mean velocity at the sheath edge >= Bohm
       speed sqrt(T_e/m_i).

Then (Week 4): time-averaged particle and energy fluxes to the wall
from WallFluxCounter -> wall heat flux -> 1D transient conduction
estimate of wall temperature rise. Discuss relevance to Hall thruster
channel walls, citing Jorns / Hara erosion papers.

Outputs:
    report/figures/sheath_potential_profile.png
    report/figures/sheath_density_profiles.png
    report/figures/bohm_criterion.png
    report/figures/wall_heat_flux.png
TODO (Weeks 3-4).
"""
