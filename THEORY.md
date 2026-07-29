# Theory: The Physics Behind This Code

This document explains the plasma physics that each component of `pic1d`
discretizes. It is written to be read alongside the source, and it forms the
basis of the theory and numerical method sections of the technical report.

---

## 1. The problem: a self-consistent, collisionless plasma

A plasma is a collection of charged particles that interact through the
electromagnetic fields they themselves create. This self-referential structure
is what makes the problem hard: the particles produce fields, and those same
fields determine how the particles move, which changes the fields again. No
external field is prescribed. The field *is* the particles, smeared out.

The exact description of a collisionless plasma is the **Vlasov-Poisson
system**. Rather than tracking individual particles, it evolves a distribution
function `f(x, v, t)` giving the density of particles at each point in phase
space. In one spatial and one velocity dimension:

```
    df/dt + v (df/dx) + (q/m) E (df/dv) = 0          (Vlasov)

    d2(phi)/dx2 = -rho / eps_0,   E = -d(phi)/dx     (Poisson)

    rho(x,t) = sum_species q integral f dv
```

The Vlasov equation states that `f` is conserved along particle trajectories in
phase space. Poisson's equation closes the system by building the field from the
charge density that `f` implies.

Solving this directly on a phase-space grid (the continuum "Vlasov" approach) is
expensive and suffers from filamentation: structure in velocity space becomes
progressively finer until no practical grid resolves it.

---

## 2. Why particle-in-cell works

PIC avoids the phase-space grid using the **method of characteristics**. Because
`f` is constant along particle trajectories, a set of markers scattered
according to `f` at `t = 0` and then advanced along Newton's equations will
continue to sample `f` correctly at all later times.

These markers are **macroparticles**. A macroparticle is not a single electron;
it stands in for a large number of real particles occupying nearly the same
point in phase space. Simulating 50,000 macroparticles can represent a plasma
containing many orders of magnitude more real electrons, because what is being
sampled is a smooth distribution, not a particle inventory.

The "in-cell" half of the name refers to the second idea: particles live at
continuous positions, but fields are far cheaper to compute on a grid. Each
timestep therefore alternates between the two representations:

```
    particles --[deposit]--> grid
                             grid: solve Poisson, differentiate for E
    particles <--[gather]---  grid
    particles: push (leapfrog)
```

Those four operations are the four components implemented in Week 1.

---

## 3. Charge deposition, and why macroparticles are clouds

`deposition.deposit_charge` maps particles at continuous positions onto grid
node densities using cloud-in-cell (linear) weighting: a particle a fraction `f`
of the way from node `j` to node `j+1` contributes `(1-f)` of its charge to `j`
and `f` to `j+1`.

Two consequences, one numerical and one physical.

**Numerically**, this weighting conserves total charge exactly. The unit test
asserting that `sum(rho) * dx` equals the total particle charge to machine
precision is not a coding nicety; it enforces a conservation law that the
Poisson solve depends on.

**Physically**, linear weighting means each macroparticle is not a point charge
but a finite-width cloud with a triangular shape function. This matters more
than it first appears. Point charges produce singular `1/r` fields and
experience violent close encounters; finite-width clouds do not. The finite
particle size smooths short-range interactions and suppresses artificial
collisions between macroparticles. **A PIC plasma is collisionless precisely
because its particles are clouds.** That is what keeps the method faithful to
the collisionless Vlasov equation it is meant to solve.

---

## 4. The field solve: screening and quasineutrality

Poisson's equation is the statement that charge separation creates potential.
In normalized units (see section 6) it reads `d2(phi)/dx2 = -rho`.

Physically it encodes the defining behavior of a plasma: **Debye screening**.
Introduce a positive test charge into a plasma and electrons swarm toward it
while ions drift away, until the charge is electrically hidden from the rest of
the plasma. The characteristic distance over which this occurs is the **Debye
length**,

```
    lambda_D = sqrt(eps_0 k_B T_e / (n_e e^2))
```

Solving Poisson each step is how the simulation knows that charge imbalance
produces a restoring field.

### The periodic (FFT) solver

In Fourier space each mode decouples:

```
    phi_k = rho_k / k^2
```

This single relation is the screening physics expressed mode by mode: short
wavelengths (large `k`) produce weak potentials, long wavelengths produce strong
ones.

The `k = 0` mode must be handled separately and set to zero. This is not a
numerical hack. It reflects two physical facts: the mean level of an
electrostatic potential is arbitrary (gauge freedom), and a periodic domain must
be globally charge-neutral for the problem to be well posed at all. A
net-charged periodic plasma has no solution, because the field it would generate
cannot satisfy periodicity.

### The Dirichlet (tridiagonal) solver

Bounded problems with walls at fixed potential require a different treatment.
Central differencing on interior nodes gives

```
    phi[i-1] - 2 phi[i] + phi[i+1] = -rho[i] dx^2
```

with the known boundary values folded into the right-hand side of the first and
last interior rows. The resulting tridiagonal system is solved by the Thomas
algorithm in `O(n)` operations. The matrix is diagonally dominant, so no
pivoting is required.

Taylor expansion of the central difference gives

```
    (phi[i-1] - 2 phi[i] + phi[i+1]) / dx^2 = phi'' + (dx^2/12) phi'''' + O(dx^4)
```

so the scheme is second-order accurate. Halving `dx` should reduce the maximum
error by a factor of four, and the convergence test asserts exactly this ratio
rather than a bare tolerance: it verifies the order of accuracy that the
truncation analysis predicts.

This solver is what makes the Week 3 sheath problem possible, since a sheath is
defined by its walls.

---

## 5. The particle push, and why leapfrog is not optional

`pusher.push` advances particles by the leapfrog scheme, with velocity and
position staggered by half a timestep:

```
    v^(n+1/2) = v^(n-1/2) + (q/m) E(x^n) dt
    x^(n+1)   = x^n + v^(n+1/2) dt
```

Leapfrog is chosen over apparently more accurate schemes (such as RK4) for a
structural reason: it is **symplectic**. It exactly preserves the geometric
structure of Hamiltonian mechanics, including phase-space volume as required by
Liouville's theorem.

The practical consequence is decisive. Non-symplectic integrators systematically
gain or lose energy over long runs, so a simulated plasma oscillation would
spuriously damp or grow. A symplectic integrator instead conserves a slightly
perturbed "shadow" Hamiltonian *exactly*, so the true energy oscillates within a
small bounded band indefinitely rather than drifting.

This is what the energy-conservation unit test measures. Energy that wobbles but
does not walk is the signature of a correct symplectic implementation. For a
method whose entire purpose is stable integration over many thousands of steps,
this is the difference between a usable code and a useless one.

### The self-force constraint

`pusher.gather_field` must interpolate the field to particle positions using the
**same** linear weighting as deposition. If the two differ, a particle feels a
force arising from its own deposited charge: an unphysical self-force that
corrupts the dynamics in ways that are difficult to diagnose. The symmetry
between `deposit_charge` and `gather_field` is a physical requirement, not a
stylistic one.

---

## 6. Normalization: the two scales that define a plasma

All quantities in this code are dimensionless. This is not merely numerical
convenience; the chosen scales are the two fundamental scales of plasma physics.

| Quantity  | Normalized by                      |
|-----------|------------------------------------|
| Length    | Debye length `lambda_D`            |
| Time      | Inverse plasma frequency `1/w_pe`  |
| Velocity  | Electron thermal speed `v_th`      |
| Potential | `k_B T_e / e`                      |
| Density   | Reference density `n_0`            |

**The plasma frequency** `w_pe = sqrt(n_e e^2 / (eps_0 m_e))` is the rate at
which electrons collectively respond to charge imbalance. Displace the electrons
from the ions and release them, and they oscillate back at `w_pe`. It is the
plasma's natural ringing frequency and the fastest electrostatic timescale in
the system.

**The Debye length** is its spatial partner, roughly the thermal speed divided
by the plasma frequency: the distance over which charge is screened.

These two scales set the stability constraints enforced in `Simulation`:

```
    dt < ~0.2 / w_pe        resolve plasma oscillations
    dx < ~0.5 lambda_D      avoid the finite-grid instability
    v_max dt < dx           particles must not skip cells
```

The first is straightforward: stepping over the oscillation the plasma is trying
to execute produces instability. The second is subtler. If the grid is too
coarse to resolve charge separation before it is screened, aliasing between the
particle and grid scales spuriously heats the plasma, an artifact known as the
**finite-grid instability** or numerical self-heating. Neither constraint is
arbitrary numerical hygiene: both are the requirement that the discretization
resolve the scales that define a plasma.

---

## 7. What this machinery makes accessible

### The two-stream instability (Week 2)

Two electron populations streaming through one another are unstable: free energy
stored in their relative motion feeds growing electrostatic waves. This is a
**kinetic** instability, meaning it depends on the shape of the velocity
distribution rather than on bulk fluid quantities. A fluid model averages away
precisely the information that drives it. PIC captures it because PIC represents
the distribution directly.

The phase-space vortices that appear at saturation are the nonlinear end state:
electrons become trapped in the potential wells of the waves they generated. The
growth rate during the linear phase can be compared against the cold two-stream
dispersion relation, providing a quantitative validation of the whole loop.

### The plasma sheath (Weeks 3-4)

Wherever a plasma contacts a wall, a thin non-neutral layer forms. Electrons,
being far lighter and faster than ions, reach the wall first and charge it
negatively; the resulting potential drop then repels further electrons until
electron and ion fluxes balance. This layer is the **sheath**, and it is
typically a few Debye lengths thick.

The sheath governs the energy with which ions strike the wall, which in turn
determines sputtering erosion rates and wall heat loading. In a Hall thruster,
these set the channel wall lifetime and the thermal design constraints, making
sheath physics a direct limit on thruster durability.

The **Bohm criterion**, which requires that ions enter the sheath at or above
the ion sound speed `sqrt(k_B T_e / m_i)`, is one of the cleaner results in
plasma physics. In a kinetic simulation it is not imposed: it emerges
self-consistently from the particle dynamics. Recovering it is strong evidence
that the simulation is capturing the correct physics.

Because ions are far heavier than electrons, sheath simulations conventionally
use a **reduced mass ratio** (for example `m_i/m_e = 100` or `400`) to keep
runtimes tractable. This is standard practice and is stated explicitly in the
results, since the sheath potential drop depends on the mass ratio used.

---

## 8. Assumptions and limits of this model

Worth stating plainly, both for the report and because these define the
boundaries of what the results can claim:

- **Electrostatic.** No magnetic field, no radiation, no retarded potentials.
  Valid when particle velocities are far below `c` and collective electric
  interactions dominate. Real Hall thrusters are strongly magnetized, so this
  code models sheath physics rather than a complete thruster.
- **Collisionless.** No Coulomb collisions, no neutral collisions, no
  ionization. Adding Monte Carlo collisions (MCC) is the natural extension.
- **One-dimensional.** Captures the essential physics of sheaths and streaming
  instabilities, which are genuinely 1D phenomena, but excludes cross-field
  transport and azimuthal instabilities.
- **Finite macroparticle count.** Statistical noise scales roughly as
  `1/sqrt(N_particles)`. Results are quoted with the particle count used.

---

## References

- C. K. Birdsall and A. B. Langdon, *Plasma Physics via Computer Simulation*
- R. W. Hockney and J. W. Eastwood, *Computer Simulation Using Particles*
- L. Brieda, *Plasma Simulations by Example*
- F. F. Chen, *Introduction to Plasma Physics and Controlled Fusion*
  (Debye screening, sheaths, the Bohm criterion)
