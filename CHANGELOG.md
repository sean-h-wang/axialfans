# Changelog

All notable changes to this multistage axial fan analytical solver are documented in this file.

The format follows Keep a Changelog and Semantic Versioning.

---

# Changelog

## [2.0.0] — 2026-02-26

### Breaking Changes
- Solver now requires a `State` object for inlet conditions — scalar `T0`, `P0`, `vax0`, `rho0` kwargs removed
- `beta` now accepts 2D `(N x M)` arrays for full spanwise twist
- `sigma` is now user-supplied externally; internal Qiu model removed
- Result access via `solver.states[n].P_avg` etc. — direct array attributes `solver.P`, `solver.T` removed

### Added
- Pointwise spanwise resolution: `M` independent radial stations per blade row
- Area-weighted averaging applied after nonlinear solve (correct Jensen ordering)
- Module-level tunable constants: `MAX_ITER`, `TOL`, `ALPHA`
- `State.summary()` and `State.__repr__()` for diagnostics
- `performance()` returns `PR`, `TR`, and `eta_isen`

### Improved
- Pressure ratio error reduced from ~9% to <3% on NASA Rotor 37 and Rotor 67
- Newton-Raphson uses analytical Jacobian throughout — no finite differences
- Variable area handled via local vax scaling, preserving continuity across stage interfaces


## [1.0.2] - 2026-02-14

### Fixed
- Corrected error in slip-factor ($\sigma$) contribution to the tangential velocity formulation.
- Ensured consistency between reduced parameter definition and integrated energy equation.
- Resolved discrepancy affecting stage work computation under certain operating conditions.

### Notes
- No structural changes to governing equations.
- No modifications to residual system or Newton solver.
- Numerical results may differ slightly from previous versions due to corrected $\sigma$ term.

---

## [1.0.1] - 2026-02-14

### Fixed
- Corrected area-scaling output behavior.
- Prevented area adjustment message from printing when no change in annular area occurred between stages.

### Notes
- No changes to governing equations or numerical formulation.
- Purely a logic/output correction.

---

## [1.0.0] - 2026-02-14

### Fixed
- Corrected handling of rotation-direction changes between sequential stages.
- Properly implemented temporary sign reversal of swirl parameters ($\Psi$, $\Xi$) during counter-rotation transitions.

### Notes
- Stabilized multistage counter-rotating configurations.
- Established first fully stable multistage implementation.

