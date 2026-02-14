# Changelog

All notable changes to this multistage axial fan analytical solver are documented in this file.

The format follows Keep a Changelog and Semantic Versioning.

---

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

