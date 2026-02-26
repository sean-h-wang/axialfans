# axialfans

**Spanwise-resolved analytical solver for multistage axial turbomachinery**

[![PyPI version](https://badge.fury.io/py/axialfans.svg)](https://pypi.org/project/axialfans/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Quasi-1D pointwise solver for arbitrary sequences of rotating and stationary axial blade rows (rotors, stators, counter-rotating fans) using Newton-Raphson iteration on coupled thermodynamic equations at each radial station. Designed for preliminary turbomachinery design and uncertainty quantification.

---

## Installation

```bash
pip install axialfans
```

**Requirements:** Python 3.9+, NumPy

---

## Quick Start

```python
from axialfans.fan_solver import MultistageFanSolver, State
import numpy as np

M = 25  # radial stations

solver = MultistageFanSolver(
    N=1, M=M,
    direction=[1],
    sigma=0.9,
    omega=1800,
    beta=45,
    rp=0.25, rm=0.18,
    eta=0.85,
    R=287.05, cp=1004.7, gamma=1.4
)

T0   = 288.15
P0   = 101325.0
rho0 = P0 / (287.05 * T0)

inlet = State(0, M, 0.18, 0.25, T0, P0, 50.0, rho0, 0)
solver.solve(inlet)

perf = solver.performance()
print(f"Pressure ratio:        {perf['PR']:.3f}")
print(f"Temperature ratio:     {perf['TR']:.3f}")
print(f"Isentropic efficiency: {perf['eta_isen']:.3f}")
```

---

## Examples

See the [`examples/`](examples/) directory:

| Example | Description |
|---------|-------------|
| `01_single_rotor.py` | Basic single-stage fan |
| `02_rotor_stator.py` | Fan with straightening vanes |
| `03_counter_rotating.py` | Swirl cancellation with opposite-spinning rotors |
| `04_variable_area.py` | Converging/diverging ducts |
| `05_multi_spool.py` | LP + stator + HP compressor (jet engine architecture) |
| `06_monte_carlo_uq.py` | Uncertainty quantification with LHS |

---

## Key Concepts

### State Object

All inlet conditions are passed as a `State`:

```python
inlet = State(
    n=0,        # stage index (0 = inlet)
    M=25,       # radial stations
    rm=0.18,    # hub radius [m]
    rp=0.25,    # tip radius [m]
    T=288.15,   # temperature [K]
    P=101325,   # pressure [Pa]
    vax=50.0,   # axial velocity [m/s]
    rho=1.225,  # density [kg/m³]
    vtheta=0    # tangential velocity [m/s]
)
```

### Spanwise Resolution

`M` radial stations are solved independently at each blade row. Spanwise variation arises from blade twist. Pass a 2D `(N x M)` beta array for full control, or `np.linspace(beta_hub, beta_tip, M)` for linear twist:

```python
beta_row = np.linspace(65, 71, M)           # single twisted row
beta_matrix = np.array([row1, row2])        # N x M for multistage
```

### Slip Factor

Supplied externally before calling `solve()`. The **modified Stodola model** developed in Wang (2026) is recommended for axial machines:
```python
sigma = 1 - C / Z   # Z = blade count, C calibrated from NASA Rotor 37
```

Unlike Qiu, Wiesner, and classical Stodola — which were derived for centrifugal machines and degenerate toward σ ≈ 1 for low-twist axial geometries — this model explicitly encodes blade count as the primary slip mechanism, which is physically correct for axial blade rows. It applies consistently to both rotors and stators.

### Direction Array

| Configuration | direction |
|---------------|-----------|
| Single rotor | `[1]` |
| Rotor + stator | `[1, -1]` |
| Counter-rotating pair | `[1, -1]` |
| LP + stator + HP | `[1, -1, 1, -1, ...]` |

Stators have `omega=0`. Direction still alternates to correctly handle the velocity triangle frame at each interface.

### Variable Area

Pass per-stage radii as arrays:

```python
rp=[0.30, 0.25, 0.20]
rm=[0.20, 0.18, 0.15]
```

### Tunable Solver Constants

```python
import axialfans.fan_solver as solver_module

solver_module.MAX_ITER = 2000   # Newton-Raphson iteration limit
solver_module.TOL      = 1e-10  # convergence tolerance
solver_module.ALPHA    = 0.5    # line search step (reduce if NR diverges)
```

---

## Performance

- ~10,000 Monte Carlo evaluations in ~20 minutes on a standard workstation
- Enables UQ studies that would be computationally prohibitive with RANS CFD
- A single RANS CFD evaluation for an equivalent geometry takes hours on a cluster

---

## Validation

Validated against NASA experimental data across two independent transonic rotors:

| Test Case | Configuration | PR Error |
|-----------|---------------|----------|
| **NASA Rotor 37** | Transonic single rotor, deterministic | 2.0% |
| **NASA Rotor 67** | Transonic single rotor, uncertain geometry, 50,000-sample UQ | 2.54% |
| **NACA TR 729** | Low-speed rotor-stator, regime of validity characterization | — |

The NACA TR 729 case identifies the solver's validity boundary: the inviscid compressible model breaks down at low rotor speed and high inlet axial velocity where viscous effects dominate.

---

## Important Note

This solver is a **preliminary design and research tool**. Outputs must be independently verified before any engineering application. See LICENSE for full liability disclaimer.

---

## Citation

If you use this software in your research, please cite:

```bibtex
@software{wang2026axialfans,
  title={axialfans: Spanwise-Resolved Analytical Solver for Multistage Axial Turbomachinery with Uncertainty Quantification},
  author={Wang, Sean H.},
  year={2026},
  url={https://github.com/sean-h-wang/axialfans},
  version={2.0.0}
}

@article{wang2026analytical,
  title={Spanwise-Resolved Analytical Solver for Multistage Axial Turbomachinery with Uncertainty Quantification},
  author={Wang, Sean H.},
  journal={AIAA Journal},
  year={2026},
  note={submitted}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contributing

Issues and pull requests welcome. For bug reports include Python version, NumPy version, and a minimal reproducible example.

---

## Author

**Sean H. Wang**
- Personal: [seanwangpiano@gmail.com](mailto:seanwangpiano@gmail.com)
- Institutional: [swang122@montgomerycollege.edu](mailto:swang122@montgomerycollege.edu)
- GitHub: [@sean-h-wang](https://github.com/sean-h-wang)