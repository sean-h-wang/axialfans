# axialfans

**Fast analytical solver for multistage axial turbomachinery**

[![PyPI version](https://badge.fury.io/py/axialfans.svg)](https://pypi.org/project/axialfans/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Generalized solver for arbitrary sequences of rotating and stationary axial stages (rotors, stators, counter-rotating fans) using Newton-Raphson iteration on coupled thermodynamic equations. Designed for preliminary turbomachinery design and rapid parameter exploration.

**Key Features:**
- ⚡ **Fast**: 333 solves/second enables real-time Monte Carlo (10,000 samples in ~30s)
- 🎯 **Validated**: 6% accuracy vs NASA experimental data (Rotor 37, Report 729)
- 🔧 **General**: Rotors, stators, counter-rotation, variable area in unified framework
- 📊 **Production-ready**: Type hints, error handling, comprehensive documentation

---

## Installation
```bash
pip install axialfans
```

**Requirements:** Python 3.9+, NumPy

---

## Quick Start
```python
from axialfans.fan_solver import MultistageFanSolver

# Single stage fan
solver = MultistageFanSolver(
    N=1,                # Number of stages
    direction=[1],      # Rotation: +1=CCW, -1=CW
    sigma=0.9,          # Slip factor
    omega=1800,         # Angular velocity [rad/s]
    beta=45,            # Blade angle [degrees]
    rp=0.25,            # Tip radius [m]
    rm=0.18,            # Hub radius [m]
    eta=0.85,           # Efficiency
    R=287,              # Gas constant [J/kg-K]
    cp=1005,            # Specific heat [J/kg-K]
    gamma=1.4           # Heat capacity ratio
)

# Solve from inlet conditions
solver.solve(
    T0=288,         # Temperature [K]
    vax0=50,        # Axial velocity [m/s]
    P0=101325,      # Pressure [Pa]
    rho0=1.225      # Density [kg/m³]
)

# Results
print(f"Pressure ratio: {solver.P[1]/solver.P[0]:.2f}")
print(f"Exit velocity: {solver.vax[1]:.1f} m/s")
```

---

## Examples

See the [`examples/`](examples/) directory for detailed usage:

| Example | Description |
|---------|-------------|
| `01_single_rotor.py` | Basic single-stage fan |
| `02_rotor_stator.py` | Fan with straightening vanes (NASA 729 config) |
| `03_counter_rotating.py` | Swirl cancellation with opposite-spinning rotors |
| `04_variable_area.py` | Converging/diverging ducts |
| `05_multi_spool.py` | LP + stator + HP compressor configuration |
| `06_monte_carlo_uq.py` | Uncertainty quantification (1000 samples) |

---

## Key Concepts

### Indexing Convention
State arrays indexed 0 → N:
- `solver.T[0]` = inlet temperature
- `solver.T[n]` = temperature exiting stage n

### Direction Array
Specifies rotation sense:
- `direction=[1]` = counter-clockwise rotation
- `direction=[-1]` = clockwise rotation / mirror frame
- `direction=[1, -1]` = rotor + stator pair

### Stage Types
Configure via parameters:
- **Rotor**: `omega > 0`, any `direction`
- **Stator**: `omega = 0`, set opposite `direction` from previous
- **Counter-rotating**: Adjacent stages with opposite `direction` values

### Variable Area
Different radii per stage:
```python
rp=[0.30, 0.25, 0.20]  # Decreasing tip radius
rm=[0.20, 0.18, 0.15]  # Decreasing hub radius
```
Area decrease → velocity increase, pressure ratio decrease

---

## Performance

- **Single solve**: ~3 ms
- **1,000 solves**: ~3 seconds
- **10,000 solves**: ~30 seconds
- **4-5 orders of magnitude faster than RANS CFD**

Enables real-time design space exploration and Monte Carlo uncertainty quantification.

---

## Validation

Validated against NASA experimental data:

| Test Case | Configuration | Error |
|-----------|---------------|-------|
| **Rotor 37** | High-speed single rotor (Re ~ 10⁶) | 6.2% on pressure ratio |
| **NASA 729** | Low-speed rotor-stator (Re ~ 10⁵) | Correctly identifies inviscid breakdown |

Uncertainty quantification with 10,000 Monte Carlo samples shows 99.8% convergence rate. See [`notebooks/rotor37_validation.ipynb`](notebooks/rotor37_validation.ipynb) for detailed analysis.

---

## Documentation

Full API documentation available in code docstrings. Quick reference:
```python
MultistageFanSolver(
    N,              # Number of stages
    direction,      # Rotation sense: +1 (CCW) or -1 (CW)
    sigma,          # Slip factor (0.85-0.95 typical)
    omega,          # Angular velocity [rad/s] (unsigned, use 0 for stators)
    beta,           # Blade angle [degrees] (unsigned)
    rp, rm,         # Tip and hub radius [m]
    eta,            # Efficiency (0-1)
    R, cp, gamma    # Gas properties
)
```

All parameters accept scalars (same for all stages) or arrays (per-stage values).

---

## Citation

If you use this software in your research, please cite:
```bibtex
@software{wang2026axialfans,
  title={axialfans: Fast Analytical Solver for Multistage Axial Turbomachinery},
  author={Wang, Sean},
  year={2026},
  url={https://github.com/sean-h-wang/axialfans},
  version={1.0.0}
}
```

Associated paper:
```bibtex
@article{wang2026analytical,
  title={Analytical Solver for Multistage Axial Turbomachinery with Variable Area and Counter-Rotation},
  author={Wang, Sean},
  journal={AIAA Journal},
  year={2026},
  note={submitted}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

## Contributing

Issues and pull requests welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

For bug reports, include:
- Python version
- NumPy version
- Minimal reproducible example
- Error message / unexpected behavior

---

## Author

**Sean H. Wang**
- Email: [seanwangpiano@gmail.com](mailto:seanwangpiano@gmail.com)
- GitHub: [@sean-h-wang](https://github.com/sean-h-wang)

Developed as part of hypersonic propulsion research.