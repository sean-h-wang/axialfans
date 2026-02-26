# ═══════════════════════════════════════════════════════════
# Example 4 — Variable area (converging duct)
# ═══════════════════════════════════════════════════════════
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axialfans.fan_solver import State, MultistageFanSolver

import numpy as np
M    = 25
T0   = 288.15; P0 = 101325.0; rho0 = P0 / (287.05 * T0)
R    = 287.05; cp = 1004.7; gamma = 1.4
vax0 = 10.0

solver = MultistageFanSolver(
    N=2, M=M, direction=[1, -1],
    sigma=0.9, omega=[3000, 0], beta=[65, 20],
    rp=[0.252, 0.100], rm=[0.177, 0.050],
    eta=0.7, R=R, cp=cp, gamma=gamma
)
inlet = State(0, M, 0.177, 0.252, T0, P0, vax0, rho0, 0)
solver.solve(inlet)
perf = solver.performance()

A1 = np.pi * (0.252**2 - 0.177**2)
A2 = np.pi * (0.100**2 - 0.050**2)
print("Example 4 — Variable area")
print("=" * 40)
print(f"Area ratio (A2/A1): {A2/A1:.3f}")
print(f"PR:                 {perf['PR']:.3f}")
print(f"TR:                 {perf['TR']:.3f}")
print(f"Exit vax:           {solver.states[-1].vax_avg:.2f} m/s")