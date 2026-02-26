# ═══════════════════════════════════════════════════════════
# Example 1 — Single rotor
# ═══════════════════════════════════════════════════════════
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axialfans.fan_solver import State, MultistageFanSolver

import numpy as np
import axialfans.fan_solver as solver_module
solver_module.MAX_ITER = 1000
solver_module.TOL      = 1e-8
solver_module.ALPHA    = 1.0

M    = 25
T0   = 288.15; P0 = 101325.0; rho0 = P0 / (287.05 * T0)
R    = 287.05; cp = 1004.7; gamma = 1.4
vax0 = 50.0

solver = MultistageFanSolver(
    N=1, M=M, direction=[1],
    sigma=0.9, omega=1800, beta=45,
    rp=0.25, rm=0.18, eta=0.85,
    R=R, cp=cp, gamma=gamma
)
inlet = State(0, M, 0.18, 0.25, T0, P0, vax0, rho0, 0)
solver.solve(inlet)
perf = solver.performance()
print("Example 1 — Single rotor")
print("=" * 40)
print(f"PR:  {perf['PR']:.3f}")
print(f"TR:  {perf['TR']:.3f}")
print(f"eta: {perf['eta_isen']:.3f}")