# ═══════════════════════════════════════════════════════════
# Example 2 — Rotor-stator pair
# ═══════════════════════════════════════════════════════════
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axialfans.fan_solver import State, MultistageFanSolver

import numpy as np
M    = 25
T0   = 288.15; P0 = 101325.0; rho0 = P0 / (287.05 * T0)
R    = 287.05; cp = 1004.7; gamma = 1.4
vax0 = 17.8

solver = MultistageFanSolver(
    N=2, M=M, direction=[1, -1],
    sigma=[0.85, 0.90], omega=[1800, 0], beta=[65, 20],
    rp=0.2667, rm=0.18415, eta=0.877,
    R=R, cp=cp, gamma=gamma
)
inlet = State(0, M, 0.18415, 0.2667, T0, P0, vax0, rho0, 0)
solver.solve(inlet)
perf = solver.performance()
print("Example 2 — Rotor-stator pair")
print("=" * 40)
print(f"Overall PR:  {perf['PR']:.3f}")
print(f"Rotor PR:    {solver.states[1].P_avg / solver.states[0].P_avg:.3f}")
print(f"Stator PR:   {solver.states[2].P_avg / solver.states[1].P_avg:.3f}")
print(f"TR:          {perf['TR']:.3f}")