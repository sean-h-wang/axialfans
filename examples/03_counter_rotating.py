# ═══════════════════════════════════════════════════════════
# Example 3 — Counter-rotating pair
# ═══════════════════════════════════════════════════════════
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axialfans.fan_solver import State, MultistageFanSolver

import numpy as np
M    = 25
T0   = 288.15; P0 = 101325.0; rho0 = P0 / (287.05 * T0)
R    = 287.05; cp = 1004.7; gamma = 1.4
vax0 = 50.0

solver = MultistageFanSolver(
    N=2, M=M, direction=[1, -1],
    sigma=0.9, omega=[1600, 1600], beta=45,
    rp=0.25, rm=0.18, eta=0.88,
    R=R, cp=cp, gamma=gamma
)
inlet = State(0, M, 0.18, 0.25, T0, P0, vax0, rho0, 0)
solver.solve(inlet)
perf = solver.performance()
print("Example 3 — Counter-rotating pair")
print("=" * 40)
print(f"Stage 1 PR:  {solver.states[1].P_avg / solver.states[0].P_avg:.3f}")
print(f"Stage 2 PR:  {solver.states[2].P_avg / solver.states[1].P_avg:.3f}")
print(f"Overall PR:  {perf['PR']:.3f}")
print(f"TR:          {perf['TR']:.3f}")