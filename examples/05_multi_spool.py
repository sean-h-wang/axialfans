# ═══════════════════════════════════════════════════════════
# Example 5 — Multi-spool compressor (LP + stator + HP)
# Mimics a jet engine axial compressor architecture:
# LP spool (3 rotor-stator pairs) → HP spool (3 rotor-stator pairs)
# ═══════════════════════════════════════════════════════════
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axialfans.fan_solver import State, MultistageFanSolver

import numpy as np
M    = 25
T0   = 288.15; P0 = 101325.0; rho0 = P0 / (287.05 * T0)
R    = 287.05; cp = 1004.7; gamma = 1.4
vax0 = 80.0

#                   R1   S1   R2   S2   R3   S3   R4   S4   R5   S5   R6   S6
solver = MultistageFanSolver(
    N=12, M=M,
    direction=       [1,  -1,   1,  -1,   1,  -1,   1,  -1,   1,  -1,   1,  -1],
    omega=   [1200,   0,1200,   0,1200,   0,2400,   0,2400,   0,2400,   0],
    sigma=   [0.86, 1.0, 0.86, 1.0, 0.86, 1.0, 0.88, 1.0, 0.88, 1.0, 0.88, 1.0],
    beta=    [  50,  20,   47,  20,   44,  20,   41,  18,   38,  18,   35,  18],
    rp=      [0.30, 0.30, 0.29, 0.29, 0.28, 0.28, 0.25, 0.25, 0.23, 0.23, 0.21, 0.21],
    rm=      [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.18, 0.18, 0.17, 0.17, 0.16, 0.16],
    eta=     [0.85, 0.98, 0.85, 0.98, 0.85, 0.98, 0.87, 0.98, 0.87, 0.98, 0.87, 0.98],
    R=R, cp=cp, gamma=gamma
)
inlet = State(0, M, 0.20, 0.30, T0, P0, vax0, rho0, 0)
solver.solve(inlet)
perf = solver.performance()

# LP = stages 1-6 (3 rotor-stator pairs), HP = stages 7-12
PR_LP = solver.states[6].P_avg  / solver.states[0].P_avg
PR_HP = solver.states[12].P_avg / solver.states[6].P_avg

print("Example 5 — Multi-spool compressor")
print("=" * 40)
print(f"LP spool PR (3 stages): {PR_LP:.3f}")
print(f"HP spool PR (3 stages): {PR_HP:.3f}")
print(f"Overall PR:             {perf['PR']:.3f}")
print(f"Overall TR:             {perf['TR']:.3f}")
print(f"Overall eta_isen:       {perf['eta_isen']:.3f}")