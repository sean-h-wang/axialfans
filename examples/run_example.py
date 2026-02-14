from pathlib import Path
import sys

# Add project root (parent of this file's parent) to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from axialfans.fan_solver import MultistageFanSolver

solver = MultistageFanSolver(
    N=2,
    direction=[1,-1],
    sigma=0.9,
    omega=[3000,0],
    beta=[65,20],
    rp=[0.252,0.1],
    rm=[0.177,0.05],
    eta=0.7,
    R=287,
    cp=1005
)
solver.solve(T0=288, vax0=10, P0=101325, rho0=1.225, verbose=True)

# Check output
print(f"Final pressure ratio: {solver.P[-1] / solver.P[0]:.2f}")