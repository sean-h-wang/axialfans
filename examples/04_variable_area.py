"""Variable area (converging duct) example."""
from axialfans.fan_solver import MultistageFanSolver

# Area decreases between stages
solver = MultistageFanSolver(
    N=2,
    direction=[1, -1],
    sigma=0.9,
    omega=[3000, 0],
    beta=[65, 20],
    rp=[0.252, 0.100],      # Decreasing tip radius
    rm=[0.177, 0.050],      # Decreasing hub radius
    eta=0.7,
    R=287,
    cp=1005
)

solver.solve(T0=288, vax0=10, P0=101325, rho0=1.225, verbose=True)

# Compare velocities
import numpy as np
A1 = np.pi * (solver.rp[1]**2 - solver.rm[1]**2)
A2 = np.pi * (solver.rp[2]**2 - solver.rm[2]**2)

print("\n" + "="*60)
print("AREA CHANGE EFFECTS:")
print(f"Area 1: {A1:.4f} m²")
print(f"Area 2: {A2:.4f} m²")
print(f"Area ratio (A2/A1): {A2/A1:.3f}")
print(f"\nStage 1 exit velocity: {solver.vax[1]:.2f} m/s")
print(f"Stage 2 exit velocity: {solver.vax[2]:.2f} m/s")
print(f"Velocity increased by: {solver.vax[2]/solver.vax[1]:.2f}x")
print(f"\nPressure ratio: {solver.P[2]/solver.P[0]:.2f}")