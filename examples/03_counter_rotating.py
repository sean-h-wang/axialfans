"""Counter-rotating fans (reduced exit swirl, high compression)."""
from axialfans.fan_solver import MultistageFanSolver
import numpy as np

# Two rotors spinning opposite directions
solver = MultistageFanSolver(
    N=2,
    direction=[1, -1],      # Opposite rotations
    sigma=0.9,
    omega=[1800, 1800],     # Same speed, opposite directions
    beta=45,                # Same blade angle
    rp=0.25,
    rm=0.18,
    eta=0.88,
    R=287,
    cp=1005
)

solver.solve(T0=288, vax0=50, P0=101325, rho0=1.225, verbose=True)

# Compute swirl in global frame
# Stage 1 adds positive swirl
# Stage 2 adds negative swirl (direction flipped)
r_avg = (solver.rp[2] + solver.rm[2]) / 2

# Local frame swirl coefficients (always positive)
print("\n" + "="*60)
print("LOCAL FRAME SWIRL COEFFICIENTS:")
print(f"Stage 1 Psi: {solver.Psi[1]:.1f} (positive in CCW frame)")
print(f"Stage 2 Psi: {solver.Psi[2]:.1f} (positive in CW mirror frame)")

# Global frame tangential velocities (estimate at mean radius)
v_theta_1 = solver.Psi[1] * r_avg + solver.Xi[1]
v_theta_2_local = solver.Psi[2] * r_avg + solver.Xi[2]
v_theta_2_global = -v_theta_2_local  # Flip sign for direction change

print("\nGLOBAL FRAME TANGENTIAL VELOCITIES (at r={:.3f}m):".format(r_avg))
print(f"After stage 1: v_θ = {v_theta_1:.2f} m/s (CCW)")
print(f"After stage 2: v_θ ≈ {v_theta_2_global:.2f} m/s (sign flipped)")
print(f"Swirl reduction: ~{abs(v_theta_1 - abs(v_theta_2_global)):.1f} m/s")

print("\nCOMPRESSION PERFORMANCE:")
print(f"Stage 1 PR: {solver.P[1]/solver.P[0]:.3f}")
print(f"Stage 2 PR: {solver.P[2]/solver.P[1]:.3f}")
print(f"Overall PR: {solver.P[2]/solver.P[0]:.3f} (very high!)")
print(f"\nBoth stages contribute to compression.")
print(f"Counter-rotation reduces exit swirl but doesn't eliminate it.")