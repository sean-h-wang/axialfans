"""Rotor-stator pair."""
from axialfans.fan_solver import MultistageFanSolver

# Fan with straightening vanes
solver = MultistageFanSolver(
    N=2,
    direction=[1, -1],      # Rotor CCW, stator in mirror frame
    sigma=[0.85, 0.90],
    omega=[1800, 0],        # Stator has zero rotation
    beta=[65, 20],          # Steep rotor, shallow stator
    rp=0.2667,
    rm=0.18415,
    eta=0.877,
    R=287.053,
    cp=1005
)

solver.solve(T0=288.15, vax0=17.8, P0=101325, rho0=1.229, verbose=True)

# Results
print("\n" + "="*60)
print("PERFORMANCE BREAKDOWN:")
print(f"Overall pressure ratio: {solver.P[2]/solver.P[0]:.3f}")
print(f"Rotor pressure ratio: {solver.P[1]/solver.P[0]:.3f}")
print(f"Stator pressure ratio: {solver.P[2]/solver.P[1]:.3f}")
print(f"\nOverall temperature ratio: {solver.T[2]/solver.T[0]:.3f}")
print(f"Exit velocity: {solver.vax[2]:.2f} m/s")