"""Single rotor (basic fan) example."""
from axialfans.fan_solver import MultistageFanSolver

# Simple axial fan
solver = MultistageFanSolver(
    N=1,
    direction=[1],      # Counter-clockwise
    sigma=0.9,
    omega=1800,         # 1800 rad/s
    beta=45,            # 45° blade angle
    rp=0.25,
    rm=0.18,
    eta=0.85,
    R=287,
    cp=1005
)

solver.solve(T0=288, vax0=50, P0=101325, rho0=1.225, verbose=True)

# Extract results
print("\n" + "="*60)
print("RESULTS:")
print(f"Pressure ratio: {solver.P[1]/solver.P[0]:.3f}")
print(f"Temperature ratio: {solver.T[1]/solver.T[0]:.3f}")
print(f"Exit velocity: {solver.vax[1]:.2f} m/s")
print(f"Exit density: {solver.rho[1]:.4f} kg/m³")