"""Multi-spool compressor (LP + stator + HP)."""
from axialfans.fan_solver import MultistageFanSolver

# Three LP stages, one stator, two HP stages
solver = MultistageFanSolver(
    N=6,
    direction=[1, 1, 1, -1, 1, 1],
    sigma=0.88,
    omega=[1200, 1200, 1200, 0, 2400, 2400],
    beta=[50, 45, 40, 25, 38, 35],
    rp=[0.30, 0.28, 0.26, 0.24, 0.22, 0.20],
    rm=[0.20, 0.19, 0.18, 0.17, 0.16, 0.15],
    eta=[0.85, 0.86, 0.87, 0.90, 0.88, 0.87],
    R=287,
    cp=1005
)

solver.solve(T0=288, vax0=80, P0=101325, rho0=1.225, verbose=True)

# Analyze spool contributions
print("\n" + "="*60)
print("MULTI-SPOOL PERFORMANCE:")
print(f"LP section (stages 1-3): PR = {solver.P[3]/solver.P[0]:.3f}")
print(f"Stator (stage 4): PR = {solver.P[4]/solver.P[3]:.3f}")
print(f"HP section (stages 5-6): PR = {solver.P[6]/solver.P[4]:.3f}")
print(f"\nOverall pressure ratio: {solver.P[6]/solver.P[0]:.3f}")
print(f"Overall temperature ratio: {solver.T[6]/solver.T[0]:.3f}")