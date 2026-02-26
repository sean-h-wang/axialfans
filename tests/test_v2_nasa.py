"""
test_solver.py — Validation tests for MultistageFanSolver
=================================================================
All slip factor models (Qiu, etc.) live here, not in the solver.
The solver only takes sigma as a numeric input.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


import numpy as np
import traceback
from axialfans.fan_solver import State, MultistageFanSolver

# ── Qiu et al. (2007) slip factor model ───────────────────────────────────────
# Fig 5C solid line, linear fit to pink endpoints:
#   sigma(phi2) = 1.0043 - 0.2143 * phi2
# where phi2 = vax_exit_avg / (omega * R2),  R2 = tip radius
# Valid for axial rotors (gamma2 = 0, radial term vanishes).
# Reference: Qiu, Mallikarachchi, Anderson, ASME GT2007-27064.

def qiu_sigma(vax_exit_avg, omega, R2):
    """
    Bulk slip factor from Qiu et al. (2007) Fig 5C.

    Parameters
    ----------
    vax_exit_avg : area-averaged axial velocity at stage exit [m/s]
    omega        : angular velocity [rad/s]
    R2           : tip radius [m]

    Returns
    -------
    sigma : float
    """
    phi2 = vax_exit_avg / (omega * R2)
    return float(np.clip(1.0043 - 0.2143 * phi2, 0.5, 1.0))


def picard_sigma(solver_factory, inlet, omega, R2,
                 sigma0=0.90, tol=1e-6, max_iter=50):
    """
    Picard iteration to find sigma consistent with Qiu model.

    Iterates:
        1. Run solver with current sigma
        2. Compute phi2 = vax_exit_avg / (omega * R2)
        3. Update sigma = qiu_sigma(vax_exit_avg, omega, R2)
        4. Repeat until |delta_sigma| < tol

    Parameters
    ----------
    solver_factory : callable(sigma) -> MultistageFanSolver
    inlet          : State at n=0
    omega          : angular velocity [rad/s]
    R2             : tip radius [m]
    sigma0         : initial sigma guess
    tol            : convergence tolerance on sigma

    Returns
    -------
    sigma_converged : float
    states          : list of State from final solve
    n_iters         : int
    """
    sigma = sigma0
    states = None
    for i in range(max_iter):
        solver = solver_factory(sigma)
        states = solver.solve(inlet, verbose=False)
        vax_exit_avg = states[1].vax_avg
        sigma_new = qiu_sigma(vax_exit_avg, omega, R2)
        delta = abs(sigma_new - sigma)
        sigma = sigma_new
        if delta < tol:
            return sigma, states, i + 1
    return sigma, states, max_iter


# ══════════════════════════════════════════════════════════════════════════════

M  = 200
rm = 0.177
rp = 0.252
omega = 1800.0
rho0  = 101325.0 / (287.05 * 288.15)

def make_inlet(M=M, rm=rm, rp=rp, vax=161.0):
    return State(n=0, M=M, rm=rm, rp=rp,
                 T=288.15, P=101325.0,
                 vax=vax, rho=rho0, vtheta=0.0)

# ──────────────────────────────────────────────────────────────────────────────
print("=" * 60)
print("TEST 1 — Rotor 37, uniform beta=38.9°, uniform sigma=0.90")
print("=" * 60)
try:
    inlet  = make_inlet()
    solver = MultistageFanSolver(
        N=1, M=M, direction=[1],
        sigma=0.90, omega=omega,
        beta=38.9,
        rp=rp, rm=rm, eta=0.877,
    )
    states = solver.solve(inlet, verbose=True)
    perf   = solver.performance()
    err    = abs(perf["PR"] - 2.106) / 2.106 * 100
    print(f"\n  PR={perf['PR']:.4f}  (NASA 2.106)  error={err:.2f}%")
    print(f"  TR={perf['TR']:.4f}  (NASA 1.270)")
except Exception:
    traceback.print_exc()

# ──────────────────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("TEST 2 — Rotor 37, blade twist β(r) = linspace(37.62°, 38.9°)")
print("         beta_hub = 38.9 - 17*(rp-rm),  uniform sigma=0.90")
print("=" * 60)
try:
    beta_hub = 38.9 - 17.0 * (rp - rm)
    beta_tip = 38.9
    beta_arr = np.linspace(beta_hub, beta_tip, M)
    print(f"  beta: hub={beta_arr[0]:.3f}°  tip={beta_arr[-1]:.3f}°")

    inlet  = make_inlet()
    solver = MultistageFanSolver(
        N=1, M=M, direction=[1],
        sigma=0.90, omega=omega,
        beta=beta_arr,
        rp=rp, rm=rm, eta=0.877,
    )
    states = solver.solve(inlet, verbose=True)
    perf   = solver.performance()
    err    = abs(perf["PR"] - 2.106) / 2.106 * 100
    print(f"\n  PR={perf['PR']:.4f}  (NASA 2.106)  error={err:.2f}%")
    print(f"  TR={perf['TR']:.4f}  (NASA 1.270)")
except Exception:
    traceback.print_exc()

# ──────────────────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("TEST 3 — vax=0 inlet (static air, area averaging must not break)")
print("=" * 60)
try:
    inlet  = make_inlet(M=10, vax=0.0)
    solver = MultistageFanSolver(
        N=1, M=10, direction=[1],
        sigma=0.90, omega=omega,
        beta=38.9, rp=rp, rm=rm, eta=0.877,
    )
    states = solver.solve(inlet, verbose=True)
    perf   = solver.performance()
    print(f"\n  PR={perf['PR']:.4f}  TR={perf['TR']:.4f}")
    print(f"  inlet P_avg={inlet.P_avg:.1f} Pa  (should be 101325)")
except Exception:
    traceback.print_exc()

# ──────────────────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("TEST 4 — 3-stage counter-rotating, scalar inputs")
print("=" * 60)
try:
    inlet  = State(n=0, M=5, rm=0.18, rp=0.25,
                   T=288.15, P=101325.0,
                   vax=100.0, rho=rho0, vtheta=0.0)
    solver = MultistageFanSolver(
        N=3, M=5, direction=[1, -1, 1],
        sigma=0.90, omega=1500.0,
        beta=40.0, rp=0.25, rm=0.18, eta=0.85,
    )
    states = solver.solve(inlet, verbose=True)
    perf   = solver.performance()
    print(f"\n  PR={perf['PR']:.4f}  TR={perf['TR']:.4f}")
except Exception:
    traceback.print_exc()

# ──────────────────────────────────────────────────────────────────────────────
print()
print("=" * 60)
print("TEST 5 — Rotor 37, blade twist + Qiu sigma via Picard iteration")
print("         sigma self-consistent with phi2 = vax_avg / (omega * R2)")
print("=" * 60)
try:
    beta_hub = 38.9 - 17.0 * (rp - rm)
    beta_tip = 38.9
    beta_arr = np.linspace(beta_hub, beta_tip, M)

    inlet = make_inlet()

    def make_solver(sigma):
        return MultistageFanSolver(
            N=1, M=M, direction=[1],
            sigma=sigma, omega=omega,
            beta=beta_arr,
            rp=rp, rm=rm, eta=0.877,
        )

    sigma_conv, states, n_iters = picard_sigma(
        make_solver, inlet, omega=omega, R2=rp, sigma0=0.90
    )

    vax_exit_avg = states[1].vax_avg
    phi2_final   = vax_exit_avg / (omega * rp)
    print(f"  Converged in {n_iters} Picard iterations")
    print(f"  sigma={sigma_conv:.4f}  phi2={phi2_final:.4f}")

    # Final solve with converged sigma for reporting
    solver = make_solver(sigma_conv)
    states = solver.solve(inlet, verbose=True)
    perf   = solver.performance()
    err    = abs(perf["PR"] - 2.106) / 2.106 * 100
    print(f"\n  PR={perf['PR']:.4f}  (NASA 2.106)  error={err:.2f}%")
    print(f"  TR={perf['TR']:.4f}  (NASA 1.270)")
    print(f"  eta_isen={perf['eta_isen']:.4f}")
except Exception:
    traceback.print_exc()