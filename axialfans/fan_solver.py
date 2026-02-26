"""
solver.py — Quasi-1D Pointwise Multistage Axial Turbomachinery Solver
======================================================================
Architecture
------------
- Each blade row is solved pointwise at M independent radial stations.
- Newton-Raphson (analytical Jacobian + line search) at each station.
- Area averaging applied AFTER the nonlinear solve (correct Jensen order).
- State averaging is area-weighted — valid even when vax = 0.
- Sequential stage coupling: stage n takes stage n-1 output as fixed input.
  No outer loop needed — the problem is causal (inlet-driven).

Inputs
------
sigma, omega, rp, rm, eta : scalar or 1-D length-N array (one per stage)
beta                       : scalar | 1-D length-N | 1-D length-M (N=1) | 2-D (N x M)

Slip factor sigma
-----------------
Pass a scalar or per-stage array. The solver does NOT compute sigma
internally — supply it from your chosen model (Qiu, Carter, Wiesner, etc.)
before calling solve(). For Qiu-consistent sigma use Picard iteration
externally (see test_solver.py).
"""

import numpy as np
import math
from typing import Union, Sequence


# ── Physical constants (defaults) ─────────────────────────────────────────────
GAMMA_DEFAULT = 1.4
R_DEFAULT     = 287.05   # J/(kg·K)
CP_DEFAULT    = GAMMA_DEFAULT * R_DEFAULT / (GAMMA_DEFAULT - 1)  # ~1004.7


# ══════════════════════════════════════════════════════════════════════════════
#  STATE
# ══════════════════════════════════════════════════════════════════════════════

class State:
    """
    Holds the full radial distribution of flow properties at one axial station.

    Parameters
    ----------
    n   : stage index (0 = inlet)
    M   : number of radial stations
    rm  : hub radius [m]
    rp  : tip radius [m]
    T, P, vax, rho, vtheta : arrays of length M (or scalars → broadcast)
    """

    def __init__(self, n, M, rm, rp, T, P, vax, rho, vtheta):
        self.n  = int(n)
        self.M  = int(M)
        self.rm = float(rm)
        self.rp = float(rp)
        self.r  = np.linspace(rm, rp, M)

        def _to_arr(x):
            x = np.asarray(x, dtype=float).ravel()
            if x.size == 1:
                return np.full(M, x[0])
            if x.size == M:
                return x.copy()
            raise ValueError(f"Expected scalar or length-{M} array, got {x.size}")

        self.T      = _to_arr(T)
        self.P      = _to_arr(P)
        self.vax    = _to_arr(vax)
        self.rho    = _to_arr(rho)
        self.vtheta = _to_arr(vtheta)

        # Precompute annular area weights  (trapezoidal, shape M-1)
        r      = self.r
        dr     = np.diff(r)
        r_mid  = 0.5 * (r[1:] + r[:-1])
        self.dA       = 2.0 * np.pi * r_mid * dr   # annular ring areas
        self.A_total  = float(np.sum(self.dA))

        self._compute_averages()

    # ── Area-weighted average (works even when vax = 0) ──────────────────────
    def _area_avg(self, var):
        """Trapezoid-rule area average of a radial array."""
        mid = 0.5 * (var[1:] + var[:-1])
        return float(np.sum(mid * self.dA) / self.A_total)

    def _compute_averages(self):
        self.T_avg      = self._area_avg(self.T)
        self.P_avg      = self._area_avg(self.P)
        self.vax_avg    = self._area_avg(self.vax)
        self.rho_avg    = self._area_avg(self.rho)
        self.vtheta_avg = self._area_avg(self.vtheta)

        # Mass flow — still useful as a diagnostic even if vax can be zero
        mid_flux      = 0.5 * ((self.rho * self.vax)[1:] + (self.rho * self.vax)[:-1])
        self.mdot     = float(np.sum(mid_flux * self.dA))

    def summary(self):
        return {
            "n":           self.n,
            "P_avg [Pa]":  self.P_avg,
            "T_avg [K]":   self.T_avg,
            "vax_avg":     self.vax_avg,
            "rho_avg":     self.rho_avg,
            "vtheta_avg":  self.vtheta_avg,
            "mdot [kg/s]": self.mdot,
        }

    def __repr__(self):
        s = self.summary()
        return (f"State(n={s['n']}, P={s['P_avg [Pa]']:.1f} Pa, "
                f"T={s['T_avg [K]']:.2f} K, vax={s['vax_avg']:.2f} m/s)")


# ══════════════════════════════════════════════════════════════════════════════
#  POINTWISE NEWTON-RAPHSON STREAMTUBE SOLVER
# ══════════════════════════════════════════════════════════════════════════════

def _solve_radial_station(r, work_local, vt_exit, vt_in,
                      T_prev, P_prev, vax_prev, rho_prev,
                      eta, Rgas, cp, gamma, max_iter=1000):
    """
    Solve 4 equations at a single radial station r using NewtonRaphsonSolver.

    Unknowns:  x = [T_n, vax_n, P_n, rho_n]

    Equations (UNINTEGRATED / pointwise):
      R1 — energy:   work + 0.5*vax_prev^2 + 0.5*vt_in^2
                     = cp*(T_n - T_prev) + 0.5*vax_n^2 + 0.5*vt_exit^2
      R2 — isentropic + eta:
                     P_n = P_prev * [1 + eta*(T_n - T_prev)/T_prev]^(g/(g-1))
      R3 — ideal gas:  rho_n = P_n / (R * T_n)
      R4 — continuity: rho_n * vax_n = rho_prev * vax_prev
    """
    from axialfans.newton_raphson import NewtonRaphsonSolver

    exp = gamma / (gamma - 1.0)
    E   = work_local + 0.5 * vax_prev**2 + 0.5 * vt_in**2

    def residual(x):
        T_n, vax_n, P_n, rho_n = x
        ratio = max(1.0 + eta * (T_n - T_prev) / T_prev, 1e-10)
        return np.array([
            E - cp * (T_n - T_prev) - 0.5 * vax_n**2 - 0.5 * vt_exit**2,        # R1 energy
            P_n - P_prev * ratio**exp,                                          # R2 isentropic
            rho_n - P_n / (Rgas * T_n),                                         # R3 ideal gas
            rho_n * vax_n - rho_prev * vax_prev,                                # R4 continuity
        ])

    def jacobian(x):
        T_n, vax_n, P_n, rho_n = x
        ratio = max(1.0 + eta * (T_n - T_prev) / T_prev, 1e-10)
        return np.array([
            [-cp,                                           -vax_n,              0.0,              0.0   ],  # dR1
            [-P_prev * exp * (eta/T_prev) * ratio**(exp-1), 0.0,                1.0,              0.0   ],  # dR2
            [ P_n / (Rgas * T_n**2),                        0.0,   -1.0/(Rgas*T_n),              1.0   ],  # dR3
            [ 0.0,                                          rho_n,               0.0,            vax_n  ],  # dR4
        ])

    # Initial guess
    T0   = T_prev * 1.05
    P0   = P_prev * 1.10
    rho0 = P0 / (Rgas * T0)
    vax0 = (rho_prev * vax_prev / rho0) if abs(rho_prev * vax_prev) > 1e-12 else vax_prev

    solver = NewtonRaphsonSolver(residual, jacobian, max_iter=max_iter, tol=1e-8)
    return solver.solve(np.array([T0, vax0, P0, rho0]))


# ══════════════════════════════════════════════════════════════════════════════
#  MULTISTAGE SOLVER
# ══════════════════════════════════════════════════════════════════════════════

class MultistageFanSolver:
    """
    Pointwise quasi-1D multistage axial fan/compressor solver.

    Stage indexing
    --------------
    Internal arrays are length N+1.  Index 0 = inlet (no blade row).
    Indices 1..N correspond to user-supplied blade rows.

    beta handling
    -------------
    If beta is supplied as a scalar or 1-D length-N array → uniform across span.
    If beta is 2-D (N × M) → full spanwise twist specified by user.
    If beta_twist (d_beta_dm in 1/m) is also supplied, it is applied on top of
    the tip beta to generate beta(r) = beta_tip + d_beta_dm * (r_tip - r).
    """

    def __init__(self,
                 N: int,
                 M: int,
                 direction: Sequence[int],
                 sigma:  Union[float, Sequence[float]],
                 omega:  Union[float, Sequence[float]],
                 beta:   Union[float, Sequence[float], Sequence[Sequence[float]]],
                 rp:     Union[float, Sequence[float]],
                 rm:     Union[float, Sequence[float]],
                 eta:    Union[float, Sequence[float]],
                 R:      float = R_DEFAULT,
                 cp:     float = CP_DEFAULT,
                 gamma:  float = GAMMA_DEFAULT,
):
        """
        Parameters
        ----------
        N         : number of blade rows
        M         : radial stations per stage
        direction : length-N array of +1 / -1
        sigma     : slip factor(s)
        omega     : angular velocity [rad/s]
        beta      : blade exit angle(s) [deg], measured from axial
        rp, rm    : tip / hub radii [m]
        eta       : isentropic efficiency
        R, cp     : gas constants
        gamma     : ratio of specific heats
        beta : scalar / 1-D length-N / 1-D length-M (N=1) / 2-D (N x M)
                Blade angle(s) in degrees, positive from axial. Pass
                np.linspace(beta_hub, beta_tip, M) for a twisted blade.
        """
        self.N     = int(N)
        self.M     = int(M)
        self.R     = float(R)
        self.cp    = float(cp)
        self.gamma = float(gamma)

        # Expand all scalar/1D inputs to length N+1
        self.sigma  = self._expand1d(sigma, 'sigma')
        self.omega  = self._expand1d(omega,  "omega")
        self.rp     = self._expand1d(rp,     "rp")
        self.rm     = self._expand1d(rm,     "rm")
        self.eta    = self._expand1d(eta,    "eta")

        # Direction
        direction = np.array(direction, dtype=int)
        if direction.size != N:
            raise ValueError(f"direction must have length N={N}")
        if not np.all(np.isin(direction, [-1, 1])):
            raise ValueError("direction entries must be +1 or -1")
        # Pad with 0 at index 0 (inlet — not used)
        self.direction = np.concatenate([[0], direction])

        # Beta — 2D array (N+1) × M  [radians, positive]
        self.beta = self._expand2d(beta, 'beta', deg2rad=True)

        # Geometry validation
        bad = np.where(self.rp[1:] <= self.rm[1:])[0] + 1
        if len(bad):
            raise ValueError(f"rp must be > rm at stage(s): {bad.tolist()}")

        self.states: list[State] = []
        self._verbose = False

    # ── Input expansion helpers ───────────────────────────────────────────────

    def _expand1d(self, x, name):
        """
        Expand into length-(N+1) array, index 0 unused.
        Accepts: scalar  |  1-D length-N
        """
        x = np.asarray(x, dtype=float).ravel()
        out = np.zeros(self.N + 1, dtype=float)
        if x.size == 1:
            out[1:] = x[0]
        elif x.size == self.N:
            out[1:] = x
        else:
            raise ValueError(
                f"{name} must be scalar or length-{self.N} array, got {x.size}"
            )
        return out

    def _expand2d(self, x, name, deg2rad=False):
        """
        Expand into (N+1) x M array, row 0 unused.
        Accepts:
          scalar            -> broadcast to all stages and radii
          1-D length N      -> one value per stage, broadcast across span
          1-D length M      -> full spanwise profile (only valid when N=1)
          2-D (N x M)       -> full per-stage spanwise profile
        If deg2rad=True, converts degrees → radians after expansion.
        """
        N, M = self.N, self.M
        x   = np.asarray(x, dtype=float)
        out = np.zeros((N + 1, M), dtype=float)

        if x.ndim == 2:
            if x.shape != (N, M):
                raise ValueError(
                    f"{name} 2D must be shape ({N},{M}), got {x.shape}"
                )
            out[1:, :] = x

        else:
            x = x.ravel()
            if x.size == 1:
                out[1:, :] = x[0]
            elif x.size == N:
                out[1:, :] = x[:, None]
            elif x.size == M and N == 1:
                out[1, :] = x
            else:
                raise ValueError(
                    f"{name} must be scalar, length-{N} (per stage), "
                    f"length-{M} (spanwise, N=1 only), or shape ({N},{M}). "
                    f"Got {x.shape}."
                )

        if deg2rad:
            out = np.deg2rad(out)
        return out

    # ── Logging ───────────────────────────────────────────────────────────────

    def _log(self, msg):
        if self._verbose:
            print(msg)

    # ── Main solve ────────────────────────────────────────────────────────────

    def solve(self, inlet_state: State, verbose: bool = False) -> list[State]:
        """
        Solve all N stages sequentially.

        Parameters
        ----------
        inlet_state : State object with n=0 and matching M, rm, rp

        Returns
        -------
        List of State objects [inlet, stage1, stage2, ...]
        """
        # Validate inlet
        if inlet_state.n != 0:
            raise ValueError("inlet_state.n must be 0")
        if inlet_state.M != self.M:
            raise ValueError("inlet_state.M must match solver M")

        self._verbose = verbose
        self.states   = [inlet_state]

        self._log("=" * 60)
        self._log(f"Pointwise solver: N={self.N} stages, M={self.M} stations")

        for n in range(1, self.N + 1):
            prev  = self.states[n - 1]
            r_n   = np.linspace(self.rm[n], self.rp[n], self.M)

            # Output arrays
            T_n   = np.zeros(self.M)
            vax_n = np.zeros(self.M)
            P_n   = np.zeros(self.M)
            rho_n = np.zeros(self.M)
            vt_n  = np.zeros(self.M)

            # Same radial grid across all stages — no interpolation needed
            T_prev   = prev.T.copy()
            P_prev   = prev.P.copy()
            vax_prev = prev.vax.copy()
            rho_prev = prev.rho.copy()
            vt_prev  = prev.vtheta.copy()

            # Variable area: temporarily scale vax using incompressible continuity
            # vax_scaled = vax_prev * (A_prev / A_curr)
            # Restored automatically — vax_prev is a local copy
            A_prev = prev.A_total
            A_curr = float(np.pi * (self.rp[n]**2 - self.rm[n]**2))
            if not np.isclose(A_prev, A_curr):
                vax_prev = vax_prev * (A_prev / A_curr)

            # Direction flip for counter-rotating stages
            if n > 1 and self.direction[n] != self.direction[n - 1]:
                vt_prev = -vt_prev

            # ── Pointwise solve at each radial station ────────────────────
            sigma_n = self.sigma[n]   # slip factor for this stage (user-supplied)

            for m in range(self.M):
                r      = r_n[m]
                beta_m = self.beta[n, m]   # [radians]

                # Exit tangential velocity from velocity triangle
                vt_exit    = sigma_n * (self.omega[n] * r - vax_prev[m] * math.tan(beta_m))
                vt_n[m]    = vt_exit
                work_local = self.omega[n] * r * (vt_exit - vt_prev[m])

                sol = _solve_radial_station(
                    r          = r,
                    work_local = work_local,
                    vt_exit    = vt_exit,
                    vt_in      = vt_prev[m],
                    T_prev     = T_prev[m],
                    P_prev     = P_prev[m],
                    vax_prev   = vax_prev[m],
                    rho_prev   = rho_prev[m],
                    eta        = self.eta[n],
                    Rgas       = self.R,
                    cp         = self.cp,
                    gamma      = self.gamma,
                )
                T_n[m], vax_n[m], P_n[m], rho_n[m] = sol

            new_state = State(n, self.M, self.rm[n], self.rp[n],
                              T_n, P_n, vax_n, rho_n, vt_n)
            self.states.append(new_state)

            s = new_state.summary()
            self._log(f"  Stage {n}: P_avg={s['P_avg [Pa]']/1e3:.3f} kPa  "
                      f"T_avg={s['T_avg [K]']:.2f} K  "
                      f"vax_avg={s['vax_avg']:.2f} m/s")

        self._log("Done.")
        return self.states

    # ── Convenience: pressure ratio and temperature ratio ─────────────────────

    def performance(self):
        """Return area-averaged PR and TR relative to inlet."""
        inlet = self.states[0]
        exit  = self.states[-1]
        PR = exit.P_avg / inlet.P_avg
        TR = exit.T_avg / inlet.T_avg
        return {"PR": PR, "TR": TR,
                "eta_isen": (PR**((self.gamma-1)/self.gamma) - 1) / (TR - 1)}