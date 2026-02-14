import numpy as np
import math
from pathlib import Path
from typing import Union, Sequence
from axialfans.newton_raphson import NewtonRaphsonSolver

class MultistageFanSolver:
    """
    Analytical solver for arbitrary turbomachinery cascades with variable area.
    
    Solves sequences of rotating/stationary axial stages (rotors, stators, 
    counter-rotating) using Newton-Raphson on coupled thermodynamic equations.
    
    Indexing Convention
    -------------------
    State arrays indexed 0 → N:
        - Index 0: inlet conditions
        - Index n: exit of stage n
    
    Stage parameters (Psi, Xi) use UNSIGNED values in local coordinates:
        Psi[n] = sigma[n] * omega[n]  (always positive)
        Xi[n] = -v_ax,inlet * tan(beta[n])
    
    Direction array handles rotation sense:
        direction[n] = +1 (CCW) or -1 (CW/mirror frame)
        Coordinate transforms occur at interfaces when direction changes
    
    Parameters
    ----------
    N : int
        Number of stages
    direction : array_like, shape (N,)
        Rotation direction: +1 (CCW) or -1 (CW). Example: [1, -1] for rotor-stator
    sigma : float or array_like, shape (N,)
        Slip factor, typically 0.85-0.95
    omega : float or array_like, shape (N,)
        Angular velocity [rad/s], UNSIGNED. Use 0 for stators
    beta : float or array_like, shape (N,)
        Blade angle [degrees], UNSIGNED, typically 20-70°
    rp, rm : float or array_like, shape (N,)
        Outer and inner radius [m]. Must have rp > rm
    eta : float or array_like, shape (N,)
        Isentropic efficiency, typically 0.75-0.92
    R : float
        Gas constant [J/(kg·K)]. Air: 287.053
    cp : float
        Specific heat [J/(kg·K)]. Air: 1005
    gamma : float, default=1.4
        Ratio of specific heats
    
    Attributes
    ----------
    T, P, vax, rho : ndarray, shape (N+1,)
        Temperature [K], pressure [Pa], axial velocity [m/s], density [kg/m³]
        at each station. Index n = exit of stage n.
    Psi, Xi : ndarray, shape (N+1,)
        Swirl parameters in local frame (always positive output convention)
    
    Examples
    --------
    Rotor-stator pair:
    >>> solver = MultistageFanSolver(
    ...     N=2, direction=[1, -1], omega=[1800, 0], beta=[65, 20],
    ...     sigma=0.9, rp=0.267, rm=0.184, eta=0.877, R=287, cp=1005)
    >>> solver.solve(T0=288, vax0=18, P0=101325, rho0=1.225)
    >>> print(f"PR: {solver.P[-1]/solver.P[0]:.2f}")
    
    Notes
    -----
    - Inviscid actuator disk model, valid for Re > 10^6
    - ~333 solves/sec enables real-time Monte Carlo (10k samples in 30s)
    - Area changes handled via mass flux conservation
    - See Wang (2026, submitted) for validation and regime boundaries
    """

    MAX_ITER = 60000
    TOL = 1e-8
    ALPHA = 0.5

    def __init__(self,
                 N: int,
                 direction: Sequence[int],
                 sigma: Union[float, Sequence[float]],
                 omega: Union[float, Sequence[float]],
                 beta: Union[float, Sequence[float]],  # degrees
                 rp: Union[float, Sequence[float]],
                 rm: Union[float, Sequence[float]],
                 eta: Union[float, Sequence[float]],
                 R: float,
                 cp: float,
                 gamma: float = 1.4):
        """
        Initialize the multistage fan solver.

        Args:
            N: Number of stages
            direction: Array of length N; 1 = CCW, -1 = CW
            sigma: Slip factor (scalar or array length N)
            omega: Blade angular speed (scalar or array length N, unsigned)
            beta: Blade angle in degrees (scalar or array length N)
            rp: Outer radius (scalar or array length N)
            rm: Inner radius (scalar or array length N)
            eta: Stage efficiency factor (scalar or array length N)
            R: Gas constant
            cp: Specific heat at constant pressure
            gamma: Ratio of specific heats
        """
        self.N = N

        # --- Helper: Expand scalar or array to 1-indexed array ---
        def expand_param(param, name):
            arr = np.zeros(N + 1)
            if np.isscalar(param):
                arr[1:] = param
            else:
                param = np.asarray(param, dtype=float)
                if len(param) != N:
                    raise ValueError(f"{name} must be scalar or length N")
                arr[1:] = param
            return arr

        # --- Direction validation ---
        direction = np.asarray(direction, dtype=int)
        if len(direction) != N:
            raise ValueError("direction array length must equal N")
        if not np.all(np.isin(direction, [-1, 1])):
            raise ValueError("direction entries must be 1 (CCW) or -1 (CW)")
        self.direction = np.zeros(N + 1, dtype=int)
        self.direction[1:] = direction

        # --- Stage parameters ---
        self.sigma = expand_param(sigma, "sigma")
        self.omega = expand_param(omega, "omega")  # UNSIGNED
        self.beta = np.deg2rad(expand_param(beta, "beta"))  # convert to radians
        self.rp = expand_param(rp, "rp")
        self.rm = expand_param(rm, "rm")
        if not np.all(self.rp[1:] > self.rm[1:]):
            raise ValueError("Each rp[i] must be greater than rm[i]")
        self.eta = expand_param(eta, "eta")

        # --- Precompute radial constants ---
        self.Y2 = self.rp**2 - self.rm**2
        self.Y3 = self.rp**3 - self.rm**3
        self.Y4 = self.rp**4 - self.rm**4
        self.A = self.Y2 * math.pi

        # --- Gas properties ---
        self.R = float(R)
        self.cp = float(cp)
        self.gamma = float(gamma)

        # --- State variables (0 → N) ---
        self.T = np.zeros(N + 1)
        self.vax = np.zeros(N + 1)
        self.P = np.zeros(N + 1)
        self.rho = np.zeros(N + 1)
        self.Psi = np.zeros(N + 1)
        self.Xi = np.zeros(N + 1)

        # Internal logger
        self._verbose = False

    def _log(self, msg: str):
        """Internal logger controlled by verbose flag."""
        if self._verbose:
            print(msg)

    def solve(self, T0: float, vax0: float, P0: float, rho0: float,
              Psi0: float = 0, Xi0: float = 0, verbose: bool = False):
        """
        Solve the multistage fan problem.

        Args:
            T0: Inlet temperature
            vax0: Inlet axial velocity
            P0: Inlet pressure
            rho0: Inlet density
            Psi0: Optional initial Psi
            Xi0: Optional initial Xi
            verbose: Print progress
        """
        self._verbose = verbose

        # --- Set inlet (stage 0) ---
        self.T[0], self.vax[0], self.P[0], self.rho[0] = T0, vax0, P0, rho0
        self.Psi[0], self.Xi[0] = Psi0, Xi0

        self._log("="*60)
        self._log(f"Starting solver for N={self.N} stages")
        self._log(f"INITIAL VALUES: T={self.T[0]:.2f} K, P={self.P[0]/1000:.2f} kPa, "
                      f"vax={self.vax[0]:.2f} m/s, rho={self.rho[0]:.4f} kg/m^3")

        for n in range(1, self.N + 1):
            
            if n > 1:
                old_vax_nm1 = self.vax[n-1]
                self.vax[n-1] = old_vax_nm1 * self.A[n-1] / self.A[n] # account for area change
                if not np.isclose(self.A[n-1], self.A[n]):
                    self._log(
                        f"Area change: {self.A[n-1]:.3f} → {self.A[n]:.3f}, "
                        f"input axial velocity for fan {n} = {self.vax[n-1]:.3f}"
                    )

            # Compute stage parameters
            self.Psi[n] = self.sigma[n] * self.omega[n]
            self.Xi[n] = -self.sigma[n] * self.vax[n-1] * math.tan(self.beta[n])

            # Determine Direction Change
            if n > 1: 
                if self.direction[n] != self.direction[n-1]:
                    self._log(f"Direction Change from {self.direction[n-1]} to {self.direction[n]} at stage {n}!") 
                    self.Psi[n-1] *= -1
                    self.Xi[n-1] *= -1
            
            # Solve stage n using Newton-Raphson
            sol = self._solve_stage(n)
            self.T[n], self.vax[n], self.P[n], self.rho[n] = sol

            if n > 1: 
                if self.direction[n] != self.direction[n-1]:
                    # Change Psi[n-1] and Xi[n-1] back to original values.
                    self.Psi[n-1] *= -1
                    self.Xi[n-1] *= -1
            
            if n > 1:
                self.vax[n-1] = old_vax_nm1 # Undo vax change.

            # Log stage
            self._log(f"Stage {n}: T={self.T[n]:.2f} K, P={self.P[n]/1000:.2f} kPa, "
                      f"vax={self.vax[n]:.2f} m/s, rho={self.rho[n]:.4f} kg/m^3")

        self._log("Solver completed.")
        self._verbose = False

    def _solve_stage(self, n: int) -> np.ndarray:
        """
        Solve a single stage using Newton-Raphson.

        Returns:
            np.ndarray: [T_n, vax_n, P_n, rho_n]
        """
        # Precompute E_CONST
        E_CONST = (- self.omega[n] * self.Psi[n-1] * self.Y4[n] / 4
                   - self.omega[n] * self.Xi[n-1] * self.Y3[n] / 3
                   + self.sigma[n] * self.omega[n]**2 * self.Y4[n] / 4
                   - self.sigma[n] * self.omega[n] * math.tan(self.beta[n]) * self.vax[n-1] * self.Y3[n] / 3
                   + 0.25 * self.vax[n-1]**2 * self.Y2[n]
                   + 0.125 * self.Psi[n-1]**2 * self.Y4[n]
                   + (1/3) * self.Psi[n-1] * self.Xi[n-1] * self.Y3[n]
                   + 0.25 * self.Xi[n-1]**2 * self.Y2[n])

        def residual(x):
            T, vax, P, rho = x
            R1 = (E_CONST
                  - 0.5 * self.cp * (T - self.T[n-1]) * self.Y2[n]
                  - 0.25 * vax**2 * self.Y2[n]
                  - 0.125 * self.Psi[n]**2 * self.Y4[n]
                  - (1/3) * self.Psi[n] * self.Xi[n] * self.Y3[n]
                  - 0.25 * self.Xi[n]**2 * self.Y2[n])
            R2 = P - self.P[n-1] * (1 + self.eta[n] * (T - self.T[n-1]) / self.T[n-1])**(self.gamma / (self.gamma - 1))
            R3 = rho - P / (self.R * T)
            R4 = rho * vax - self.rho[n-1] * self.vax[n-1]
            return np.array([R1, R2, R3, R4])

        def jacobian(x):
            T, vax, P, rho = x
            J = np.zeros((4,4))
            J[0,0] = -0.5 * self.cp * self.Y2[n]
            J[0,1] = -0.5 * vax * self.Y2[n]
            J[0,2] = 0
            J[0,3] = 0

            J[1,0] = -self.P[n-1] * self.eta[n] * self.gamma / ((self.gamma - 1) * self.T[n-1]) \
                     * (1 + self.eta[n] * (T - self.T[n-1]) / self.T[n-1])**(1/(self.gamma - 1))
            J[1,1] = 0
            J[1,2] = 1
            J[1,3] = 0

            J[2,0] = P / (self.R * T**2)
            J[2,1] = 0
            J[2,2] = -1 / (self.R * T)
            J[2,3] = 1

            J[3,0] = 0
            J[3,1] = rho
            J[3,2] = 0
            J[3,3] = vax
            return J

        solver = NewtonRaphsonSolver(residual, jacobian,
                                     max_iter=self.MAX_ITER,
                                     tol=self.TOL,
                                     alpha=self.ALPHA,
                                     verbose=False)
        x0 = np.array([self.T[n-1], self.vax[n-1], self.P[n-1], self.rho[n-1]])
        sol = solver.solve(x0)
        return sol

    def reset(self):
        """Reset solver state to initial conditions."""
        self.T[:] = 0
        self.vax[:] = 0
        self.P[:] = 0
        self.rho[:] = 0
        self.Psi[:] = 0
        self.Xi[:] = 0

    def __repr__(self):
        return f"MultistageFanSolver(N={self.N}, stages={self.direction[1:].tolist()})"