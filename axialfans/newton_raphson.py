import numpy as np
from typing import Callable, Tuple

class NewtonRaphsonSolver:
    """
    Newton-Raphson solver for systems of nonlinear equations.
    
    Parameters
    ----------
    residual_func : Callable[[np.ndarray], np.ndarray]
        Function returning the residual vector f(x).
    jacobian_func : Callable[[np.ndarray], np.ndarray]
        Function returning the Jacobian matrix J(x).
    max_iter : int, default=100
        Maximum number of iterations.
    tol : float, default=1e-12
        Convergence tolerance on residual norm.
    alpha : float, default=0.5
        Step damping factor (0 < alpha <= 1).
    verbose : bool, default=False
        If True, prints iteration info.
    """
    def __init__(
        self,
        residual_func: Callable[[np.ndarray], np.ndarray],
        jacobian_func: Callable[[np.ndarray], np.ndarray],
        max_iter: int = 100,
        tol: float = 1e-12,
        alpha: float = 0.5,
        verbose: bool = False
    ):
        self.residual_func = residual_func
        self.jacobian_func = jacobian_func
        self.max_iter = max_iter
        self.tol = tol
        self.alpha = alpha
        self.verbose = verbose

    def solve(self, x0: np.ndarray) -> np.ndarray:
        """
        Solve the nonlinear system f(x) = 0 using Newton-Raphson.
        
        Parameters
        ----------
        x0 : np.ndarray
            Initial guess for the solution.
        
        Returns
        -------
        np.ndarray
            Approximate solution x.
        """
        x = np.array(x0, dtype=float)
        
        for i in range(self.max_iter):
            f = self.residual_func(x)
            norm_f = np.linalg.norm(f)
            if norm_f < self.tol:
                if self.verbose:
                    print(f"[Newton-Raphson] Converged in {i} iterations (residual={norm_f:e})")
                return x
            
            J = self.jacobian_func(x)
            try:
                delta_x = np.linalg.solve(J, -f)
            except np.linalg.LinAlgError:
                raise RuntimeError(f"Jacobian is singular at iteration {i}.")
            
            x += self.alpha * delta_x
            
            if self.verbose:
                print(f"[Iteration {i}] residual norm = {norm_f:e}, |Δx| = {np.linalg.norm(delta_x):e}")
        
        raise RuntimeError(f"Newton-Raphson did not converge after {self.max_iter} iterations.")
