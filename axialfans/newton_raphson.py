import numpy as np
from typing import Callable, Tuple

class NewtonRaphsonSolver:
    """
    Newton-Raphson solver for nonlinear systems f(x) = 0.
    
    Uses iteration: x_{k+1} = x_k - α * J(x_k)^{-1} * f(x_k)
    where J is the Jacobian and α is a damping factor.
    
    Parameters
    ----------
    residual_func : callable
        f(x) returning residual vector
    jacobian_func : callable
        J(x) returning Jacobian matrix where J[i,j] = ∂f_i/∂x_j
    max_iter : int, default=100
        Maximum iterations
    tol : float, default=1e-12
        Convergence tolerance on ||f(x)||
    alpha : float, default=0.5
        Damping factor (0 < alpha <= 1). Smaller = more stable, slower
    verbose : bool, default=False
        Print iteration info
    
    Methods
    -------
    solve(x0)
        Solve from initial guess x0
    
    Examples
    --------
    >>> def f(x): return np.array([x[0]**2 + x[1]**2 - 1, x[0] - x[1]])
    >>> def J(x): return np.array([[2*x[0], 2*x[1]], [1, -1]])
    >>> solver = NewtonRaphsonSolver(f, J, tol=1e-10)
    >>> x = solver.solve(np.array([1.0, 0.0]))
    
    Notes
    -----
    - Analytical Jacobians give 10-20× speedup vs finite differences
    - Typical convergence: 4-6 iterations for well-conditioned systems
    - Raises RuntimeError if singular Jacobian or no convergence
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
