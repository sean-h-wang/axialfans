"""Monte Carlo uncertainty quantification example."""
import numpy as np
from axialfans.fan_solver import MultistageFanSolver
import time


# Nominal configuration
N_samples = 1000

# Parameter distributions
beta_samples = np.random.normal(45, 5, N_samples)  # Mean 45°, std 5°
sigma_samples = np.random.uniform(0.85, 0.95, N_samples)

results = []
failed = 0

print(f"Running {N_samples} Monte Carlo samples...")
start_time = time.time()

for i in range(N_samples):
    try:
        solver = MultistageFanSolver(
            N=1,
            direction=[1],
            sigma=sigma_samples[i],
            omega=1800,
            beta=beta_samples[i],
            rp=0.25,
            rm=0.18,
            eta=0.85,
            R=287,
            cp=1005
        )
        solver.solve(T0=288, vax0=50, P0=101325, rho0=1.225, verbose=False)
        results.append(solver.P[1] / solver.P[0])
    except:
        failed += 1

elapsed = time.time() - start_time

# Statistics
results = np.array(results)
print(f"\n{'='*60}")
print(f"Completed in {elapsed:.2f} seconds")
print(f"Solve rate: {len(results)/elapsed:.0f} solves/second")
print(f"Success rate: {len(results)}/{N_samples} ({100*len(results)/N_samples:.1f}%)")
print(f"\nPRESSURE RATIO STATISTICS:")
print(f"Mean: {np.mean(results):.3f}")
print(f"Std Dev: {np.std(results):.3f}")
print(f"95% CI: [{np.percentile(results, 2.5):.3f}, {np.percentile(results, 97.5):.3f}]")
print(f"Min: {np.min(results):.3f}")
print(f"Max: {np.max(results):.3f}")