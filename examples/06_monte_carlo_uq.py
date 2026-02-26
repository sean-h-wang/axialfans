# ═══════════════════════════════════════════════════════════
# Example 6 — Monte Carlo UQ
# ═══════════════════════════════════════════════════════════
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from axialfans.fan_solver import State, MultistageFanSolver

import numpy as np
import time
from scipy.stats import qmc

M    = 25
T0   = 288.15; P0 = 101325.0; rho0 = P0 / (287.05 * T0)
R    = 287.05; cp = 1004.7; gamma = 1.4
vax0 = 50.0
N_samples = 1000

sampler = qmc.LatinHypercube(d=2, seed=42)
sample  = qmc.scale(sampler.random(N_samples), [35, 0.85], [55, 0.95])
beta_samples  = sample[:, 0]
sigma_samples = sample[:, 1]

results = []; failed = 0
start = time.time()

for i in range(N_samples):
    try:
        s = MultistageFanSolver(
            N=1, M=M, direction=[1],
            sigma=sigma_samples[i], omega=1800, beta=beta_samples[i],
            rp=0.25, rm=0.18, eta=0.85,
            R=R, cp=cp, gamma=gamma
        )
        inlet = State(0, M, 0.18, 0.25, T0, P0, vax0, rho0, 0)
        s.solve(inlet)
        results.append(s.performance()['PR'])
    except Exception:
        failed += 1

elapsed = time.time() - start
results = np.array(results)
print("Example 6 — Monte Carlo UQ")
print("=" * 40)
print(f"Samples:  {len(results)} / {N_samples}  ({failed} failed)")
print(f"Time:     {elapsed:.2f} s  ({len(results)/elapsed:.0f} solves/s)")
print(f"PR mean:  {np.mean(results):.3f}")
print(f"PR std:   {np.std(results):.3f}")
print(f"95% CI:   [{np.percentile(results, 2.5):.3f}, {np.percentile(results, 97.5):.3f}]")