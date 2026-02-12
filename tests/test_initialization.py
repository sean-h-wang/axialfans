import numpy as np
from axialfans import MultistageFanSolver

def test_initialization():
    solver = MultistageFanSolver(
        N=2,
        direction=[1, -1],
        sigma=0.9,
        omega=[3000,3000],
        beta=[60, 60],
        rp=[0.2,0.25],
        rm=[0.05,0.1],
        eta=0.7,
        R=287,
        cp=1005
    )
    assert solver.N == 2
    assert solver.rp[1] == 0.2
