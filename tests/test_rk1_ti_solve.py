import numpy as np
import sdepack


def drift_func(x):
    return 1.0


def diffusion_func(x):
    return 1.0


def _run_case():
    print("   Drift function:      f(x) = 1.0")
    print("   Diffusion function:  g(x) = 1.0")
    print()

    print("-" * 40)

    x0 = 0.0
    t0 = 0.0
    tn = 1.0
    n = 10
    q = 1.0
    seed = 123456789

    x = np.zeros(n + 1, dtype=np.float64)

    print(f"   Initial value x0 = {x0}")
    print(f"   Time: {t0} -> {tn}, steps = {n}")
    print(f"   Noise strength Q = {q}")
    print()

    sdepack.rk1_ti_solve(drift_func, diffusion_func, x, t0, tn, x0, n, q, seed)

    print()

    return x, x0, n


def test_rk1_ti_solve_with_seed_output():
    x, x0, n = _run_case()

    assert x.shape == (n + 1,)
    assert np.isclose(x[0], x0)
    assert np.isfinite(x).all()
    expected = np.array(
        [
            0.0,
            0.630959,
            0.581457,
            0.502453,
            0.529365,
            1.01293,
            1.28212,
            1.78354,
            2.21543,
            1.78857,
            1.29873,
        ],
        dtype=np.float64,
    )
    assert np.allclose(x, expected)
