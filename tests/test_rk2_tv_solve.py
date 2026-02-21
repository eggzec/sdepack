import numpy as np
import sdepack


def drift_func(x, t):
    return -0.5 * x + 0.1 * np.cos(2.0 * t)


def diffusion_func(x, t):
    return 1.0 + 0.25 * np.sin(t)


def _run_case():
    print("   Drift function:      f(x,t) = -0.5*x + 0.1*cos(2*t)")
    print("   Diffusion function:  g(x,t) = 1.0 + 0.25*sin(t)")
    print()

    print("-" * 40)

    x0 = 1.0
    t0 = 0.0
    tn = 2.0
    n = 10
    q = 0.05
    seed = 99999

    x = np.zeros(n + 1, dtype=np.float64)

    print(f"   Initial value x0 = {x0}")
    print(f"   Time: {t0} -> {tn}, steps = {n}")
    print(f"   Noise strength Q = {q}")
    print()

    sdepack.rk2_tv_solve(drift_func, diffusion_func, x, t0, tn, x0, n, q, seed)

    print()

    return x, x0, n


def test_rk2_tv_solve_with_seed_output():
    x, x0, n = _run_case()

    assert x.shape == (n + 1,)
    assert np.isclose(x[0], x0)
    assert np.isfinite(x).all()
    expected = np.array(
        [
            1.0,
            0.84872418,
            0.71193827,
            0.46195585,
            0.3306327,
            0.29421564,
            0.32318298,
            0.36999102,
            0.16215647,
            0.04937996,
            0.24627925,
        ],
        dtype=np.float64,
    )
    assert np.allclose(x, expected)
