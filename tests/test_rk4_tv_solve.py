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

    sdepack.rk4_tv_solve(drift_func, diffusion_func, x, t0, tn, x0, n, q, seed)

    print()

    return x, x0, n


def test_rk4_tv_solve_with_seed_output():
    x, x0, n = _run_case()

    assert x.shape == (n + 1,)
    assert np.isclose(x[0], x0)
    assert np.isfinite(x).all()
    expected = np.array(
        [
            1.0,
            0.84965124,
            0.58111185,
            0.54310931,
            0.58651215,
            0.40478869,
            0.15210296,
            0.15575955,
            0.17411699,
            0.17575902,
            0.29481523,
        ],
        dtype=np.float64,
    )
    assert np.allclose(x, expected)
