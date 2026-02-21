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

    sdepack.rk3_ti_solve(drift_func, diffusion_func, x, t0, tn, x0, n, q, seed)

    print()

    return x, x0, n


def test_rk3_ti_solve_with_seed_output():
    x, x0, n = _run_case()

    assert x.shape == (n + 1,)
    assert np.isclose(x[0], x0)
    assert np.isfinite(x).all()
    expected = np.array(
        [
            0.0,
            0.32969837,
            0.64810164,
            0.9667549,
            0.3845353,
            0.37553716,
            0.48784343,
            1.16036777,
            1.49431369,
            1.59940185,
            2.1654792,
        ],
        dtype=np.float64,
    )
    assert np.allclose(x, expected)
