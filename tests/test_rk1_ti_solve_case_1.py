import numpy as np
import sdepack


def drift_func(x):
    return -0.5 * x


def diffusion_func(x):
    return 1.0


def _run_case():
    print("   Drift function:      f(x) = -0.5*x")
    print("   Diffusion function:  g(x) = 1.0")
    print()

    print("4. Solving SDE with rk1_ti_solve:")
    print("-" * 40)

    x0 = 1.0
    t0 = 0.0
    tn = 2.0
    n = 20
    q = 0.05

    x = np.zeros(n + 1, dtype=np.float64)

    print(f"   Initial value x0 = {x0}")
    print(f"   Time: {t0} -> {tn}, steps = {n}")
    print(f"   Noise strength Q = {q}")
    print()

    sdepack.rk1_ti_solve(drift_func, diffusion_func, x, t0, tn, x0, n, q, 99999)

    print("   Solver completed!")
    print()

    print("   Result summary:")
    print(f"   x[0] = {x[0]:.8f} (initial value at t=0)")
    print(f"   x[n] = {x[-1]:.8f} (final value at t={tn})")
    print()

    return x, x0, n


def test_rk1_ti_solve_output_and_shape(capsys):
    x, x0, n = _run_case()
    captured = capsys.readouterr().out

    assert "Drift function:      f(x) = -0.5*x" in captured
    assert "Diffusion function:  g(x) = 1.0" in captured
    assert "Solving SDE with rk1_ti_solve:" in captured
    assert "Solver completed!" in captured
    assert "Result summary:" in captured
    assert "x[0] =" in captured
    assert "x[n] =" in captured

    assert x.shape == (n + 1,)
    assert np.isclose(x[0], x0)
    assert np.isfinite(x).all()