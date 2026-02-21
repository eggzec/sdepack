import numpy as np
import sdepack


def drift_func(x):
    """Mean-reverting drift: f(x) = -0.5*x"""
    return -0.5 * x


def diffusion_func(x):
    """Constant diffusion: g(x) = 1.0"""
    return 1.0


print("   Drift function:      f(x) = -0.5*x")
print("   Diffusion function:  g(x) = 1.0")
print()

print("4. Solving SDE with rk1_ti_solve:")
print("-" * 40)

x0 = 1.0  # Initial value
t0 = 0.0  # Initial time
tn = 2.0  # Final time
n = 20  # Number of steps
q = 0.05  # Noise strength

x = np.zeros(n + 1, dtype=np.float64)

print(f"   Initial value x0 = {x0}")
print(f"   Time: {t0} -> {tn}, steps = {n}")
print(f"   Noise strength Q = {q}")
print()

sdepack.rk1_ti_solve(drift_func, diffusion_func, x, t0, tn, x0, n, q, 99999)

print("   Solver completed!")
print()

t_values = np.linspace(t0, tn, n + 1)
print("   Result summary:")
print(f"   x[0] = {x[0]:.8f} (initial value at t=0)")
print(f"   x[n] = {x[-1]:.8f} (final value at t={tn})")
print()
