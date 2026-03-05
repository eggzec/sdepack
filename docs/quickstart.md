# Quickstart

This page walks through progressively richer examples of using `sdepack`.
Each example is self-contained and produces deterministic output thanks to the
seeded internal RNG.

---

## Concepts

Before diving in, here are the key ideas:

1. **Allocate the output array** — create a NumPy array of shape `(N+1,)` to
   hold the trajectory (including the initial condition at index 0).
2. **Define drift and diffusion** — write plain Python functions matching the
   solver's signature (`f(x)` for time-invariant, `f(x, t)` for time-variant).
3. **Call the solver** — the solver fills the array in place.

```python
import numpy as np
import sdepack

x = np.zeros(N + 1, dtype=np.float64)
sdepack.rkX_YY_solve(drift, diffusion, x, t0, tn, x0, N, Q, seed)
# x now contains the full trajectory
```

---

## Example 1: Pure Brownian motion (RK1, time-invariant)

The simplest possible SDE — pure additive noise with zero drift:

$$
dX(t) = 0\,dt + 1\cdot dW(t),\quad X(0)=0.
$$

This produces a standard Wiener process trajectory.

```python
import numpy as np
import sdepack

def drift(x):
    return 0.0

def diffusion(x):
    return 1.0

x = np.zeros(11, dtype=np.float64)
sdepack.rk1_ti_solve(drift, diffusion, x, 0.0, 1.0, 0.0, 10, 1.0, 123456789)
print(x)
```

## Example 2: Constant drift with noise (RK1, time-invariant)

Solve

$$
dX(t)=1\,dt + 1\cdot dW(t),\quad X(0)=0,
$$

with $N=10$, $T_0=0$, $T_N=1$, $Q=1$, seed $=123456789$.

This is a Brownian motion with constant upward drift — the solution is a
random walk biased to increase by 1 per unit time on average.

```python
import numpy as np
import sdepack

def drift(x):
    return 1.0

def diffusion(x):
    return 1.0

x = np.zeros(11, dtype=np.float64)
sdepack.rk1_ti_solve(drift, diffusion, x, 0.0, 1.0, 0.0, 10, 1.0, 123456789)
print(x)
```

Expected trajectory:

```text
[0.0, 0.630959, 0.581457, 0.502453, 0.529365,
 1.01293, 1.28212, 1.78354, 2.21543, 1.78857, 1.29873]
```

## Example 3: Ornstein-Uhlenbeck process (RK2, time-invariant)

The **Ornstein-Uhlenbeck (OU) process** is a fundamental mean-reverting SDE:

$$
dX(t) = -\theta X(t)\,dt + \sigma\, dW(t),\quad X(0) = X_0.
$$

It models systems that tend to drift back toward zero (or a long-run mean)
while being perturbed by noise. Used extensively in physics (velocity of a
Brownian particle under friction) and finance (Vasicek interest rate model).

Here we solve with $\theta = 1$ and $\sigma = 1$:

```python
import numpy as np
import sdepack

def drift(x):
    return -x  # mean-reverting drift

def diffusion(x):
    return 1.0  # constant diffusion

N = 100
x = np.zeros(N + 1, dtype=np.float64)
sdepack.rk2_ti_solve(drift, diffusion, x, 0.0, 5.0, 2.0, N, 1.0, 42)
print(f"Start: {x[0]:.4f}, End: {x[-1]:.4f}")
```

## Example 4: Geometric Brownian motion (RK4, time-invariant)

**Geometric Brownian motion (GBM)** is the foundation of the Black-Scholes
options pricing model. The SDE is:

$$
dX(t) = \mu X(t)\,dt + \sigma X(t)\,dW(t),\quad X(0) = X_0.
$$

Note that both drift and diffusion are proportional to $X$, making this a
**multiplicative noise** problem. The exact solution is:

$$
X(t) = X_0 \exp\!\left[\left(\mu - \tfrac{1}{2}\sigma^2\right)t + \sigma W(t)\right].
$$

```python
import numpy as np
import sdepack

mu = 0.05     # drift rate (5% annual return)
sigma = 0.2   # volatility (20%)

def drift(x):
    return mu * x

def diffusion(x):
    return sigma * x

N = 252  # one trading year of daily steps
x = np.zeros(N + 1, dtype=np.float64)
sdepack.rk4_ti_solve(drift, diffusion, x, 0.0, 1.0, 100.0, N, 1.0, 12345)
print(f"Initial price: ${x[0]:.2f}")
print(f"Final price:   ${x[-1]:.2f}")
```

## Example 5: Time-variant SDE (RK4)

Solve a system with explicitly time-dependent coefficients:

$$
dX(t)=\left(-0.5X + 0.1\cos(2t)\right)dt + 0.05\left(1+0.25\sin t\right)dW(t),\quad X(0)=1.
$$

This models a mean-reverting process with a periodic forcing term in the drift
and a seasonally modulated noise amplitude — useful for modeling climate
variables, cyclical economic indicators, etc.

```python
import numpy as np
import sdepack

def drift(x, t):
    return -0.5 * x + 0.1 * np.cos(2.0 * t)

def diffusion(x, t):
    return 1.0 + 0.25 * np.sin(t)

x = np.zeros(11, dtype=np.float64)
sdepack.rk4_tv_solve(drift, diffusion, x, 0.0, 2.0, 1.0, 10, 0.05, 99999)
print(x)
```

Expected trajectory:

```text
[1.0, 0.84965124, 0.58111185, 0.54310931, 0.58651215,
 0.40478869, 0.15210296, 0.15575955, 0.17411699, 0.17575902, 0.29481523]
```

## Example 6: Comparing solvers on the same problem

A powerful feature of `sdepack` is that all solvers use the same seeded RNG,
making it easy to compare methods on the same noise realization:

```python
import numpy as np
import sdepack

def drift(x):
    return -x

def diffusion(x):
    return 1.0

N = 50
seed = 777

solvers = [
    ("RK1", sdepack.rk1_ti_solve),
    ("RK2", sdepack.rk2_ti_solve),
    ("RK4", sdepack.rk4_ti_solve),
]

for name, solver in solvers:
    x = np.zeros(N + 1, dtype=np.float64)
    solver(drift, diffusion, x, 0.0, 5.0, 1.0, N, 1.0, seed)
    print(f"{name}: final value = {x[-1]:.6f}")
```

!!! note "Different RNG consumption"
    Each solver order consumes a different number of random variates per step
    (1 for RK1, 2 for RK2, 4 for RK4), so trajectories will diverge even with
    the same seed. The comparison is still useful for assessing smoothness and
    convergence behavior.

## Example 7: Monte Carlo ensemble

Run multiple trajectories to estimate the expected value and variance of the
solution:

```python
import numpy as np
import sdepack

def drift(x):
    return -x

def diffusion(x):
    return 0.5

N = 200
M = 1000  # number of sample paths

final_values = np.empty(M)
for i in range(M):
    x = np.zeros(N + 1, dtype=np.float64)
    sdepack.rk4_ti_solve(drift, diffusion, x, 0.0, 5.0, 1.0, N, 1.0, i + 1)
    final_values[i] = x[-1]

print(f"E[X(5)]  = {final_values.mean():.4f}")
print(f"Var[X(5)] = {final_values.var():.4f}")
```

For the OU process $dX = -X\,dt + 0.5\,dW$, the stationary distribution has
mean 0 and variance $\sigma^2/(2\theta) = 0.125$.

## Tips

!!! tip "Choosing step count $N$"
    Start with a moderate $N$ (100–1000) and double it to check convergence.
    If results change significantly, increase $N$. Higher-order solvers (RK4)
    converge faster, so they can use fewer steps for the same accuracy.

!!! tip "Plotting trajectories"
    Use `matplotlib` with the time grid for visualization:

    ```python
    import matplotlib.pyplot as plt

    t = np.linspace(0.0, tn, N + 1)
    plt.plot(t, x)
    plt.xlabel("Time")
    plt.ylabel("X(t)")
    plt.title("SDE trajectory")
    plt.show()
    ```

!!! tip "Seed selection"
    The seed must be a positive integer. Different seeds produce entirely
    different noise realizations. For Monte Carlo studies, use sequential seeds
    (`1, 2, 3, ...`) to generate independent sample paths.
