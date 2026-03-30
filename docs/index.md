# sdepack

**Runge-Kutta numerical integration of scalar stochastic differential equations (SDEs) for Python.**

---

## Overview

`sdepack` is a Python library for solving scalar
[stochastic differential equations](https://en.wikipedia.org/wiki/Stochastic_differential_equation)
(SDEs) using stochastic Runge-Kutta methods. It allows you to **easily** and **efficiently**
integrate SDEs driven by a Wiener process — from molecular dynamics and population genetics
to option pricing and control systems.

A **stochastic differential equation** is a differential equation in which one or more terms
are stochastic processes, producing a solution that is itself a random process.

`sdepack` targets **scalar Itô SDEs** of the general form:

$$
dX(t) = a(X, t)\,dt + b(X, t)\,dW(t),
$$

where $W(t)$ is a [Wiener process](https://en.wikipedia.org/wiki/Wiener_process) (standard
Brownian motion). Within the solvers, the SDE is parameterized as:

$$
dX(t) = F(X, t)\,dt + Q\,G(X, t)\,dW(t),
$$

where $F$ is the user-supplied **drift** callback, $G$ is the user-supplied **diffusion**
callback, and $Q$ scales the noise intensity. For **time-invariant** routines, $F$ and $G$
do not depend explicitly on time.

## Requirements

- [NumPy](http://www.numpy.org/)

## Example Usage

```python
import numpy as np
import sdepack

# Ornstein-Uhlenbeck process: dX = -X dt + dW
x = np.zeros(1001, dtype=np.float64)
sdepack.rk4_ti_solve(
    lambda x: -x,    # drift F(X) = -X
    lambda x: 1.0,   # diffusion G(X) = 1
    x,
    0.0,             # t0
    10.0,            # tn
    1.0,             # x0
    1000,            # n steps
    1.0,             # Q (noise intensity)
    42,              # seed
)
print(x[:6])
```

## Main Features

1. Seven stochastic Runge-Kutta solvers (orders 1–4).
2. Both time-invariant and time-variant SDE forms.
3. Plain Python callables for drift and diffusion — lambdas work too.
4. Seeded, fully reproducible trajectories via a built-in PRNG.
5. NumPy array output — integrates directly with matplotlib, scipy, and the scientific Python stack.

**Available solvers:**

| Solver | Stages | Order | Time-dependence | Based on |
|---|---|---|---|---|
| `rk1_ti_solve` | 1 | 1 | invariant | Euler-Maruyama |
| `rk1_tv_solve` | 1 | 1 | variant | Euler-Maruyama |
| `rk2_ti_solve` | 2 | 2 | invariant | Kasdin (1995) |
| `rk2_tv_solve` | 2 | 2 | variant | Kasdin (1995) |
| `rk3_ti_solve` | 3 | 3 | invariant | Kasdin (1995) |
| `rk4_ti_solve` | 4 | 4 | invariant | Kasdin (1995) |
| `rk4_tv_solve` | 4 | 4 | variant | Kasdin (1995) |

## References

- N. J. Kasdin, 1995, *Runge-Kutta algorithm for the numerical integration of stochastic differential equations*, Journal of Guidance, Control, and Dynamics, Vol. 18, No. 1, pp. 114–120.
