# sdepack

**Runge-Kutta numerical integration of scalar stochastic differential equations (SDEs) for Python.**

---

## What is sdepack?

`sdepack` is a Python library for solving scalar
[stochastic differential equations](https://en.wikipedia.org/wiki/Stochastic_differential_equation)
(SDEs) using stochastic Runge-Kutta methods. The numerical core is written in
Fortran and compiled via `f2py`, giving near-native performance while
providing a clean, NumPy-based Python API.

A **stochastic differential equation** is a differential equation in which one
or more terms are stochastic processes, producing a solution that is itself a
random process. SDEs are used to model phenomena across science and engineering
where systems are subject to random influences — from molecular dynamics and
population genetics to option pricing and signal processing.

`sdepack` targets **scalar Itô SDEs** of the general form:

$$
dX(t) = \underbrace{a(X, t)}_{\text{drift}}\,dt
      + \underbrace{b(X, t)}_{\text{diffusion}}\,dW(t),
$$

where $W(t)$ is a [Wiener process](https://en.wikipedia.org/wiki/Wiener_process)
(standard Brownian motion).

## Mathematical model

Within the solvers, the SDE is parameterized as:

**Time-invariant** routines integrate

$$
dX(t) = F(X)\,dt + Q\,G(X)\,dW(t),
$$

**Time-variant** routines integrate

$$
dX(t) = F(X,t)\,dt + Q\,G(X,t)\,dW(t),
$$

where $Q$ scales the noise intensity (spectral density of the driving white
noise), $F$ is the user-supplied **drift** callback and $G$ is the
user-supplied **diffusion** callback.

Integration is performed with a uniform time step:

$$
H = \frac{T_N - T_0}{N}.
$$

Stage noise for the Runge-Kutta stages is drawn from:

$$
W_i = Z_i\sqrt{\frac{Q_i Q}{H}},\quad Z_i \sim \mathcal{N}(0,1),
$$

where $Q_i$ are method-specific noise coefficients that ensure the correct
stochastic scaling at each Runge-Kutta stage.

## Available solvers

| Solver | Stages | Order | Time-dependence | Based on |
|---|---|---|---|---|
| `rk1_ti_solve` | 1 | 1 | invariant | Euler-Maruyama |
| `rk1_tv_solve` | 1 | 1 | variant | Euler-Maruyama |
| `rk2_ti_solve` | 2 | 2 | invariant | Kasdin (1995) |
| `rk2_tv_solve` | 2 | 2 | variant | Kasdin (1995) |
| `rk3_ti_solve` | 3 | 3 | invariant | Kasdin (1995) |
| `rk4_ti_solve` | 4 | 4 | invariant | Kasdin (1995) |
| `rk4_tv_solve` | 4 | 4 | variant | Kasdin (1995) |

!!! tip "Choosing a solver"
    For exploratory work, start with `rk1_ti_solve` (Euler-Maruyama).
    For production accuracy, prefer `rk4_ti_solve` or `rk4_tv_solve` — higher-order
    methods converge faster and yield smoother trajectories for a given step count.

## Quick example

```python
import numpy as np
import sdepack

# Ornstein-Uhlenbeck process: dX = -X dt + dW
x = np.zeros(1001, dtype=np.float64)
sdepack.rk4_ti_solve(
    lambda x: -x,    # drift F(X) = -X
    lambda x: 1.0,   # diffusion G(X) = 1
    x,
    0.0,              # t0
    10.0,             # tn
    1.0,              # x0
    1000,             # n steps
    1.0,              # Q (noise intensity)
    42,               # seed
)
print(x[:6])  # first few values of the trajectory
```
