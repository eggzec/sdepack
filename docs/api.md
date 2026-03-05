# API Reference

All public routines are exported via the `sdepack` Python module. Only
subroutines declared in the `src/sdepack.pyf` interface file are accessible
from Python.

---

## Overview

| Function | Stages | Time-dep. | SDE form |
|---|---|---|---|
| `rk1_ti_solve` | 1 | No | $dX = F(X)dt + QG(X)dW$ |
| `rk1_tv_solve` | 1 | Yes | $dX = F(X,t)dt + QG(X,t)dW$ |
| `rk2_ti_solve` | 2 | No | $dX = F(X)dt + QG(X)dW$ |
| `rk2_tv_solve` | 2 | Yes | $dX = F(X,t)dt + QG(X,t)dW$ |
| `rk3_ti_solve` | 3 | No | $dX = F(X)dt + QG(X)dW$ |
| `rk4_ti_solve` | 4 | No | $dX = F(X)dt + QG(X)dW$ |
| `rk4_tv_solve` | 4 | Yes | $dX = F(X,t)dt + QG(X,t)dW$ |

---

## Common interface

All solvers share the same calling convention and write the solution trajectory
into a pre-allocated NumPy array `X(0:N)` in place. The uniform step size is:

$$
H = \frac{T_N - T_0}{N}.
$$

### Arguments

| Argument | Type | Direction | Description |
|---|---|---|---|
| `F` | callable | input | Drift function. Signature: `f(x) -> float` (TI) or `f(x, t) -> float` (TV) |
| `G` | callable | input | Diffusion function. Same signature convention as `F` |
| `X` | `ndarray(N+1)` | in/out | Solution array. Must be pre-allocated as `np.zeros(N+1, dtype=np.float64)` |
| `T0` | `float` | input | Start time of the integration interval |
| `TN` | `float` | input | End time of the integration interval |
| `X0` | `float` | input | Initial condition: $X(T_0) = X_0$ |
| `N` | `int` | input | Number of time steps (determines array size and step width $H$) |
| `Q` | `float` | input | Noise intensity — spectral density of the driving white noise |
| `SEED` | `int` | input | Seed for the internal Park-Miller PRNG (positive integer) |

!!! warning "Array size"
    The array `X` must have exactly `N + 1` elements. Passing an array of the
    wrong size will raise a Fortran runtime error.

!!! warning "Dtype"
    The array `X` must be `np.float64` (double precision). Other dtypes are
    not supported.

### Callback signatures

**Time-invariant** solvers (`_ti_`) expect:

```python
def drift(x: float) -> float:
    ...

def diffusion(x: float) -> float:
    ...
```

**Time-variant** solvers (`_tv_`) expect:

```python
def drift(x: float, t: float) -> float:
    ...

def diffusion(x: float, t: float) -> float:
    ...
```

!!! note "Lambda functions"
    Lambda functions work as callbacks:
    `sdepack.rk1_ti_solve(lambda x: -x, lambda x: 1.0, ...)`

---

## Time-invariant solvers

These solvers integrate SDEs where the drift and diffusion functions do **not**
depend explicitly on time:

$$
dX(t)=F(X)\,dt+Q\,G(X)\,dW(t).
$$

### `rk1_ti_solve(F, G, X, T0, TN, X0, N, Q, SEED)`

**Euler-Maruyama / RK1** — the simplest first-order stochastic integrator.

Update rule:

$$
X_{k+1}=X_k+H F(X_k)+\sqrt{HQ}\,G(X_k)\,Z_k,\qquad Z_k\sim\mathcal{N}(0,1).
$$

Butcher-like coefficients: $A_{21}=1$, $Q_1=1$.

**Convergence:** strong order ½, weak order 1.

---

### `rk2_ti_solve(F, G, X, T0, TN, X0, N, Q, SEED)`

**Two-stage stochastic Runge-Kutta.**

$$
\begin{aligned}
K_1 &= H F(X_1)+H G(X_1)W_1,\\
X_2 &= X_1 + A_{21}K_1,\\
K_2 &= H F(X_2)+H G(X_2)W_2,\\
X_{k+1} &= X_1 + A_{31}K_1 + A_{32}K_2,
\end{aligned}
$$

Coefficients:

$$
A_{21}=1,\quad A_{31}=0.5,\quad A_{32}=0.5,\quad Q_1=Q_2=2.
$$

---

### `rk3_ti_solve(F, G, X, T0, TN, X0, N, Q, SEED)`

**Three-stage stochastic Runge-Kutta.**

$$
X_{k+1} = X_1 + A_{41}K_1 + A_{42}K_2 + A_{43}K_3,
$$

Kasdin coefficients:

| | Value |
|---|---|
| $A_{21}$ | 1.52880952525675 |
| $A_{31}$ | 0.0 |
| $A_{32}$ | 0.51578733443615 |
| $A_{41}$ | 0.53289582961739 |
| $A_{42}$ | 0.25574324768195 |
| $A_{43}$ | 0.21136092270067 |

Noise coefficients:

| | Value |
|---|---|
| $Q_1$ | 1.87653936176981 |
| $Q_2$ | 3.91017166264989 |
| $Q_3$ | 4.73124353935667 |

---

### `rk4_ti_solve(F, G, X, T0, TN, X0, N, Q, SEED)`

**Four-stage stochastic Runge-Kutta** — highest-order time-invariant solver.

$$
X_{k+1}=X_1 + A_{51}K_1 + A_{52}K_2 + A_{53}K_3 + A_{54}K_4.
$$

Kasdin time-invariant coefficients:

| Stage weights | Value |
|---|---|
| $A_{21}$ | 2.71644396264860 |
| $A_{31}$ | −6.95653259006152 |
| $A_{32}$ | 0.78313689457981 |
| $A_{41}$ | 0.0 |
| $A_{42}$ | 0.48257353309214 |
| $A_{43}$ | 0.26171080165848 |
| $A_{51}$ | 0.47012396888046 |
| $A_{52}$ | 0.36597075368373 |
| $A_{53}$ | 0.08906615686702 |
| $A_{54}$ | 0.07483912056879 |

Noise coefficients:

| | Value |
|---|---|
| $Q_1$ | 2.12709852335625 |
| $Q_2$ | 2.73245878238737 |
| $Q_3$ | 11.22760917474960 |
| $Q_4$ | 13.36199560336697 |

!!! note "Weight sum"
    The final-stage weights sum to 1:
    $A_{51}+A_{52}+A_{53}+A_{54} = 1.0$, ensuring consistency of the method.

---

## Time-variant solvers

These solvers integrate SDEs where the drift and diffusion depend explicitly
on time $t$:

$$
dX(t)=F(X,t)\,dt+Q\,G(X,t)\,dW(t).
$$

### `rk1_tv_solve(F, G, X, T0, TN, X0, N, Q, SEED)`

**Euler-Maruyama with time-dependent coefficients.**

$$
X_{k+1}=X_k+H F(X_k,T_k)+\sqrt{HQ}\,G(X_k,T_k)\,Z_k,
$$

with $A_{21}=1$, $Q_1=1$.

---

### `rk2_tv_solve(F, G, X, T0, TN, X0, N, Q, SEED)`

**Two-stage time-dependent stochastic Runge-Kutta.**

Same weight structure as the time-invariant RK2:

$$
A_{21}=1,\quad A_{31}=0.5,\quad A_{32}=0.5,\quad Q_1=Q_2=2.
$$

The time argument $T$ is advanced to $T_2 = T_1 + A_{21}H$ at the second stage.

---

### `rk4_tv_solve(F, G, X, T0, TN, X0, N, Q, SEED)`

**Four-stage time-dependent stochastic Runge-Kutta** — highest-order
time-variant solver. Uses a **different** Kasdin coefficient set than the
time-invariant method.

Kasdin time-variant coefficients:

| Stage weights | Value |
|---|---|
| $A_{21}$ | 0.66667754298442 |
| $A_{31}$ | 0.63493935027993 |
| $A_{32}$ | 0.00342761715422 |
| $A_{41}$ | −2.32428921184321 |
| $A_{42}$ | 2.69723745129487 |
| $A_{43}$ | 0.29093673271592 |
| $A_{51}$ | 0.25001351164789 |
| $A_{52}$ | 0.67428574806272 |
| $A_{53}$ | −0.00831795169360 |
| $A_{54}$ | 0.08401868181222 |

Noise coefficients:

| | Value |
|---|---|
| $Q_1$ | 3.99956364361748 |
| $Q_2$ | 1.64524970733585 |
| $Q_3$ | 1.59330355118722 |
| $Q_4$ | 0.26330006501868 |

---

## References

See [References](references.md) for full citations of the Kasdin papers and
other source literature.
