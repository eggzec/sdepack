# sdepack

Runge-Kutta numerical integration of scalar stochastic differential equations (SDEs) for Python.

## Features

- Stochastic Runge-Kutta solvers (order 1 through 4) for time-invariant and time-variant scalar SDEs.
- Deterministic seeded random number generation via `R8_UNIFORM` and `R8_NORMAL` (internal to solvers).

## Mathematical model

Time-invariant routines integrate

$$
dX(t) = F(X)\,dt + Q\,G(X)\,dW(t),
$$

and time-variant routines integrate

$$
dX(t) = F(X,t)\,dt + Q\,G(X,t)\,dW(t),
$$

with uniform step size

$$
H = \frac{T_N - T_0}{N}.
$$

Stage noise variates follow

$$
W_i = Z_i\sqrt{\frac{Q_i Q}{H}},\quad Z_i \sim \mathcal{N}(0,1),
$$

with method-specific $Q_i$ coefficients.

## Public API

- `rk1_ti_solve`
- `rk1_tv_solve`
- `rk2_ti_solve`
- `rk2_tv_solve`
- `rk3_ti_solve`
- `rk4_ti_solve`
- `rk4_tv_solve`

`R8_UNIFORM` and `R8_NORMAL` are internal Fortran helpers, not exported to Python.

## Documentation

- [Theory](theory.md) — numerical background and coefficients
- [Quickstart](quickstart.md) — runnable examples
- [API Reference](api.md) — Fortran routine signatures
- [References](references.md) — source literature
