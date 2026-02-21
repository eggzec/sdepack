# API Reference

Only routines declared in `src/sdepack.pyf` are exported to Python.

## Exported routines

- `rk1_ti_solve`
- `rk1_tv_solve`
- `rk2_ti_solve`
- `rk2_tv_solve`
- `rk3_ti_solve`
- `rk4_ti_solve`
- `rk4_tv_solve`

## Common interface

All solvers write to `X(0:N)` in place with step size:

$$
H = \frac{T_N - T_0}{N}.
$$

### Common arguments

- `F`, `G`: external drift and diffusion callbacks
- `X(0:N)`: solution array
- `T0`, `TN`: interval endpoints
- `X0`: initial condition
- `N`: number of steps
- `Q`: noise intensity parameter
- `SEED`: integer RNG seed

## Time-invariant solvers

Solve:

$$
dX(t)=F(X)dt+QG(X)dW(t).
$$

### `RK1_TI_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)`

Euler-Maruyama / RK1:

$$
X_{k+1}=X_k+H F(X_k)+H G(X_k)W_1,
$$

with $A_{21}=1$, $Q_1=1$.

### `RK2_TI_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)`

Two-stage SRK:

$$
\begin{aligned}
K_1 &= H F(X_1)+H G(X_1)W_1,\\
X_2 &= X_1 + A_{21}K_1,\\
K_2 &= H F(X_2)+H G(X_2)W_2,\\
X_{k+1} &= X_1 + A_{31}K_1 + A_{32}K_2,
\end{aligned}
$$

with $A_{21}=1$, $A_{31}=0.5$, $A_{32}=0.5$, $Q_1=Q_2=2$.

### `RK3_TI_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)`

Three-stage SRK:

$$
X_{k+1} = X_1 + A_{41}K_1 + A_{42}K_2 + A_{43}K_3,
$$

coefficients:

$$
\begin{aligned}
A_{21}&=1.52880952525675,\\
A_{31}&=0.0,\quad A_{32}=0.51578733443615,\\
A_{41}&=0.53289582961739,\quad A_{42}=0.25574324768195,\quad A_{43}=0.21136092270067,
\end{aligned}
$$

noise coefficients:

$$
Q_1=1.87653936176981,\quad Q_2=3.91017166264989,\quad Q_3=4.73124353935667.
$$

### `RK4_TI_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)`

Four-stage SRK:

$$
X_{k+1}=X_1 + A_{51}K_1 + A_{52}K_2 + A_{53}K_3 + A_{54}K_4.
$$

Kasdin time-invariant coefficients in code:

$$
\begin{aligned}
&A_{21}=2.71644396264860,\;A_{31}=-6.95653259006152,\;A_{32}=0.78313689457981,\\
&A_{41}=0.0,\;A_{42}=0.48257353309214,\;A_{43}=0.26171080165848,\\
&A_{51}=0.47012396888046,\;A_{52}=0.36597075368373,\;A_{53}=0.08906615686702,\;A_{54}=0.07483912056879,
\end{aligned}
$$

$$
Q_1=2.12709852335625,\;Q_2=2.73245878238737,\;Q_3=11.22760917474960,\;Q_4=13.36199560336697.
$$

## Time-variant solvers

Solve:

$$
dX(t)=F(X,t)dt+QG(X,t)dW(t).
$$

### `RK1_TV_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)`

Euler-Maruyama with time-dependent coefficients:

$$
X_{k+1}=X_k+H F(X_k,T_k)+H G(X_k,T_k)W_1,
$$

with $A_{21}=1$, $Q_1=1$.

### `RK2_TV_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)`

Two-stage time-dependent SRK with the same RK2 weights as the time-invariant method:

$$
A_{21}=1,\quad A_{31}=0.5,\quad A_{32}=0.5,\quad Q_1=Q_2=2.
$$

### `RK4_TV_SOLVE(F, G, X, T0, TN, X0, N, Q, SEED)`

Four-stage time-dependent SRK with Kasdin coefficients for the time-variant method:

$$
\begin{aligned}
&A_{21}=0.66667754298442,\;A_{31}=0.63493935027993,\;A_{32}=0.00342761715422,\\
&A_{41}=-2.32428921184321,\;A_{42}=2.69723745129487,\;A_{43}=0.29093673271592,\\
&A_{51}=0.25001351164789,\;A_{52}=0.67428574806272,\;A_{53}=-0.00831795169360,\;A_{54}=0.08401868181222,
\end{aligned}
$$

$$
Q_1=3.99956364361748,\;Q_2=1.64524970733585,\;Q_3=1.59330355118722,\;Q_4=0.26330006501868.
$$

## Internal helper routines (not exported in `sdepack.pyf`)

### `R8_UNIFORM(SEED)`

Returns a uniform random value in $(0,1)$ using Park-Miller LCG with Schrage overflow protection:

$$
\mathrm{seed}_{n+1} = 16807\,\mathrm{seed}_n\;\bmod(2^{31}-1),\qquad
u_n = \frac{\mathrm{seed}_n}{2^{31}-1}.
$$

### `R8_NORMAL(SEED)`

Returns a standard normal variate using Box-Muller transform with cached second value:

$$
v_1 = \sqrt{-2\ln u_1}\cos(2\pi u_2),\qquad
v_2 = \sqrt{-2\ln u_1}\sin(2\pi u_2).
$$

## References

See [References](references.md) for full citations.
