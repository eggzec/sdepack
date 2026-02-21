# Theory

## 1) Scalar Itô SDE form

All solvers target scalar Itô SDEs:

$$
dX(t) = a(X,t)\,dt + b(X,t)\,dW(t).
$$

Within the solver routines the equivalent form is:

$$
dX(t) = F(X,t)\,dt + Q\,G(X,t)\,dW(t),
$$

where $Q$ scales the stochastic forcing strength.

For time-invariant routines, $F(X,t) \to F(X)$ and $G(X,t) \to G(X)$.

## 2) Time discretization and noise scaling

For $N$ steps over $[T_0, T_N]$:

$$
H = \frac{T_N - T_0}{N}, \qquad T_k = T_0 + kH.
$$

Stage noise uses:

$$
W_i = Z_i\sqrt{\frac{Q_i Q}{H}},\qquad Z_i \sim \mathcal{N}(0,1),
$$

giving $\sqrt{H}$ stochastic scaling once multiplied by $H$ in the stage equations.

## 3) RK1 (Euler-Maruyama)

For both time-invariant and time-variant variants:

$$
K_1 = H\,F(X_1,T_1) + H\,G(X_1,T_1)W_1,
$$

$$
X_{k+1} = X_1 + A_{21}K_1,\quad A_{21}=1,\;Q_1=1.
$$

For time-invariant routines, omit the explicit time argument in $F$ and $G$.

## 4) RK2

The two-stage update is:

$$
\begin{aligned}
K_1 &= H\,F(X_1,T_1) + H\,G(X_1,T_1)W_1,\\
X_2 &= X_1 + A_{21}K_1,\quad T_2 = T_1 + A_{21}H,\\
K_2 &= H\,F(X_2,T_2) + H\,G(X_2,T_2)W_2,\\
X_{k+1} &= X_1 + A_{31}K_1 + A_{32}K_2,
\end{aligned}
$$

with

$$
A_{21}=1,\;A_{31}=\tfrac{1}{2},\;A_{32}=\tfrac{1}{2},\;Q_1=Q_2=2.
$$

## 5) RK3 (time-invariant)

The time-invariant 3-stage method uses:

$$
X_{k+1}=X_1 + A_{41}K_1 + A_{42}K_2 + A_{43}K_3,
$$

with coefficients:

$$
\begin{aligned}
A_{21}&=1.52880952525675,\\
A_{31}&=0.0,\quad A_{32}=0.51578733443615,\\
A_{41}&=0.53289582961739,\quad A_{42}=0.25574324768195,\quad A_{43}=0.21136092270067,
\end{aligned}
$$

and stochastic stage parameters:

$$
Q_1=1.87653936176981,\;Q_2=3.91017166264989,\;Q_3=4.73124353935667.
$$

## 6) RK4 (time-invariant and time-variant)

The 4-stage update is:

$$
X_{k+1} = X_1 + A_{51}K_1 + A_{52}K_2 + A_{53}K_3 + A_{54}K_4.
$$

The time-invariant and time-variant methods use different Kasdin coefficient sets (see [API Reference](api.md)).

## 7) Random number generation

`R8_UNIFORM` and `R8_NORMAL` are internal Fortran routines; they are not exported to Python.

### `R8_UNIFORM`

Park-Miller Lehmer generator:

$$
\mathrm{seed}_{n+1} = 16807\,\mathrm{seed}_n \bmod (2^{31}-1),
$$

and returns

$$
u_n = \frac{\mathrm{seed}_n}{2^{31}-1} \in (0,1).
$$

Schrage's method avoids integer overflow.

### `R8_NORMAL`

Box-Muller transform:

$$
v_1 = \sqrt{-2\ln u_1}\cos(2\pi u_2),\qquad
v_2 = \sqrt{-2\ln u_1}\sin(2\pi u_2),
$$

with one value cached for the next call.
