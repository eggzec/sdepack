# Theory

This page covers the mathematical foundations behind `sdepack`: what stochastic
differential equations are, how Itô calculus underpins them, and how the
stochastic Runge-Kutta methods implemented in this library discretize and solve
them numerically.

---

## 1) Background: Stochastic Differential Equations

A **stochastic differential equation (SDE)** is a differential equation in which
one or more terms are stochastic processes, producing a solution that is itself a
random process. SDEs originated in the theory of Brownian motion through the work
of Albert Einstein and Marian Smoluchowski in 1905, with the mathematical
foundation established by Kiyosi Itô in the 1940s.

SDEs are used to model systems subject to random perturbations across many
disciplines:

- **Physics** — molecular dynamics, thermal fluctuations, Langevin equations
- **Finance** — option pricing (Black-Scholes), interest rate models (Vasicek, CIR)
- **Biology** — population dynamics, gene expression, neural firing
- **Engineering** — noise in control systems, signal processing, spacecraft attitude dynamics

### The Wiener process

The randomness in an SDE is driven by a **Wiener process** $W(t)$ (standard
Brownian motion), characterized by:

1. $W(0) = 0$
2. $W(t)$ has **independent increments**: $W(t) - W(s)$ is independent of
   $\{W(u) : u \le s\}$ for $s < t$
3. $W(t) - W(s) \sim \mathcal{N}(0, t - s)$ — increments are normally
   distributed with variance equal to the elapsed time
4. $W(t)$ is continuous in $t$ (with probability 1)

The Wiener process is almost surely **nowhere differentiable**, which is why
standard calculus cannot be directly applied to SDEs. This motivates the
development of stochastic calculus.

### Itô vs. Stratonovich interpretation

There are two principal ways to interpret the stochastic integral that appears
in an SDE:

| | **Itô** | **Stratonovich** |
|---|---|---|
| Evaluation point | Left endpoint of each sub-interval | Midpoint of each sub-interval |
| Martingale property | Preserved | Not preserved in general |
| Chain rule | Modified (Itô's lemma) — has an extra $\frac{1}{2}\sigma^2 f''$ term | Standard chain rule applies |
| Typical use | Mathematical finance, filtering theory | Physics (Langevin equations), manifolds |

`sdepack` uses the **Itô interpretation** exclusively. For additive noise
(constant diffusion coefficient $G$), both interpretations coincide, so this
distinction matters primarily for multiplicative noise.

!!! info "Itô's lemma"
    For a twice-differentiable function $f$ and Itô process
    $dX = \mu\,dt + \sigma\,dW$:

    $$df(X_t) = f'(X_t)\,\mu_t\,dt + \frac{1}{2}f''(X_t)\,\sigma_t^2\,dt + f'(X_t)\,\sigma_t\,dW_t.$$

    The additional $\frac{1}{2}f''\sigma^2$ term (absent in ordinary calculus)
    arises because Brownian motion has non-zero quadratic variation:
    $[W, W]_t = t$.

## 2) Scalar Itô SDE form

All solvers in `sdepack` target **scalar** Itô SDEs:

$$
dX(t) = a(X,t)\,dt + b(X,t)\,dW(t).
$$

Within the solver routines, this is parameterized as:

$$
dX(t) = F(X,t)\,dt + Q\,G(X,t)\,dW(t),
$$

where:

- $F$ is the user-supplied **drift** function
- $G$ is the user-supplied **diffusion** function
- $Q$ is the **noise intensity** (spectral density of the driving white noise)
- $dW(t)$ denotes a Wiener process increment

For **time-invariant** routines, $F(X,t) \to F(X)$ and $G(X,t) \to G(X)$ —
the drift and diffusion do not depend explicitly on time.

!!! note "Additive vs. multiplicative noise"
    If $G(X) = \text{const}$, the noise is **additive** and the Itô/Stratonovich
    distinction vanishes. If $G$ depends on $X$, the noise is **multiplicative**
    and the choice of calculus matters.

## 3) Time discretization and noise scaling

Given the interval $[T_0, T_N]$ and $N$ time steps:

$$
H = \frac{T_N - T_0}{N}, \qquad T_k = T_0 + kH.
$$

At each Runge-Kutta stage $i$, a noise variate is generated:

$$
W_i = Z_i\sqrt{\frac{Q_i Q}{H}},\qquad Z_i \sim \mathcal{N}(0,1),
$$

where $Q_i$ is a method-specific coefficient that controls the variance of the
noise at stage $i$. When multiplied by $H$ in the stage equations, this
produces the correct $\sqrt{H}$ stochastic scaling required for Wiener
increments.

## 4) RK1 — Euler-Maruyama method

The simplest and most widely known method for numerically solving SDEs. Named
after Leonhard Euler and Gisiro Maruyama, it is the stochastic analog of the
forward Euler method for ODEs.

For both time-invariant and time-variant variants:

$$
K_1 = H\,F(X_1,T_1) + H\,G(X_1,T_1)W_1,
$$

$$
X_{k+1} = X_1 + A_{21}K_1,\quad A_{21}=1,\;Q_1=1.
$$

This simplifies to:

$$
X_{k+1} = X_k + H\,F(X_k, T_k) + \sqrt{HQ}\,G(X_k, T_k)\,Z_k,\quad Z_k \sim \mathcal{N}(0,1).
$$

For time-invariant routines, omit the explicit time argument in $F$ and $G$.

!!! info "Convergence of Euler-Maruyama"
    The Euler-Maruyama method has **strong order** $\gamma_s = \tfrac{1}{2}$
    and **weak order** $\gamma_w = 1$, provided the drift $\mu$ and diffusion
    $\sigma$ satisfy Lipschitz continuity and linear growth conditions.
    Strong convergence means:
    $E[|X_T - Y_N|] = O(\Delta t^{1/2})$,
    while weak convergence means:
    $|E[g(X_T)] - E[g(Y_N)]| = O(\Delta t)$
    for any sufficiently smooth test function $g$.

## 5) RK2 — Two-stage stochastic Runge-Kutta

A two-stage method that provides improved accuracy over Euler-Maruyama:

$$
\begin{aligned}
K_1 &= H\,F(X_1,T_1) + H\,G(X_1,T_1)W_1,\\
X_2 &= X_1 + A_{21}K_1,\quad T_2 = T_1 + A_{21}H,\\
K_2 &= H\,F(X_2,T_2) + H\,G(X_2,T_2)W_2,\\
X_{k+1} &= X_1 + A_{31}K_1 + A_{32}K_2,
\end{aligned}
$$

with Kasdin coefficients:

$$
A_{21}=1,\;A_{31}=\tfrac{1}{2},\;A_{32}=\tfrac{1}{2},\;Q_1=Q_2=2.
$$

This is analogous to Heun's method (improved Euler) for deterministic ODEs,
with the additional stochastic noise terms at each stage.

## 6) RK3 — Three-stage stochastic Runge-Kutta (time-invariant only)

A three-stage method providing higher accuracy for time-invariant SDEs:

$$
X_{k+1}=X_1 + A_{41}K_1 + A_{42}K_2 + A_{43}K_3,
$$

with Kasdin coefficients:

$$
\begin{aligned}
A_{21}&=1.52880952525675,\\
A_{31}&=0.0,\quad A_{32}=0.51578733443615,\\
A_{41}&=0.53289582961739,\quad A_{42}=0.25574324768195,\quad A_{43}=0.21136092270067,
\end{aligned}
$$

and stochastic stage noise parameters:

$$
Q_1=1.87653936176981,\;Q_2=3.91017166264989,\;Q_3=4.73124353935667.
$$

!!! note
    The RK3 solver is only available for time-invariant systems. There is no
    `rk3_tv_solve` in sdepack.

## 7) RK4 — Four-stage stochastic Runge-Kutta

The highest-order methods in `sdepack`, using four stages for maximum accuracy.

The general update formula is:

$$
X_{k+1} = X_1 + A_{51}K_1 + A_{52}K_2 + A_{53}K_3 + A_{54}K_4.
$$

The time-invariant and time-variant variants use **different** Kasdin coefficient
sets, optimized for their respective problem classes. See the
[API Reference](api.md) for the complete coefficient tables.

## 8) Convergence and accuracy considerations

### Strong vs. weak convergence

When assessing the quality of a numerical SDE solver, two notions of convergence
are distinguished:

- **Strong convergence** — accuracy of individual sample paths. A method has
  strong order $\gamma$ if:
  $$E[|X_T - Y_N|] \le C \cdot (\Delta t)^\gamma.$$

- **Weak convergence** — accuracy of statistical quantities (expectations). A
  method has weak order $\gamma$ if for all sufficiently smooth $g$:
  $$|E[g(X_T)] - E[g(Y_N)]| \le C \cdot (\Delta t)^\gamma.$$

Weak convergence is often sufficient for Monte Carlo applications (computing
expectations, pricing derivatives), while strong convergence matters when
tracking individual trajectories or coupling paths.

### Relationship between solvers and order

The Kasdin Runge-Kutta schemes in `sdepack` are designed to achieve higher
weak-order convergence. The Euler-Maruyama method (RK1) achieves weak order 1
and strong order ½. Higher-order stochastic RK methods generally improve the
accuracy-per-step, allowing larger step sizes for the same error tolerance.

!!! tip "Step size guidance"
    Unlike deterministic ODEs, halving the step $H$ in an SDE solver only
    reduces the strong error by a factor of $\sqrt{2}$ (for Euler-Maruyama).
    Higher-order methods like RK4 benefit more from step refinement.

## 9) Random number generation

The solvers use a self-contained random number generation pipeline with two
components:

### Park-Miller LCG

A multiplicative linear congruential generator (LCG) based on the Park-Miller
"minimal standard" PRNG:

$$
\mathrm{seed}_{n+1} = 16807\,\mathrm{seed}_n \bmod (2^{31}-1),
$$

returning:

$$
u_n = \frac{\mathrm{seed}_n}{2^{31}-1} \in (0,1).
$$

**Schrage's method** is used to compute the modular arithmetic without
intermediate overflow on 32-bit integers. The modulus $2^{31}-1 = 2\,147\,483\,647$
is a Mersenne prime, ensuring the generator has maximal period $2^{31}-2$.

### Box-Muller transform

Converts pairs of uniform variates into standard normal variates using the
Box-Muller transform:

$$
v_1 = \sqrt{-2\ln u_1}\cos(2\pi u_2),\qquad
v_2 = \sqrt{-2\ln u_1}\sin(2\pi u_2),
$$

where $u_1, u_2 \sim \mathrm{Uniform}(0,1)$. One value is returned immediately,
and the other is cached for the next call, amortizing the cost of the
transcendental functions.

### Why a built-in RNG?

The use of a self-contained, deterministic PRNG (rather than the system or
NumPy RNG) ensures **exact bitwise reproducibility** of trajectories given the
same seed — regardless of platform, Python version, or NumPy version. This is
critical for regression testing and for comparing solver outputs across
different configurations.
