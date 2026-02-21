# Quickstart

## Example 1: RK1 time-invariant solver

Solve

$$
dX(t)=1\,dt + 1\cdot dW(t),\quad X(0)=0,
$$

with $N=10$, $T_0=0$, $T_N=1$, $Q=1$, seed $=123456789$.

```py
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

## Example 2: RK4 time-variant solver

Solve

$$
dX(t)=\left(-0.5X + 0.1\cos(2t)\right)dt + 0.05\left(1+0.25\sin t\right)dW(t),\quad X(0)=1.
$$

with $N=10$, $T_0=0$, $T_N=2$, seed $=99999$.

```py
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



