# sdepack

**Runge-Kutta Numerical Integration of Stochastic Differential Equations for Python**

[![Tests](https://github.com/eggzec/sdepack/actions/workflows/test.yml/badge.svg)](https://github.com/eggzec/sdepack/actions/workflows/test.yml)
[![Documentation](https://github.com/eggzec/sdepack/actions/workflows/docs.yml/badge.svg)](https://github.com/eggzec/sdepack/actions/workflows/docs.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![codecov](https://codecov.io/github/eggzec/sdepack/graph/badge.svg)](https://codecov.io/github/eggzec/sdepack)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[![PyPI Downloads](https://img.shields.io/pypi/dm/sdepack.svg?label=PyPI%20downloads)](https://pypi.org/project/sdepack/)
[![Python versions](https://img.shields.io/pypi/pyversions/sdepack.svg)](https://pypi.org/project/sdepack/)

`sdepack` is a high-performance Python library for numerically solving scalar
[stochastic differential equations (SDEs)](https://en.wikipedia.org/wiki/Stochastic_differential_equation)
using stochastic Runge-Kutta methods. It offers deterministic, seed-controlled integration of Itô SDEs of the form:

$$dX(t) = F(X, t)\,dt + Q\,G(X, t)\,dW(t)$$

Solvers range from the first-order **Euler-Maruyama** scheme to fourth-order
stochastic Runge-Kutta methods using Kasdin coefficients.

## Quick example

```python
import numpy as np
import sdepack

x = np.zeros(101, dtype=np.float64)
sdepack.rk4_ti_solve(
    lambda x: -0.5 * x,     # drift
    lambda x: 1.0,           # diffusion
    x, 0.0, 10.0, 1.0, 100, 0.1, 42
)
```

## Installation

```bash
pip install sdepack
```

Requires Python 3.10+ and NumPy. See the
[full installation guide](https://eggzec.github.io/sdepack/installation/) for
conda, poetry, and source builds.

## Documentation

- [Theory](https://eggzec.github.io/sdepack/theory/) — mathematical background, Itô SDEs, convergence
- [Quickstart](https://eggzec.github.io/sdepack/quickstart/) — runnable examples
- [API Reference](https://eggzec.github.io/sdepack/api/) — solver signatures and parameters
- [References](https://eggzec.github.io/sdepack/references/) — literature citations

## License

MIT — see [LICENSE](LICENSE).
