# Installation

`sdepack` is distributed as a compiled wheel on PyPI and can also be installed
from source via GitHub.

---

## Prerequisites

- **Python 3.10+**
- **NumPy** (installed automatically as a dependency)

For source builds you additionally need:

- A Fortran compiler (`gfortran` recommended)
- `meson` and `meson-python` build system
- `numpy` (for `f2py` compilation)

## PyPI (recommended)

### pip

```bash
pip install --upgrade sdepack
```

### pyproject.toml dependency

```toml
[project]
dependencies = [
    "sdepack"
]
```

### requirements.txt

```text
sdepack
```

## Package managers

### uv

```bash
# Add to a uv project
uv add sdepack

# Or install into the current environment
uv pip install sdepack
```

### pipenv

```bash
pipenv install sdepack
```

### poetry

```bash
poetry add sdepack
```

### pdm

```bash
pdm add sdepack
```

### hatch

```bash
hatch add sdepack
```

## Installing from source (GitHub)

Install the latest development version directly from the repository:

```bash
pip install --upgrade "git+https://github.com/eggzec/sdepack.git#egg=sdepack"
```

### Building locally

Clone and build from source if you want to modify the Fortran code or test
local changes:

```bash
git clone https://github.com/eggzec/sdepack.git
cd sdepack
pip install -e .
```

This invokes the `meson` build system to compile the Fortran sources via
`f2py` and install the resulting extension module in development mode.

!!! warning "Fortran compiler required"
    Source builds require a working Fortran compiler. On most Linux
    distributions install `gfortran`:

```bash
# Debian/Ubuntu
sudo apt install gfortran

# Fedora
sudo dnf install gcc-gfortran

# macOS (Homebrew)
brew install gcc
```

On Windows, install MinGW-w64 with gfortran or use MSYS2.

## Verifying the installation

After installation, verify that the package loads correctly:

```python
import sdepack
import numpy as np

x = np.zeros(11, dtype=np.float64)
sdepack.rk1_ti_solve(lambda x: 0.0, lambda x: 1.0, x, 0.0, 1.0, 0.0, 10, 1.0, 42)
print("sdepack is working! Trajectory:", x)
```

## Dependencies

| Package | Purpose |
|---|---|
| `numpy` | Array handling, `f2py` integration |

No other runtime dependencies are required.
