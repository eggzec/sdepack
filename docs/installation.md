# Installation

`sdepack` can be installed from PyPI, conda-forge, or directly from git.

## [PyPI](https://pypi.org/project/sdepack)

For using the PyPI package in your project, you can update your configuration file by adding
the following snippet.

=== "pyproject.toml"

    ```toml
    [project.dependencies]
    sdepack = "*" # (1)!
    ```

    1. Specifying a version is recommended

=== "requirements.txt"

    ```
    sdepack>=0.1.0
    ```

### pip

=== "Installation for user"

    ```bash
    pip install --upgrade --user sdepack # (1)!
    ```

    1. You may need to use `pip3` instead of `pip` depending on your Python installation.

=== "Installation in virtual environment"

    ```bash
    python -m venv .venv
    source .venv/bin/activate
    pip install --require-virtualenv --upgrade sdepack # (1)!
    ```

    1. You may need to use `pip3` instead of `pip` depending on your Python installation.

    !!! note
        Command to activate the virtual env depends on your platform and shell. [More info](https://docs.python.org/3/library/venv.html#how-venvs-work)

### pipenv

    pipenv install sdepack

### uv

=== "Adding to uv project"

    ```bash
    uv add sdepack
    uv sync
    ```

=== "Installing to uv environment"

    ```bash
    uv venv
    uv pip install sdepack
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

## [conda-forge](https://anaconda.org/conda-forge/sdepack)

You can update your environment spec file by adding the following snippet.

```yaml title="environment.yml"
    channels:
    - conda-forge
    dependencies:
    - pip
    - pip:
        - sdepack # (1)!
```

1. Specifying a version is recommended

Installation can be done using the updated environment spec file.

=== "conda"

    ```bash
    conda env update --file environment.yml
    ```

=== "micromamba"

    ```bash
    micromamba env update --file environment.yml
    ```

!!! note
    Replace `environment.yml` with your actual environment spec file name if it's different.

## [git](https://github.com/eggzec/sdepack)

Install the latest development version directly from the repository:

```bash
pip install --upgrade "git+https://github.com/eggzec/sdepack.git#egg=sdepack"
```

### Building from source

Clone and build from source if you want to modify the Fortran code or test local changes:

```bash
git clone https://github.com/eggzec/sdepack.git
cd sdepack
pip install -e .
```

!!! warning "Fortran compiler required"
    Source builds require a working Fortran compiler (`gfortran` recommended) as well as
    `meson` and `meson-python`.

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

- Python >= 3.10
- [numpy](https://pypi.org/project/numpy)
