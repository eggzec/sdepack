# Installation

`sdepack` can be installed from `pypi`, `conda-forge`, and `git`.

## PyPI

### `pyproject.toml`

```toml
[project]
dependencies = [
	"sdepack"
]
```

### `requirements.txt`

```text
sdepack
```

## pip

```bash
pip install --upgrade sdepack
```

## pipenv

```bash
pipenv install sdepack
```

## uv

### Add to a uv project

```bash
uv add sdepack
```

### Install into the current uv environment

```bash
uv pip install sdepack
```

## poetry

```bash
poetry add sdepack
```

## pdm

```bash
pdm add sdepack
```

## hatch

```bash
hatch add sdepack
```

## GitHub

Install the latest code from the repository:

```bash
pip install --upgrade "git+https://github.com/eggzec/sdepack.git#egg=sdepack"
```

## Dependencies

- Python 3.10+
- numpy
