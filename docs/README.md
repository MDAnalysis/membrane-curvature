# Building MembraneCurvature's Documentation

The MembraneCurvature documentation is built with [Sphinx] using the
[mdanalysis-sphinx-theme]. The published docs are available on
[Read the Docs][documentation].

To build the docs locally, install the documentation dependencies and run
`sphinx-build` from the root of the repository.

## Installing dependencies

Documentation dependencies are defined in the `docs` dependency group in
`pyproject.toml`. Note that by installing the `dev` group, the `docs` dependency
group is also installed.

### With uv (recommended)

We recommend using [uv] as a development environment for MembraneCurvature. From
the root of the repository, run:

```bash
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
uv sync --group docs
```

To install all development dependencies (including `docs` and `tests`), use
`uv sync --group dev` instead.

### With pip

Install MembraneCurvature in editable mode and the `docs` dependency group:

```bash
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
python -m pip install -e .
pip install --group docs
```

Alternatively, install the dependencies listed in `docs/requirements.txt`:

```bash
python -m pip install -e .
pip install -r docs/requirements.txt
```

> [!NOTE]
>
> The command `pip install --group docs` requires `pip >= 25.1`.

### With conda

To create a conda environment with the documentation dependencies, run:

```bash
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
conda env create -f docs/requirements.yaml
conda activate mc_docs
python -m pip install -e .
```

Some notebooks and examples in the documentation use test data from
[MDAnalysisTests] and [MDAnalysisData]. These packages are included in
`docs/requirements.txt` and `docs/requirements.yaml`.

## Building the documentation

From the root of the repository, run:

```bash
uv run sphinx-build -b html docs/ build/html
```

If you are not using [uv], run `sphinx-build -b html docs/ build/html` instead.

The compiled docs will be in `build/html/` and can be viewed by opening
`build/html/index.html`.

Alternatively, move to the `docs/` directory and use the `Makefile`:

```bash
cd docs
make html
```

The compiled docs will be in `docs/_build/html/` and can be viewed by opening
`docs/_build/html/index.html`.

On Windows, use `make.bat` instead of `make`.

> [!TIP]
>
> Some documentation pages include Jupyter notebooks built with `myst-nb`. If
> the build fails on a notebook, check that the `docs` dependencies are
> installed and that MembraneCurvature is available in editable mode.

For more information on contributing to MembraneCurvature, see the
[Contributing] page.


[Sphinx]: https://www.sphinx-doc.org/
[mdanalysis-sphinx-theme]: https://github.com/MDAnalysis/mdanalysis-sphinx-theme
[documentation]: https://membrane-curvature.readthedocs.io/en/latest/
[uv]: https://docs.astral.sh/uv/
[MDAnalysisTests]: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
[MDAnalysisData]: https://www.mdanalysis.org/MDAnalysisData/
[Contributing]: https://github.com/MDAnalysis/membrane-curvature/blob/main/.github/CONTRIBUTING.md
