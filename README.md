Membrane Curvature
==============================
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![GitHub Actions Status](https://github.com/MDAnalysis/membrane-curvature/workflows/CI/badge.svg)](https://github.com/MDAnalysis/membrane-curvature/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/MDAnalysis/membrane-curvature/branch/main/graph/badge.svg)](https://codecov.io/gh/MDAnalysis/membrane-curvature/branch/main)
[![docs](https://readthedocs.org/projects/membrane-curvature/badge/?version=latest)](https://membrane-curvature.readthedocs.io/en/latest/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![ty](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ty/main/assets/badge/v0.json)](https://github.com/astral-sh/ty)
[![License](https://img.shields.io/github/license/MDAnalysis/membrane-curvature?color=yellow)](https://github.com/MDAnalysis/membrane-curvature/blob/main/LICENSE)

[![Python versions](https://img.shields.io/pypi/pyversions/membrane-curvature)](https://pypi.org/project/membrane-curvature/)
[![PyPI](https://img.shields.io/pypi/v/membrane-curvature)](https://pypi.org/project/membrane-curvature/)
[![Conda version](https://anaconda.org/conda-forge/membrane-curvature/badges/version.svg)](https://anaconda.org/conda-forge/membrane-curvature)

![](https://github.com/MDAnalysis/membrane-curvature/blob/main/docs/source/_static/PM_Membrane_EBO.png?raw=true)

MembraneCurvature is an [MDAnalysis] tool to calculate membrane curvature from
Molecular Dynamics simulations.

> **Interested in becoming a maintainer?** We welcome your passion and expertise to help shape and grow this open-source project! Please contact estefania@ojeda-e.com for more details.

Features
--------------

With MembraneCurvature you can:

- Calculate mean and Gaussian curvature from MD simulations.
- Derive 2D curvature profiles from atoms of reference with two different methods: binning or Fourier-based.
- Get per-frame or averaged results for surface, mean and Gaussian curvature.
- Live a happier life.

Installation
--------------

The main dependency in MembraneCurvature is [MDAnalysis]. You can find
instructions to install the latest stable version of MDAnalysis in the [UserGuide].

MembraneCurvature is available via `pip`:

```bash
pip install membrane-curvature
```

To install from source:

```bash
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
python -m pip install -e .
```

Some of the examples included in the MembraneCurvature documentation use test
data from [MDAnalysisTests] and [MDAnalysisData]. To install these dependencies
with conda, run:

```bash
conda install -c conda-forge MDAnalysisTests MDAnalysisData
```

or via `pip`:

```bash
pip install --upgrade MDAnalysisTests MDAnalysisData
```

Usage
--------------

This is a quick example on how to run MembraneCurvature with the default
surface method (Fourier):

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# run with the default surface_method - Fourier
curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4'
                                            ).run()

# extract average mean curvature
mean_upper_leaflet = curvature_upper_leaflet.results.z_surface

# extract average mean curvature
mean_upper_leaflet = curvature_upper_leaflet.results.mean

# extract average Gaussian curvature
gaussian_upper_leaflet = curvature_upper_leaflet.results.gaussian
```

In this example, we use the PO4 beads in the upper leaflet as reference to
derive a surface and calculate its respective mean and Gaussian curvature.
If you want per-frame arrays instead, use `results.z_surface`, `results.mean`,
and `results.gaussian`.

The same example run with the binning surface method looks like:

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# run with the binning surface_method
curvature_upper_leaflet_binning = MembraneCurvature(universe,
                                                    select='resid 1-225 and name PO4',
                                                    surface_method='binning',
                                                    n_x_bins=8,
                                                    n_y_bins=8,
                                                    wrap=True).run()

# extract average mean curvature
mean_upper_leaflet_binning = curvature_upper_leaflet_binning.results.z_surface

# extract average mean curvature
mean_upper_leaflet_binning = curvature_upper_leaflet_binning.results.mean

# extract average Gaussian curvature
gaussian_upper_leaflet_binning = curvature_upper_leaflet_binning.results.gaussian
```

You can find more examples on how to run MembraneCurvature in the [Usage] page.
To plot results from MembraneCurvature please check the [Visualization] page.

Documentation
---------------

To help you get the most out of MembraneCurvature, we have [documentation]
available where you can find:

- The standard [API] documentation.
- Quick examples of how to run MembraneCurvature in the [Usage] page.
- Detailed explanation of the [Algorithm] implemented in MembraneCurvature.
- Examples on how to plot the results obtained from MembraneCurvature in the [Visualization] page.
- Detailed [Tutorials] to run MembraneCurvature in membrane-only and membrane-protein systems.


Contributing
---------------

Contributions are very welcome!

MembraneCurvature is compatible with [uv] (recommended for development):

```bash
# create an environment and install the project + dev tools
uv sync --extra dev

# add test dependencies and run the test suite
uv sync --extra dev --extra tests
uv run pytest
```

This repository uses [pre-commit] hooks to run quick checks
before commits such as whitespace cleanup, TOML/YAML validation, and [Ruff] linting/formatting.
Using these hooks is highly encouraged because it helps catch common issues early and
keeps pull requests easier to review.

<<<<<<< Updated upstream
To set up hooks locally, with pip:

```bash
pip install -e ".[dev]"
pre-commit install
```

Or with [uv]:
=======
To set up hooks locally, with with [uv]:
>>>>>>> Stashed changes

```bash
uv sync --extra dev
uv run pre-commit install
```

Or with pip:

```bash
pip install -e ".[dev]"
pre-commit install
```


License
---------------

Source code included in this project is available in the GitHub repository
https://github.com/MDAnalysis/membrane-curvature under the GNU General Public
License v3 (see [LICENSE]).

MembraneCurvature was developed as a [Google Summer of Code 2021][GSoC]
project with [MDAnalysis] and it is linked to a [Code of Conduct][code_of_conduct].


[GSoC]: https://summerofcode.withgoogle.com/
[MDAnalysis]: https://www.mdanalysis.org
[NumPy]: https://numpy.org
[SciPy]: https://www.scipy.org
[code_of_conduct]: https://www.mdanalysis.org/conduct/
[Usage]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Usage.html
[License]: https://github.com/MDAnalysis/membrane-curvature/blob/main/LICENSE
[documentation]: https://membrane-curvature.readthedocs.io/en/latest/index.html#
[Visualization]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Visualization.html
[Algorithm]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Algorithm.html
[API]: https://membrane-curvature.readthedocs.io/en/latest/api/membrane_curvature.html
[Tutorials]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Tutorials.html
[MDAnalysisTests]: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
[MDAnalysisData]: https://www.mdanalysis.org/MDAnalysisData/
[UserGuide]: https://userguide.mdanalysis.org/stable/installation.html
[pre-commit]: https://pre-commit.com
[ruff]: https://docs.astral.sh/ruff/
[uv]: https://docs.astral.sh/uv/
