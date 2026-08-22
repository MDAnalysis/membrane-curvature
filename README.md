# MembraneCurvature

[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![GitHub Actions Status](https://github.com/MDAnalysis/membrane-curvature/workflows/CI/badge.svg)](https://github.com/MDAnalysis/membrane-curvature/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/MDAnalysis/membrane-curvature/branch/main/graph/badge.svg)](https://codecov.io/gh/MDAnalysis/membrane-curvature/branch/main)
[![docs](https://readthedocs.org/projects/membrane-curvature/badge/?version=latest)](https://mdanalysismembrane-curvature.readthedocs.io/en/latest/index.html)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![ty](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ty/main/assets/badge/v0.json)](https://github.com/astral-sh/ty)
[![License](https://img.shields.io/github/license/MDAnalysis/membrane-curvature?color=yellow)](https://github.com/MDAnalysis/membrane-curvature/blob/main/LICENSE)

[![Python versions](https://img.shields.io/pypi/pyversions/membrane-curvature)](https://pypi.org/project/membrane-curvature/)
[![PyPI](https://img.shields.io/pypi/v/membrane-curvature)](https://pypi.org/project/membrane-curvature/)
[![Conda version](https://anaconda.org/conda-forge/membrane-curvature/badges/version.svg)](https://anaconda.org/conda-forge/membrane-curvature)

![](https://github.com/MDAnalysis/membrane-curvature/blob/main/docs/_static/mc-logo-whitebg-EBO.png?raw=true)

MembraneCurvature is an [MDAnalysis] MDAKit to calculate membrane curvature from
Molecular Dynamics simulations.

## Features

With MembraneCurvature you can:

- Derive 2D surface profiles from MD simulations using an atom selection as reference with three
  different methods via the `surface_method` parameter: `'fourier'` (default),
  `'binning'`, or `'binning_nearest'`.
- Calculate mean and Gaussian curvature from the derived surface.
- Choose where curvature is evaluated with the `curvature_on` parameter:
  `'per_frame'`, to get per-frame curvature maps, or
  `'average_surface'`, to get curvature from the time-averaged surface.
- Control the surface calculation for the binning methods with two optional parameters:
   - `padding`: periodic edge padding to reduce finite difference artifacts.
   - `fft_filter`: a brick-wall FFT filter to remove high-frequency noise from the averaged surface.

## Installation

MembraneCurvature is available via pip and conda. Please refer to the [Installation] section in the Getting Started Documentation page for
detailed installation instructions.

### With pip

MembraneCurvature is available via `pip`:

```bash
pip install membrane-curvature
```

Or to install from source:

```bash
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
python -m pip install -e .
```

### With conda

MembraneCurvature is available via `conda`:

```bash
conda install -c conda-forge membrane-curvature
```

Or to install from source:

```bash
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
conda env create -f devtools/conda-envs/environment.yaml
conda activate membrane-curvature
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

## Usage

### Fourier method (default)

We can calculate membrane curvature using the Fourier surface method by
omitting `surface_method`, or by setting `surface_method='fourier'`.
In this example, we use the PO4 beads in the upper leaflet as reference:

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# run with all default values
curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4',
                                            ).run()

average_surface = curvature_upper_leaflet.results.average_z_surface

mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
```

To access the per-frame arrays for the example above, use
`results.z_surface[<frame_id>]`, `results.mean[<frame_id>]`, and
`results.gaussian[<frame_id>]`:

```python
# to access the surface for the first frame
surface_first_frame = curvature_upper_leaflet.results.z_surface[0]

# access the mean curvature for the last frame
mean_last_frame = curvature_upper_leaflet.results.mean[-1]

# access the Gaussian curvature for the frame 10
gaussian_frame_10 = curvature_upper_leaflet.results.gaussian[10]
```

The sections below cover the binning surface methods and optional parameters.
Expand a section to view a usage example code.

### Binning methods

Binning methods calculate curvature using finite differences on the derived
surfaces. Two binning methods are available: `'binning'` and
`'binning_nearest'`.

<details>
<summary><code>surface_method='binning'</code></summary>

Set `surface_method='binning'` to run membrane curvature with raw binning by
providing the values for the grid `n_x_bins` and `n_y_bins`.

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# run with binning and coordinate wrapping
curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4',
                                            surface_method='binning',
                                            n_x_bins=8,
                                            n_y_bins=8,
                                            wrap=True,
                                            ).run()

mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
```

</details>

<details>
<summary><code>surface_method='binning_nearest'</code></summary>

Set `surface_method='binning_nearest'`. Wrapping is not valid with this method.
Omit the `wrap` parameter, or set it to `False`.

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# run with binning nearest and no coordinate wrapping
curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4',
                                            surface_method='binning_nearest',
                                            n_x_bins=8,
                                            n_y_bins=8,
                                            ).run()

mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
```

</details>

### Binning methods with padding

Padding is optional and is set to `False` by default. To apply edge periodic
padding, run with `padding=True`.

> [!WARNING]
> Padding is only available with `surface_method='binning'` or
> `surface_method='binning_nearest'` and **requires orthorhombic boxes**.

<details>
<summary><code>surface_method='binning_nearest'</code> with padding</summary>

The example below uses
`surface_method='binning_nearest'`. The same option works with
`surface_method='binning'`.

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# run with binning nearest and padding
curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4',
                                            surface_method='binning_nearest',
                                            n_x_bins=8,
                                            n_y_bins=8,
                                            padding=True,
                                            ).run()

mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
```

</details>

### Binning methods with FFT filtering

FFT filtering is available for the binning methods
`surface_method='binning'` and `surface_method='binning_nearest'` via the
parameter `fft_filter`.

FFT filtering is disabled by default (`fft_filter=None`). To enable automatic
filtering, pass `fft_filter='auto'`.

> [!NOTE]
>
> Per-frame `results.z_surface`, `results.mean`, and `results.gaussian` are not
> FFT-filtered. With filtering enabled, `results.average_z_surface` contains the
> filtered average height. `results.average_mean` and `results.average_gaussian`
> are calculated from that filtered height.

<details>
<summary><code>surface_method='binning'</code> with FFT filtering</summary>

The example below uses
`surface_method='binning'`. The same option works with
`surface_method='binning_nearest'`.

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# run with binning and automatic FFT filtering
curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4',
                                            surface_method='binning',
                                            n_x_bins=8,
                                            n_y_bins=8,
                                            wrap=True,
                                            fft_filter='auto',
                                            ).run()

mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
```

</details>

### Average curvature mode (`curvature_on`)

The parameter `curvature_on` controls how `results.average_mean` and
`results.average_gaussian` are obtained.

By default, MembraneCurvature runs with `curvature_on='per_frame'` (equivalent to `curvature_on=None`)and calculates calculates mean and Gaussian curvature in every frame, and the average maps are
the time averages of those per-frame arrays.

To calculate curvature once from the time-averaged surface, set `curvature_on='average_surface'` explicitly:

```python
import MDAnalysis as mda
from membrane_curvature import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

# calculate curvature from the average surface
curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4',
                                            curvature_on='average_surface',
                                            ).run()

mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
```

You can find more examples on how to run MembraneCurvature in the [Usage] page.
To plot results from MembraneCurvature please check the [Visualization] page.

## Documentation

To help you get the most out of MembraneCurvature, we have [documentation]
available where you can find:

- The standard [API] documentation.
- Quick examples of how to run MembraneCurvature in the [Usage] page.
- Detailed explanation of the [Algorithm] implemented in MembraneCurvature.
- Examples on how to plot the results obtained from MembraneCurvature in the [Visualization] page.
- Detailed [Tutorials] to run MembraneCurvature in membrane-only and membrane-protein systems.

## Contributing

Contributions are very welcome!

MembraneCurvature follows the [MDAnalysis AI Policy][AI_policy]. As a scientific software package, contributions should show real understanding of the code and the science. Please read the [Contributing] page before opening a pull request. **Pull requests that are evidently LLM-generated may be closed at the maintainers' discretion.**

If you are interested in contributing to MembraneCurvature, installation of the development dependencies is required.

There are three dependency groups defined in `pyproject.toml`:

- `dev`: development tools (includes `tests` and `docs`).
- `tests`: testing dependencies.
- `docs`: documentation build dependencies (e.g. Sphinx, themes).

Note that by installing the `dev` group, the `tests` and `docs` dependency groups are also installed.

We recommend using [uv] as a development environment for MembraneCurvature. To set up a `uv` dev environment, clone the repository and move to the root directory:

```bash
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
```

To create a new `uv` environment and install the dependencies included in the `dev` group, run:

```bash
uv sync --group dev
uv run pre-commit install
```

By syncing the `dev` group with [uv], all the development dependencies are installed. MembraneCurvature uses [pre-commit] hooks to run quick checks
before commits such as whitespace cleanup, TOML/YAML validation, and [Ruff] linting/formatting.
Using these hooks is highly encouraged because it helps catch common issues early and keeps pull requests easier to review.

To run the hooks manually without committing:

```bash
uv run pre-commit run --all-files
```

For instructions on how to create a development environment with pip, see the [Installing development dependencies with pip] page.

For more information on how to contribute to MembraneCurvature, see the [Contributing] page.

> **Interested in becoming a maintainer?** We welcome your passion and expertise to help shape and grow this open-source project! Please contact estefania@ojeda-e.com for more details.

## License

Source code included in this project is available in the GitHub repository
https://github.com/MDAnalysis/membrane-curvature under the GNU General Public
License v3 (see [LICENSE]).

MembraneCurvature was developed as a [Google Summer of Code 2021][GSoC]
project with [MDAnalysis] and it is linked to a [Code of Conduct][code_of_conduct].


[GSoC]: https://summerofcode.withgoogle.com/
[MDAnalysis]: https://www.mdanalysis.org
[code_of_conduct]: https://www.mdanalysis.org/conduct/
[Installation]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/getting_started.html#installation
[Installing development dependencies with pip]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/getting_started.html#installing-development-dependencies-with-pip
[Usage]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/source/pages/Usage.html
[License]: https://github.com/MDAnalysis/membrane-curvature/blob/main/LICENSE
[documentation]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/index.html#
[Visualization]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/source/pages/Visualization.html
[Algorithm]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/source/pages/Algorithm.html
[API]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/api/membrane_curvature.html
[Tutorials]: https://mdanalysismembrane-curvature.readthedocs.io/en/latest/source/pages/Tutorials.html
[MDAnalysisTests]: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
[MDAnalysisData]: https://www.mdanalysis.org/MDAnalysisData/
[pre-commit]: https://pre-commit.com
[Ruff]: https://docs.astral.sh/ruff/
[uv]: https://docs.astral.sh/uv/
[Contributing]: https://github.com/MDAnalysis/membrane-curvature/blob/main/.github/CONTRIBUTING.md
[AI_policy]: https://github.com/MDAnalysis/mdanalysis/blob/develop/AI_POLICY.md
