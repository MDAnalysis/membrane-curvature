Membrane Curvature
==============================
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![GitHub Actions Status](https://github.com/MDAnalysis/membrane-curvature/workflows/CI/badge.svg)](https://github.com/MDAnalysis/membrane-curvature/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/MDAnalysis/membrane-curvature/branch/main/graph/badge.svg)](https://codecov.io/gh/MDAnalysis/membrane-curvature/branch/main)
[![docs](https://readthedocs.org/projects/membrane-curvature/badge/?version=latest)](https://membrane-curvature.readthedocs.io/en/latest/)
![PyPI](https://img.shields.io/pypi/v/membrane-curvature?color=lightgray)

![](https://github.com/MDAnalysis/membrane-curvature/blob/main/docs/source/_static/PM_Membrane_EBO.png?raw=true)

MembraneCurvature is an [MDAnalysis] tool to calculate membrane curvature from 
Molecular Dynamics simulations. 

> **Interested in becoming a maintainer?** We welcome your passion and expertise to help shape and grow this open-source project! Please contact estefania@ojeda-e.com for more details.

Features
--------------

With MembraneCurvature you can:

- Calculate mean and Gaussian curvature from MD simulations.
- Derive 2D curvature profiles.
- Live a happier life.


Installation
--------------

The main dependency in MembraneCurvature is [MDAnalysis]. You can find
instructions to install the latest stable version of MDAnalysis via `conda` in the [UserGuide].

MembraneCurvature is available via `pip`:

```
pip install membrane-curvature
```

To install from source:

```
git clone https://github.com/MDAnalysis/membrane-curvature.git
cd membrane-curvature
conda env create -f devtools/conda-envs/environment.yaml
conda activate membrane-curvature
python setup.py install
```

Some of the examples included in the MembraneCurvature documentation use test
cases from [MDAnalysisTests]. To install the unit tests via `conda`:

```
conda install -c conda-forge MDAnalysisTests
```

or via `pip`:

```
pip install --upgrade MDAnalysisTests
```

> ⚠️ In comparison to the previous version, `membrane-curvature==0.0.3` shows a significant improvement in performance, particularly notable for membrane-protein systems. Installing the last available version is highly encouraged. 

Usage
--------------

This is a quick example on how to run MembraneCurvature:

```Python
import MDAnalysis as mda
from membrane_curvature.base import MembraneCurvature
from MDAnalysis.tests.datafiles import Martini_membrane_gro

universe = mda.Universe(Martini_membrane_gro)

curvature_upper_leaflet = MembraneCurvature(universe,
                                            select='resid 1-225 and name PO4',
                                            n_x_bins=8,
                                            n_y_bins=8,
                                            wrap=True).run()

# extract mean curvature
mean_upper_leaflet = curvature_upper_leaflet.results.z_surface

# extract mean curvature
mean_upper_leaflet = curvature_upper_leaflet.results.mean

# extract Gaussian
gaussian_upper_leaflet = curvature_upper_leaflet.results.gaussian
```

In this example, we use the PO4 beads in the upper leaflet as reference to
derive a surface and calculate its respective mean and Gaussian curvature.

You can find more examples on how to run MembraneCurvature in the [Usage] page.
To plot results from MembraneCurvature please check the [Visualization] page.

Documentation
---------------

To help you get the most out MembraneCurvature, we have [documentation] available 
where you can find:

- The standard [API] documentation.
- Quick examples of how to run Membrane Curvature in the [Usage] page.
- Detailed explanation of the [Algorithm] implemented in MembraneCurvature.
- Examples on how to plot the results obtained from MembraneCurvature in the [Visualization] page.


License
---------------

Source code included in this project is available in the GitHub repository
https://github.com/MDAnalysis/membrane-curvature under the GNU Public License
v3 , version 3 (see [LICENSE]).

MembraneCurvature was developed as a [Google Summer of Code 2021][GSoC] 
project with [MDAnalysis] and it is linked to a [Code of Conduct][code_of_conduct].

Copyright (c) 2021-2022, Estefania Barreto-Ojeda


[GSoC]: https://summerofcode.withgoogle.com/
[MDAnalysis]: https://www.mdanalysis.org
[NumPy]: https://numpy.org
[SciPy]: https://www.scipy.org
[code_of_conduct]: https://www.mdanalysis.org/pages/conduct/
[Usage]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Usage.html
[License]: https://github.com/MDAnalysis/membrane-curvature/blob/main/LICENSE
[documentation]: https://membrane-curvature.readthedocs.io/en/latest/index.html#
[Visualization]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Visualization.html
[Algorithm]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Algorithm.html
[API]: https://membrane-curvature.readthedocs.io/en/latest/api/membrane_curvature.html
[MDAnalysisTests]: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
[UserGuide]: https://userguide.mdanalysis.org/2.0.0-dev0/installation.html