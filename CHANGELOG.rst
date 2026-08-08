Membranecurvature CHANGELOG
=============================

---

* Update Algorithm doc page with three surface methods and optional padding/FFT filtering (PR #266)
* Calculate analytic derivatives for curvature calculations in Fourier surface method for
  ``curvature_on='average_surface'`` (PR #272)
* Improve API documentation: MembraneCurvature class docstring (PR #258), move optional params to subsection
  in surface methods and simplify Fourier docstrings (PR #241)
* Expose and document ``nyquist_q`` to set manual FFT filter bounds (PR #241)
* Enabled ``curvature_on`` kwarg to calculate curvature as ``per-frame`` average or as
  ``'average_surface'`` (PR #250)
* Added ``surface_method='binning_nearest'`` as a new surface method (PR #246)
* Feature: Added periodic edge padding for the binning surface method (PR #238)
* Refactored binning surface method to use ``np.add.at`` for improved performance (PR #240)
* Set ``fft_filter`` to ``None`` as default (PR #239)
* Refact: Calculate Fourier coeffs in SVD with ``np.dot`` instead of matmul (PR #223)
* Remove side effect - MDAnalysis logging at import (PR #221)
* Fixed stale project metadata and updated README files (PR #218)
* Replaced action version with specific commit hash in github actions and added dependabot cfg file (PR #224 #234)
* Include tests in coverage report; set Codecov threshold to 100% (PR #219)
* Added brick-wall FFT filtering to the binning surface method (PR #188)
* Added `uv.lock` file and step in CI to validate uv environment with dev dependency group (PR #202)
* Improve API documentation: layout, content, and navigation menu (PR #193)
* Modernized project metadata: comply with PEP 639 (PR #185) and PEP 735 (PR #194)
* Updated installation instructions (conda, pip, uv) and conda environment files (PR #186)
* Refactored Fourier surface method to use singular-value decomposition (SVD) instead of least-squares (PR #177)
* Apply coordinate wrapping only for ``surface_method='binning'``, set ``wrap`` to ``None`` by default,
  raise ValueError when ``surface_method='fourier'`` (PR #170).
* Set `'fourier'` as default surface method. (PR #168)
* Updated documentation to reflect Fourier surface method (Algorithm page PR #160, Usage PR #163, Tutorials PR #164 #173, Visualization #171)
* Added surface_method='fourier'; derives surface from atoms of reference using a truncated Fourier series (PR #146)
* Added maintenance tools (ruff PR #147, pre-commit hooks PR #151, ty PR #155)

1.1.2 (24-04-2026)
------------------
* Fixed deployment workflow with TestPyPI and PyPI (PR #137)
* Dropped support for Python 3.9, 3.10; extended support to 3.13, 3.14 (PR #141)
* Refactored curvature calculations with Monge gauge helpers (PR #140)

1.1.1 (09-01-2024)
------------------
* Add automated deployment workflow (PR #124)

1.1.0 (04-01-2024)
------------------
* Switch to mdanalysis-sphinx-theme and update docs (PR #122).
* Comply with PEP518 (PR #119).
* Drop support for Python 3.8 (PR #112).

1.0.0 (03-11-2022)
-------------------

* Dropped support for Python 3.6 and 3.7 (PR #103).
* Extended support to Python 3.10 (PR #94) and 3.11 (PR #103).
* Fixed plots in documentation pages and updated plots in tutorials (PR #89).
* Fixed warning messages that affected performance (PR #95).


0.0.2 (31-10-2021)
-------------------

* Fixed bug sign mean curvature (PR #75).
* Updated tutorial membrane only system (PR #78).
* Added ipynb tutorial membrane-protein (all-atom). (PR #69)


0.0.1 (21-09-2021)
-------------------

* First release on PyPI.
