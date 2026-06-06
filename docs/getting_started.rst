Getting Started
===============

**MembraneCurvature** is an MDAKit to calculate membrane curvature from
molecular dynamics simulations. This tool enables the user to derive surfaces from atoms of 
reference (contained in an ``AtomGroup``) and calculate its associated mean and Gaussian curvature.


Installation
------------

There are three main ways to install MembraneCurvature:

:ref:`via_pip`

:ref:`via_conda`

:ref:`with_uv`

.. _via_pip:

1. Via pip
^^^^^^^^^^

The following command will install or upgrade the latest stable version of MembraneCurvature with the core dependencies
(`MDAnalysis`_ and `NumPy`_):

.. code-block:: bash

   pip install membrane-curvature

Some of the examples included in the MembraneCurvature documentation use test
cases from `MDAnalysisTests`_ or `MDAnalysisData`_. To install the unit tests via pip:

.. code-block:: bash

   pip install --upgrade MDAnalysisTests MDAnalysisData


The pip installation of `MembraneCurvature`, `MDAnalysisTests` and `MDAnalysisData` in most of the cases will be enough,
but if you want to install the development dependencies of MembraneCurvature,
you can do so with the following commands:


.. code-block:: bash

   pip install -e .
   pip install --group dev


Note that MembraneCurvature defines development dependency groups in ``pyproject.toml``, which are 
not included in the published package, so installation of these dependencies is required.
These dependency groups are:

- ``dev``: development tools (includes ``tests`` and ``docs`` via ``include-group``).
- ``tests``: testing dependencies.
- ``docs``: documentation build dependencies (e.g. Sphinx, themes).

.. warning::
    The following command requires `pip >= 25.1`

Installation from source is also available with pip by cloning the repository and running the following commands:

.. code-block:: bash

   git clone https://github.com/MDAnalysis/membrane-curvature.git
   cd membrane-curvature
   python -m pip install -e .
   pip install --group dev



.. _via_conda:

2. Via conda
^^^^^^^^^^^^

MembraneCurvature is also available via conda with the ``conda-forge`` channel:

.. code-block:: bash

   conda install -c conda-forge membrane-curvature

To install with conda from source, there is an environment file available to
create a new conda environment with the required dependencies:

.. code-block:: bash

   git clone https://github.com/MDAnalysis/membrane-curvature.git
   cd membrane-curvature
   conda env create -f devtools/conda-envs/environment.yaml
   conda activate membrane-curvature
   python -m pip install -e .

This will create a new environment named ``membrane-curvature`` with the core dependencies
and install MembraneCurvature in development mode.

To install the `MDAnalysisTests`_ and `MDAnalysisData`_ dependencies needed to run some of the examples in the documentation 
via conda, run:

.. code-block:: bash

   conda install -c conda-forge MDAnalysisTests MDAnalysisData

.. _with_uv:

3. With uv (recommended for development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install the latest stable version of MembraneCurvature with `uv`_, run the following command:

.. code-block:: bash

   uv pip install membrane-curvature

For development, install the development dependency group with the following commands:

.. code-block:: bash

   uv sync --group dev



.. _MDAnalysis: https://www.mdanalysis.org
.. _NumPy: https://numpy.org
.. _`github.com/MDAnalysis/membrane_curvature`: https://github.com/MDAnalysis/membrane-curvature
.. _`MDAnalysisTests`: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
.. _`MDAnalysisData`: https://www.mdanalysis.org/MDAnalysisData/
.. _`Installation Quick Start`: https://www.mdanalysis.org/pages/installation_quick_start/#installation-quick-start
.. _`conda`: https://conda.io/en/latest/
.. _`uv`: https://docs.astral.sh/uv/