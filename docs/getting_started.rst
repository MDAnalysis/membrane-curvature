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


MembraneCurvature has optional dependency groups that you can install with pip to enable additional features:

- ``dev``: optional dependencies for development.
- ``tests``: optional dependencies for testing.
- ``docs``: install documentation build dependencies (e.g. Sphinx, themes).

To install all optional dependencies, run:

.. code-block:: bash

   pip install -e ".[dev,tests,docs]"


Installation from source is also available with pip by cloning the repository and running:

.. code-block:: bash

   git clone https://github.com/MDAnalysis/membrane-curvature.git
   cd membrane-curvature
   python -m pip install -e .

Note that some of the examples included in the MembraneCurvature documentation use test
cases from `MDAnalysisTests`_ or `MDAnalysisData`_. To install the unit tests via conda:

.. code-block:: bash

   pip install --upgrade MDAnalysisTests MDAnalysisData


.. _via_conda:

2. Via conda
^^^^^^^^^^^^

MembraneCurvature is available via conda with the ``conda-forge`` channel:

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

For development, you can install the optional dependencies with the following commands:

.. code-block:: bash

   uv sync --extra dev # install only the development dependencies
   uv sync --extra tests --extra docs # install the testing and documentation dependencies

or to install all optional dependencies with uv:

.. code-block:: bash

   uv sync --all-extras


.. _MDAnalysis: https://www.mdanalysis.org
.. _NumPy: https://numpy.org
.. _`github.com/MDAnalysis/membrane_curvature`: https://github.com/MDAnalysis/membrane-curvature
.. _`MDAnalysisTests`: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
.. _`MDAnalysisData`: https://www.mdanalysis.org/MDAnalysisData/
.. _`Installation Quick Start`: https://www.mdanalysis.org/pages/installation_quick_start/#installation-quick-start
.. _`conda`: https://conda.io/en/latest/
.. _`uv`: https://docs.astral.sh/uv/