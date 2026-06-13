Getting Started
===============

**MembraneCurvature** is an MDAKit to calculate membrane curvature from
molecular dynamics simulations. This tool enables the user to derive surfaces from atoms of 
reference (contained in an ``AtomGroup``) and calculate its associated mean and Gaussian curvature.


Installation
------------

There are three main ways to install MembraneCurvature:

- :ref:`with_pip`

- :ref:`with_conda`

- :ref:`with_uv`

.. _with_pip:

With pip
^^^^^^^^

The following command will install or upgrade the latest stable version of MembraneCurvature with the core dependencies
(`MDAnalysis`_ and `NumPy`_):

.. code-block:: bash

   pip install membrane-curvature

Some of the examples included in the MembraneCurvature documentation use test
cases from `MDAnalysisTests`_ or `MDAnalysisData`_. To install these dependencies with pip, run:

.. code-block:: bash

   pip install --upgrade MDAnalysisTests MDAnalysisData


Installing MembraneCurvature from source
+++++++++++++++++++++++++++++++++++++++++

Installation from source is also available with pip by cloning the repository and installing the package
in editable mode:

.. code-block:: bash

   git clone https://github.com/MDAnalysis/membrane-curvature.git
   cd membrane-curvature
   python -m pip install -e .
   pip install --group dev

.. _installing_development_dependencies_with_pip:

Installing development dependencies with pip
+++++++++++++++++++++++++++++++++++++++++++++

For most users, installing MembraneCurvature (optionally with `MDAnalysisTests`_ and `MDAnalysisData`_) is sufficient.

If you plan to develop MembraneCurvature, you should install the package in editable mode along with its development
dependencies. Note that these dependencies are not included in the published package MembraneCurvature, so installation
of these dependencies is required **only for development**. 

Development dependencies are defined in ``pyproject.toml``, and includes the following dependency groups:

- ``dev``: development tools (includes ``tests`` and ``docs`` via ``include-group``).
- ``tests``: testing dependencies.
- ``docs``: documentation build dependencies (e.g. Sphinx, themes).

.. note::
   
   These dependencies are **only required for development** and not needed for users of MembraneCurvature.

The development dependencies can be installed with:

.. code-block:: bash

   pip install -e .
   pip install --group dev

.. warning::
    The command ``pip install --group dev`` requires ``pip >= 25.1``


.. _with_conda:

With conda
^^^^^^^^^^

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

With uv (recommended for development)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install the latest stable version of MembraneCurvature with `uv`_, run the following command:

.. code-block:: bash

   uv pip install membrane-curvature

Installing development dependencies with uv
++++++++++++++++++++++++++++++++++++++++++++

For development, install the development dependency group with:

.. code-block:: bash

   uv sync --group dev

This command creates a virtual environment, installs MembraneCurvature in editable mode, and resolves all ``dev``
dependency groups (including ``tests`` and ``docs``) in one step.

.. note::
   **Why is uv recommended for development?**

   We recommend working with a `uv`_ environment because it is the simplest way to set up a reproducible and
   isolated environment for the development of MembraneCurvature. Development dependencies can be installed
   with a single command ``uv sync --group dev``.

   Additionally, in comparison to :ref:`installing_development_dependencies_with_pip` there is no need to install
   MembraneCurvature in editable mode.
    

Installation guides :ref:`with_pip` and :ref:`with_conda` are also available if you prefer them.

.. _MDAnalysis: https://www.mdanalysis.org
.. _NumPy: https://numpy.org
.. _`github.com/MDAnalysis/membrane_curvature`: https://github.com/MDAnalysis/membrane-curvature
.. _`MDAnalysisTests`: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
.. _`MDAnalysisData`: https://www.mdanalysis.org/MDAnalysisData/
.. _`Installation Quick Start`: https://www.mdanalysis.org/pages/installation_quick_start/#installation-quick-start
.. _`conda`: https://conda.io/en/latest/
.. _`uv`: https://docs.astral.sh/uv/