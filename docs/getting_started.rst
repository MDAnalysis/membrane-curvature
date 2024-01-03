Getting Started
===============

The `MDAnalysis`_ **MembraneCurvature** is a Python toolkit to calculate membrane curvature from
molecular dynamics simulations. This tool enables the user to derive surfaces from `atoms` of 
reference (contained in an `AtomGroup`) and calculate its associated mean and Gaussian curvature.

.. _MDAnalysis: https://www.mdanalysis.org


Installation
--------------

The main dependency in MembraneCurvature is `MDAnalysis`_. You can find
instructions to install the latest stable version of MDAnalysis via 
`conda`_ in the `Installation Quick Start`_ page.

MembraneCurvature is available via pip:

.. code-block:: bash

   pip install membrane-curvature

To install from source:

.. code-block:: bash

   git clone https://github.com/MDAnalysis/membrane-curvature.git
   cd membrane-curvature
   conda env create -f devtools/conda-envs/environment.yaml
   conda activate membrane-curvature
   python setup.py install

Some of the examples included in the MembraneCurvature documentation use test
cases from `MDAnalysis Tests`_. To install the unit tests via conda:

.. code-block:: bash

   conda install -c conda-forge MDAnalysisTests

or via pip:

.. code-block:: bash

   pip install --upgrade MDAnalysisTests



.. _MDAnalysis: https://www.mdanalysis.org
.. _`github.com/MDAnalysis/membrane_curvature`: https://github.com/MDAnalysis/membrane-curvature
.. _`MDAnalysis Tests`: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests
.. _`Installation Quick Start`: https://www.mdanalysis.org/pages/installation_quick_start/#installation-quick-start
.. _`conda`: https://conda.io/en/latest/
.. |MDAnalysis_version| replace:: 1.1.1