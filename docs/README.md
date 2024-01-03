# Compiling Membrane Curvature's Documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/) using the [mdanalysis-sphinx-theme](https://github.com/MDAnalysis/mdanalysis-sphinx-theme).
To compile the docs, first ensure that:

1. All necessary dependencies are installed.
1. The editable version of this repo is installed.
1. You are working in the `/docs` directory.

## Installing dependencies

### Via conda

To install a new environment with conda to compile the documentation, run:

```bash
conda install sphinx mdanalysis-sphinx-theme 
```

### Via pip

Alternatively, to install all doc dependencies in your environment run:

```bash
pip install -r requirements.txt
```

## Installing the editable version

Please souble check you have installed the editabe version of the package. **From the root of the project**, run:

```bash
pip install -e .
```

## Compiling MembraneCurvature Sphinx docs

Once all dependencies are installed, you can use the `Makefile` in this directory to compile static HTML pages with

```bash
make html
```

The compiled docs will be in the `_build` directory and can be viewed by opening `index.html` (which may itself be inside a directory called `html/` depending on what version of Sphinx is installed).

A configuration file for [Read The Docs](https://readthedocs.org/) (readthedocs.yaml) is included in the top level of the repository. To use Read the Docs to host your documentation, go to https://readthedocs.org/ and connect this repository. You may need to change your default branch to `main` under Advanced Settings for the project.
