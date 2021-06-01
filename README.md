MDA Membrane Curvature
==============================
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://www.numfocus.org/)
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![Build Status](https://travis-ci.com/MDAnalysis/membrane-curvature.svg?branch=main)](https://travis-ci.com/MDAnalysis/membrane-curvature)
[![codecov](https://github.com/MDAnalysis/membrane-curvature/actions/workflows/CI.yaml/badge.svg)](https://codecov.io/gh/MDAnalysis/membrane-curvature)
[![docs](https://img.shields.io/badge/docs-unknown-lightgray)](https://www.mdanalysis.org/membrane-curvature/)

This is an [MDAnalysis][mdanalysis] module to calculate membrane curvature from 
molecular dynamics simulations. 

The MDAkit for membrane curvature analysis is part of the [Google Summer of Code][GSoC] program 
and it is linked to a [Code of Conduct][code_of_conduct].

## Abstract
Elements of differential geometry enable us to quantify the curvature of a surface. 
The core elements of biological membranes, phospholipids, provide tridimensional 
configurations from which a surface can be derived to calculate curvature descriptors. 
We would like to integrate to [MDAnalysis][mdanalysis] an analysis module to calculate 
average mean and Gaussian curvature from Molecular Dynamics simulations. By integrating 
a membrane curvature analysis module in [MDAnalysis][mdanalysis], users will benefit from 
a tool that enables rapid extraction of relevant properties of lipid bilayers in 
biomolecular systems. Our approach extracts key elements of the membrane using 
phospholipid head groups to define a surface, followed by a transformation into 
[NumPy] arrays. In this way, the functionality offered by [NumPy][numpy] and 
[SciPy] can be used to derive values of mean and gaussian curvature of biological 
membranes. Since [MDAnalysis][mdanalysis] already works very well to explore data 
interactively, visualization of the membrane curvature analysis outputs in 2D-maps can be
easily performed.

[GSoC]: https://summerofcode.withgoogle.com/
[mdanalysis]: https://www.mdanalysis.org
[NumPy]: https://numpy.org
[SciPy]: https://www.scipy.org
[code_of_conduct]: https://www.mdanalysis.org/pages/conduct/
