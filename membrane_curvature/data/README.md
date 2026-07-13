# MembraneCurvature Package Data

This directory contains GROMACS structure (`.gro`) and trajectory (`.xtc`) files
shipped with MembraneCurvature for unit tests and documentation examples.

## Accessing the data

Paths are exposed in `membrane_curvature.tests.datafiles`.

For example, to load the test data `GRO_PO4_SMALL`, `MEMB_GRO`, and `MEMB_XTC`, we do:

```python
from membrane_curvature.tests.datafiles import GRO_PO4_SMALL, MEMB_GRO, MEMB_XTC

import MDAnalysis as mda

universe = mda.Universe(MEMB_GRO, MEMB_XTC)
```

## Files

File names and test constants are defined in [`membrane_curvature/tests/datafiles.py`](../tests/datafiles.py).

- `test_po4_small.{gro,xtc}`: small POPC/POPE membrane (10 lipids), PO4 beads.
- `test_curvature_po4_only.{gro,xtc}`: POPC/POPE membrane (914 lipids), PO4 beads.
- `MEMB_traj_short.{gro,xtc}`: large membrane patch in a squared cell (~240 Å).
- `test_curvature_abca1.{gro,xtc}`: POPC/POPE/CHOL membrane–protein system.
- `Membrane_protein_fit.{gro,xtc}`: membrane–protein trajectory with rotational and translational fit.
- `test_po4_inverted_indexes.gro`: small PO4 membrane structure (not used in the test suite).

These files are included in the package via `MANIFEST.in` with the ``graft membrane_curvature`` line.
