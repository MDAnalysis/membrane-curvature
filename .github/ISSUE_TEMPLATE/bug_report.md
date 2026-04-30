---
name: Bug report
about: Report an error or unexpected behavior
title: "[Bug]: "
labels: ["bug"]
---

<!--
Thanks for taking the time to file a bug report.
Please fill out the sections below as completely as you can.
-->

## Actual behavior

<!-- Describe what is actually happening. Include full error output/stack traces. -->

## Expected behavior

<!-- Describe what you expected to happen. -->

## Steps / code to reproduce

<!--
Provide minimal steps or code that reproduces the issue. If possible, use MDAnalysis test data files.
Replace `<your_atoms_of_reference>` with the atoms of reference you are using.
-->

```python
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import GRO, XTC

# from membrane_curvature import MembraneCurvature  # adjust import if different

u = mda.Universe(GRO, XTC)

curvature_upper_leaflet = MembraneCurvature(
    u,
    select="<your_atoms_of_reference>",
    n_x_bins=<n>,
    n_y_bins=<n>,
    wrap=True,
).run()
```

## Additional context

<!-- Anything else that might help (input data details, constraints, etc.). -->

## Version information

- membrane-curvature: <!-- e.g. `0.3.1` or commit `5dfaedd` -->
- MDAnalysis: <!-- `python -c "import MDAnalysis as mda; print(mda.__version__)"` -->
- Python: <!-- `python -V` -->
- OS: <!-- Select one: macOS / Windows / Linux -->
