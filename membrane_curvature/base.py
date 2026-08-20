"""
MembraneCurvature
=======================================

:Author: Estefania Barreto-Ojeda
:Year: 2021
:Copyright: GNU Public License v3

:class:`~membrane_curvature.base.MembraneCurvature` is the main analysis class
for calculating mean and Gaussian curvature from Molecular Dynamics (MD) simulations.

Users must provide at least two inputs: an MDAnalysis :class:`~MDAnalysis.core.universe.Universe`,
and an :class:`~MDAnalysis.core.groups.AtomGroup` defining the reference atoms used to derive a height surface
from which curvature is calculated.

The height surface is derived from the selected reference atoms using a method specified by
:attr:`~membrane_curvature.base.MembraneCurvature.surface_method`. This can be either the ``'fourier'`` method
(default) or one of the binning methods (``'binning'`` or ``'binning_nearest'``).

Depending on the selected ``surface_method`` more parameters are required:

- **Fourier method**:

    Uses ``fourier_m`` and ``fourier_n`` to define the number of Fourier modes. Fourier modes default to 2.

    - Optional (advanced) ``fourier_rcond`` uses a singular-value cutoff for the Fourier fit.


- **Binning methods**:

    Required parameters are the number of bins in the $x$ and $y$ directions ``n_x_bins`` and ``n_y_bins``.

    - Optional ``padding`` applies a periodic buffer before finite differences to avoid edge artifacts.
      Available for orthorhombic boxes only.
    - Optional ``wrap`` to control whether to wrap the coordinates exceeding the simulation box dimensions.
    - Optional (adavanced) ``fft_filter`` applies a brick-wall FFT filter to reduce noisy signals.
      Manual and auto modes are available.

    This applies for the two binning methods: ``'binning'`` and ``'binning_nearest'``.

Once the height surface is derived, mean and Gaussian curvature are computed according to the mode specified by
:attr:`~membrane_curvature.base.MembraneCurvature.curvature_on`. With ``'per_frame'``, curvature is calculated
as the mean of the per-frame curvature values. With ``'average_surface'``, curvature is computed directly
from the average height surface.

Mean curvature is calculated in units of Å :sup:`-1` and Gaussian curvature in units of Å :sup:`-2`.
"""

import numpy as np
import warnings
from .binning_surface import get_z_surface
from .binning_nearest_surface import get_z_surface_nearest
from .curvature import (
    mean_curvature,
    gaussian_curvature,
    fourier_curvature,
    fourier_curvature_from_theta,
    curvature_from_primary_with_edge_pad,
    curvature_with_edge_pad,
)
from .padding import (
    get_z_surface_padded,
    padded_grid_spec,
)
from .padding_validators import validate_edge_pad_bins, validate_orthorhombic
from .fft_filtering import apply_fft_filter, resolve_fft_filter
from .fourier_surface import (
    n_fourier_parameters,
    fourier_mode_list,
)
from .fourier_validators import validate_positive_bin_counts

from MDAnalysis.analysis.base import AnalysisBase

import logging
from enum import StrEnum

logger = logging.getLogger('MDAnalysis.MDAKit.membrane_curvature')


class MembraneCurvature(AnalysisBase):
    r"""

    Parameters
    ----------
    universe : Universe or AtomGroup
        An MDAnalysis Universe object.
    select : str or iterable of str, optional.
        The selection string of an atom selection to use as a
        reference to derive a surface.
    surface_method : {'binning', 'binning_nearest', 'fourier'}, optional, [``'fourier'``]
        Method to derive the height surface.
        See :class:`~membrane_curvature.base.MembraneCurvature.SurfaceMethod`
        for more details.
    fourier_m : int, optional, [``2``]
        Maximum Fourier mode index in $x$ when ``surface_method='fourier'``.
    fourier_n : int, optional, [``2``]
        Maximum Fourier mode index in $y$ when ``surface_method='fourier'``.
    wrap : bool or None, optional, [``True``]
        Wrap $x$ and $y$ of the atom selection into the unit cell.
    n_x_bins : int, optional, [``100``]
        Number of bins in grid in the $x$ dimension.
    n_y_bins : int, optional, [``100``]
        Number of bins in grid in the $y$ dimension.
    curvature_on : {None, 'per_frame', 'average_surface'}, optional, [``None``]
        Mode to calculate curvature on. ``None`` runs as 'per_frame'.
        See :class:`~membrane_curvature.base.MembraneCurvature.CurvatureOn` for more details.
    padding : bool, optional, [``False``]
        Apply periodic edge padding for binning surface methods.
    edge_pad_bins : int, optional, [``2``]
        Buffer width in bins on each side when ``padding=True``.
    fft_filter : None, str, or dict, optional, [``None``]
        Mode to apply brick-wall filter on.
    fourier_rcond : float, optional, [``None``]
        Singular-value cutoff for the Fourier fit via truncated SVD.
    grid_origin : {'box', 'lipid_bbox'}, optional, [``'box'``]
        Grid origin for binning nearest surface method.
    x_range : tuple of (float, float), optional, [``(0, universe.dimensions[0])``]
        Range of coordinates (min, max) in the $x$ dimension.
    y_range : tuple of (float, float), optional, [``(0, `universe.dimensions[1]`)``]
        Range of coordinates (min, max) in the $y$ dimension.

    Attributes
    ----------
    results.z_surface : ndarray
        Surface derived from the atom selection in each frame.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`).
    results.mean : ndarray
        Mean curvature associated with the surface. Derived on a per-frame basis.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.gaussian : ndarray
        Gaussian curvature associated with the surface. Derived on a per-frame basis.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.average_z_surface : ndarray
        Average of the array elements in `z_surface`.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_mean : ndarray
        Average of the array elements in `mean_curvature`.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_gaussian: ndarray
        Average of the array elements in `gaussian_curvature`.
        Each array has shape (`n_x_bins`, `n_y_bins`)

    Raises
    ------
    ValueError
        If ``n_x_bins`` or ``n_y_bins`` is not a positive integer.
        If the atom selection is empty.
        If ``surface_method`` is not one of ``'binning'``, ``'binning_nearest'``,
        or ``'fourier'``.
        If ``curvature_on`` is not ``None``, ``'per_frame'``, or ``'average_surface'``.
        If ``surface_method='fourier'``, when: ``fourier_m``, ``fourier_n`` are negative
        or the selection has fewer atoms than Fourier parameters, ``wrap=True``, or
        ``fft_filter`` is passed.
        If ``padding=True`` with ``surface_method='fourier'``, or passed with a
        non-orthorhombic box, or with ``edge_pad_bins`` when it is not an integer ``>= 1``.

    See also
    --------
    :class:`~MembraneCurvature.padding`
        Apply periodic edge padding for binning surface methods.

    :class:`~MDAnalysis.transformations.wrap.wrap`
        Wrap/unwrap the atoms of a given AtomGroup in the unit cell.

    :class:`~MembraneCurvature.fft_filtering`
        Apply a brick-wall filter to the arrays.

    Notes
    -----

    **Fourier surface method (default)**

    ``surface_method='fourier'`` uses ``fourier_m = fourier_n = 2`` as default.
    Since the method performs periodic boundary conditions by itself, ``wrap`` defaults
    to ``False`` and it is not required. Curvature is calculated using a Fourier sum fit to the
    atom positions and the derivatives are computed analytically.

    **Binning methods**

    MembraneCurvature  supports two different binning methods: ``'binning'`` and ``'binning_nearest'``.

    .. note::

        The two binning methods differ in how they assign the $z$ coordinate to the bins:

        ``'binning'`` assigns atoms to bins and stores the mean $z$ coordinate at each
        bin center. Empty bins are ``NaN``.

        ``'binning_nearest'`` assigns each bin corner the $z$ coordinate of the
        nearest lipid in the $xy$ plane. For the same bin counts, its maps are offset
        by half a bin relative to ``'binning'``.

    Both methods run with ``wrap`` set to ``True`` by default. Curvature is calculated using finite
    differences, and the derivatives are computed via numerical differentiation.

    Two optional functionalities are supported by both binning methods:

    - **Padding** (``padding``, ``edge_pad_bins``): Only available for orthorhombic boxes. Creates a padded height
      field built from periodic atom images, then the buffer of ``edge_pad_bins`` on each side is clipped so stored
      arrays keep the primary grid shape.

    - **Brick-wall filter** (``fft_filter``): Applies a brick-wall filter to the height field to reduce noisy signals.
      Available in manual and auto modes. The auto mode uses a heuristic to determine the filter parameters based on the
      box dimensions and the number of bins. The manual mode allows the user to specify the filter values.


    Example
    -----------
    You can pass a universe containing your selection of reference::

        import MDAnalysis as mda
        from membrane_curvature.base import MembraneCurvature

        u = mda.Universe(coordinates, trajectory)
        mc = MembraneCurvature(u, select='<atomgroup>').run()

        surface =  mc.results.average_z_surface
        mean_curvature =  mc.results.average_mean
        gaussian_curvature = mc.results.average_gaussian

    **MembraneCurvature.results**

    For any method of choice, the derived surface and calculated curvatures are available
    in the :attr:`results` attributes. The attributes :attr:`~MembraneCurvature.results.average_mean` and
    :attr:`~MembraneCurvature.results.average_gaussian` contain mean and
    Gaussian curvature outputs.

    .. warning::

        With ``curvature_on='per_frame'``, the attributes ``results.average_mean`` and
        ``results.average_gaussian`` contain averages of curvature calculated in every frame.

        With ``curvature_on='average_surface'``,  the attributes ``results.average_mean`` and
        ``results.average_gaussian`` contain the mean and Gaussian curvatures calculated from
        the average surface.

    The respective 2D curvature plots can be obtained using the :mod:`matplotlib`
    package for data visualization via :func:`~matplotlib.pyplot.contourf` or
    :func:`~matplotlib.pyplot.imshow`.

    For specific examples visit the :ref:`usage` page.
    Check the :ref:`visualization` page for examples on how to plot
    MembraneCurvature results.

    """

    class SurfaceMethod(StrEnum):
        """Supported values for ``surface_method``.

        - ``'binning'``:
            Derives the surface by creating a grid and assigning atoms to bins.
            Calculates curvature using finite differences.
            See :mod:`~membrane_curvature.binning_surface` for more details.


        - ``'binning_nearest'``:
            Assigns each grid corner the $z$ coordinate of the nearest lipid in a cell.
            Calculates curvature using finite differences.
            See :mod:`~membrane_curvature.binning_nearest_surface` for more details.

        - ``'fourier'``:
            Fits a periodic Fourier sum to atom positions at each frame and builds
            a design matrix for the Fourier coefficients. Calculates curvature using
            analytic derivatives.
            See :mod:`~membrane_curvature.fourier_surface` for more details.

        """

        BINNING = 'binning'
        BINNING_NEAREST = 'binning_nearest'
        FOURIER = 'fourier'

    class CurvatureOn(StrEnum):
        r"""Supported values for ``curvature_on``. Where to evaluate curvature.

        ``'per_frame'``
            Computes curvature as the mean of the per-frame curvature arrays:
            :math:`\langle H \rangle`, :math:`\langle K \rangle`.
        ``'average_surface'``
            Computes curvature from the average surface: :math:`H(\langle S \rangle)`,
            :math:`K(\langle S \rangle)`.
        ``None`` (default)
            Equivalent to ``'per_frame'``. Keeps behaviour from v1.1.2 and prior.
        """

        PER_FRAME = 'per_frame'
        AVERAGE_SURFACE = 'average_surface'

    def __init__(
        self,
        universe,
        select='all',
        n_x_bins=100,
        n_y_bins=100,
        x_range=None,
        y_range=None,
        wrap=None,
        surface_method='fourier',
        fourier_m=2,
        fourier_n=2,
        fourier_rcond=None,
        fft_filter=None,
        curvature_on=None,
        padding=False,
        edge_pad_bins=2,
        grid_origin='box',
        **kwargs,
    ):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.ag = universe.select_atoms(select)
        if len(self.ag) == 0:
            raise ValueError('Invalid selection. Empty AtomGroup.')
        validate_positive_bin_counts(n_x_bins, n_y_bins)
        self.n_x_bins = n_x_bins
        self.n_y_bins = n_y_bins
        self._x_range_user = x_range is not None
        self._y_range_user = y_range is not None
        self.x_range = x_range if x_range else (0, universe.dimensions[0])
        self.y_range = y_range if y_range else (0, universe.dimensions[1])
        self.dx = (self.x_range[1] - self.x_range[0]) / n_x_bins
        self.dy = (self.y_range[1] - self.y_range[0]) / n_y_bins
        if not isinstance(padding, (bool, np.bool_)):
            raise ValueError(f'padding must be True or False, got {padding!r}')
        self.padding = padding
        if self.padding:
            self.edge_pad_bins = validate_edge_pad_bins(edge_pad_bins, n_x_bins, n_y_bins)
        else:
            self.edge_pad_bins = None
        self._pad_spec = None

        try:
            self.surface_method = self.SurfaceMethod(surface_method)
        except ValueError as err:
            allowed = ', '.join(repr(method.value) for method in self.SurfaceMethod)
            raise ValueError(f'surface_method must be one of ({allowed}), got {surface_method!r}') from err
        self.grid_origin = grid_origin
        if grid_origin not in {'box', 'lipid_bbox'}:
            raise ValueError("grid_origin must be 'box' or 'lipid_bbox'")
        self.fourier_m = int(fourier_m)
        self.fourier_n = int(fourier_n)
        self.fourier_rcond = fourier_rcond

        if self.surface_method == self.SurfaceMethod.FOURIER:
            if wrap is True:
                raise ValueError("wrap=True is only valid when surface_method='binning'")
            if fft_filter is not None:
                raise ValueError("fft_filter is only valid for surface_method='binning' or 'binning_nearest'")
            if self.padding:
                raise ValueError(
                    "padding=True is only valid for surface_method='binning' or surface_method='binning_nearest'"
                )
            self.wrap = False
            self.fft_filter = None
            self._fft_q_bounds = None
            if self.fourier_m < 0 or self.fourier_n < 0:
                raise ValueError('fourier_m and fourier_n must be non-negative')
            n_param = n_fourier_parameters(self.fourier_m, self.fourier_n)
            n_atom = len(self.ag)
            if n_atom < n_param:
                max_square_mode = max(int((np.sqrt(n_atom) - 1) // 2), 0)
                # If too few atoms for the passed modes fourier_m / fourier_n, raise error.
                # The message suggests the maximum fourier_m / fourier_n given the number
                # of atoms in the selection n_atom = len(self.ag)
                raise ValueError(
                    f"surface_method='fourier' needs at least {n_param} atoms in the selection "
                    f'(fourier_m={self.fourier_m}, fourier_n={self.fourier_n}), got {n_atom}. '
                    f'Suggested max modes for {n_atom} atoms is '
                    f'fourier_m = fourier_n = {max_square_mode}.'
                )
        elif self.surface_method == self.SurfaceMethod.BINNING_NEAREST:
            if wrap is True:
                raise ValueError("wrap=True is not valid when surface_method='binning_nearest'")
            if grid_origin == 'lipid_bbox' and fft_filter is not None:
                raise ValueError("grid_origin='lipid_bbox' cannot be combined with fft_filter")
            if grid_origin == 'lipid_bbox' and self.padding:
                raise ValueError("grid_origin='lipid_bbox' cannot be combined with padding")
            if self.padding:
                validate_orthorhombic(universe.dimensions)
            self.wrap = False
            self.fft_filter = fft_filter
            self._fft_q_bounds = None
            if self.fft_filter is not None:
                self._fft_q_bounds = resolve_fft_filter(self.fft_filter, self.dx, self.dy)
                logger.info(
                    f'fft_filter={self.fft_filter!r} resolved to '
                    f'q_low={self._fft_q_bounds[0]:.6g}, q_high={self._fft_q_bounds[1]:.6g} (rad/A)',
                )
        else:
            self.wrap = True if wrap is None else wrap
            self.fft_filter = fft_filter
            self._fft_q_bounds = None
            if self.padding:
                # padding only makes sense for orthorhombic boxes
                validate_orthorhombic(universe.dimensions)
                self._pad_spec = padded_grid_spec(
                    self.n_x_bins,
                    self.n_y_bins,
                    self.x_range,
                    self.y_range,
                    self.edge_pad_bins,
                )
            if self.fft_filter is not None:
                self._fft_q_bounds = resolve_fft_filter(self.fft_filter, self.dx, self.dy)
                logger.info(
                    f'fft_filter={self.fft_filter!r} resolved to '
                    f'q_low={self._fft_q_bounds[0]:.6g}, q_high={self._fft_q_bounds[1]:.6g} (rad/A)',
                )

        self.curvature_on = self._resolve_curvature_on(curvature_on)

        # Only checks the first frame. NPT simulations not properly covered here.
        # Warning message if range doesn't cover entire dimensions of simulation box
        for dim_string, dim_range, num in [('x', self.x_range, 0), ('y', self.y_range, 1)]:
            if dim_range[1] < universe.dimensions[num]:
                msg = (
                    f'Grid range in {dim_string} does not cover entire '
                    'dimensions of simulation box.\n Minimum dimensions '
                    'must be equal to simulation box.'
                )
                warnings.warn(msg)
                logger.warning(msg)

        # Warn when binning without wrapping coordinates
        if not self.wrap and self.surface_method == self.SurfaceMethod.BINNING and not self.padding:
            # Warning
            msg = (
                ' `wrap == False` may result in inaccurate calculation '
                "of membrane curvature with `surface_method='binning'`. "
                'Surfaces will be derived from a reduced number of atoms. \n '
                ' Ignore this warning if your trajectory has '
                ' rotational/translational fit rotations! '
            )
            warnings.warn(msg)
            logger.warning(msg)

    def _resolve_curvature_on(self, curvature_on):
        """Resolve curvature calculation mode depending on the value of `curvature_on`."""
        if curvature_on is None:
            if self._fft_q_bounds is not None:
                return self.CurvatureOn.AVERAGE_SURFACE
            return self.CurvatureOn.PER_FRAME
        try:
            return self.CurvatureOn(curvature_on)
        except ValueError as err:
            allowed = ', '.join(repr(method.value) for method in self.CurvatureOn)
            raise ValueError(f'curvature_on must be None or one of ({allowed}), got {curvature_on!r}') from err

    def _wrap_xy(self):
        """Wrap x,y into the unit cell and leave z unchanged."""
        z = self.ag.positions[:, 2].copy()
        self.ag.wrap()
        positions = self.ag.positions
        positions[:, 2] = z
        self.ag.positions = positions

    def _prepare(self):
        # Initialize empty np.array with results
        self.results.z_surface = np.full((self.n_frames, self.n_x_bins, self.n_y_bins), np.nan)
        self.results.mean = np.full((self.n_frames, self.n_x_bins, self.n_y_bins), np.nan)
        self.results.gaussian = np.full((self.n_frames, self.n_x_bins, self.n_y_bins), np.nan)
        self._fourier_theta_sum = None
        self._fourier_modes = None
        if self.surface_method == self.SurfaceMethod.FOURIER and self.curvature_on == self.CurvatureOn.AVERAGE_SURFACE:
            self._fourier_modes = fourier_mode_list(self.fourier_m, self.fourier_n)
            self._fourier_theta_sum = np.zeros(
                n_fourier_parameters(self.fourier_m, self.fourier_n),
                dtype=np.float64,
            )
            self._fourier_theta_n_finite = 0

    def _single_frame(self):
        if self.surface_method == self.SurfaceMethod.BINNING:
            if self.wrap:
                self._wrap_xy()
            if not self.padding:
                z_surface = get_z_surface(
                    self.ag.positions,
                    n_x_bins=self.n_x_bins,
                    n_y_bins=self.n_y_bins,
                    x_range=self.x_range,
                    y_range=self.y_range,
                )
                self.results.z_surface[self._frame_index] = z_surface
                self.results.mean[self._frame_index] = mean_curvature(z_surface, self.dx, self.dy)
                self.results.gaussian[self._frame_index] = gaussian_curvature(z_surface, self.dx, self.dy)
            else:
                z_padded = get_z_surface_padded(
                    self.ag.positions,
                    self.n_x_bins,
                    self.n_y_bins,
                    self.x_range,
                    self.y_range,
                    self.edge_pad_bins,
                )
                z_surface, mean_c, gauss_c = curvature_with_edge_pad(z_padded, self.dx, self.dy, self.edge_pad_bins)
                self.results.z_surface[self._frame_index] = z_surface
                self.results.mean[self._frame_index] = mean_c
                self.results.gaussian[self._frame_index] = gauss_c
        elif self.surface_method == self.SurfaceMethod.BINNING_NEAREST:
            box = self._ts.dimensions.copy()
            if self._x_range_user or self.grid_origin == 'box':
                x_range = self.x_range
            else:
                x_range = None
            if self._y_range_user or self.grid_origin == 'box':
                y_range = self.y_range
            else:
                y_range = None
            z_surface, grid_axes = get_z_surface_nearest(
                self.ag,
                n_x_bins=self.n_x_bins,
                n_y_bins=self.n_y_bins,
                box=box,
                grid_origin=self.grid_origin,
                x_range=x_range,
                y_range=y_range,
            )
            frame_dx = grid_axes['dx']
            frame_dy = grid_axes['dy']
            self.results.z_surface[self._frame_index] = z_surface
            if self.padding:
                _, mean_c, gauss_c = curvature_from_primary_with_edge_pad(
                    z_surface, frame_dx, frame_dy, self.edge_pad_bins
                )
                self.results.mean[self._frame_index] = mean_c
                self.results.gaussian[self._frame_index] = gauss_c
            else:
                self.results.mean[self._frame_index] = mean_curvature(z_surface, frame_dx, frame_dy)
                self.results.gaussian[self._frame_index] = gaussian_curvature(z_surface, frame_dx, frame_dy)
        else:
            z_f, mean_f, gauss_f, theta = fourier_curvature(
                self.ag.positions,
                self.x_range,
                self.y_range,
                self.n_x_bins,
                self.n_y_bins,
                self.fourier_m,
                self.fourier_n,
                rcond=self.fourier_rcond,
                return_theta=True,
            )
            if self.curvature_on == self.CurvatureOn.AVERAGE_SURFACE and np.all(np.isfinite(theta)):
                self._fourier_theta_sum += theta
                self._fourier_theta_n_finite += 1
            self.results.z_surface[self._frame_index] = z_f
            self.results.mean[self._frame_index] = mean_f
            self.results.gaussian[self._frame_index] = gauss_f

    def _conclude(self):
        if self.surface_method == self.SurfaceMethod.FOURIER and self.curvature_on == self.CurvatureOn.AVERAGE_SURFACE:
            n_kept = self._fourier_theta_n_finite
            n_dropped = self.n_frames - n_kept
            if n_kept == 0:
                warnings.warn(
                    f'All {self.n_frames} frames were excluded from the average '
                    'surface because their Fourier fit did not produce valid '
                    'coefficients (NaN or infinite values). The average surface '
                    'and curvature are undefined (NaN).',
                    UserWarning,
                    stacklevel=2,
                )
                avg_theta = np.full_like(self._fourier_theta_sum, np.nan)
            else:
                assert self._fourier_theta_sum is not None
                if n_dropped:
                    warnings.warn(
                        f'{n_dropped} of {self.n_frames} frames were excluded from the '
                        'average surface because their Fourier fit did not produce valid '
                        f'coefficients (NaN or infinite values). The average curvature is '
                        f'computed from the remaining {n_kept} frames.',
                        UserWarning,
                        stacklevel=2,
                    )
                avg_theta = self._fourier_theta_sum / n_kept
            z_avg, avg_mean, avg_gauss = fourier_curvature_from_theta(
                avg_theta,
                self._fourier_modes,
                self.x_range,
                self.y_range,
                self.n_x_bins,
                self.n_y_bins,
            )
            # Same averaged coefficients for height and analytic curvature.
            self.results.average_z_surface = z_avg
            self.results.average_mean = avg_mean
            self.results.average_gaussian = avg_gauss
            return

        z_average = np.nanmean(self.results.z_surface, axis=0)
        if self._fft_q_bounds is not None:
            self.results.average_z_surface = apply_fft_filter(z_average, self.dx, self.dy, self._fft_q_bounds)
        else:
            self.results.average_z_surface = z_average

        if self.curvature_on == self.CurvatureOn.PER_FRAME:
            self.results.average_mean = np.nanmean(self.results.mean, axis=0)
            self.results.average_gaussian = np.nanmean(self.results.gaussian, axis=0)
            return

        z_for_curvature = self.results.average_z_surface
        if self.padding:
            _, avg_mean, avg_gauss = curvature_from_primary_with_edge_pad(
                z_for_curvature, self.dx, self.dy, self.edge_pad_bins
            )
            self.results.average_mean = avg_mean
            self.results.average_gaussian = avg_gauss
        else:
            self.results.average_mean = mean_curvature(z_for_curvature, self.dx, self.dy)
            self.results.average_gaussian = gaussian_curvature(z_for_curvature, self.dx, self.dy)
