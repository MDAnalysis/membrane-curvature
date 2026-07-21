"""
MembraneCurvature
=======================================

:Author: Estefania Barreto-Ojeda
:Year: 2021
:Copyright: GNU Public License v3

:class:`~membrane_curvature.base.MembraneCurvature` is the main analysis class
for calculating mean and Gaussian curvature from atom selections.

It derives a height surface from the selected reference atoms using either the
``"fourier"`` method (default) or binning methods (``"binning"`` or ``"binning_nearest"``),
method, then evaluates curvature on the resulting surface. The specific operations used to derive
the surface depend on the method selected by the user in the parameter
:attr:`~membrane_curvature.base.MembraneCurvature.surface_method`.

The set of required parameters to run :class:`~membrane_curvature.base.MembraneCurvature` varies depending
on the selected ``surface_method``:

- **Fourier method**:

    Required parameters are the maximum Fourier mode indices ``fourier_m`` and ``fourier_n`` (Default: ``2``).
    Optional parameters include a singular-value cutoff for the Fourier fit via truncated SVD with ``fourier_rcond``.

- **Binning method**:

    Required parameters are the number of bins in the x and y directions ``n_x_bins`` and ``n_y_bins``.
    Optional ``wrap`` parameter to control whether to wrap the coordinates exceeding the simulation box dimensions.
    Optional ``padding`` applies a periodic buffer before finite differences (in orthorhombic boxes only);
    ``edge_pad_bins`` sets the buffer width in bins (default ``2``). Also available for
    ``surface_method='binning_nearest'`` (wrap-pad of the nearest-lipid height field).

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
    curvature_from_primary_with_edge_pad,
    curvature_with_edge_pad,
)
from .padding import (
    get_z_surface_padded,
    padded_grid_spec,
)
from .padding_validators import validate_edge_pad_bins, validate_orthorhombic
from .fft_filtering import apply_fft_filter, resolve_fft_filter
from .fourier_surface import n_fourier_parameters
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
    wrap : bool or None, optional
        Wrap ``x``/``y`` into the unit cell for ``surface_method='binning'``
        but leave ``z`` unchanged. Defaults to ``True`` for binning and must be
        omitted or explicitly set to ``False`` for ``surface_method='fourier'``.
        With ``padding=True``, lateral PBC are also handled by image tiling.
    n_x_bins : int, optional, default: 100
        Number of bins in grid in the x dimension.
    n_y_bins : int, optional, default: 100
        Number of bins in grid in the y dimension.
    x_range : tuple of (float, float), optional, default: (0, `universe.dimensions[0]`)
        Range of coordinates (min, max) in the x dimension.
    y_range : tuple of (float, float), optional, default: (0, `universe.dimensions[1]`)
        Range of coordinates (min, max) in the y dimension.
    surface_method : {'binning', 'binning_nearest', 'fourier'}, optional
        ``fourier`` (default) fits a periodic Fourier sum to atom positions at
        each frame and evaluates Monge-gauge curvature from analytic derivatives
        on the same bin grid (bin centers). ``binning`` derives the surface by creating
        a grid and assigning atoms to bins. It uses :func:`numpy.gradient` for derivatives.
        ``binning_nearest`` assigns each grid corner the ``z`` of the nearest lipid in
        ``xy`` (minimum-image).
    grid_origin : {'box', 'lipid_bbox'}, optional
        xy grid extent for ``surface_method='binning_nearest'``. Default ``'box'``
        (``[0, lx]``, ``[0, ly]``).
    fourier_m : int, optional
        Maximum Fourier mode index in ``x`` when ``surface_method='fourier'``.
        Default ``2``.
    fourier_n : int, optional
        Maximum Fourier mode index in ``y`` when ``surface_method='fourier'``.
        Default ``2``.
    fourier_rcond : float, optional
        Singular-value cutoff for the Fourier fit via truncated SVD
        with :func:`~membrane_curvature.fourier_surface._solve_design_least_squares_svd`
        when ``surface_method='fourier'``. The cutoff is interpreted as a
        relative threshold on singular values.
    fft_filter : None, ``'auto'``, or dict, optional
        Brick-wall filter on the binned height field for ``surface_method='binning'`` and
        ``surface_method='binning_nearest'``. Default is ``None``. Pass ``'auto'`` to enable
        automatic filtering with low-pass set to ``(0, 0.5 * q_Nyq)``. For custom bounds in
        rad/Å, pass ``{'q': (q_low, q_high)}``. For **average** maps: time-average of
        ``z_surface``, filter once, then curvature on that filtered height when
        ``curvature_on`` resolves to ``'average_surface'``.
        Per-frame arrays are not FFT-filtered. Ignored for ``surface_method='fourier'``.
        Compatible with ``padding``. Not allowed with ``grid_origin='lipid_bbox'``.
    curvature_on : str, optional, default: None
        Where to evaluate curvature for ``results.average_mean`` and ``results.average_gaussian``.
        See :class:`CurvatureOn` for more details.
    padding : bool, optional
        Apply periodic edge padding for binning surface methods. Default ``False``.
        For ``surface_method='binning'``, pads the simulation box using periodic images
        in ``x`` and ``y``, computes surface and curvature on the padded grid, then
        clips back to ``(n_x_bins, n_y_bins)``. For ``surface_method='binning_nearest'``,
        builds the nearest-lipid height field on the primary grid and evaluates curvature
        with a wrap-pad buffer of ``edge_pad_bins``. Requires an orthorhombic box.
        Invalid with ``surface_method='fourier'`` and with ``grid_origin='lipid_bbox'``.
    edge_pad_bins : int, optional
        Buffer width in bins on each side when ``padding=True`` (default ``2``).
        Ignored when ``padding=False``.

    Attributes
    ----------
    results.z_surface : ndarray
        Per-frame height field from the atom selection (unfiltered when
        ``fft_filter`` is used with ``binning`` or ``binning_nearest``).
        Shape (`n_frames`, `n_x_bins`, `n_y_bins`).
    results.mean : ndarray
        Per-frame mean curvature associated with the surface.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.gaussian : ndarray
        Per-frame Gaussian curvature associated with the surface.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.average_z_surface : ndarray
        Average of the array elements in `z_surface`. With ``binning`` or
        ``binning_nearest`` and ``fft_filter`` enabled, this is the FFT-filtered
        temporal mean of ``z_surface``, not the mean of per-frame filtered surfaces.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_mean : ndarray
        Mean curvature map for analysis and plotting.
        With ``curvature_on='per_frame'``, the time average of ``results.mean``.
        With ``curvature_on='average_surface'``, mean curvature of
        ``results.average_z_surface``.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_gaussian : ndarray
        Gaussian curvature map for analysis and plotting.
        With ``curvature_on='per_frame'``, the time average of ``results.gaussian``.
        With ``curvature_on='average_surface'``, Gaussian curvature of
        ``results.average_z_surface``.
        Each array has shape (`n_x_bins`, `n_y_bins`)

    Raises
    ------
    ValueError
        If ``n_x_bins`` or ``n_y_bins`` is not a positive integer
        (see :func:`~membrane_curvature.fourier_validators.validate_positive_bin_counts`),
        if the selection is empty, if ``surface_method`` is not one of
        ``'binning'``, ``'binning_nearest'``, or ``'fourier'``, or, when
        ``surface_method='fourier'``,
        if ``fourier_m`` / ``fourier_n`` are negative or the selection has fewer
        atoms than Fourier parameters, if ``wrap=True`` with
        ``surface_method='fourier'``, if a manual ``fft_filter`` dict is passed
        with ``surface_method='fourier'``, if ``padding=True`` with
        ``surface_method='fourier'``, if ``padding=True`` on a
        non-orthorhombic box, if ``edge_pad_bins`` is not an integer
        ``>= 1`` when ``padding=True``, or if ``curvature_on`` is not
        ``None``, ``'per_frame'``, or ``'average_surface'``.

    See also
    --------
    :class:`~MDAnalysis.transformations.wrap.wrap`
        Wrap/unwrap the atoms of a given AtomGroup in the unit cell.

    Notes
    -----

    **Fourier surface method (default)**

    ``surface_method='fourier'`` uses ``fourier_m = fourier_n = 2`` as default.
    Do not modify the default values unless you need shorter wavelengths.
    Since the method performs periodic boundary conditions by itself, ``wrap`` defaults
    to ``False`` and is not required.

    **Binning mode**

    ``surface_method='binning'`` runs with ``wrap`` set to ``True`` by default:
    only ``x`` and ``y`` are wrapped into the unit cell while ``z`` remains unchanged.
    Omit ``wrap`` or pass ``wrap=True`` for raw trajectories so atoms outside the unit
    cell in ``x``/``y`` are included in the bins. Run with ``wrap=False`` for trajectories
    already wrapped into the primary cell, or after rotational / translational fitting.
    For membrane-protein systems without position restraints, preprocessing should include
    rotational and translational fitting around the protein.

    Padding is only available for orthorhombic boxes. With ``padding=True`` and
    ``surface_method='binning'``, finite differences run on a padded height field built
    from periodic atom images, then the buffer of ``edge_pad_bins`` on each side is
    clipped so stored arrays keep the primary grid shape. With
    ``surface_method='binning_nearest'``, the nearest-lipid height is built on the
    primary grid and curvature uses a wrap-pad buffer of the same width.
    The default ``edge_pad_bins=2``, which uses two bins on each side for padding, is enough to
    reduce finite difference artifacts at edges and corners.

    For any method of choice, the derived surface and calculated curvatures are available
    in the :attr:`results` attributes.

    The attribute :attr:`~MembraneCurvature.results.average_z_surface` contains
    the time-averaged derived surface. When ``fft_filter`` is set, the brick-wall filter
    is applied to that averaged surface.

    The attributes :attr:`~MembraneCurvature.results.average_mean` and
    :attr:`~MembraneCurvature.results.average_gaussian` contain mean and
    Gaussian curvature maps for analysis and plotting. With
    ``curvature_on='per_frame'`` the arrays are time averages of the per-frame
    curvature arrays. With ``curvature_on='average_surface'`` the arrays are curvatures
    of ``average_z_surface`` (including after ``fft_filter`` when that is set).
    Leaving ``curvature_on=None`` preserves the previous defaults:
    ``'per_frame'`` when ``fft_filter`` is ``None``, and ``'average_surface'`` when
    ``fft_filter`` is set.



    Example
    -----------
    You can pass a universe containing your selection of reference::

        import MDAnalysis as mda
        from membrane_curvature.base import MembraneCurvature

        u = mda.Universe(coordinates, trajectory)
        mc = MembraneCurvature(u).run()

        surface =  mc.results.average_z_surface
        mean_curvature =  mc.results.average_mean
        gaussian_curvature = mc.results.average_gaussian

    The respective 2D curvature plots can be obtained using the `matplotlib`
    package for data visualization via :func:`~matplotlib.pyplot.contourf` or
    :func:`~matplotlib.pyplot.imshow`.

    For specific examples visit the :ref:`usage` page.
    Check the :ref:`visualization` page for more examples to plot
    MembraneCurvature results using :func:`~matplotlib.pyplot.contourf`
    and :func:`~matplotlib.pyplot.imshow`.


    """

    class SurfaceMethod(StrEnum):
        """Supported values for ``surface_method``."""

        BINNING = 'binning'
        BINNING_NEAREST = 'binning_nearest'
        FOURIER = 'fourier'

    class CurvatureOn(StrEnum):
        """Supported values for ``curvature_on``.

        Where to evaluate curvature for ``results.average_mean`` and ``results.average_gaussian``.

        -  ``'per_frame'`` time-averages the per-frame curvature arrays
            (:math:`\langle H \rangle`, :math:`\langle K \rangle`).

        - ``'average_surface'`` evaluates curvature on ``results.average_z_surface``
            (:math:`H(\langle S \rangle)`, :math:`K(\langle S \rangle)`).

        - ``None`` (default) keeps prior behavior: ``'per_frame'`` when ``fft_filter`` is
            ``None``, and ``'average_surface'`` when ``fft_filter`` is set.
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
            if isinstance(fft_filter, dict):
                raise ValueError('fft_filter dict is only allowed for binning surface methods')
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
            z_f, mean_f, gauss_f = fourier_curvature(
                self.ag.positions,
                self.x_range,
                self.y_range,
                self.n_x_bins,
                self.n_y_bins,
                self.fourier_m,
                self.fourier_n,
                rcond=self.fourier_rcond,
            )
            self.results.z_surface[self._frame_index] = z_f
            self.results.mean[self._frame_index] = mean_f
            self.results.gaussian[self._frame_index] = gauss_f

    def _conclude(self):
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
