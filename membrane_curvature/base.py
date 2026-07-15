"""
MembraneCurvature
=======================================

:Author: Estefania Barreto-Ojeda
:Year: 2021
:Copyright: GNU Public License v3

:class:`~membrane_curvature.base.MembraneCurvature` is the main analysis class
for calculating mean and Gaussian curvature from atom selections.

It derives a height surface from the selected reference atoms using either the
``"fourier"`` method (default) or the ``"binning"`` method, then evaluates
curvature on the resulting surface. The specific operations used to derive
the surface depend on the method selected by the user
in the parameter :attr:`~membrane_curvature.base.MembraneCurvature.surface_method`.

The set of required parameters to run :class:`~membrane_curvature.base.MembraneCurvature` varies depending
on the selected ``surface_method``:

- **Fourier method**:

    Required parameters are the maximum Fourier mode indices ``fourier_m`` and ``fourier_n`` (Default: ``2``).
    Optional parameters include a singular-value cutoff for the Fourier fit via truncated SVD with ``fourier_rcond``.

- **Binning method**:

    Required parameters are the number of bins in the x and y directions ``n_x_bins`` and ``n_y_bins``.
    Optional ``wrap`` parameter to control whether to wrap the coordinates exceeding the simulation box dimensions.
    Optional ``padding`` applies a periodic buffer before finite differences (in orthorhombic boxes only);
    ``edge_pad_bins`` sets the buffer width in bins (default ``2``).

Mean curvature is calculated in units of Å :sup:`-1` and Gaussian curvature in units of Å :sup:`-2`.
"""

import numpy as np
import warnings
from .binning_surface import get_z_surface
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
    surface_method : {'binning', 'fourier'}, optional
        ``fourier`` (default) fits a periodic Fourier sum to atom positions at
        each frame and evaluates Monge-gauge curvature from analytic derivatives
        on the same bin grid (bin centers). ``binning`` derives the surface by creating
        a grid and assigning atoms to bins. It uses :func:`numpy.gradient` for derivatives.
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
        Brick-wall filter on the binned height field when ``surface_method='binning'``.
        Default is ``None``. Pass ``'auto'`` to enable automatic filtering with low-pass set
        as ``(0, 0.5 * q_Nyq)``. For custom bounds in rad/Å, pass ``{'q': (q_low, q_high)}``.
        For **average** maps: time-average of ``z_surface``, filter once, then curvature
        on that filtered height. Per-frame arrays are not FFT-filtered. Ignored for
        ``surface_method='fourier'``. Compatible with ``padding``.
    padding : bool, optional
        Apply periodic edge padding for ``surface_method='binning'``. Default ``False``.
        Pads the simulation box using periodic images in ``x`` and ``y``, compute surface
        and curvature on the padded grid, then clip back to ``(n_x_bins, n_y_bins)``.
        Requires an orthorhombic box. Invalid with ``surface_method='fourier'``.
    edge_pad_bins : int, optional
        Buffer width in bins on each side when ``padding=True`` (default ``2``).
        Ignored when ``padding=False``.

    Attributes
    ----------
    results.z_surface : ndarray
        Per-frame height field from the atom selection (unfiltered when
        ``fft_filter`` is used with `surface_method='binning'`).
        Shape (`n_frames`, `n_x_bins`, `n_y_bins`).
    results.mean : ndarray
        Per-frame mean curvature associated with the surface.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.gaussian : ndarray
        Per-frame Gaussian curvature associated with the surface.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.average_z_surface : ndarray
        Average of the array elements in `z_surface`. With binning and
        ``fft_filter`` enabled, this is the FFT-filtered temporal mean of ``z_surface``,
        not the mean of per-frame filtered surfaces.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_mean : ndarray
        Average of the array elements in `mean_curvature`. With binning
        and ``fft_filter`` enabled, curvature of the filtered time-averaged height
        (not the time average of per-frame ``results.mean``).
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_gaussian : ndarray
        Average of the array elements in `gaussian_curvature`. With
        binning and ``fft_filter`` enabled, curvature of the filtered time-averaged
        height (not the time average of per-frame ``results.gaussian``).
        Each array has shape (`n_x_bins`, `n_y_bins`)

    Raises
    ------
    ValueError
        If ``n_x_bins`` or ``n_y_bins`` is not a positive integer
        (see :func:`~membrane_curvature.fourier_validators.validate_positive_bin_counts`),
        if the selection is empty, if ``surface_method`` is not one of
        ``'binning'`` or ``'fourier'``, or, when ``surface_method='fourier'``,
        if ``fourier_m`` / ``fourier_n`` are negative or the selection has fewer
        atoms than Fourier parameters, if ``wrap=True`` with
        ``surface_method='fourier'``, if a manual ``fft_filter`` dict is passed
        with ``surface_method='fourier'``, if ``padding=True`` with
        ``surface_method='fourier'``, if ``padding=True`` on a
        non-orthorhombic box, or if ``edge_pad_bins`` is not an integer
        ``>= 1`` when ``padding=True``.

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

    Padding is only available for orthorhombic boxes. With ``padding=True``, finite differences
    run on a padded height field built from periodic images, then the buffer of ``edge_pad_bins``
    on each side is clipped so stored arrays keep the primary grid shape.
    The default ``edge_pad_bins=2``, which uses two bins on each side for padding, is enough to
    reduce finite difference artifacts at edges and corners.

    For any method of choice, the derived surface and calculated curvatures are available
    in the :attr:`results` attributes.

    The attribute :attr:`~MembraneCurvature.results.average_z_surface` contains
    the time-averaged derived surface. When ``fft_filter`` is set, the brick-wall filter
    is applied to that averaged surface.

    The attributes :attr:`~MembraneCurvature.results.average_mean` and
    :attr:`~MembraneCurvature.results.average_gaussian` contain mean and
    Gaussian curvature maps for analysis and plotting; with ``fft_filter`` on
    binning surfaces they are curvatures of the filtered average height, not the
    time average of per-frame curvatures.

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
        padding=False,
        edge_pad_bins=2,
        **kwargs,
    ):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.ag = universe.select_atoms(select)
        # Validate selection up front so an empty AtomGroup produces a clear
        # message before any downstream check (bin counts, surface method,
        # Fourier-parameter sizing) can mask it.
        if len(self.ag) == 0:
            raise ValueError('Invalid selection. Empty AtomGroup.')
        validate_positive_bin_counts(n_x_bins, n_y_bins)
        self.n_x_bins = n_x_bins
        self.n_y_bins = n_y_bins
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

        valid_methods = ('binning', 'fourier')
        if surface_method not in valid_methods:
            raise ValueError(f'surface_method must be one of {valid_methods}, got {surface_method!r}')
        self.surface_method = surface_method
        self.fourier_m = int(fourier_m)
        self.fourier_n = int(fourier_n)
        self.fourier_rcond = fourier_rcond

        if self.surface_method == 'fourier':
            if wrap is True:
                raise ValueError("wrap=True is only valid when surface_method='binning'")
            if isinstance(fft_filter, dict):
                raise ValueError("fft_filter dict is only allowed when surface_method='binning'")
            if self.padding:
                raise ValueError("padding=True is only valid when surface_method='binning'")
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

        if not self.wrap and self.surface_method == 'binning' and not self.padding:
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
        if self.surface_method == 'binning':
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
        # resolve bounds and apply filter ONLY when binning + filtering is enabled
        if self.surface_method == 'binning' and self._fft_q_bounds is not None:
            self.results.average_z_surface = apply_fft_filter(z_average, self.dx, self.dy, self._fft_q_bounds)
            z_filtered = self.results.average_z_surface
            if self.padding:
                _, avg_mean, avg_gauss = curvature_from_primary_with_edge_pad(
                    z_filtered, self.dx, self.dy, self.edge_pad_bins
                )
                self.results.average_mean = avg_mean
                self.results.average_gaussian = avg_gauss
            else:
                self.results.average_mean = mean_curvature(z_filtered, self.dx, self.dy)
                self.results.average_gaussian = gaussian_curvature(z_filtered, self.dx, self.dy)
        else:
            self.results.average_z_surface = z_average
            self.results.average_mean = np.nanmean(self.results.mean, axis=0)
            self.results.average_gaussian = np.nanmean(self.results.gaussian, axis=0)
