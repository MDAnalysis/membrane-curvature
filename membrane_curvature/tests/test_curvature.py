"""Unit tests for Fourier curvature helpers in membrane_curvature.curvature."""

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from membrane_curvature.curvature import (
    fourier_curvature,
    fourier_curvature_from_coefficients,
    fourier_curvature_from_theta,
)
from membrane_curvature.fourier_surface import fourier_mode_list

# Dummy domain for testing
DOMAIN = dict(x_range=(0.0, 2.0), y_range=(0.0, 2.0), n_x_bins=4, n_y_bins=4)


def test_flat_surface_has_zero_curvature():
    z, H, K = fourier_curvature_from_coefficients(A00=12.5, coeffs={}, **DOMAIN)
    assert_almost_equal(z, np.full_like(z, 12.5))
    assert_almost_equal(H, np.zeros_like(z))
    assert_almost_equal(K, np.zeros_like(z))


@pytest.fixture()
def single_cosine_fit():
    xs = np.linspace(0.0, 2.0, 4, endpoint=False)
    xx, yy = np.meshgrid(xs, xs, indexing='ij')
    zz = 9.0 + 0.5 * np.cos(np.pi * xx)
    positions = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    z, H, K, theta = fourier_curvature(positions, **DOMAIN, M=1, N=0, return_theta=True)
    return theta, (z, H, K)


def test_from_theta_matches_fourier_curvature_fit(single_cosine_fit):
    theta, expected = single_cosine_fit
    result = fourier_curvature_from_theta(theta, fourier_mode_list(1, 0), **DOMAIN)
    for got, exp in zip(result, expected):
        assert_almost_equal(got, exp)
