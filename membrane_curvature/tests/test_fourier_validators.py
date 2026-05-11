"""
Unit tests for the Fourier-surface input validators.
"""

import numpy as np
import pytest

from membrane_curvature.fourier_validators import (
    _coerce_positions,
    validate_positions,
    validate_positive_bin_counts,
    validate_positive_domain_widths,
)


class TestCoercePositions:
    def test_accepts_valid_input(self):
        positions = np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
        coerced = _coerce_positions(positions)
        assert coerced.dtype == np.float64
        assert coerced.shape == (2, 3)

    def test_coerces_integer_to_float64(self):
        positions = np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int64)
        coerced = _coerce_positions(positions)
        assert coerced.dtype == np.float64
        assert coerced.shape == (2, 3)

    def test_accepts_more_than_three_columns(self):
        positions = np.array([[0.0, 1.0, 2.0, 99.0], [3.0, 4.0, 5.0, 99.0]])
        coerced = _coerce_positions(positions)
        assert coerced.shape == (2, 4)

    def test_invalid_ragged_input(self):
        with pytest.raises(ValueError, match=r'positions could not be converted') as exc_info:
            _coerce_positions([[1.0, 2.0, 3.0], [4.0, 5.0]])
        assert isinstance(exc_info.value.__cause__, (TypeError, ValueError))

    def test_invalid_non_numeric_input(self):
        with pytest.raises(ValueError, match=r'positions could not be converted') as exc_info:
            _coerce_positions([['a', 'b', 'c']])
        assert isinstance(exc_info.value.__cause__, (TypeError, ValueError))

    def test_handles_none(self):
        # NumPy's handling of ``None`` with ``dtype=float64`` has changed across
        # versions: it may yield a 0-d ``nan`` array (then rejected by
        # ``validate_positions``) or raise (then wrapped by ``_coerce_positions``).
        # Either way the public validator pair must reject it.
        with pytest.raises(
            ValueError,
            match=r'could not be converted|must be a 2D array|must have at least 3 columns',
        ):
            validate_positions(_coerce_positions(None))


class TestValidatePositions:
    def test_accepts_valid_input(self):
        validate_positions(np.zeros((2, 3)))

    def test_invalid_1d_array(self):
        with pytest.raises(ValueError, match=r'positions must be a 2D array; got ndim=1'):
            validate_positions(np.array([1.0, 2.0, 3.0]))

    def test_invalid_0d_array(self):
        with pytest.raises(ValueError, match=r'positions must be a 2D array; got ndim=0'):
            validate_positions(np.asarray(1.0))

    def test_invalid_fewer_than_three_columns(self):
        with pytest.raises(ValueError, match=r'positions must have at least 3 columns'):
            validate_positions(np.zeros((5, 2)))

    def test_invalid_empty(self):
        with pytest.raises(ValueError, match=r'positions must contain at least one row'):
            validate_positions(np.empty((0, 3)))


class TestValidatePositiveBinCounts:
    @pytest.mark.parametrize('n_x_bins,n_y_bins', [(0, 1), (1, 0), (-1, 3), (0, 0), (5, -1)])
    def test_invalid_nonpositive(self, n_x_bins, n_y_bins):
        with pytest.raises(ValueError, match=r'n_x_bins and n_y_bins must be positive'):
            validate_positive_bin_counts(n_x_bins, n_y_bins)

    def test_accepts_positive(self):
        validate_positive_bin_counts(1, 1)
        validate_positive_bin_counts(100, 200)

    def test_accepts_numpy_integers(self):
        validate_positive_bin_counts(np.int64(10), np.int32(20))

    @pytest.mark.parametrize(
        'n_x_bins,n_y_bins',
        [
            (1.5, 2),
            (2, 1.5),
            (1.5, 1.5),
            (True, 2),
            (2, False),
            ('10', 10),
            (None, 10),
            (np.bool_(True), 2),
            (2, np.bool_(False)),
            (np.float64(2.0), 2),
        ],
    )
    def test_invalid_non_integer(self, n_x_bins, n_y_bins):
        with pytest.raises(ValueError, match=r'must be integers'):
            validate_positive_bin_counts(n_x_bins, n_y_bins)

    def test_non_integer_message_format(self):
        """
        When both inputs are non-integer, the error message must contain both
        name=value pairs together with their type, so a user can see exactly
        what triggered the error.
        """
        with pytest.raises(ValueError) as exc_info:
            validate_positive_bin_counts(1.5, 2.5)
        msg = str(exc_info.value)
        assert 'n_x_bins=1.5' in msg
        assert 'n_y_bins=2.5' in msg
        assert 'float' in msg


class TestValidatePositiveDomainWidths:
    @pytest.mark.parametrize(
        'x_range,y_range',
        [
            ((0.0, 10.0), (0.0, 10.0)),
            ((-5.0, 5.0), (0.0, 1.0)),
            ((0.0, 1e-3), (0.0, 1e-3)),
        ],
    )
    def test_accepts_positive(self, x_range, y_range):
        validate_positive_domain_widths(x_range, y_range)

    @pytest.mark.parametrize(
        'x_range,y_range',
        [
            ((0.0, 0.0), (0.0, 10.0)),
            ((10.0, 10.0), (0.0, 10.0)),
            ((0.0, 300.0), (5.0, 3.0)),
            ((20.0, 10.0), (0.0, 10.0)),
            ((0.0, 1.0), (5.0, 5.0)),
        ],
    )
    def test_error_nonpositive(self, x_range, y_range):
        with pytest.raises(ValueError, match=r'x_range and y_range must have positive width'):
            validate_positive_domain_widths(x_range, y_range)
