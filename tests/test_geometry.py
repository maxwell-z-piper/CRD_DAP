import numpy as np

from crd_utils.geometry import inclination_from_axis_ratio
from crd_utils.disk_model import make_ring_grid, interpolate_rotation_curve


def test_face_on_axis_ratio_is_small_inclination():
    inc = inclination_from_axis_ratio(1.0, 0.2)
    assert np.isclose(inc, 0.0)


def test_ring_grid_snaps_to_nearest_node():
    grid = make_ring_grid(0.8, 4.55)
    assert np.isclose(grid.r_start, 0.8)
    assert np.isclose(grid.r_final, 4.8)
    assert np.isclose(grid.delta, 0.4)
    assert np.allclose(grid.radii, [0.8, 1.6, 2.4, 3.2, 4.0, 4.8])


def test_rotation_curve_inner_interpolation_and_outer_nan():
    r = np.array([0.0, 0.4, 0.8, 1.2, 2.0])
    values = interpolate_rotation_curve(r, np.array([0.8, 1.6]), np.array([80.0, 120.0]))
    assert np.allclose(values[:4], [0.0, 40.0, 80.0, 100.0])
    assert np.isnan(values[-1])
