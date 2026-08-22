import numpy as np

from crd_utils import plotting


def test_sn_plots_accept_undefined_bins(tmp_path):
    bin_map = np.array([[0, 0, 1], [2, 2, -1]], dtype=int)
    values = np.array([5.0, np.nan, 10.0])
    out1 = tmp_path / "single.png"
    plotting.plot_bin_value_map(
        bin_map,
        values,
        out1,
        title="test",
        colorbar_label="S/N",
        vmin=0.0,
        vmax=10.0,
    )
    assert out1.is_file()

    out2 = tmp_path / "comparison.png"
    plotting.plot_bl_rh3_sn_comparison(
        bin_map,
        np.array([5.0, 6.0, 10.0]),
        values,
        out2,
    )
    assert out2.is_file()
