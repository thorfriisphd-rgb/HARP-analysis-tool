import numpy as np

from harp.plotting import plot_panel_null


def test_money_shot_renderer_writes_png(tmp_path):
    out = tmp_path / "panel_null.png"
    null = np.linspace(-0.2, 0.2, 101)
    plot_panel_null(
        null,
        0.5,
        p_value=0.001,
        null_mean=float(null.mean()),
        null_sd=float(null.std(ddof=1)),
        null_q95=float(np.quantile(null, 0.95)),
        n_taxa=26,
        n_permutations=100,
        title="HARP v4.1 panel null — shared phase alignment",
        path=out,
    )
    assert out.is_file()
    assert out.stat().st_size > 0
