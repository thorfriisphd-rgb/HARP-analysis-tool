from datetime import datetime, timezone

from harp.provenance import panel_output_name


def test_panel_output_name_contains_taxon_count_and_timestamp():
    dt = datetime(2026, 8, 13, 20, 0, 0, tzinfo=timezone.utc)
    name = panel_output_name(26, dt)
    assert name.startswith("HARP_v4.1_panel_n26_")
    assert "20260813T200000+0000" in name
