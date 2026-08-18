from pathlib import Path

from harp.deepcoil import parse_deepcoil2


def test_canonical_register(tmp_path: Path):
    path = tmp_path / "dc.out"
    rows = ["aa cc raw_cc prob_a prob_d"]
    # Confident a anchors at positions 1,8; d anchors at 4,11.
    for i in range(1, 15):
        pa = 0.9 if i in {1, 8} else 0.0
        pd = 0.9 if i in {4, 11} else 0.0
        rows.append(f"A 0.9 0.9 {pa} {pd}")
    path.write_text("\n".join(rows) + "\n")
    reg = parse_deepcoil2(path)
    assert reg.d_offset == 3
    assert reg.a_origin == 1
    assert reg.cc_start == 1
    assert reg.cc_end == 14
