from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

PHASES = tuple("abcdefg")


@dataclass(frozen=True)
class DeepCoilRegister:
    table: pd.DataFrame
    a_origin: int
    d_offset: int
    cc_start: int
    cc_end: int
    anchor_threshold: float
    cc_threshold: float
    a_purity: float
    d_purity: float
    n_a_calls: int
    n_d_calls: int
    a_tie: bool
    d_tie: bool

    @property
    def register_purity(self) -> float:
        """Fraction of confident a/d calls supporting the assigned dominant mod-7 frame."""
        total = self.n_a_calls + self.n_d_calls
        if total == 0:
            return 0.0
        return (
            self.a_purity * self.n_a_calls + self.d_purity * self.n_d_calls
        ) / total

    @property
    def sequence(self) -> str:
        return "".join(self.table["aa"].astype(str))

    @property
    def length(self) -> int:
        return int(len(self.table))

    def write_tsv(self, path: str | Path) -> None:
        self.table.to_csv(path, sep="\t", index=False)


def _modal_modulo(values: Iterable[int], modulo: int = 7) -> tuple[int, float, bool]:
    """Return modal residue class, its purity, and whether the mode is tied.

    Purity is the fraction of calls in the returned class. Values below 1.0
    diagnose disagreement with a single mod-``modulo`` frame; they do not, by
    themselves, identify the biological or prediction-level cause.
    """
    values = np.asarray(list(values), dtype=int)
    if values.size == 0:
        raise ValueError("No anchor positions passed to modulo inference")
    residues = values % modulo
    counts = np.bincount(residues, minlength=modulo)
    top = int(counts.max())
    winners = np.flatnonzero(counts == top)
    modal = int(winners[0])
    purity = float(top) / float(values.size)
    return modal, purity, bool(winners.size > 1)


def parse_deepcoil2(
    path: str | Path,
    *,
    cc_threshold: float = 0.5,
    anchor_threshold: float = 0.5,
) -> DeepCoilRegister:
    """Parse DeepCoil2 tabular output and infer a fixed a–g register.

    DeepCoil2 supplies probabilities for a and d positions. HARP uses the
    dominant modulo-7 a anchor as phase zero and checks that d is a+3.
    The register is assigned only within the contiguous CC-positive region.
    """
    path = Path(path)
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    required = {"aa", "cc", "raw_cc", "prob_a", "prob_d"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"DeepCoil2 file lacks columns: {sorted(missing)}")

    df = df.copy()
    df.insert(0, "seq_pos", np.arange(1, len(df) + 1, dtype=int))
    cc_mask = df["cc"].to_numpy(float) >= cc_threshold
    if not np.any(cc_mask):
        raise ValueError(f"No residues exceed cc threshold {cc_threshold}")

    cc_positions = df.loc[cc_mask, "seq_pos"].to_numpy(int)
    # Use the largest contiguous CC segment, rather than silently spanning gaps.
    splits = np.where(np.diff(cc_positions) > 1)[0] + 1
    segments = np.split(cc_positions, splits)
    segment = max(segments, key=len)
    cc_start, cc_end = int(segment[0]), int(segment[-1])
    in_segment = df["seq_pos"].between(cc_start, cc_end)

    a_calls = df.loc[in_segment & (df["prob_a"] >= anchor_threshold), "seq_pos"].to_numpy(int)
    d_calls = df.loc[in_segment & (df["prob_d"] >= anchor_threshold), "seq_pos"].to_numpy(int)
    if a_calls.size < 2 or d_calls.size < 2:
        raise ValueError(
            "Too few confident DeepCoil a/d calls to establish a register; "
            f"a={a_calls.size}, d={d_calls.size}, threshold={anchor_threshold}"
        )

    a_mod, a_purity, a_tie = _modal_modulo(a_calls)
    d_mod, d_purity, d_tie = _modal_modulo(d_calls)
    d_offset = (d_mod - a_mod) % 7
    if d_offset != 3:
        raise ValueError(
            "DeepCoil anchors are not internally canonical: dominant d phase is "
            f"a+{d_offset}, expected a+3. Inspect register before analysis."
        )

    # Pick the earliest confident a call in the dominant modulo class as a human-readable origin.
    a_origin = int(a_calls[(a_calls % 7) == a_mod].min())
    phases: list[str | None] = []
    phase_index: list[float] = []
    for pos in df["seq_pos"].to_numpy(int):
        if cc_start <= pos <= cc_end:
            idx = int((pos - a_origin) % 7)
            phases.append(PHASES[idx])
            phase_index.append(float(idx))
        else:
            phases.append(None)
            phase_index.append(np.nan)

    df["phase"] = phases
    df["phase_index"] = phase_index
    df["is_cc_segment"] = in_segment
    df["is_a_anchor"] = df["prob_a"] >= anchor_threshold
    df["is_d_anchor"] = df["prob_d"] >= anchor_threshold

    # Explicit consistency diagnostics are retained in the table.
    df["register_support"] = np.where(
        df["phase"].eq("a"), df["prob_a"],
        np.where(df["phase"].eq("d"), df["prob_d"], np.nan),
    )

    return DeepCoilRegister(
        table=df,
        a_origin=a_origin,
        d_offset=d_offset,
        cc_start=cc_start,
        cc_end=cc_end,
        anchor_threshold=anchor_threshold,
        cc_threshold=cc_threshold,
        a_purity=a_purity,
        d_purity=d_purity,
        n_a_calls=int(a_calls.size),
        n_d_calls=int(d_calls.size),
        a_tie=a_tie,
        d_tie=d_tie,
    )
