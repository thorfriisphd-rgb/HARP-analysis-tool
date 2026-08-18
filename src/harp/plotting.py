from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PHASES = tuple("abcdefg")


def _taxon_title(taxon: str | None, description: str) -> str:
    """Build a consistent plot title without requiring taxon metadata."""
    if taxon:
        return f"HARP v4.1 — {taxon} — {description}"
    return f"HARP v4.1 — {description}"


def plot_contact_profile(df: pd.DataFrame, path: str | Path, taxon: str | None = None) -> None:
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(df["seq_pos"], df["contact_occupancy"], linewidth=1.5)
    cc = df["is_cc_segment"].astype(bool)
    if cc.any():
        ax.axvspan(df.loc[cc, "seq_pos"].min(), df.loc[cc, "seq_pos"].max(), alpha=0.12)
    label_rows = df.loc[cc & df["phase"].isin(["a", "d"])]
    for _, row in label_rows.iterrows():
        ax.text(row["seq_pos"], -0.035, row["phase"], ha="center", va="top", fontsize=7)
    ax.set_xlabel("MyhT sequence position")
    ax.set_ylabel("MG contact occupancy")
    ax.set_ylim(bottom=-0.08)
    ax.set_title(_taxon_title(taxon, "trajectory-derived MyhT contact profile"))
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_phase_summary(phase_df: pd.DataFrame, path: str | Path, taxon: str | None = None) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(phase_df["phase"], phase_df["contact_fraction"])
    ax.plot(phase_df["phase"], phase_df["baseline_fraction"], marker="o", label="Residue-count baseline")
    ax.set_xlabel("DeepCoil-anchored heptad phase")
    ax.set_ylabel("Fraction")
    ax.set_title(_taxon_title(taxon, "contact occupancy by heptad phase"))
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_time_blocks(block_vectors: np.ndarray, path: str | Path, taxon: str | None = None) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(block_vectors, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(7), PHASES)
    ax.set_xlabel("Heptad phase")
    ax.set_ylabel("Trajectory block")
    ax.set_title(_taxon_title(taxon, "temporal persistence of phase occupancy"))
    fig.colorbar(im, ax=ax, label="Contact fraction")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_null(null: np.ndarray, observed: float, title: str, path: str | Path, taxon: str | None = None) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(null, bins=50)
    ax.axvline(observed, linestyle="--", linewidth=2, label=f"Observed = {observed:.4g}")
    ax.set_xlabel("Null statistic")
    ax.set_ylabel("Count")
    ax.set_title(_taxon_title(taxon, title))
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

def plot_panel_null(
    null: np.ndarray,
    observed: float,
    *,
    p_value: float,
    null_mean: float,
    null_sd: float,
    null_q95: float,
    n_taxa: int,
    n_permutations: int,
    title: str,
    path: str | Path,
) -> None:
    """Publication-facing panel-null figure with inferential annotation.

    Presentation only: all statistics are calculated upstream.
    """
    null = np.asarray(null, dtype=float)

    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.hist(null, bins=50)

    ax.axvline(
        observed,
        linestyle="--",
        linewidth=2.2,
        label=f"Observed = {observed:.4f}",
    )
    ax.axvline(
        null_q95,
        linestyle=":",
        linewidth=2.0,
        label=f"Null 95th percentile = {null_q95:.4f}",
    )

    p_text = f"{p_value:.1e}" if p_value < 0.001 else f"{p_value:.4f}"
    stats_text = (
        f"n = {int(n_taxa)} taxa\n"
        f"{int(n_permutations):,} permutations\n"
        f"Permutation p = {p_text}\n"
        f"Null mean = {null_mean:.4f}\n"
        f"Null SD = {null_sd:.4f}"
    )

    ax.text(
        0.985,
        0.965,
        stats_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10.5,
        bbox=dict(
            boxstyle="round,pad=0.45",
            facecolor="white",
            edgecolor="0.65",
            alpha=0.92,
        ),
    )

    ax.set_xlabel("Mean cross-taxon phase similarity")
    ax.set_ylabel("Permutation count")
    ax.set_title(title)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)

