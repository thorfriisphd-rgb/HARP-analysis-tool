#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from collections import Counter
from dataclasses import dataclass, asdict
import csv
import json
from datetime import datetime
from typing import Dict, Tuple, List, Any, Literal, Optional, Sequence

import numpy as np
import pandas as pd
import argparse
import subprocess
import sys


from harp_run_scaffold import load_prco_table


# ---------------------------------------------------------------------
# Default paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DEFAULT_PRCO_DIR = DATA_DIR / "prco"
DEFAULT_SAMPLES_PATH = DATA_DIR / "samples.csv"
BASE_OUT_DIR = PROJECT_ROOT / "HARP-output"
OUT_CSV_PREFIX = "harp_panel_summary"

# ---------------------------------------------------------------------
# Null-model configuration
# ---------------------------------------------------------------------

NULL_MODE: Literal["linear_shift", "circular_shift", "shuffle_weights"] = "linear_shift"
LINEAR_SHIFT_ON_OCCUPIED_BLOCK = True
SHUFFLE_ITER = 10000
RNG_SEED = 123


# ---------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------

NULL_MODELS = [
    "shuffle_full_span",
    "shuffle_nonzero_support",
    "shuffle_occupied_span",
]

def get_permutable_indices(weights, null_model):
    """
    Return the indices whose values are allowed to be permuted.

    null_model:
      - shuffle_full_span: whole vector
      - shuffle_nonzero_support: only nonzero-contact residues
      - shuffle_occupied_span: contiguous tract from first to last nonzero
    """
    weights = np.asarray(weights)
    nonzero = np.flatnonzero(weights > 0)

    if null_model == "shuffle_full_span":
        idx = np.arange(len(weights), dtype=int)

    elif null_model == "shuffle_nonzero_support":
        idx = nonzero.astype(int)

    elif null_model == "shuffle_occupied_span":
        if len(nonzero) == 0:
            idx = np.array([], dtype=int)
        else:
            idx = np.arange(nonzero[0], nonzero[-1] + 1, dtype=int)

    else:
        raise ValueError(f"Unknown null model: {null_model}")

    return idx


def permute_profile(weights, null_model, rng):
    """
    Return a permuted copy of weights according to the chosen null model.
    """
    weights = np.asarray(weights, dtype=float)
    out = weights.copy()
    idx = get_permutable_indices(weights, null_model)

    if len(idx) > 1:
        out[idx] = rng.permutation(out[idx])

    return out


def get_support_metadata(weights):
    """
    Handy metadata for later interpretation and CSV output.
    """
    weights = np.asarray(weights)
    nonzero = np.flatnonzero(weights > 0)

    if len(nonzero) == 0:
        return {
            "nonzero_count": 0,
            "occupied_start": -1,
            "occupied_end": -1,
            "occupied_span_len": 0,
        }

    return {
        "nonzero_count": int(len(nonzero)),
        "occupied_start": int(nonzero[0]),
        "occupied_end": int(nonzero[-1]),
        "occupied_span_len": int(nonzero[-1] - nonzero[0] + 1),
    }

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run HARP panel analysis across PRCO tables"
    )

    parser.add_argument(
        "--samples",
        default=None,
        help=(
            "CSV or TSV sample metadata file. Expected header: "
            "SAMPLE/ #SAMPLE, WORKDIR, C12_START, C12_END, MYH_START, MYH_END. "
            "Default: data/samples.csv."
        ),
    )

    parser.add_argument(
        "--prco-dir",
        default=str(DEFAULT_PRCO_DIR),
        help="Directory containing *_prco.csv files",
    )

    parser.add_argument(
        "--out-dir",
        default=str(BASE_OUT_DIR),
        help="Base directory for timestamped HARP output folders",
    )

    parser.add_argument(
        "--null-model",
        default="shuffle_occupied_span",
        choices=[
            "shuffle_full_span",
            "shuffle_nonzero_support",
            "shuffle_occupied_span",
            "all",
        ],
        help="Permutation null model to use",
    )

    parser.add_argument(
        "--shuffle-iter",
        type=int,
        default=10000,
        help="Number of permutations",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Random seed",
    )

    parser.add_argument(
        "--out-prefix",
        default=OUT_CSV_PREFIX,
        help="Prefix for timestamped output files",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and exit without running HARP",
    )

    return parser.parse_args()


def load_sample_map(samples_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Load the samples metadata file.

    Accepted formats:
    - CSV: SAMPLE,WORKDIR,C12_START,C12_END,MYH_START,MYH_END
    - TSV: SAMPLE\tWORKDIR\tC12_START\tC12_END\tMYH_START\tMYH_END

    A legacy leading #SAMPLE header is also accepted.
    """
    sample_map: Dict[str, Dict[str, Any]] = {}

    with samples_path.open("r", newline="") as fh:
        first_line = fh.readline()
        fh.seek(0)

        if "\t" in first_line:
            delimiter = "\t"
            delimiter_name = "TSV"
        elif "," in first_line:
            delimiter = ","
            delimiter_name = "CSV"
        else:
            raise ValueError(
                f"Could not detect delimiter in {samples_path}.\n"
                "Expected either comma-separated CSV or tab-separated TSV metadata."
            )

        reader = csv.reader(fh, delimiter=delimiter)
        rows = list(reader)

    if not rows:
        raise ValueError(f"Samples file is empty: {samples_path}")

    header = [h.strip() for h in rows[0]]
    expected = ["SAMPLE", "WORKDIR", "C12_START", "C12_END", "MYH_START", "MYH_END"]
    normalized = ["SAMPLE" if h == "#SAMPLE" else h for h in header[:6]]

    if normalized != expected:
        raise ValueError(
            f"Unexpected samples header in {samples_path}: {header}\n"
            f"Detected format: {delimiter_name}\n"
            f"Expected header: {expected}\n"
            "Accepted formats: CSV or TSV."
        )

    for line_no, row in enumerate(rows[1:], start=2):
        if not row or all(not cell.strip() for cell in row):
            continue
        if len(row) < 6:
            raise ValueError(
                f"Malformed row {line_no} in {samples_path}: {row}\n"
                f"Detected format: {delimiter_name}. "
                "Expected six fields: SAMPLE, WORKDIR, C12_START, C12_END, MYH_START, MYH_END."
            )

        taxon = row[0].strip()
        workdir = row[1].strip()

        try:
            c12_start = int(row[2])
            c12_end = int(row[3])
            myh_start = int(row[4])
            myh_end = int(row[5])
        except ValueError as exc:
            raise ValueError(
                f"Residue ranges must be integers at row {line_no} in {samples_path}: {row}"
            ) from exc

        if not taxon:
            raise ValueError(f"Empty SAMPLE value at row {line_no} in {samples_path}")
        if c12_start > c12_end or myh_start > myh_end:
            raise ValueError(f"Invalid residue range at row {line_no} in {samples_path}: {row}")

        sample_map[taxon] = {
            "ibam_chain": "A",
            "myht_chain": "B",
            "ibam_range": (c12_start, c12_end),
            "myht_range": (myh_start, myh_end),
            "workdir": workdir,
        }

    if not sample_map:
        raise ValueError(f"No sample rows found in {samples_path}")

    print(f"Loaded {len(sample_map)} samples from {samples_path} ({delimiter_name})")

    return sample_map

def choose_default_samples_path() -> Path:
    """Return the default samples metadata path."""
    return DEFAULT_SAMPLES_PATH


def validate_inputs(samples_path: Path, prco_dir: Path) -> Dict[str, Dict[str, Any]]:
    """Validate metadata and PRCO file availability before running the expensive analysis."""
    if not samples_path.exists():
        raise FileNotFoundError(f"Samples file not found: {samples_path}")
    if not prco_dir.exists():
        raise FileNotFoundError(f"PRCO directory not found: {prco_dir}")
    if not prco_dir.is_dir():
        raise NotADirectoryError(f"PRCO path is not a directory: {prco_dir}")

    sample_map = load_sample_map(samples_path)

    missing = []
    for taxon in sample_map:
        expected_prco = prco_dir / f"{taxon}_prco.csv"
        if not expected_prco.exists():
            missing.append(str(expected_prco))

    if missing:
        preview = "\n".join(missing[:10])
        more = "" if len(missing) <= 10 else f"\n...and {len(missing) - 10} more"
        raise FileNotFoundError(
            "Missing PRCO file(s) for samples listed in metadata:\n"
            f"{preview}{more}"
        )

    return sample_map


def make_rankings(df: pd.DataFrame) -> pd.DataFrame:
    """Create a compact ranking table from the main HARP summary."""
    if df.empty or "taxon" not in df.columns:
        return pd.DataFrame()

    ranking_cols = [
        "taxon",
        "full_best_phase",
        "full_second_phase",
        "full_best_score",
        "full_margin",
        "full_best_minus_flat",
        "full_p_best",
        "full_q_best",
        "full_p_margin",
        "full_q_margin",
        "full_null_stability_rate",
        "top10_best_phase",
        "top10_best_score",
        "phase_stability",
        "interpretation",
    ]
    cols = [c for c in ranking_cols if c in df.columns]
    out = df[cols].copy()

    if "full_best_score" in out.columns:
        out["rank_full_best_score"] = pd.to_numeric(out["full_best_score"], errors="coerce").rank(
            ascending=False, method="min"
        ).astype("Int64")
    if "full_best_minus_flat" in out.columns:
        out["rank_full_best_minus_flat"] = pd.to_numeric(out["full_best_minus_flat"], errors="coerce").rank(
            ascending=False, method="min"
        ).astype("Int64")
    if "full_margin" in out.columns:
        out["rank_full_margin"] = pd.to_numeric(out["full_margin"], errors="coerce").rank(
            ascending=False, method="min"
        ).astype("Int64")

    sort_col = "full_best_minus_flat" if "full_best_minus_flat" in out.columns else "full_best_score"
    if sort_col in out.columns:
        out = out.sort_values(sort_col, ascending=False, na_position="last")

    return out


# ---------------------------------------------------------------------
# Core HARP scoring
# ---------------------------------------------------------------------


def heptad_letter(resid: int, phase: int) -> str:
    """Map a residue index to a heptad letter under a given phase."""
    letters = "abcdefg"
    return letters[(resid - phase) % 7]



def build_weighted_profile(
    rows: List[Dict[str, Any]],
    myht_range: Tuple[int, int],
    partner_occ_cutoff: float = 0.0,
) -> Counter:
    """
    Build weighted MyhT residue profile from PRCO-decoded rows.

    Weight = IBAM occupancy * partner occupancy.
    """
    weighted: Counter = Counter()
    myht_start, myht_end = myht_range

    for r in rows:
        ibam_occ = r.get("occupancy")
        if ibam_occ is None:
            ibam_occ = 1.0

        for p in r["partners"]:
            if not p:
                continue

            resid = p.get("resid")
            p_occ = p.get("occ")

            if resid is None:
                continue
            if p_occ is None or p_occ < partner_occ_cutoff:
                continue
            if not (myht_start <= resid <= myht_end):
                continue

            weighted[int(resid)] += float(ibam_occ) * float(p_occ)

    return weighted



def phase_scan(profile: Counter) -> List[Tuple[int, float, Counter]]:
    """
    For each of the 7 phases, compute the weighted a/d fraction.

    Returns a ranked list:
      [(phase, ad_fraction, letter_profile), ...]
    sorted by descending ad_fraction.
    """
    phase_scores: List[Tuple[int, float, Counter]] = []

    for phase in range(7):
        total = 0.0
        ad = 0.0
        letter_profile: Counter = Counter()

        for resid, score in profile.items():
            letter = heptad_letter(int(resid), phase)
            letter_profile[letter] += score
            total += float(score)
            if letter in {"a", "d"}:
                ad += float(score)

        frac = ad / total if total else 0.0
        phase_scores.append((phase, frac, letter_profile))

    ranked = sorted(phase_scores, key=lambda x: x[1], reverse=True)
    return ranked



def summarize_profile(profile: Counter) -> Dict[str, Any]:
    """Summarize one weighted profile."""
    ranked = phase_scan(profile)

    best_phase, best_score, _ = ranked[0]
    second_phase, second_score, _ = ranked[1]
    margin = best_score - second_score
    flat_expectation = 2 / 7  # a/d out of 7 positions
    best_minus_flat = best_score - flat_expectation

    return {
        "n_nonzero_myht": len(profile),
        "best_phase": best_phase,
        "second_phase": second_phase,
        "best_score": round(best_score, 3),
        "margin": round(margin, 3),
        "best_minus_flat": round(best_minus_flat, 3),
        "phase_ranking": ",".join(f"{p}:{score:.3f}" for p, score, _ in ranked),
    }


# ---------------------------------------------------------------------
# Empirical null utilities
# ---------------------------------------------------------------------

NullMode = Literal["linear_shift", "circular_shift", "shuffle_weights"]


@dataclass
class HarpObservedResult:
    phase_scores: List[float]
    best_phase: int
    best_score: float
    second_best_score: float
    margin: float


@dataclass
class HarpNullSummary:
    null_mode: str
    n_null: int
    null_mean_best: float
    null_sd_best: float
    null_p95_best: float
    null_mean_margin: float
    null_sd_margin: float
    null_p95_margin: float
    p_best: float
    p_margin: float
    null_stability_rate: float



def profile_counter_to_dense(profile: Counter, myht_range: Tuple[int, int]) -> np.ndarray:
    """Convert sparse residue-keyed Counter to dense vector over the MyhT span."""
    myht_start, myht_end = myht_range
    positions = list(range(myht_start, myht_end + 1))
    return np.array([float(profile.get(resid, 0.0)) for resid in positions], dtype=float)



def dense_to_profile_counter(vec: Sequence[float], myht_range: Tuple[int, int]) -> Counter:
    """Convert dense vector back to residue-keyed Counter, dropping zero entries."""
    myht_start, myht_end = myht_range
    positions = range(myht_start, myht_end + 1)
    out: Counter = Counter()
    for resid, weight in zip(positions, vec):
        if float(weight) > 0.0:
            out[int(resid)] = float(weight)
    return out



def summarize_observed_dense(vec: Sequence[float], myht_range: Tuple[int, int]) -> HarpObservedResult:
    profile = dense_to_profile_counter(vec, myht_range)
    ranked = phase_scan(profile)
    best_phase, best_score, _ = ranked[0]
    second_best_phase, second_best_score, _ = ranked[1]
    return HarpObservedResult(
        phase_scores=[float(score) for _, score, _ in sorted(ranked, key=lambda x: x[0])],
        best_phase=int(best_phase),
        best_score=float(best_score),
        second_best_score=float(second_best_score),
        margin=float(best_score - second_best_score),
    )



def _occupied_span(vec: np.ndarray) -> Tuple[int, int]:
    nz = np.flatnonzero(vec > 0)
    if len(nz) == 0:
        return 0, max(0, len(vec) - 1)
    return int(nz[0]), int(nz[-1])



def generate_linear_shift_profiles(
    vec: Sequence[float],
    shift_on_occupied_block: bool = True,
    include_zero_shift: bool = False,
) -> List[np.ndarray]:
    arr = np.asarray(vec, dtype=float)
    length = len(arr)
    out: List[np.ndarray] = []

    if length == 0:
        return out

    if shift_on_occupied_block:
        start, end = _occupied_span(arr)
        block = arr[start : end + 1].copy()
        block_len = len(block)

        for new_start in range(0, length - block_len + 1):
            shifted = np.zeros(length, dtype=float)
            shifted[new_start : new_start + block_len] = block
            if not include_zero_shift and new_start == start:
                continue
            out.append(shifted)
    else:
        for shift in range(-(length - 1), length):
            if shift == 0 and not include_zero_shift:
                continue
            shifted = np.zeros(length, dtype=float)
            if shift > 0:
                shifted[shift:] = arr[: length - shift]
            else:
                k = -shift
                shifted[: length - k] = arr[k:]
            out.append(shifted)

    return out



def generate_circular_shift_profiles(vec: Sequence[float], include_zero_shift: bool = False) -> List[np.ndarray]:
    arr = np.asarray(vec, dtype=float)
    out: List[np.ndarray] = []
    for shift in range(len(arr)):
        if shift == 0 and not include_zero_shift:
            continue
        out.append(np.roll(arr, shift))
    return out



def generate_shuffle_weight_profiles(
    vec: Sequence[float],
    n_iter: int = 10000,
    seed: Optional[int] = None,
) -> List[np.ndarray]:
    arr = np.asarray(vec, dtype=float)
    rng = np.random.default_rng(seed)
    return [rng.permutation(arr) for _ in range(n_iter)]



def empirical_pvalue(null_values: Sequence[float], observed: float) -> float:
    arr = np.asarray(null_values, dtype=float)
    ge = int(np.sum(arr >= observed))
    return (ge + 1.0) / (len(arr) + 1.0)


# ---------------------------------------------------------------------
# Null phase stability test
# ---------------------------------------------------------------------

# Weight-percentile cutoffs applied to each null profile.
# These mirror the progressive contact removal of the real
# occupancy-cutoff stability test (0.00, 0.01, 0.05, 0.10, 0.20).
# Percentile 0 = no filtering (analogous to cutoff 0.00);
# higher percentiles progressively strip weaker weights.
NULL_STABILITY_PERCENTILES = [0, 10, 25, 50, 75]


def null_phase_stability_dense(
    vec: np.ndarray,
    myht_range: Tuple[int, int],
) -> bool:
    """
    Apply 5 weight-percentile cutoffs to a dense profile and check
    whether the best heptad phase is invariant across all 5.

    Returns True if phase is stable (same best phase at all cutoffs).
    """
    vec = np.asarray(vec, dtype=float)
    nonzero_vals = vec[vec > 0]

    if len(nonzero_vals) == 0:
        return True  # degenerate case, no contacts

    best_phases: List[int] = []

    for pctl in NULL_STABILITY_PERCENTILES:
        if pctl == 0:
            threshold = 0.0
        else:
            threshold = float(np.percentile(nonzero_vals, pctl))

        filtered = vec.copy()
        filtered[filtered < threshold] = 0.0

        profile = dense_to_profile_counter(filtered, myht_range)

        if not profile:
            # All weights zeroed at this cutoff — skip
            # (use the last valid phase to avoid penalizing extreme cutoffs)
            if best_phases:
                best_phases.append(best_phases[-1])
            continue

        ranked = phase_scan(profile)
        best_phases.append(ranked[0][0])

    if len(best_phases) < 2:
        return True  # not enough valid cutoffs to test

    return len(set(best_phases)) == 1


def evaluate_null(
    profile: Counter,
    myht_range: Tuple[int, int],
    null_mode: str = "shuffle_occupied_span",
    linear_shift_on_occupied_block: bool = True,
    shuffle_iter: int = 10000,
    seed: Optional[int] = None,
) -> Tuple[HarpObservedResult, HarpNullSummary]:
    dense = profile_counter_to_dense(profile, myht_range)
    obs = summarize_observed_dense(dense, myht_range)

    if null_mode == "linear_shift":
        null_profiles = generate_linear_shift_profiles(
            dense,
            shift_on_occupied_block=linear_shift_on_occupied_block,
            include_zero_shift=False,
        )
    elif null_mode == "circular_shift":
        null_profiles = generate_circular_shift_profiles(
            dense,
            include_zero_shift=False
        )
    elif null_mode == "shuffle_weights":
        null_profiles = generate_shuffle_weight_profiles(
            dense,
            n_iter=shuffle_iter,
            seed=seed
        )

    # --- new null family ---
    elif null_mode in {
        "shuffle_full_span",
        "shuffle_nonzero_support",
        "shuffle_occupied_span",
    }:
        rng = np.random.default_rng(seed)
        null_profiles = [
            permute_profile(dense, null_mode, rng)
            for _ in range(shuffle_iter)
        ]

    else:
        raise ValueError(f"Unknown null_mode: {null_mode}")

    if not null_profiles:
        raise ValueError("No null profiles were generated. Check profile length and shift settings.")

    null_best: List[float] = []
    null_margin: List[float] = []
    null_stable_count: int = 0

    for dense_prof in null_profiles:
        res = summarize_observed_dense(dense_prof, myht_range)
        null_best.append(res.best_score)
        null_margin.append(res.margin)

        # Phase stability test on this null replicate
        if null_phase_stability_dense(dense_prof, myht_range):
            null_stable_count += 1

    nb = np.asarray(null_best, dtype=float)
    nm = np.asarray(null_margin, dtype=float)

    null_stability_rate = float(null_stable_count) / len(null_profiles)

    summary = HarpNullSummary(
        null_mode=null_mode,
        n_null=len(null_profiles),
        null_mean_best=float(np.mean(nb)),
        null_sd_best=float(np.std(nb, ddof=1)) if len(nb) > 1 else 0.0,
        null_p95_best=float(np.quantile(nb, 0.95)),
        null_mean_margin=float(np.mean(nm)),
        null_sd_margin=float(np.std(nm, ddof=1)) if len(nm) > 1 else 0.0,
        null_p95_margin=float(np.quantile(nm, 0.95)),
        p_best=float(empirical_pvalue(nb, obs.best_score)),
        p_margin=float(empirical_pvalue(nm, obs.margin)),
        null_stability_rate=null_stability_rate,
    )
    return obs, summary




def benjamini_hochberg(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]

    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev

    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0.0, 1.0)
    return out



def phase_stability_from_cutoff_results(best_phases: Sequence[int]) -> float:
    vals = list(best_phases)
    if not vals:
        return float("nan")
    mode = max(set(vals), key=vals.count)
    return float(sum(1 for x in vals if x == mode) / len(vals))


# ---------------------------------------------------------------------
# Taxon summarization
# ---------------------------------------------------------------------


def classify_taxon(
    full_best_score: float,
    full_margin: float,
    cutoff_best_phases: List[int],
    top10_margin: float,
) -> str:
    """Crude first-pass classification for triage."""
    cutoff_stable = len(set(cutoff_best_phases)) == 1

    if full_best_score < 0.30:
        return "weak_or_flat"
    if cutoff_stable and full_margin >= 0.05 and top10_margin >= 0.03:
        return "stable_phase_biased"
    if full_margin > 0.0:
        return "competing_phases"
    return "ambiguous"



def summarize_taxon(
    prco_path: Path,
    meta: Dict[str, Any],
    null_mode: str,
    shuffle_iter: int,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Run HARP summary analyses for one taxon, including empirical null testing."""
    taxon = prco_path.name.replace("_prco.csv", "")

    rows = load_prco_table(
        prco_path,
        ibam_chain=meta["ibam_chain"],
        myht_chain=meta["myht_chain"],
    )

    # Full profile
    full_profile = build_weighted_profile(
        rows,
        myht_range=meta["myht_range"],
        partner_occ_cutoff=0.0,
    )
    full_summary = summarize_profile(full_profile)
    _, full_null = evaluate_null(
        full_profile,
        meta["myht_range"],
        null_mode=null_mode,
        linear_shift_on_occupied_block=LINEAR_SHIFT_ON_OCCUPIED_BLOCK,
        shuffle_iter=shuffle_iter,
        seed=seed,
    )

    # Top-10 profile
    top10_profile = Counter(dict(full_profile.most_common(10)))
    top10_summary = summarize_profile(top10_profile)

    # Occupancy-cutoff stability
    cutoffs = [0.0, 0.01, 0.05, 0.10, 0.20]
    cutoff_best_phases: List[int] = []

    for cutoff in cutoffs:
        profile = build_weighted_profile(
            rows,
            myht_range=meta["myht_range"],
            partner_occ_cutoff=cutoff,
        )
        ranked = phase_scan(profile)
        cutoff_best_phases.append(ranked[0][0] if ranked else -1)

    phase_stability = phase_stability_from_cutoff_results(cutoff_best_phases)

    interpretation = classify_taxon(
        full_best_score=full_summary["best_score"],
        full_margin=full_summary["margin"],
        cutoff_best_phases=cutoff_best_phases,
        top10_margin=top10_summary["margin"],
    )

    return {
        "taxon": taxon,
        "ibam_chain": meta["ibam_chain"],
        "ibam_range": f"{meta['ibam_range'][0]}-{meta['ibam_range'][1]}",
        "myht_chain": meta["myht_chain"],
        "myht_range": f"{meta['myht_range'][0]}-{meta['myht_range'][1]}",
        "n_nonzero_myht": full_summary["n_nonzero_myht"],
        "full_best_phase": full_summary["best_phase"],
        "full_second_phase": full_summary["second_phase"],
        "full_best_score": full_summary["best_score"],
        "full_margin": full_summary["margin"],
        "full_best_minus_flat": full_summary["best_minus_flat"],
        "full_phase_ranking": full_summary["phase_ranking"],
        "full_null_mode": full_null.null_mode,
        "full_n_null": full_null.n_null,
        "full_null_mean_best": round(full_null.null_mean_best, 3),
        "full_null_sd_best": round(full_null.null_sd_best, 3),
        "full_null_p95_best": round(full_null.null_p95_best, 3),
        "full_p_best": round(full_null.p_best, 6),
        "full_null_mean_margin": round(full_null.null_mean_margin, 3),
        "full_null_sd_margin": round(full_null.null_sd_margin, 3),
        "full_null_p95_margin": round(full_null.null_p95_margin, 3),
        "full_p_margin": round(full_null.p_margin, 6),
        "full_null_stability_rate": round(full_null.null_stability_rate, 6),
        "top10_best_phase": top10_summary["best_phase"],
        "top10_second_phase": top10_summary["second_phase"],
        "top10_best_score": top10_summary["best_score"],
        "top10_margin": top10_summary["margin"],
        "top10_best_minus_flat": top10_summary["best_minus_flat"],
        "top10_phase_ranking": top10_summary["phase_ranking"],
        "cutoff_best_phases": ",".join(map(str, cutoff_best_phases)),
        "cutoff_stable": len(set(cutoff_best_phases)) == 1,
        "phase_stability": round(phase_stability, 3),
        "interpretation": interpretation,
    }


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main(args) -> None:

    samples_path = Path(args.samples) if args.samples else choose_default_samples_path()
    prco_dir = Path(args.prco_dir)

    base_out_dir = Path(args.out_dir)
    run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out_dir = base_out_dir / run_timestamp
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.shuffle_iter < 1:
        raise ValueError("--shuffle-iter must be >= 1")

    if args.null_model == "all":
        null_models_to_run = NULL_MODELS
    else:
        null_models_to_run = [args.null_model]

    sample_map = validate_inputs(samples_path, prco_dir)
    prco_files = sorted(prco_dir.glob("*_prco.csv"))

    if not prco_files:
        raise FileNotFoundError(f"No PRCO files found in {prco_dir}")

    # Preserve the historical behavior of scanning PRCO_DIR, but warn about extras.
    prco_taxa = {p.name.replace("_prco.csv", "") for p in prco_files}
    extra_prco = sorted(prco_taxa - set(sample_map))

    print(f"Using samples file: {samples_path}")
    print(f"Scanning PRCO files in: {prco_dir}")
    print(f"Writing outputs to: {out_dir}")
    print(f"Requested null model(s): {null_models_to_run}")
    print(f"Samples in metadata: {len(sample_map)}")
    print(f"PRCO files found: {len(prco_files)}")
    if extra_prco:
        print(f"Warning: {len(extra_prco)} PRCO file(s) are not listed in samples metadata and will be skipped.")
    print()

    if args.dry_run:
        print("Dry run passed: samples metadata and required PRCO files are present.")
        return

    run_config = {
        "samples": str(samples_path),
        "prco_dir": str(prco_dir),
        "out_dir": str(out_dir),
        "null_model": args.null_model,
        "null_models_to_run": null_models_to_run,
        "shuffle_iter": args.shuffle_iter,
        "seed": args.seed,
        "out_prefix": args.out_prefix,
    }

    with (out_dir / "harp_run_config.json").open("w") as fh:
        json.dump(run_config, fh, indent=2)

    summaries: List[Dict[str, Any]] = []
    log_lines: List[str] = []

    for null_model in null_models_to_run:

        msg = f"=== Running null model: {null_model} ==="
        print(f"\n{msg}")
        log_lines.append(msg)

        for idx, prco_path in enumerate(prco_files):

            taxon = prco_path.name.replace("_prco.csv", "")

            if taxon not in sample_map:
                msg = f"{taxon}: SKIPPED - not found in samples metadata"
                print(msg)
                log_lines.append(msg)
                continue

            try:

                result = summarize_taxon(
                    prco_path,
                    sample_map[taxon],
                    null_mode=null_model,
                    shuffle_iter=args.shuffle_iter,
                    seed=args.seed + idx,
                )

                summaries.append(result)

                msg = (
                    f"{result['taxon']}: "
                    f"phase={result['full_best_phase']} "
                    f"(score {result['full_best_score']}, margin {result['full_margin']}, "
                    f"p={result['full_p_best']}, "
                    f"null_stab={result['full_null_stability_rate']:.3f})"
                )
                print(msg)
                log_lines.append(msg)

            except Exception as e:

                msg = f"{taxon}: ERROR - {e}"
                print(msg)
                log_lines.append(msg)

                summaries.append({
                    "taxon": taxon,
                    "interpretation": f"ERROR: {e}",
                })

    df = pd.DataFrame(summaries)

    if not df.empty and "full_p_best" in df.columns:

        valid_best = pd.to_numeric(df["full_p_best"], errors="coerce")
        valid_margin = pd.to_numeric(df["full_p_margin"], errors="coerce")

        df["full_q_best"] = ""
        df["full_q_margin"] = ""

        if valid_best.notna().any():
            df.loc[valid_best.notna(), "full_q_best"] = np.round(
                benjamini_hochberg(valid_best.dropna().values), 6
            )

        if valid_margin.notna().any():
            df.loc[valid_margin.notna(), "full_q_margin"] = np.round(
                benjamini_hochberg(valid_margin.dropna().values), 6
            )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_path = out_dir / f"{args.out_prefix}_{timestamp}.tsv"
    stable_summary_path = out_dir / "harp_panel_summary.tsv"
    rankings_path = out_dir / "harp_rankings.tsv"
    log_path = out_dir / "harp_run.log"

    df.to_csv(summary_path, sep="\t", index=False)
    df.to_csv(stable_summary_path, sep="\t", index=False)

    rankings = make_rankings(df)
    rankings.to_csv(rankings_path, sep="\t", index=False)

    with log_path.open("w") as fh:
        fh.write("\n".join(log_lines) + "\n")

    print(f"\nWrote timestamped summary TSV to: {summary_path}")
    print(f"Wrote stable summary TSV to:      {stable_summary_path}")
    print(f"Wrote rankings TSV to:            {rankings_path}")
    print(f"Wrote run config JSON to:         {out_dir / 'harp_run_config.json'}")
    print(f"Wrote run log to:                 {log_path}")

    # ---------------------------------------------------------
    # Generate manuscript plots
    # ---------------------------------------------------------
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    plot_script = PROJECT_ROOT / "scripts" / "generate_harp_plots.py"

    print(f"[DEBUG] Summary file: {stable_summary_path}")
    print("\n[INFO] Generating HARP plots...")

    subprocess.run(
        [
            sys.executable,
            str(plot_script),
            "--summary", str(stable_summary_path),
            "--out-dir", str(figures_dir),
            "--formats", "svg,png,pdf",
        ],
        check=True,
    )

    print("[INFO] Plot generation complete.\n")


if __name__ == "__main__":
    args = parse_args()
    main(args)
