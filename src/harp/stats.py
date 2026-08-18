from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Sequence

import numpy as np

PHASES = tuple("abcdefg")


@dataclass(frozen=True)
class PermutationResult:
    observed: float
    p_value: float
    null_mean: float
    null_sd: float
    null_q95: float
    n_permutations: int

    def to_dict(self) -> dict:
        return asdict(self)


def _plus_one_p(observed: float, null: np.ndarray, *, upper: bool = True) -> float:
    null = np.asarray(null, dtype=float)
    exceed = np.sum(null >= observed) if upper else np.sum(null <= observed)
    return float((exceed + 1) / (null.size + 1))


def phase_sums(signal: np.ndarray, phase_index: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    signal = np.asarray(signal, dtype=float)
    phase_index = np.asarray(phase_index, dtype=float)
    valid = np.isfinite(phase_index) & np.isfinite(signal)
    phase_int = np.full(phase_index.shape, -1, dtype=int)
    phase_int[valid] = phase_index[valid].astype(int)
    sums = np.zeros(7, dtype=float)
    counts = np.zeros(7, dtype=float)
    for k in range(7):
        mask = valid & (phase_int == k)
        sums[k] = signal[mask].sum()
        counts[k] = mask.sum()
    return sums, counts


def phase_fraction(signal: np.ndarray, phase_index: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sums, counts = phase_sums(signal, phase_index)
    total = sums.sum()
    fractions = sums / total if total > 0 else np.zeros(7, dtype=float)
    baseline = counts / counts.sum() if counts.sum() > 0 else np.full(7, 1 / 7)
    return fractions, baseline


def phase_enrichment(signal: np.ndarray, phase_index: np.ndarray) -> np.ndarray:
    fractions, baseline = phase_fraction(signal, phase_index)
    with np.errstate(divide="ignore", invalid="ignore"):
        enrichment = np.divide(fractions, baseline, out=np.zeros_like(fractions), where=baseline > 0)
    return enrichment


def phase_nonuniformity(signal: np.ndarray, phase_index: np.ndarray) -> float:
    """Pearson-style concentration statistic against residue-count baseline."""
    sums, counts = phase_sums(signal, phase_index)
    if sums.sum() <= 0 or counts.sum() <= 0:
        return 0.0
    expected = sums.sum() * counts / counts.sum()
    return float(np.sum((sums - expected) ** 2 / np.maximum(expected, 1e-12)))


def target_enrichment(
    signal: np.ndarray,
    phase_index: np.ndarray,
    target_phases: Sequence[str] = ("a", "d"),
) -> float:
    fractions, baseline = phase_fraction(signal, phase_index)
    idx = np.array([PHASES.index(p) for p in target_phases], dtype=int)
    obs = fractions[idx].sum()
    exp = baseline[idx].sum()
    return float(obs / exp) if exp > 0 else 0.0


def period7_amplitude(signal: np.ndarray, seq_pos: np.ndarray, valid_mask: np.ndarray | None = None) -> float:
    """Normalized complex amplitude at exactly one cycle per seven residues."""
    signal = np.asarray(signal, dtype=float)
    seq_pos = np.asarray(seq_pos, dtype=float)
    valid = np.isfinite(signal) & np.isfinite(seq_pos)
    if valid_mask is not None:
        valid &= np.asarray(valid_mask, dtype=bool)
    x = signal[valid]
    p = seq_pos[valid]
    if x.size == 0 or x.sum() <= 0:
        return 0.0
    z = np.sum(x * np.exp(-2j * np.pi * p / 7.0))
    return float(np.abs(z) / np.sum(np.abs(x)))


def seven_bin_modes(signal: np.ndarray, phase_index: np.ndarray) -> np.ndarray:
    """DFT amplitudes m=1,2,3 of the seven phase bins, normalized to total signal."""
    sums, _ = phase_sums(signal, phase_index)
    if sums.sum() <= 0:
        return np.zeros(3)
    k = np.arange(7)
    return np.array([
        np.abs(np.sum(sums * np.exp(-2j * np.pi * m * k / 7.0))) / sums.sum()
        for m in (1, 2, 3)
    ])



def max_seven_bin_mode(signal: np.ndarray, phase_index: np.ndarray) -> float:
    """Largest of heptad-bin Fourier modes 1–3.

    The maximization is performed identically for every null permutation, so
    the p-value is corrected for selecting the strongest heptad harmonic.
    """
    return float(np.max(seven_bin_modes(signal, phase_index)))

def centered_phase_signature(signal: np.ndarray, phase_index: np.ndarray) -> np.ndarray:
    enrich = phase_enrichment(signal, phase_index)
    vec = enrich - np.mean(enrich)
    norm = np.linalg.norm(vec)
    return vec / norm if norm > 0 else vec


def template_similarity(signal: np.ndarray, phase_index: np.ndarray, template: np.ndarray) -> float:
    vec = centered_phase_signature(signal, phase_index)
    template = np.asarray(template, dtype=float)
    template = template - template.mean()
    norm = np.linalg.norm(template)
    if norm <= 0 or np.linalg.norm(vec) <= 0:
        return 0.0
    return float(np.dot(vec, template / norm))


def block_shuffle(values: np.ndarray, block_size: int, rng: np.random.Generator) -> np.ndarray:
    values = np.asarray(values)
    if block_size < 1:
        raise ValueError("block_size must be >= 1")
    blocks = [values[i : i + block_size] for i in range(0, len(values), block_size)]
    order = rng.permutation(len(blocks))
    return np.concatenate([blocks[i] for i in order])


def permutation_test(
    signal: np.ndarray,
    scorer,
    *,
    n_permutations: int = 9999,
    block_size: int = 4,
    seed: int = 20260801,
    shuffle_mask: np.ndarray | None = None,
) -> tuple[PermutationResult, np.ndarray]:
    signal = np.asarray(signal, dtype=float)
    observed = float(scorer(signal))
    rng = np.random.default_rng(seed)
    if shuffle_mask is None:
        shuffle_mask = np.ones(signal.shape, dtype=bool)
    else:
        shuffle_mask = np.asarray(shuffle_mask, dtype=bool)
        if shuffle_mask.shape != signal.shape:
            raise ValueError("shuffle_mask must match signal shape")
    null = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        permuted = signal.copy()
        permuted[shuffle_mask] = block_shuffle(signal[shuffle_mask], block_size, rng)
        null[i] = scorer(permuted)
    result = PermutationResult(
        observed=observed,
        p_value=_plus_one_p(observed, null),
        null_mean=float(np.mean(null)),
        null_sd=float(np.std(null, ddof=1)) if len(null) > 1 else 0.0,
        null_q95=float(np.quantile(null, 0.95)),
        n_permutations=int(n_permutations),
    )
    return result, null


def time_block_phase_vectors(
    frame_signal: np.ndarray,
    phase_index: np.ndarray,
    n_blocks: int = 20,
) -> tuple[np.ndarray, np.ndarray]:
    frame_signal = np.asarray(frame_signal, dtype=float)
    if frame_signal.ndim != 2:
        raise ValueError("frame_signal must have shape (n_frames, n_residues)")
    n_blocks = max(1, min(int(n_blocks), frame_signal.shape[0]))
    frame_indices = np.array_split(np.arange(frame_signal.shape[0]), n_blocks)
    vectors = []
    signatures = []
    for idx in frame_indices:
        profile = frame_signal[idx].mean(axis=0)
        fractions, _ = phase_fraction(profile, phase_index)
        vectors.append(fractions)
        signatures.append(centered_phase_signature(profile, phase_index))
    return np.asarray(vectors), np.asarray(signatures)


def temporal_signature_consistency(signatures: np.ndarray) -> float:
    signatures = np.asarray(signatures, dtype=float)
    if signatures.ndim != 2 or signatures.shape[0] < 2:
        return 1.0
    sims = signatures @ signatures.T
    upper = sims[np.triu_indices_from(sims, k=1)]
    return float(np.mean(upper)) if upper.size else 1.0


def panel_leave_one_out_test(
    signatures: np.ndarray,
    *,
    n_permutations: int = 9999,
    seed: int = 20260801,
) -> tuple[PermutationResult, np.ndarray, np.ndarray]:
    """Test shared a–g phase alignment across taxa.

    Each row is a centered, unit-length seven-phase signature. The null
    independently rotates each taxon's signature, preserving its strength and
    shape while breaking common phase registration.
    """
    signatures = np.asarray(signatures, dtype=float)
    if signatures.ndim != 2 or signatures.shape[1] != 7 or signatures.shape[0] < 3:
        raise ValueError("Panel test requires at least 3 taxa with seven-phase signatures")

    def loo_scores(arr: np.ndarray) -> np.ndarray:
        out = np.zeros(arr.shape[0], dtype=float)
        for i in range(arr.shape[0]):
            template = np.mean(np.delete(arr, i, axis=0), axis=0)
            tnorm = np.linalg.norm(template)
            vnorm = np.linalg.norm(arr[i])
            out[i] = np.dot(arr[i], template) / (vnorm * tnorm) if tnorm > 0 and vnorm > 0 else 0.0
        return out

    observed_taxon = loo_scores(signatures)
    observed = float(np.mean(observed_taxon))
    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for j in range(n_permutations):
        rotated = np.vstack([np.roll(row, int(rng.integers(0, 7))) for row in signatures])
        null[j] = np.mean(loo_scores(rotated))
    result = PermutationResult(
        observed=observed,
        p_value=_plus_one_p(observed, null),
        null_mean=float(np.mean(null)),
        null_sd=float(np.std(null, ddof=1)),
        null_q95=float(np.quantile(null, 0.95)),
        n_permutations=int(n_permutations),
    )
    return result, observed_taxon, null
