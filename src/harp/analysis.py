from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .deepcoil import PHASES, parse_deepcoil2
from .plotting import plot_contact_profile, plot_null, plot_panel_null, plot_phase_summary, plot_time_blocks
from .stats import (
    centered_phase_signature,
    panel_leave_one_out_test,
    period7_amplitude,
    max_seven_bin_mode,
    permutation_test,
    phase_fraction,
    phase_nonuniformity,
    seven_bin_modes,
    target_enrichment,
    template_similarity,
    temporal_signature_consistency,
    time_block_phase_vectors,
)
from .trajectory import compute_contact_trajectory
from .validation import validate_config
from .provenance import build_provenance, write_provenance, panel_output_name, timestamp_now

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","HSD":"H","HSE":"H","HSP":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M",
    "PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def _jsonable(obj: Any):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(type(obj).__name__)


def analyze_profile(
    *,
    deepcoil_path: str | Path,
    profile: np.ndarray,
    outdir: str | Path,
    frame_signal: np.ndarray | None = None,
    topology_resids: np.ndarray | None = None,
    topology_resnames: np.ndarray | None = None,
    cc_threshold: float = 0.5,
    anchor_threshold: float = 0.5,
    n_permutations: int = 9999,
    block_size: int = 4,
    seed: int = 20260801,
    n_time_blocks: int = 20,
    target_phases: tuple[str, ...] = ("a", "d"),
    locked_template: np.ndarray | None = None,
    primary_test: str = "period7_max_mode",
    trajectory_params: dict | None = None,
    taxon: str | None = None,
) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    reg = parse_deepcoil2(
        deepcoil_path,
        cc_threshold=cc_threshold,
        anchor_threshold=anchor_threshold,
    )
    profile = np.asarray(profile, dtype=float)
    if len(profile) != reg.length:
        raise ValueError(f"Profile length {len(profile)} != DeepCoil sequence length {reg.length}")

    table = reg.table.copy()
    table["contact_occupancy"] = profile
    if topology_resids is not None:
        if len(topology_resids) != reg.length:
            raise ValueError(
                f"topology_resids has {len(topology_resids)} entries but "
                f"DeepCoil sequence has {reg.length}"
            )
        table["topology_resid"] = np.asarray(topology_resids)
    if topology_resnames is not None:
        if len(topology_resnames) != reg.length:
            raise ValueError(
                f"topology_resnames has {len(topology_resnames)} entries but "
                f"DeepCoil sequence has {reg.length}"
            )
        table["topology_resname"] = np.asarray(topology_resnames)
        seq_top = "".join(AA3_TO_1.get(str(x).upper(), "X") for x in topology_resnames)
        mismatch = [i + 1 for i, (a, b) in enumerate(zip(reg.sequence, seq_top)) if b != "X" and a != b]
        if mismatch:
            raise ValueError(
                f"DeepCoil/topology sequence mismatch at {len(mismatch)} positions; "
                f"first positions: {mismatch[:10]}"
            )
    table.to_csv(outdir / "per_residue_contacts.tsv", sep="\t", index=False)
    reg.write_tsv(outdir / "deepcoil_register.tsv")

    phase_idx = table["phase_index"].to_numpy(float)
    seq_pos = table["seq_pos"].to_numpy(int)
    cc_mask = table["is_cc_segment"].to_numpy(bool)
    fractions, baseline = phase_fraction(profile, phase_idx)
    sums = np.array([profile[np.isfinite(phase_idx) & (phase_idx == k)].sum() for k in range(7)])
    counts = np.array([np.sum(np.isfinite(phase_idx) & (phase_idx == k)) for k in range(7)])
    phase_df = pd.DataFrame({
        "phase": PHASES,
        "n_residues": counts,
        "contact_sum": sums,
        "contact_fraction": fractions,
        "baseline_fraction": baseline,
        "enrichment": np.divide(fractions, baseline, out=np.zeros(7), where=baseline > 0),
    })
    phase_df.to_csv(outdir / "phase_summary.tsv", sep="\t", index=False)

    tests = {}
    nulls = {}

    tests["period7_max_mode"], nulls["period7_max_mode"] = permutation_test(
        profile,
        lambda x: max_seven_bin_mode(x, phase_idx),
        n_permutations=n_permutations,
        block_size=block_size,
        seed=seed,
        shuffle_mask=cc_mask,
    )
    tests["phase_nonuniformity"], nulls["phase_nonuniformity"] = permutation_test(
        profile,
        lambda x: phase_nonuniformity(x, phase_idx),
        n_permutations=n_permutations,
        block_size=block_size,
        seed=seed + 1,
        shuffle_mask=cc_mask,
    )
    tests["target_phase_enrichment"], nulls["target_phase_enrichment"] = permutation_test(
        profile,
        lambda x: target_enrichment(x, phase_idx, target_phases),
        n_permutations=n_permutations,
        block_size=block_size,
        seed=seed + 2,
        shuffle_mask=cc_mask,
    )
    if locked_template is not None:
        locked_template = np.asarray(locked_template, dtype=float)
        if locked_template.shape != (7,):
            raise ValueError("locked_template must contain exactly seven a–g values")
        tests["locked_template_similarity"], nulls["locked_template_similarity"] = permutation_test(
            profile,
            lambda x: template_similarity(x, phase_idx, locked_template),
            n_permutations=n_permutations,
            block_size=block_size,
            seed=seed + 3,
            shuffle_mask=cc_mask,
        )

    allowed_primary = set(tests)
    if primary_test not in allowed_primary:
        raise ValueError(f"primary_test {primary_test!r} is unavailable; choose one of {sorted(allowed_primary)}")

    block_vectors = None
    temporal_consistency = None
    if frame_signal is not None:
        frame_signal = np.asarray(frame_signal, dtype=float)
        if frame_signal.shape[1] != reg.length:
            raise ValueError("frame_signal residue dimension does not match DeepCoil length")
        block_vectors, block_signatures = time_block_phase_vectors(frame_signal, phase_idx, n_time_blocks)
        temporal_consistency = temporal_signature_consistency(block_signatures)
        pd.DataFrame(block_vectors, columns=PHASES).to_csv(
            outdir / "time_block_phase_fractions.tsv", sep="\t", index_label="block"
        )

    signature = centered_phase_signature(profile, phase_idx)
    summary = {
        "harp_version": "4.1",
        "taxon": taxon,
        "status": "HARP v4.1 validated execution; interpretation remains bounded by the documented statistical model",
        "deepcoil": {
            "sequence_length": reg.length,
            "cc_segment": [reg.cc_start, reg.cc_end],
            "a_origin": reg.a_origin,
            "d_offset": reg.d_offset,
            "cc_threshold": cc_threshold,
            "anchor_threshold": anchor_threshold,
        },
        "register_quality": {
            "a_purity": reg.a_purity,
            "d_purity": reg.d_purity,
            "combined_purity": reg.register_purity,
            "n_a_calls": reg.n_a_calls,
            "n_d_calls": reg.n_d_calls,
            "a_modal_tie": reg.a_tie,
            "d_modal_tie": reg.d_tie,
            "cc_segment_length": int(reg.cc_end - reg.cc_start + 1),
            "note": (
                "Purity below 1.0 indicates that confident DeepCoil a/d calls "
                "do not all support a single mod-7 frame across the selected "
                "CC segment. This is a register-quality diagnostic and does "
                "not by itself identify the cause of the disagreement."
            ),
        },
        "trajectory": trajectory_params or {"recorded": False},
        "primary_test": primary_test,
        "primary_result": tests[primary_test].to_dict(),
        "statistics": {
            name: result.to_dict() for name, result in tests.items()
        },
        "descriptive": {
            "phase_fractions": dict(zip(PHASES, fractions.tolist())),
            "phase_enrichment": dict(zip(PHASES, phase_df["enrichment"].tolist())),
            "fundamental_sequence_amplitude_m1": period7_amplitude(profile, seq_pos, cc_mask),
            "seven_bin_fourier_modes_m1_m2_m3": seven_bin_modes(profile, phase_idx).tolist(),
            "period7_max_mode": max_seven_bin_mode(profile, phase_idx),
            "phase_signature": signature.tolist(),
            "temporal_signature_consistency": temporal_consistency,
            "target_phases": list(target_phases),
            "locked_template": locked_template.tolist() if locked_template is not None else None,
        },
        "null_model": {
            "type": "contiguous block shuffle of the residue contact profile",
            "block_size": block_size,
            "n_permutations": n_permutations,
            "plus_one_p": True,
        },
    }
    with open(outdir / "harp_v4_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, default=_jsonable)

    np.savez_compressed(outdir / "null_distributions.npz", **nulls)
    plot_contact_profile(table, outdir / "01_contact_profile.png", taxon=taxon)
    plot_phase_summary(phase_df, outdir / "02_phase_summary.png", taxon=taxon)
    if block_vectors is not None:
        plot_time_blocks(block_vectors, outdir / "03_time_block_phase_heatmap.png", taxon=taxon)
    null_order = ["period7_max_mode", "phase_nonuniformity", "target_phase_enrichment"]
    if "locked_template_similarity" in tests:
        null_order.append("locked_template_similarity")
    for i, name in enumerate(null_order, start=4):
        plot_null(nulls[name], tests[name].observed, f"null — {name}", outdir / f"{i:02d}_null_{name}.png", taxon=taxon)
    return summary


def run_from_config(config_path: str | Path) -> dict:
    config_path = Path(config_path).resolve()
    validation = validate_config(config_path)
    if not validation.passed:
        codes = ", ".join(issue.code for issue in validation.issues)
        raise ValueError(
            "HARP v4.1 validation failed; analysis was not started. "
            f"Failure codes: {codes}"
        )
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    base = config_path.parent
    resolve = lambda p: str((base / p).resolve()) if p and not Path(p).is_absolute() else p

    dc = cfg["deepcoil"]
    tr = cfg["trajectory"]
    st = cfg.get("statistics", {})
    taxon = cfg.get("taxon")
    outdir = resolve(cfg.get("output_dir", "harp_v4_output"))
    locked_template = st.get("locked_template")
    if locked_template is None and st.get("locked_template_file"):
        template_path = Path(resolve(st["locked_template_file"]))
        with open(template_path, "r", encoding="utf-8") as fh:
            template_obj = json.load(fh)
        locked_template = template_obj.get("phase_signature", template_obj.get("template", template_obj))
    deepcoil_path = resolve(dc["file"])
    register = parse_deepcoil2(
        deepcoil_path,
        cc_threshold=float(dc.get("cc_threshold", 0.5)),
        anchor_threshold=float(dc.get("anchor_threshold", 0.5)),
    )

    result = compute_contact_trajectory(
        topology=resolve(tr["topology"]),
        trajectory=resolve(tr["trajectory"]),
        mg_selection=tr["mg_selection"],
        myht_selection=tr["myht_selection"],
        cutoff_angstrom=float(tr.get("cutoff_angstrom", 4.5)),
        contact_mode=tr.get("contact_mode", "hard"),
        smooth_power=int(tr.get("smooth_power", 6)),
        start=tr.get("start"),
        stop=tr.get("stop"),
        step=int(tr.get("step", 1)),
    )
    if len(result.myht_resids) != register.length:
        raise ValueError(
            f"MyhT selection has {len(result.myht_resids)} residues but DeepCoil has {register.length}. "
            "Use a selection containing exactly the sequence submitted to DeepCoil."
        )

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    validation.write_json(out / "validation_report.json")
    provenance_files = [
        config_path,
        deepcoil_path,
        resolve(tr["topology"]),
        resolve(tr["trajectory"]),
    ]
    fasta_value = cfg.get("myht_fasta") or cfg.get("fasta") or dc.get("fasta")
    if fasta_value:
        provenance_files.append(resolve(fasta_value))
    write_provenance(
        out / "provenance.json",
        build_provenance(
            files=provenance_files,
            repo=config_path.parent,
            extra={
                "analysis_type": "per_taxon",
                "config": str(config_path),
                "validation_status": validation.status,
            },
        ),
    )
    np.savez_compressed(
        out / "trajectory_contacts.npz",
        frame_times_ps=result.frame_times_ps,
        myht_frame_signal=result.myht_frame_signal,
        pair_occupancy=result.pair_occupancy,
        myht_resids=result.myht_resids,
        myht_resnames=result.myht_resnames,
        mg_resids=result.mg_resids,
        mg_resnames=result.mg_resnames,
    )
    pair_df = pd.DataFrame(
        result.pair_occupancy,
        index=[f"MG_{resid}_{resname}" for resid, resname in zip(result.mg_resids, result.mg_resnames)],
        columns=[f"MYHT_{resid}_{resname}" for resid, resname in zip(result.myht_resids, result.myht_resnames)],
    )
    pair_df.to_csv(out / "mg_myht_pair_occupancy.tsv", sep="\t", index_label="mg_residue")
    with open(out / "resolved_config.yaml", "w", encoding="utf-8") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    trajectory_params = {
        "recorded": True,
        "cutoff_angstrom": float(tr.get("cutoff_angstrom", 4.5)),
        "contact_mode": tr.get("contact_mode", "hard"),
        "smooth_power": int(tr.get("smooth_power", 6)),
        "start": tr.get("start"),
        "stop": tr.get("stop"),
        "step": int(tr.get("step", 1)),
        "n_frames": int(result.myht_frame_signal.shape[0]),
        "first_frame_ps": float(result.frame_times_ps[0]),
        "last_frame_ps": float(result.frame_times_ps[-1]),
        "n_myht_residues": int(len(result.myht_resids)),
        "n_mg_residues": int(len(result.mg_resids)),
        "myht_sequence": "".join(
            AA3_TO_1.get(str(x).upper(), "X") for x in result.myht_resnames
        ),
    }
    return analyze_profile(
        deepcoil_path=deepcoil_path,
        profile=result.myht_frame_signal.mean(axis=0),
        frame_signal=result.myht_frame_signal,
        topology_resids=result.myht_resids,
        topology_resnames=result.myht_resnames,
        outdir=out,
        cc_threshold=float(dc.get("cc_threshold", 0.5)),
        anchor_threshold=float(dc.get("anchor_threshold", 0.5)),
        n_permutations=int(st.get("n_permutations", 9999)),
        block_size=int(st.get("block_size", 4)),
        seed=int(st.get("seed", 20260801)),
        n_time_blocks=int(st.get("n_time_blocks", 20)),
        target_phases=tuple(st.get("target_phases", ["a", "d"])),
        locked_template=np.asarray(locked_template, dtype=float) if locked_template is not None else None,
        primary_test=st.get("primary_test", "period7_max_mode"),
        trajectory_params=trajectory_params,
        taxon=taxon,
    )


def run_panel(
    manifest_path: str | Path,
    outdir: str | Path,
    n_permutations: int,
    seed: int,
    *,
    direct_outdir: bool = False,
    run_name_override: str | None = None,
) -> dict:



    manifest_path = Path(manifest_path).resolve()
    manifest = pd.read_csv(manifest_path)
    if not {"taxon", "summary_json"}.issubset(manifest.columns):
        raise ValueError("Panel manifest must contain columns: taxon,summary_json")
    if manifest.empty:
        raise ValueError("Panel manifest contains no taxa")
    if manifest["taxon"].duplicated().any():
        dup = manifest.loc[manifest["taxon"].duplicated(), "taxon"].astype(str).tolist()
        raise ValueError(f"Panel manifest contains duplicate taxa: {dup}")

    signatures = []
    resolved_summaries = []
    versions = []
    param_rows = []
    quality_rows = []
    for _, row in manifest.iterrows():
        p = Path(row["summary_json"])
        if not p.is_absolute():
            p = manifest_path.parent / p
        p = p.resolve()
        if not p.is_file():
            raise ValueError(f"Missing panel summary for taxon {row['taxon']!r}: {p}")
        with open(p, "r", encoding="utf-8") as fh:
            summary = json.load(fh)
        signature = np.asarray(summary["descriptive"]["phase_signature"], dtype=float)
        if signature.shape != (7,) or not np.all(np.isfinite(signature)):
            raise ValueError(
                f"Taxon {row['taxon']!r} has invalid phase_signature; "
                "expected seven finite values"
            )
        norm = float(np.linalg.norm(signature))
        if not np.isclose(norm, 1.0, atol=1e-8):
            raise ValueError(
                f"Taxon {row['taxon']!r} has a degenerate phase_signature "
                f"(norm={norm:.3e}, expected 1.0). A zero-norm signature is "
                "finite and would silently dilute the panel statistic. This "
                "usually means the contact profile is empty or constant."
            )
        signatures.append(signature)
        resolved_summaries.append(p)
        version = str(summary.get("harp_version", "unknown"))
        versions.append(version)
        dc_block = summary.get("deepcoil", {})
        tr_block = summary.get("trajectory", {})
        param_rows.append({
            "taxon": str(row["taxon"]),
            "harp_version": version,
            "cc_threshold": dc_block.get("cc_threshold"),
            "anchor_threshold": dc_block.get("anchor_threshold"),
            "cutoff_angstrom": tr_block.get("cutoff_angstrom"),
            "contact_mode": tr_block.get("contact_mode"),
            "smooth_power": tr_block.get("smooth_power"),
            "step": tr_block.get("step"),
            "start": tr_block.get("start"),
            "stop": tr_block.get("stop"),
            "n_frames": tr_block.get("n_frames"),
            "first_frame_ps": tr_block.get("first_frame_ps"),
            "last_frame_ps": tr_block.get("last_frame_ps"),
            "trajectory_recorded": bool(tr_block.get("recorded", False)),
        })
        quality_rows.append({
            "taxon": str(row["taxon"]),
            **{k: summary.get("register_quality", {}).get(k)
               for k in ("a_purity", "d_purity", "combined_purity",
                         "n_a_calls", "n_d_calls", "a_modal_tie", "d_modal_tie",
                         "cc_segment_length")},
        })

    param_df = pd.DataFrame(param_rows)
    quality_df = pd.DataFrame(quality_rows)

    distinct_versions = sorted(set(versions))
    if len(distinct_versions) > 1:
        raise ValueError(
            f"Panel mixes HARP versions: {distinct_versions}. "
            "Regenerate all taxa under a single version before running the panel."
        )

    homogeneity_fields = [
        "cc_threshold", "anchor_threshold", "cutoff_angstrom",
        "contact_mode", "smooth_power", "step"
    ]
    homogeneity_verified = bool(param_df["trajectory_recorded"].all())
    disagreements = {}
    for field in homogeneity_fields:
        values = param_df[field].dropna().unique().tolist()
        if len(values) > 1:
            disagreements[field] = values
    if disagreements:
        raise ValueError(
            "Panel taxa were analysed under different core parameters: "
            f"{disagreements}. Contact and register parameters must be "
            "homogeneous for a frozen panel analysis."
        )

    signatures = np.asarray(signatures, dtype=float)
    result, taxon_scores, null = panel_leave_one_out_test(
        signatures, n_permutations=n_permutations, seed=seed
    )

    results_root = Path(outdir).resolve()
    results_root.mkdir(parents=True, exist_ok=True)
    run_time = timestamp_now()
    n_taxa = int(len(manifest))

    if direct_outdir:
        run_dir = results_root
        run_name = run_name_override or run_dir.name
    else:
        run_name = panel_output_name(n_taxa, run_time)
        run_dir = results_root / run_name
        suffix = 0
        while run_dir.exists():
            suffix += 1
            run_dir = results_root / f"{run_name}_{suffix:02d}"
        run_dir.mkdir(parents=True, exist_ok=False)
        run_name = run_dir.name

    manifest_used = manifest.copy()
    manifest_used["resolved_summary_json"] = [str(p) for p in resolved_summaries]
    manifest_used.to_csv(run_dir / "panel_manifest_used.tsv", sep="\t", index=False)

    score_df = manifest[["taxon"]].copy()
    score_df["leave_one_out_similarity"] = taxon_scores
    score_df.to_csv(run_dir / "panel_taxon_scores.tsv", sep="\t", index=False)
    param_df.to_csv(run_dir / "panel_parameter_audit.tsv", sep="\t", index=False)
    quality_df.to_csv(run_dir / "panel_register_quality.tsv", sep="\t", index=False)

    summary = {
        "harp_version": "4.1",
        "harp_version_panel": distinct_versions[0],
        "homogeneity_verified": homogeneity_verified,
        "homogeneity_fields_checked": homogeneity_fields,
        "run_timestamp": run_time.isoformat(),
        "run_name": run_name,
        "test": "panel shared phase alignment by leave-one-out consensus",
        "null": "independent circular rotation of each taxon's seven-phase signature",
        "result": result.to_dict(),
        "n_taxa": n_taxa,
        "n_permutations": int(n_permutations),
        "seed": int(seed),
        "manifest": str(manifest_path),
        "output_dir": str(run_dir),
    }

    np.save(run_dir / "panel_null.npy", null)
    plot_panel_null(
        null,
        result.observed,
        p_value=result.p_value,
        null_mean=result.null_mean,
        null_sd=result.null_sd,
        null_q95=result.null_q95,
        n_taxa=n_taxa,
        n_permutations=n_permutations,
        title="HARP v4.1 panel null — shared phase alignment",
        path=run_dir / "panel_null.png",
    )

    write_provenance(
        run_dir / "provenance.json",
        build_provenance(
            files=[manifest_path, *resolved_summaries],
            repo=manifest_path.parent,
            extra={
                "analysis_type": "panel",
                "run_name": run_name,
                "n_taxa": n_taxa,
                "n_permutations": int(n_permutations),
                "seed": int(seed),
            },
        ),
    )

    primary_summary = run_dir / f"{run_name}_summary.json"
    summary["primary_output"] = str(primary_summary)
    with open(primary_summary, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    # Stable convenience name retained for scripts and backward compatibility.
    with open(run_dir / "panel_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    return summary
