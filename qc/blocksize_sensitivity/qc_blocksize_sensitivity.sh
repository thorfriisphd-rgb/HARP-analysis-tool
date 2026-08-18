#!/usr/bin/env bash
set -Eeuo pipefail

# =============================================================================
# HARP v4.1 — block-size sensitivity QC
# =============================================================================
#
# Purpose
# -------
# Release-QC sweep for the taxon-level contiguous-block null model.
# Re-runs the frozen 26-taxon HARP analysis de novo at block sizes 1..7
# (configurable), while keeping every QC config/output inside this QC folder.
#
# This script does NOT modify HARP statistical code and does NOT add a new
# production feature. It exercises the released CLI against the frozen reference
# corpus and reports how period7_max_mode inference changes with block_size.
# The two secondary taxon statistics are harvested as well because HARP already
# calculates them during every analysis.
#
# Expected location
# -----------------
#   <release-root>/qc/blocksize_sensitivity/qc_blocksize_sensitivity.sh
#
# Expected release files
# ----------------------
#   <release-root>/run_reference_26taxa.sh
#   <release-root>/src/harp/...
#   <release-root>/reference/26taxa/taxa.txt
#   <release-root>/reference/configs/26taxa/*.yaml
#
# Usage
# -----
#   chmod +x qc_blocksize_sensitivity.sh
#   ./qc_blocksize_sensitivity.sh all
#
# Optional environment variables
# ------------------------------
#   HARP_QC_BLOCK_SIZES="1 2 3 4 5 6 7"   # default
#   HARP_QC_N_PERMUTATIONS=9999             # default
#   HARP_QC_ANALYSIS_JOBS=1                  # default; increase cautiously
#   HARP_QC_ALPHA=0.05                       # reporting threshold only
#
# Output
# ------
#   runs/HARP_v4.1_blocksize_sensitivity_n26_<timestamp>/
#       block_01/ ... block_07/
#       sensitivity_results.tsv
#       period7_max_mode_summary.tsv
#       sensitivity_report.md
#       period7_max_mode_sensitivity.png
#       panel_invariance_signatures.tsv
#       panel_invariance.tsv
#       panel_invariance_report.md
#       run_metadata.txt
#       source_hashes.tsv
#       qc.log
#
# Notes
# -----
# block_size affects the taxon-level permutation null. It does not enter the
# phase_signature construction or panel circular-rotation statistic, so the
# panel is intentionally NOT rerun seven times here.
# =============================================================================

IFS=$'\n\t'
umask 002

MODE="${1:-all}"
case "${MODE}" in
    all) ;;
    -h|--help|help)
        sed -n '3,62p' "$0" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
    *)
        echo "ERROR: unknown command: ${MODE}" >&2
        exit 2
        ;;
esac

QC_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${QC_DIR}/../.." && pwd)"
REFERENCE_RUNNER="${ROOT}/run_reference_26taxa.sh"
REFERENCE_ROOT="${ROOT}/reference"
REFERENCE_RUN_CONFIG_ROOT="${REFERENCE_ROOT}/run_configs/26taxa"
TAXA_FILE="${REFERENCE_ROOT}/26taxa/taxa.txt"
RUNS_ROOT="${QC_DIR}/runs"

BLOCK_SIZES_RAW="${HARP_QC_BLOCK_SIZES:-1 2 3 4 5 6 7}"
N_PERMUTATIONS="${HARP_QC_N_PERMUTATIONS:-9999}"
ANALYSIS_JOBS="${HARP_QC_ANALYSIS_JOBS:-1}"
ALPHA="${HARP_QC_ALPHA:-0.05}"
REFERENCE_BLOCK_SIZE=4
EXPECTED_TAXA=26

# Normalize the user-supplied block-size list while preserving order.
IFS=' ' read -r -a BLOCK_SIZES <<< "${BLOCK_SIZES_RAW}"

stamp="$(date '+%Y%m%dT%H%M%S%z')"
RUN_DIR="${RUNS_ROOT}/HARP_v4.1_blocksize_sensitivity_n26_${stamp}"
mkdir -p "${RUN_DIR}"

exec > >(tee -a "${RUN_DIR}/qc.log") 2>&1

now() { date '+%Y-%m-%dT%H:%M:%S%z'; }
log() { printf '[%s] %s\n' "$(now)" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }

harp_cli() {
    if command -v harp >/dev/null 2>&1; then
        harp "$@"
    else
        PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
            python -m harp.cli "$@"
    fi
}

validate_parameters() {
    [[ -x "${REFERENCE_RUNNER}" ]] || die "Reference runner missing/not executable: ${REFERENCE_RUNNER}"
    [[ -f "${TAXA_FILE}" ]] || die "Missing taxon roster: ${TAXA_FILE}"
    [[ "${N_PERMUTATIONS}" =~ ^[0-9]+$ ]] && (( N_PERMUTATIONS >= 1 )) || \
        die "HARP_QC_N_PERMUTATIONS must be a positive integer"
    [[ "${ANALYSIS_JOBS}" =~ ^[0-9]+$ ]] && (( ANALYSIS_JOBS >= 1 )) || \
        die "HARP_QC_ANALYSIS_JOBS must be a positive integer"

    local b
    for b in "${BLOCK_SIZES[@]}"; do
        [[ "${b}" =~ ^[0-9]+$ ]] && (( b >= 1 )) || die "Invalid block size: ${b}"
    done

    python - "${ALPHA}" <<'PY'
import sys
x = float(sys.argv[1])
if not 0.0 < x < 1.0:
    raise SystemExit("HARP_QC_ALPHA must be between 0 and 1")
PY
}

load_taxa() {
    mapfile -t TAXA < <(
        sed -e 's/\r$//' -e '/^[[:space:]]*#/d' -e '/^[[:space:]]*$/d' "${TAXA_FILE}"
    )
    [[ ${#TAXA[@]} -eq ${EXPECTED_TAXA} ]] || \
        die "Expected ${EXPECTED_TAXA} taxa, found ${#TAXA[@]}"
}

record_provenance() {
    local out="${RUN_DIR}/source_hashes.tsv"
    printf 'sha256\tpath\n' > "${out}"
    local files=(
        "${ROOT}/run_reference_26taxa.sh"
        "${ROOT}/src/harp/analysis.py"
        "${ROOT}/src/harp/stats.py"
        "${ROOT}/src/harp/deepcoil.py"
        "${ROOT}/src/harp/trajectory.py"
    )
    local f h
    for f in "${files[@]}"; do
        if [[ -f "${f}" ]]; then
            h="$(sha256sum "${f}" | awk '{print $1}')"
            printf '%s\t%s\n' "${h}" "${f}" >> "${out}"
        fi
    done

    {
        echo "run_timestamp=${stamp}"
        echo "release_root=${ROOT}"
        echo "block_sizes=${BLOCK_SIZES[*]}"
        echo "reference_block_size=${REFERENCE_BLOCK_SIZE}"
        echo "n_permutations=${N_PERMUTATIONS}"
        echo "analysis_jobs=${ANALYSIS_JOBS}"
        echo "alpha_reporting_threshold=${ALPHA}"
        echo "python=$(python --version 2>&1)"
        if git -C "${ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
            echo "git_commit=$(git -C "${ROOT}" rev-parse HEAD)"
            echo "git_dirty=$(if [[ -n "$(git -C "${ROOT}" status --porcelain)" ]]; then echo yes; else echo no; fi)"
        else
            echo "git_commit=not-a-git-worktree"
        fi
    } > "${RUN_DIR}/run_metadata.txt"
}

prepare_reference_configs() {
    log "Running frozen reference CONFIGS gate (includes preflight/inventory audit)..."
    "${REFERENCE_RUNNER}" configs
    [[ -d "${REFERENCE_RUN_CONFIG_ROOT}" ]] || die "Reference run configs not generated"
}

make_qc_configs() {
    local block_size="$1"
    local block_dir="$2"
    local config_dir="${block_dir}/configs"
    local taxon_root="${block_dir}/taxa"
    mkdir -p "${config_dir}" "${taxon_root}" "${block_dir}/logs/validation" "${block_dir}/logs/analysis"

    python - \
        "${REFERENCE_RUN_CONFIG_ROOT}" \
        "${config_dir}" \
        "${taxon_root}" \
        "${TAXA_FILE}" \
        "${block_size}" \
        "${N_PERMUTATIONS}" <<'PY'
import copy
import sys
from pathlib import Path
import yaml

source_root = Path(sys.argv[1]).resolve()
out_root = Path(sys.argv[2]).resolve()
taxon_root = Path(sys.argv[3]).resolve()
taxa_file = Path(sys.argv[4])
block_size = int(sys.argv[5])
n_permutations = int(sys.argv[6])

taxa = [
    line.strip() for line in taxa_file.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]

for taxon in taxa:
    src = source_root / f"{taxon}.yaml"
    cfg = yaml.safe_load(src.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise SystemExit(f"{taxon}: generated reference config is not a mapping")

    out = copy.deepcopy(cfg)
    stats = out.setdefault("statistics", {})
    stats["block_size"] = block_size
    stats["n_permutations"] = n_permutations
    out["output_dir"] = str((taxon_root / taxon).resolve())

    dst = out_root / f"{taxon}.yaml"
    dst.write_text(yaml.safe_dump(out, sort_keys=False), encoding="utf-8")

print(f"QC configs: PASS ({len(taxa)} taxa; block_size={block_size})")
PY
}

validate_block() {
    local block_size="$1"
    local block_dir="$2"
    local config_dir="${block_dir}/configs"
    local log_dir="${block_dir}/logs/validation"
    local failures=0
    local taxon cfg logf report

    log "Validating block_size=${block_size} configs..."
    for taxon in "${TAXA[@]}"; do
        cfg="${config_dir}/${taxon}.yaml"
        logf="${log_dir}/${taxon}.log"
        report="${log_dir}/${taxon}_validation_report.json"
        if harp_cli validate --config "${cfg}" --report "${report}" > "${logf}" 2>&1; then
            echo "  PASS  ${taxon}"
        else
            echo "  FAIL  ${taxon}: see ${logf}"
            failures=$((failures + 1))
        fi
    done
    (( failures == 0 )) || die "block_size=${block_size}: ${failures} validation failure(s)"
}

analyze_block() {
    local block_size="$1"
    local block_dir="$2"
    local config_dir="${block_dir}/configs"
    local taxon_root="${block_dir}/taxa"
    local log_dir="${block_dir}/logs/analysis"
    local status_dir="${log_dir}/.status"
    rm -rf -- "${status_dir}"
    mkdir -p "${status_dir}"

    analyze_one_qc() {
        local taxon="$1"
        local cfg="${config_dir}/${taxon}.yaml"
        local logf="${log_dir}/${taxon}.log"
        local summary="${taxon_root}/${taxon}/harp_v4_summary.json"
        local status="${status_dir}/${taxon}.status"

        if harp_cli analyze --config "${cfg}" > "${logf}" 2>&1 && [[ -s "${summary}" ]]; then
            printf 'PASS\t%s\n' "${taxon}" > "${status}"
            return 0
        fi
        printf 'FAIL\t%s\t%s\n' "${taxon}" "${logf}" > "${status}"
        return 1
    }

    export -f analyze_one_qc harp_cli
    export config_dir log_dir taxon_root status_dir ROOT

    log "Analysing block_size=${block_size} with ${ANALYSIS_JOBS} worker(s)..."
    local running=0 taxon
    for taxon in "${TAXA[@]}"; do
        analyze_one_qc "${taxon}" &
        running=$((running + 1))
        if (( running >= ANALYSIS_JOBS )); then
            wait -n || true
            running=$((running - 1))
        fi
    done
    while (( running > 0 )); do
        wait -n || true
        running=$((running - 1))
    done

    local failures=0 status
    for taxon in "${TAXA[@]}"; do
        status="${status_dir}/${taxon}.status"
        if [[ -f "${status}" ]] && grep -q '^PASS' "${status}"; then
            echo "  PASS  ${taxon}"
        else
            echo "  FAIL  ${taxon}: see ${log_dir}/${taxon}.log"
            failures=$((failures + 1))
        fi
    done
    (( failures == 0 )) || die "block_size=${block_size}: ${failures} analysis failure(s)"
}

aggregate_report() {
    log "Aggregating sensitivity results and generating report..."

    python - \
        "${RUN_DIR}" \
        "${TAXA_FILE}" \
        "${ALPHA}" \
        "${REFERENCE_BLOCK_SIZE}" \
        "${BLOCK_SIZES[*]}" <<'PY'
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

run_dir = Path(sys.argv[1]).resolve()
taxa_file = Path(sys.argv[2])
alpha = float(sys.argv[3])
reference_block = int(sys.argv[4])
block_sizes = [int(x) for x in sys.argv[5].split()]

taxa = [
    line.strip() for line in taxa_file.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]

stats_order = ["period7_max_mode", "phase_nonuniformity", "target_phase_enrichment"]
rows = []
by_primary = {taxon: {} for taxon in taxa}

for b in block_sizes:
    bdir = run_dir / f"block_{b:02d}" / "taxa"
    for taxon in taxa:
        summary = bdir / taxon / "harp_v4_summary.json"
        if not summary.is_file():
            raise SystemExit(f"Missing summary: {summary}")
        obj = json.loads(summary.read_text(encoding="utf-8"))
        stat_block = obj.get("statistics", {})
        for name in stats_order:
            if name not in stat_block:
                raise SystemExit(f"{taxon} block={b}: missing statistic {name}")
            r = stat_block[name]
            row = {
                "taxon": taxon,
                "block_size": b,
                "statistic": name,
                "observed": float(r["observed"]),
                "p_value": float(r["p_value"]),
                "null_mean": float(r["null_mean"]),
                "null_sd": float(r["null_sd"]),
                "null_q95": float(r["null_q95"]),
                "n_permutations": int(r["n_permutations"]),
            }
            rows.append(row)
            if name == "period7_max_mode":
                by_primary[taxon][b] = row

long_path = run_dir / "sensitivity_results.tsv"
with long_path.open("w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# Observed statistic must be invariant to block_size; only the null changes.
observed_errors = []
summary_rows = []
flagged = []
for taxon in taxa:
    rmap = by_primary[taxon]
    missing = [b for b in block_sizes if b not in rmap]
    if missing:
        raise SystemExit(f"{taxon}: missing primary results for block sizes {missing}")

    obs = np.array([rmap[b]["observed"] for b in block_sizes], dtype=float)
    if not np.allclose(obs, obs[0], rtol=0.0, atol=1e-12):
        observed_errors.append((taxon, obs.min(), obs.max()))

    pvals = {b: rmap[b]["p_value"] for b in block_sizes}
    classes = {b: (pvals[b] <= alpha) for b in block_sizes}
    stable = len(set(classes.values())) == 1
    ref_p = pvals.get(reference_block, float("nan"))
    ref_class = "SIGNIFICANT" if ref_p <= alpha else "NOT_SIGNIFICANT"
    stability = "STABLE" if stable else "CROSSES_ALPHA"
    if not stable:
        flagged.append(taxon)

    row = {
        "taxon": taxon,
        "observed": obs[0],
        "p_reference_b4": ref_p,
        "reference_class": ref_class,
        "p_min": min(pvals.values()),
        "p_max": max(pvals.values()),
        "p_range": max(pvals.values()) - min(pvals.values()),
        "stability": stability,
    }
    for b in block_sizes:
        row[f"p_b{b}"] = pvals[b]
    summary_rows.append(row)

if observed_errors:
    lines = [f"{t}: {lo:.17g}..{hi:.17g}" for t, lo, hi in observed_errors]
    raise SystemExit("Observed period7_max_mode changed with block_size:\n" + "\n".join(lines))

fieldnames = [
    "taxon", "observed", "p_reference_b4", "reference_class",
    *[f"p_b{b}" for b in block_sizes],
    "p_min", "p_max", "p_range", "stability",
]
summary_path = run_dir / "period7_max_mode_summary.tsv"
with summary_path.open("w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(summary_rows)

overall = "ROBUST" if not flagged else "REVIEW"

# Markdown report.
report = []
report.append("# HARP v4.1 block-size sensitivity report")
report.append("")
report.append(f"**Overall QC verdict: {overall}**")
report.append("")
report.append("## Question")
report.append("")
report.append(
    "Does the primary taxon-level inference (`period7_max_mode`) materially depend "
    f"on the production contiguous-block size of {reference_block} residues?"
)
report.append("")
report.append("## Design")
report.append("")
report.append(
    f"The frozen {len(taxa)}-taxon reference analysis was rerun de novo at block sizes "
    + ", ".join(map(str, block_sizes)) + "."
)
report.append(
    "For each block size, HARP recalculated the full taxon analysis and its existing "
    "permutation statistics. No HARP statistical source code was modified."
)
report.append(
    "The QC verdict is descriptive: a taxon is `STABLE` when every tested block size "
    f"lies on the same side of the reporting threshold p = {alpha:g}; otherwise it is "
    "flagged `CROSSES_ALPHA` for review. This threshold is not a new HARP statistic."
)
report.append("")
report.append("The panel analysis was not rerun for every block size because `block_size` is used by the taxon-level block-shuffle permutation test, whereas the taxon `phase_signature` and panel circular-rotation statistic do not depend on the taxon block-shuffle size.")
report.append("")
report.append("## Summary")
report.append("")
report.append(f"- Taxa analysed: **{len(taxa)}**")
report.append(f"- Block sizes: **{', '.join(map(str, block_sizes))}**")
report.append(f"- Production reference block size: **{reference_block}**")
report.append(f"- Stable taxa: **{len(taxa) - len(flagged)} / {len(taxa)}**")
report.append(f"- Taxa crossing p = {alpha:g}: **{len(flagged)}**")
if flagged:
    report.append("- Flagged taxa: **" + ", ".join(flagged) + "**")
else:
    report.append("- No taxon changed inferential class across the tested block-size range.")
report.append("")
report.append("## Primary-statistic sensitivity")
report.append("")
header = ["Taxon", "Observed", "p(b=4)", "p min", "p max", "Status"]
report.append("| " + " | ".join(header) + " |")
report.append("|" + "|".join(["---"] * len(header)) + "|")
for row in summary_rows:
    report.append(
        f"| {row['taxon']} | {row['observed']:.6g} | {row['p_reference_b4']:.6g} | "
        f"{row['p_min']:.6g} | {row['p_max']:.6g} | {row['stability']} |"
    )
report.append("")
report.append("## Interpretation")
report.append("")
if overall == "ROBUST":
    report.append(
        "Across the tested range, the qualitative `period7_max_mode` inference is not "
        "dependent on the production choice of block size 4 at the stated reporting threshold."
    )
else:
    report.append(
        "At least one taxon changes inferential class across the tested block-size range. "
        "This does not invalidate the HARP null model, but the flagged taxon/taxa should be "
        "examined before claiming block-size robustness."
    )
report.append("")
report.append("Exact p-values and null summaries for all three existing taxon statistics are provided in `sensitivity_results.tsv`.")
report.append("")
report.append("## QC boundaries")
report.append("")
report.append(
    "This exercise evaluates sensitivity to one discretionary parameter of the existing taxon null. "
    "It does not establish that any particular block size is theoretically optimal, and it does not "
    "test phylogenetic independence or alter the panel null model."
)
(run_dir / "sensitivity_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

# Heatmap of primary p-values. Use -log10(p) so small p-values remain visible.
mat = np.array([[by_primary[t][b]["p_value"] for b in block_sizes] for t in taxa], dtype=float)
plot_mat = -np.log10(np.clip(mat, np.finfo(float).tiny, None))
fig_h = max(7.0, 0.34 * len(taxa) + 2.2)
fig, ax = plt.subplots(figsize=(9.0, fig_h))
im = ax.imshow(plot_mat, aspect="auto", interpolation="nearest")
ax.set_xticks(range(len(block_sizes)), [f"{b}" + (" (ref)" if b == reference_block else "") for b in block_sizes])
ax.set_yticks(range(len(taxa)), taxa)
ax.set_xlabel("Contiguous block size (residues)")
ax.set_ylabel("Taxon")
ax.set_title("HARP v4.1 period7_max_mode block-size sensitivity")
cb = fig.colorbar(im, ax=ax)
cb.set_label("-log10(permutation p)")
fig.tight_layout()
fig.savefig(run_dir / "period7_max_mode_sensitivity.png", dpi=220)
plt.close(fig)

print(f"Overall QC verdict: {overall}")
print(f"Stable taxa: {len(taxa) - len(flagged)}/{len(taxa)}")
print(f"Flagged taxa: {', '.join(flagged) if flagged else 'none'}")
print(f"Report: {run_dir / 'sensitivity_report.md'}")
PY
}

main() {
    validate_parameters
    load_taxa

    log "============================================================"
    log "HARP v4.1 — block-size sensitivity QC"
    log "Release root:       ${ROOT}"
    log "QC run:             ${RUN_DIR}"
    log "Block sizes:        ${BLOCK_SIZES[*]}"
    log "Permutations/test:  ${N_PERMUTATIONS}"
    log "Analysis workers:   ${ANALYSIS_JOBS}"
    log "Reference block:    ${REFERENCE_BLOCK_SIZE}"
    log "Reporting alpha:    ${ALPHA}"
    log "============================================================"

    record_provenance
    prepare_reference_configs

    local b block_dir
    for b in "${BLOCK_SIZES[@]}"; do
        block_dir="${RUN_DIR}/block_$(printf '%02d' "${b}")"
        log "------------------------------------------------------------"
        log "QC block_size=${b}"
        log "Output: ${block_dir}"
        make_qc_configs "${b}" "${block_dir}"
        validate_block "${b}" "${block_dir}"
        analyze_block "${b}" "${block_dir}"
    done

    aggregate_report
    log "------------------------------------------------------------"
    log "Running panel invariance QC"

    PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
    python "${ROOT}/qc/blocksize_sensitivity/check_panel_invariance.py" \
        "${RUN_DIR}"

    log "============================================================"
    log "HARP v4.1 block-size sensitivity QC: COMPLETE"
    log "Open: ${RUN_DIR}/sensitivity_report.md"
    log "Plot: ${RUN_DIR}/period7_max_mode_sensitivity.png"
    log "Panel invariance: ${RUN_DIR}/panel_invariance_report.md"
    log "============================================================"

}

main "$@"
