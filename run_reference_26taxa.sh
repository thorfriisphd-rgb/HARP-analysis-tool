#!/usr/bin/env bash
set -Eeuo pipefail

# =============================================================================
# HARP v4.1 — GitHub 26-taxon reference orchestrator
# =============================================================================
#
# Purpose
# -------
# Reproduce the frozen 26-taxon HARP v4.1 reference analysis while keeping
# generic HARP free of reference-corpus assumptions.
#
# This runner:
#   1. verifies the 26-taxon scientific input corpus against SHA-256 inventory
#   2. audits each reference trajectory with MDAnalysis
#      (1001 frames, 0–10000 ps, 10 ps spacing)
#   3. creates portable run configs from authoritative per-taxon YAML configs
#      while preserving scientific selections/parameters
#   4. validates all 26 taxa
#   5. analyses all 26 taxa
#   6. auto-generates reference/panel_manifest_26taxa.csv
#   7. runs the 26-taxon panel and generates the panel money-shot figure
#   8. re-verifies the frozen SHA-256 inventory after the run
#
# Reference-only constraints such as 1001 frames / 10 ns belong HERE, not in
# generic HARP validation. Generic HARP remains able to analyse mixed-duration
# trajectories.
#
# Expected release/data layout
# ----------------------------
#   <release-root>/
#   ├── src/harp/...
#   ├── reference/
#   │   ├── 26taxa/
#   │   │   ├── taxa.txt
#   │   │   └── HARP_v4.1_26taxa_input_sha256.tsv
#   │   ├── configs/
#   │   │   └── 26taxa/<Taxon>.yaml       # authoritative scientific configs
#   │   ├── run_configs/
#   │   │   └── 26taxa/<Taxon>.yaml       # generated; paths rewritten only
#   │   ├── results/
#   │   │   └── 26taxa/
#   │   │       ├── taxa/<Taxon>/...
#   │   │       └── panels/HARP_v4.1_panel_n26_<timestamp>/...
#   │   └── panel_manifest_26taxa.csv      # generated automatically
#   └── run_reference_26taxa.sh
#
#   <external-data-root>/26taxa/
#       └── <Taxon>/{md.tpr,md.xtc,myht.fa,myht.out}
#
# Set HARP_REFERENCE_DATA_ROOT to the external 26taxa/ directory.
#
# Usage
# -----
#   ./run_reference_26taxa.sh preflight
#   ./run_reference_26taxa.sh configs
#   ./run_reference_26taxa.sh validate
#   ./run_reference_26taxa.sh analyze
#   ./run_reference_26taxa.sh panel
#   ./run_reference_26taxa.sh all
#
# =============================================================================

IFS=$'\n\t'
umask 002

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REFERENCE_ROOT="${ROOT}/reference"
CORPUS_ROOT="${HARP_REFERENCE_DATA_ROOT:-${REFERENCE_ROOT}/26taxa}"
SOURCE_CONFIG_ROOT="${REFERENCE_ROOT}/configs/26taxa"
RUN_CONFIG_ROOT="${REFERENCE_ROOT}/run_configs/26taxa"

RESULTS_ROOT="${ROOT}/results"
LEGACY_RESULT_ROOT="${REFERENCE_ROOT}/results/26taxa"

LOG_ROOT="${REFERENCE_ROOT}/logs/26taxa"
VALIDATION_LOG_ROOT="${LOG_ROOT}/validation"
ANALYSIS_LOG_ROOT="${LOG_ROOT}/analysis"
PANEL_LOG_ROOT="${LOG_ROOT}/panel"

# Canonical release metadata lives with the GitHub package.
# CORPUS_ROOT may point to the separately distributed Zenodo 26-taxon dataset.
TAXA_FILE="${REFERENCE_ROOT}/26taxa/taxa.txt"
INVENTORY="${REFERENCE_ROOT}/26taxa/HARP_v4.1_26taxa_input_sha256.tsv"
# Assigned after MODE is known.

N_PERMUTATIONS="${HARP_N_PERMUTATIONS:-9999}"
PANEL_SEED="${HARP_PANEL_SEED:-20260801}"
HARP_PREFLIGHT_JOBS="${HARP_PREFLIGHT_JOBS:-4}"
HARP_ANALYSIS_JOBS="${HARP_ANALYSIS_JOBS:-1}"

# Frozen HARP v4.1 n=26 regression benchmark.
# These values are release metadata only; they do not feed the analysis.
FROZEN_PANEL_OBSERVED=0.645128594598065
FROZEN_PANEL_P_VALUE=0.0001
FROZEN_PANEL_NULL_MEAN=-0.000521984561169
FROZEN_PANEL_NULL_SD=0.151023083942446
FROZEN_PANEL_NULL_Q95=0.225049488444680
FROZEN_PANEL_N_PERMUTATIONS=9999
FROZEN_PANEL_SEED=20260801
FROZEN_PANEL_ABS_TOL=1e-12

EXPECTED_TAXA=26
EXPECTED_FILES_PER_TAXON=4
EXPECTED_INVENTORY_ROWS=$((EXPECTED_TAXA * EXPECTED_FILES_PER_TAXON))
EXPECTED_FRAMES=1001
EXPECTED_FIRST_PS=0.0
EXPECTED_LAST_PS=10000.0
EXPECTED_DT_PS=10.0
TIME_TOL_PS=0.001

REQUIRED_INPUTS=(md.tpr md.xtc myht.fa myht.out)

MODE="${1:-}"
case "${MODE}" in
    preflight|configs|validate|analyze|panel|all) ;;
    -h|--help|"")
        sed -n '3,62p' "$0" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
    *)
        echo "ERROR: unknown command: ${MODE}" >&2
        exit 2
        ;;
esac

RUN_TIMESTAMP=""
RUN_NAME=""
RUN_ROOT=""

if [[ "${MODE}" == "all" ]]; then
    RUN_TIMESTAMP="$(date '+%Y%m%dT%H%M%S%z')"
    RUN_NAME="HARP_v4.1_panel_n${EXPECTED_TAXA}_${RUN_TIMESTAMP}"
    RUN_ROOT="${RESULTS_ROOT}/${RUN_NAME}"

    RESULT_ROOT="${RUN_ROOT}"
    TAXON_RESULT_ROOT="${RUN_ROOT}/taxa"
    PANEL_RESULT_ROOT="${RUN_ROOT}/panel"
    PANEL_MANIFEST="${PANEL_RESULT_ROOT}/panel_manifest_26taxa.csv"
else
    RESULT_ROOT="${LEGACY_RESULT_ROOT}"
    TAXON_RESULT_ROOT="${RESULT_ROOT}/taxa"
    PANEL_RESULT_ROOT="${RESULT_ROOT}/panels"
    PANEL_MANIFEST="${REFERENCE_ROOT}/panel_manifest_26taxa.csv"
fi



timestamp() {
    date '+%Y-%m-%dT%H:%M:%S%z'
}

log_line() {
    printf '[%s] %s\n' "$(timestamp)" "$*"
}

die() {
    log_line "ERROR: $*" >&2
    exit 1
}

ensure_dirs() {
    mkdir -p \
        "${SOURCE_CONFIG_ROOT}" \
        "${RUN_CONFIG_ROOT}" \
        "${TAXON_RESULT_ROOT}" \
        "${PANEL_RESULT_ROOT}" \
        "${VALIDATION_LOG_ROOT}" \
        "${ANALYSIS_LOG_ROOT}" \
        "${PANEL_LOG_ROOT}"
}

harp_cli() {
    if command -v harp >/dev/null 2>&1; then
        harp "$@"
    else
        PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
            python -m harp.cli "$@"
    fi
}

load_taxa() {
    [[ -f "${TAXA_FILE}" ]] || die "Missing canonical taxon roster: ${TAXA_FILE}"

    mapfile -t TAXA < <(
        sed \
            -e 's/\r$//' \
            -e '/^[[:space:]]*#/d' \
            -e '/^[[:space:]]*$/d' \
            "${TAXA_FILE}"
    )

    [[ ${#TAXA[@]} -eq ${EXPECTED_TAXA} ]] || die \
        "taxa.txt contains ${#TAXA[@]} taxa; expected ${EXPECTED_TAXA}"

    local unique
    unique="$(printf '%s\n' "${TAXA[@]}" | sort -u | wc -l)"
    [[ "${unique}" -eq ${EXPECTED_TAXA} ]] || die \
        "taxa.txt contains duplicate taxon names"
}

check_environment() {
    log_line "Checking Python/HARP environment..."

    command -v python >/dev/null 2>&1 || die "python not found"

    python - <<'PY'
import sys
if sys.version_info < (3, 11):
    raise SystemExit(f"Python >=3.11 required; found {sys.version.split()[0]}")

mods = ["numpy", "pandas", "yaml", "MDAnalysis", "matplotlib"]
for name in mods:
    mod = __import__(name)
    print(f"{name}: {getattr(mod, '__version__', 'unknown')}")
print(f"Python: {sys.version.split()[0]}")
PY

    if command -v harp >/dev/null 2>&1; then
        log_line "HARP CLI: $(command -v harp)"
    elif [[ -f "${ROOT}/src/harp/cli.py" ]]; then
        log_line "HARP CLI: using local src/ package"
    else
        die "Neither installed 'harp' command nor ${ROOT}/src/harp/cli.py found"
    fi
}

check_corpus_structure() {
    log_line "Checking 26-taxon corpus structure..."
    [[ -d "${CORPUS_ROOT}" ]] || die "Missing corpus root: ${CORPUS_ROOT}"

    local actual_dirs=()
    mapfile -t actual_dirs < <(
        find "${CORPUS_ROOT}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort
    )

    [[ ${#actual_dirs[@]} -eq ${EXPECTED_TAXA} ]] || die \
        "Corpus has ${#actual_dirs[@]} taxon directories; expected ${EXPECTED_TAXA}"

    local expected_sorted actual_sorted
    expected_sorted="$(printf '%s\n' "${TAXA[@]}" | sort)"
    actual_sorted="$(printf '%s\n' "${actual_dirs[@]}" | sort)"
    [[ "${actual_sorted}" == "${expected_sorted}" ]] || {
        echo "Expected taxa:" >&2
        printf '  %s\n' "${TAXA[@]}" >&2
        echo "Actual corpus directories:" >&2
        printf '  %s\n' "${actual_dirs[@]}" >&2
        die "Corpus taxon roster differs from taxa.txt"
    }

    local failures=0
    for taxon in "${TAXA[@]}"; do
        for f in "${REQUIRED_INPUTS[@]}"; do
            if [[ ! -f "${CORPUS_ROOT}/${taxon}/${f}" ]]; then
                echo "  MISSING ${taxon}/${f}"
                failures=$((failures + 1))
            fi
        done
    done

    [[ ${failures} -eq 0 ]] || die \
        "${failures} required scientific input file(s) missing"

    log_line "Corpus structure: PASS (${EXPECTED_TAXA} taxa, ${EXPECTED_INVENTORY_ROWS} required files)"
}

verify_inventory() {
    log_line "Verifying SHA-256 inventory..."
    [[ -f "${INVENTORY}" ]] || die "Missing SHA-256 inventory: ${INVENTORY}"

    python - "${INVENTORY}" "${CORPUS_ROOT}" "${TAXA_FILE}" <<'PY'
import csv
import hashlib
import sys
from pathlib import Path

inventory = Path(sys.argv[1])
corpus = Path(sys.argv[2])
taxa_file = Path(sys.argv[3])

taxa = [
    line.strip()
    for line in taxa_file.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]
required = ("md.tpr", "md.xtc", "myht.fa", "myht.out")

with inventory.open(encoding="utf-8", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

expected_n = len(taxa) * len(required)
if len(rows) != expected_n:
    raise SystemExit(f"Inventory has {len(rows)} rows; expected {expected_n}")

expected_pairs = {(t, f) for t in taxa for f in required}
seen_pairs = {(r["taxon"], r["file"]) for r in rows}

missing = sorted(expected_pairs - seen_pairs)
extra = sorted(seen_pairs - expected_pairs)
if missing or extra:
    raise SystemExit(f"Inventory roster mismatch; missing={missing}, extra={extra}")

failures = []
for row in rows:
    p = corpus / row["taxon"] / row["file"]
    if not p.is_file():
        failures.append(f"MISSING {p}")
        continue

    expected_size = int(row["size_bytes"])
    got_size = p.stat().st_size
    if got_size != expected_size:
        failures.append(f"SIZE {p}: {got_size} != {expected_size}")
        continue

    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    got_hash = h.hexdigest()
    if got_hash != row["sha256"]:
        failures.append(f"HASH {p}: {got_hash} != {row['sha256']}")

if failures:
    print("\n".join(failures), file=sys.stderr)
    raise SystemExit(1)

print(f"SHA-256 inventory: PASS ({len(rows)} files)")
PY
}

audit_reference_trajectories() {
    log_line "Auditing reference trajectories with MDAnalysis..."
    python - \
        "${CORPUS_ROOT}" \
        "${TAXA_FILE}" \
        "${EXPECTED_FRAMES}" \
        "${EXPECTED_FIRST_PS}" \
        "${EXPECTED_LAST_PS}" \
        "${EXPECTED_DT_PS}" \
        "${TIME_TOL_PS}" <<'PY'
import sys
from pathlib import Path

import numpy as np
import MDAnalysis as mda

corpus = Path(sys.argv[1])
taxa_file = Path(sys.argv[2])
expected_frames = int(sys.argv[3])
expected_first = float(sys.argv[4])
expected_last = float(sys.argv[5])
expected_dt = float(sys.argv[6])
tol = float(sys.argv[7])

taxa = [
    line.strip()
    for line in taxa_file.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]

failures = []

for taxon in taxa:
    top = corpus / taxon / "md.tpr"
    xtc = corpus / taxon / "md.xtc"
    try:
        u = mda.Universe(str(top), str(xtc))
        n = len(u.trajectory)
        times = np.array([float(ts.time) for ts in u.trajectory], dtype=float)

        problems = []
        if n != expected_frames:
            problems.append(f"frames={n}, expected={expected_frames}")
        if times.size != n:
            problems.append(f"times={times.size}, frames={n}")
        if times.size:
            if not np.all(np.isfinite(times)):
                problems.append("non-finite frame time(s)")
            if abs(times[0] - expected_first) > tol:
                problems.append(f"first_ps={times[0]:.6f}, expected={expected_first:.6f}")
            if abs(times[-1] - expected_last) > tol:
                problems.append(f"last_ps={times[-1]:.6f}, expected={expected_last:.6f}")
            if times.size > 1:
                diffs = np.diff(times)
                if not np.all(diffs > 0):
                    problems.append("frame times are not strictly increasing")
                if not np.allclose(diffs, expected_dt, rtol=0.0, atol=tol):
                    problems.append(
                        f"frame spacing differs from {expected_dt:g} ps "
                        f"(min={diffs.min():.6f}, max={diffs.max():.6f})"
                    )

        if problems:
            failures.append(f"{taxon}: " + "; ".join(problems))
            print(f"  FAIL  {taxon}: " + "; ".join(problems))
        else:
            print(
                f"  PASS  {taxon}: {n} frames, "
                f"{times[0]:.1f}–{times[-1]:.1f} ps, dt={expected_dt:g} ps"
            )

    except Exception as exc:
        failures.append(f"{taxon}: {type(exc).__name__}: {exc}")
        print(f"  FAIL  {taxon}: {type(exc).__name__}: {exc}")

if failures:
    raise SystemExit(
        f"\nReference trajectory audit FAILED for {len(failures)} taxon/taxa:\n"
        + "\n".join(failures)
    )

print(f"Reference trajectory audit: PASS ({len(taxa)} taxa)")
PY
}

check_source_configs() {
    log_line "Checking authoritative reference configs..."
    local failures=0

    for taxon in "${TAXA[@]}"; do
        if [[ ! -f "${SOURCE_CONFIG_ROOT}/${taxon}.yaml" ]]; then
            echo "  MISSING ${SOURCE_CONFIG_ROOT}/${taxon}.yaml"
            failures=$((failures + 1))
        fi
    done

    [[ ${failures} -eq 0 ]] || die \
        "${failures} authoritative taxon config(s) missing"

    log_line "Authoritative reference configs: PASS (${EXPECTED_TAXA})"
}

create_run_configs() {
    log_line "Creating portable 26-taxon run configs..."
    rm -rf -- "${RUN_CONFIG_ROOT}"
    mkdir -p "${RUN_CONFIG_ROOT}" "${TAXON_RESULT_ROOT}"

    python - \
        "${SOURCE_CONFIG_ROOT}" \
        "${RUN_CONFIG_ROOT}" \
        "${CORPUS_ROOT}" \
        "${TAXON_RESULT_ROOT}" \
        "${TAXA_FILE}" <<'PY'
import copy
import sys
from pathlib import Path

import yaml

source_root = Path(sys.argv[1]).resolve()
run_root = Path(sys.argv[2]).resolve()
corpus_root = Path(sys.argv[3]).resolve()
result_root = Path(sys.argv[4]).resolve()
taxa_file = Path(sys.argv[5])

taxa = [
    line.strip()
    for line in taxa_file.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]

trajectory_scientific_fields = (
    "mg_selection", "myht_selection", "cutoff_angstrom",
    "contact_mode", "smooth_power", "start", "stop", "step",
)
deepcoil_scientific_fields = ("cc_threshold", "anchor_threshold")

for taxon in taxa:
    source_path = source_root / f"{taxon}.yaml"
    cfg = yaml.safe_load(source_path.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise SystemExit(f"{taxon}: source config is not a YAML mapping")

    deepcoil = copy.deepcopy(cfg.get("deepcoil", {}))
    trajectory = copy.deepcopy(cfg.get("trajectory", {}))
    statistics = copy.deepcopy(cfg.get("statistics", {}))

    for required in ("mg_selection", "myht_selection"):
        if not trajectory.get(required):
            raise SystemExit(
                f"{taxon}: authoritative config lacks trajectory.{required}"
            )

    # Filesystem paths are the ONLY fields rewritten here.
    deepcoil.pop("path", None)
    deepcoil["file"] = str((corpus_root / taxon / "myht.out").resolve())
    trajectory["topology"] = str((corpus_root / taxon / "md.tpr").resolve())
    trajectory["trajectory"] = str((corpus_root / taxon / "md.xtc").resolve())

    generated = {
        "taxon": taxon,
        "myht_fasta": str((corpus_root / taxon / "myht.fa").resolve()),
        "deepcoil": deepcoil,
        "trajectory": trajectory,
        "statistics": statistics,
        "output_dir": str((result_root / taxon).resolve()),
    }

    out_path = run_root / f"{taxon}.yaml"
    out_path.write_text(
        yaml.safe_dump(generated, sort_keys=False),
        encoding="utf-8",
    )

    # Audit that scientific settings were preserved.
    old_dc = cfg.get("deepcoil", {})
    old_tr = cfg.get("trajectory", {})
    for field in deepcoil_scientific_fields:
        if old_dc.get(field) != generated["deepcoil"].get(field):
            raise SystemExit(
                f"{taxon}: deepcoil.{field} changed during path rewrite"
            )
    for field in trajectory_scientific_fields:
        if old_tr.get(field) != generated["trajectory"].get(field):
            raise SystemExit(
                f"{taxon}: trajectory.{field} changed during path rewrite"
            )
    if cfg.get("statistics", {}) != generated["statistics"]:
        raise SystemExit(f"{taxon}: statistics block changed during path rewrite")

    print(f"  PASS  {taxon}")

print(f"Generated and audited {len(taxa)} run configs")
PY
}

validate_all() {
    log_line "Validating all 26 taxa with ${HARP_PREFLIGHT_JOBS} worker(s)..."

    local status_dir="${VALIDATION_LOG_ROOT}/.status"
    rm -rf -- "${status_dir}"
    mkdir -p "${status_dir}"

    validate_one() {
        local taxon="$1"
        local cfg="${RUN_CONFIG_ROOT}/${taxon}.yaml"
        local log="${VALIDATION_LOG_ROOT}/${taxon}.log"
        local report="${VALIDATION_LOG_ROOT}/${taxon}_validation_report.json"
        local status="${status_dir}/${taxon}.status"

        if [[ ! -f "${cfg}" ]]; then
            printf 'FAIL\t%s\tmissing config\n' "${taxon}" > "${status}"
            return 1
        fi

        log_line "VALIDATE ${taxon}"

        if harp_cli validate \
            --config "${cfg}" \
            --report "${report}" \
            > "${log}" 2>&1
        then
            printf 'PASS\t%s\n' "${taxon}" > "${status}"
            return 0
        else
            printf 'FAIL\t%s\t%s\n' "${taxon}" "${log}" > "${status}"
            return 1
        fi
    }

    export -f validate_one harp_cli log_line
    export RUN_CONFIG_ROOT VALIDATION_LOG_ROOT status_dir ROOT

    local running=0
    local taxon

    for taxon in "${TAXA[@]}"; do
        validate_one "${taxon}" &
        running=$((running + 1))

        if (( running >= HARP_PREFLIGHT_JOBS )); then
            wait -n || true
            running=$((running - 1))
        fi
    done

    while (( running > 0 )); do
        wait -n || true
        running=$((running - 1))
    done

    local failures=0
    for taxon in "${TAXA[@]}"; do
        local status="${status_dir}/${taxon}.status"

        if [[ ! -f "${status}" ]]; then
            echo "  FAIL  ${taxon}: validation status missing"
            failures=$((failures + 1))
            continue
        fi

        if grep -q '^PASS' "${status}"; then
            echo "  PASS  ${taxon}"
        else
            echo "  FAIL  ${taxon}: see ${VALIDATION_LOG_ROOT}/${taxon}.log"
            failures=$((failures + 1))
        fi
    done

    [[ ${failures} -eq 0 ]] || die \
        "${failures} validation(s) failed; analysis will not start"

    log_line "All 26 taxon validations: PASS"
}

clean_taxon_results() {
    log_line "Cleaning generated taxon results only..."
    rm -rf -- "${TAXON_RESULT_ROOT}"
    mkdir -p "${TAXON_RESULT_ROOT}"
}

analyze_all() {
    log_line "Running HARP v4.1 analysis for all 26 taxa..."
    clean_taxon_results

    local failures=0
    for taxon in "${TAXA[@]}"; do
        local cfg="${RUN_CONFIG_ROOT}/${taxon}.yaml"
        local log="${ANALYSIS_LOG_ROOT}/${taxon}.log"

        log_line "ANALYZE ${taxon}"
        if harp_cli analyze --config "${cfg}" > "${log}" 2>&1; then
            local summary="${TAXON_RESULT_ROOT}/${taxon}/harp_v4_summary.json"
            if [[ -s "${summary}" ]]; then
                echo "  PASS  ${taxon}"
            else
                echo "  FAIL  ${taxon}: summary missing after successful command"
                failures=$((failures + 1))
            fi
        else
            echo "  FAIL  ${taxon}: see ${log}"
            failures=$((failures + 1))
        fi
    done

    [[ ${failures} -eq 0 ]] || die \
        "${failures} taxon analysis/analyses failed"

    log_line "All 26 taxon analyses: PASS"
}

audit_analysis_trajectory_metadata() {
    log_line "Auditing trajectory metadata recorded by HARP outputs..."

    python - \
        "${TAXON_RESULT_ROOT}" \
        "${TAXA_FILE}" \
        "${EXPECTED_FRAMES}" \
        "${EXPECTED_FIRST_PS}" \
        "${EXPECTED_LAST_PS}" \
        "${TIME_TOL_PS}" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
taxa_file = Path(sys.argv[2])
expected_frames = int(sys.argv[3])
expected_first = float(sys.argv[4])
expected_last = float(sys.argv[5])
tol = float(sys.argv[6])

taxa = [
    line.strip()
    for line in taxa_file.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]

errors = []
for taxon in taxa:
    p = root / taxon / "harp_v4_summary.json"
    if not p.is_file():
        errors.append(f"{taxon}: missing {p}")
        continue

    obj = json.loads(p.read_text(encoding="utf-8"))
    tr = obj.get("trajectory", {})
    n = tr.get("n_frames")
    first = tr.get("first_frame_ps")
    last = tr.get("last_frame_ps")

    problems = []
    if n != expected_frames:
        problems.append(f"n_frames={n}, expected={expected_frames}")
    if first is None or abs(float(first) - expected_first) > tol:
        problems.append(f"first_frame_ps={first}, expected={expected_first}")
    if last is None or abs(float(last) - expected_last) > tol:
        problems.append(f"last_frame_ps={last}, expected={expected_last}")

    if problems:
        errors.append(f"{taxon}: " + "; ".join(problems))

if errors:
    raise SystemExit(
        "HARP output trajectory metadata audit FAILED:\n" + "\n".join(errors)
    )

print(f"HARP output trajectory metadata audit: PASS ({len(taxa)} taxa)")
PY
}

build_panel_manifest() {
    log_line "Generating panel_manifest_26taxa.csv from fresh summaries..."

    python - \
        "${PANEL_MANIFEST}" \
        "${TAXON_RESULT_ROOT}" \
        "${TAXA_FILE}" <<'PY'


import csv
import os
import sys
from pathlib import Path

manifest = Path(sys.argv[1]).resolve()
result_root = Path(sys.argv[2]).resolve()
taxa_file = Path(sys.argv[3])



taxa = [
    line.strip()
    for line in taxa_file.read_text(encoding="utf-8").splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]

rows = []
for taxon in taxa:
    summary = result_root / taxon / "harp_v4_summary.json"
    if not summary.is_file() or summary.stat().st_size == 0:
        raise SystemExit(f"{taxon}: missing fresh summary {summary}")

    rel = Path(os.path.relpath(summary.resolve(), start=manifest.parent))

    rows.append((taxon, rel.as_posix()))

manifest.parent.mkdir(parents=True, exist_ok=True)
with manifest.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["taxon", "summary_json"])
    writer.writerows(rows)

print(f"Wrote {manifest}")
print(f"Manifest taxa: {len(rows)}")
PY

    local lines
    lines="$(wc -l < "${PANEL_MANIFEST}")"
    [[ "${lines}" -eq 27 ]] || die \
        "Generated manifest has ${lines} lines; expected 27"

    log_line "Panel manifest: PASS (${PANEL_MANIFEST})"
}

run_panel_test() {
    build_panel_manifest

    log_line "Running 26-taxon HARP panel..."
    local log="${PANEL_LOG_ROOT}/panel_26taxa.log"
        local panel_extra=()

    if [[ "${MODE}" == "all" ]]; then
        panel_extra+=(--direct-outdir --run-name "${RUN_NAME}")
    fi



    harp_cli panel \
        --manifest "${PANEL_MANIFEST}" \
        --outdir "${PANEL_RESULT_ROOT}" \
        --n-permutations "${N_PERMUTATIONS}" \
        --seed "${PANEL_SEED}" \
        "${panel_extra[@]}" \
        > "${log}" 2>&1 || die "Panel command failed; see ${log}"

    local newest

    if [[ "${MODE}" == "all" ]]; then
        newest="${PANEL_RESULT_ROOT}"
    else
        newest="$(
            find "${PANEL_RESULT_ROOT}" -mindepth 1 -maxdepth 1 -type d \
                -name 'HARP_v4.1_panel_n26_*' \
                -printf '%T@ %p\n' |
            sort -nr |
            head -n1 |
            cut -d' ' -f2-
        )"
    fi

    [[ -n "${newest}" && -d "${newest}" ]] || die \
        "Panel completed but its output directory was not found"

    local summary="${newest}/panel_summary.json"
    local money_shot="${newest}/panel_null.png"

    [[ -s "${summary}" ]] || die "Panel summary missing: ${summary}"
    [[ -s "${money_shot}" ]] || die "Money-shot figure missing: ${money_shot}"

    python - \
        "${summary}" \
        "${N_PERMUTATIONS}" \
        "${PANEL_SEED}" \
        "${FROZEN_PANEL_OBSERVED}" \
        "${FROZEN_PANEL_P_VALUE}" \
        "${FROZEN_PANEL_NULL_MEAN}" \
        "${FROZEN_PANEL_NULL_SD}" \
        "${FROZEN_PANEL_NULL_Q95}" \
        "${FROZEN_PANEL_N_PERMUTATIONS}" \
        "${FROZEN_PANEL_SEED}" \
        "${FROZEN_PANEL_ABS_TOL}" <<'PY'
import json
import math
import sys
from pathlib import Path

(
    summary_path,
    run_n_permutations,
    run_seed,
    expected_observed,
    expected_p_value,
    expected_null_mean,
    expected_null_sd,
    expected_null_q95,
    expected_n_permutations,
    expected_seed,
    abs_tol,
) = sys.argv[1:]

p = Path(summary_path)
obj = json.loads(p.read_text(encoding="utf-8"))

if int(obj.get("n_taxa", -1)) != 26:
    raise SystemExit(f"Expected n_taxa=26, found {obj.get('n_taxa')}")

run_n_permutations = int(run_n_permutations)
run_seed = int(run_seed)
expected_n_permutations = int(expected_n_permutations)
expected_seed = int(expected_seed)
abs_tol = float(abs_tol)

parameter_failures = []
if run_n_permutations != expected_n_permutations:
    parameter_failures.append(
        f"n_permutations={run_n_permutations}, expected={expected_n_permutations}"
    )
if run_seed != expected_seed:
    parameter_failures.append(f"seed={run_seed}, expected={expected_seed}")

if parameter_failures:
    raise SystemExit(
        "Frozen n=26 regression gate requires canonical panel parameters:\n  "
        + "\n  ".join(parameter_failures)
    )

result = obj["result"]

expected = {
    "observed": float(expected_observed),
    "p_value": float(expected_p_value),
    "null_mean": float(expected_null_mean),
    "null_sd": float(expected_null_sd),
    "null_q95": float(expected_null_q95),
}
actual = {
    "observed": float(result["observed"]),
    "p_value": float(result["p_value"]),
    "null_mean": float(result["null_mean"]),
    "null_sd": float(result["null_sd"]),
    "null_q95": float(result["null_q95"]),
}

failures = []
for key, expected_value in expected.items():
    actual_value = actual[key]
    if not math.isfinite(actual_value):
        failures.append(f"{key} is non-finite: {actual_value!r}")
    elif not math.isclose(
        actual_value,
        expected_value,
        rel_tol=0.0,
        abs_tol=abs_tol,
    ):
        failures.append(
            f"{key}={actual_value:.15g}, "
            f"expected={expected_value:.15g}, "
            f"|delta|={abs(actual_value - expected_value):.3g}"
        )

actual_n_permutations = int(result["n_permutations"])
if actual_n_permutations != expected_n_permutations:
    failures.append(
        f"result.n_permutations={actual_n_permutations}, "
        f"expected={expected_n_permutations}"
    )

print("26-taxon panel result")
print(f"  observed:      {actual['observed']:.15f}")
print(f"  p_value:       {actual['p_value']:.8g}")
print(f"  null_mean:     {actual['null_mean']:.15f}")
print(f"  null_sd:       {actual['null_sd']:.15f}")
print(f"  null_q95:      {actual['null_q95']:.15f}")
print(f"  permutations:  {actual_n_permutations}")
print(f"  seed:          {run_seed}")
print(f"  run_name:      {obj.get('run_name')}")

if failures:
    raise SystemExit(
        "Frozen HARP v4.1 n=26 regression gate FAILED:\n  "
        + "\n  ".join(failures)
    )

print(
    "Frozen HARP v4.1 n=26 regression gate: PASS "
    f"(absolute tolerance {abs_tol:g})"
)
PY

    log_line "Panel regression gate: PASS"
    log_line "Panel PASS: ${newest}"
    log_line "Money shot: ${money_shot}"
}

run_tests_if_present() {
    if [[ ! -d "${ROOT}/tests" ]]; then
        log_line "NOTE: no tests/ directory found; skipping pytest."
        return 0
    fi

    log_line "Running packaged pytest suite..."

    set +e
    PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
        python -m pytest -q
    local rc=$?
    set -e

    case "${rc}" in
        0)
            log_line "Packaged pytest suite: PASS"
            ;;
        5)
            log_line "NOTE: pytest collected no tests; skipping packaged-test gate."
            ;;
        *)
            die "Packaged pytest suite failed with exit code ${rc}"
            ;;
    esac
}

preflight() {
    ensure_dirs
    load_taxa
    check_environment
    check_corpus_structure
    verify_inventory
    audit_reference_trajectories
    check_source_configs
    run_tests_if_present
    log_line "26-taxon reference preflight: PASS"
}

master_log="${LOG_ROOT}/HARP_v4.1_reference26_$(date '+%Y%m%dT%H%M%S%z').log"
ensure_dirs
exec > >(tee -a "${master_log}") 2>&1

log_line "============================================================"
log_line "HARP v4.1 — 26-taxon GitHub reference orchestrator"
log_line "Command: ${MODE}"
log_line "Release root: ${ROOT}"
log_line "Corpus:       ${CORPUS_ROOT}"
if [[ "${MODE}" == "all" ]]; then
    log_line "Run root:     ${RUN_ROOT}"
fi
log_line "Master log:   ${master_log}"
log_line "============================================================"

case "${MODE}" in
    preflight)
        preflight
        ;;
    configs)
        preflight
        create_run_configs
        verify_inventory
        ;;
    validate)
        preflight
        create_run_configs
        validate_all
        verify_inventory
        ;;
    analyze)
        preflight
        create_run_configs
        validate_all
        analyze_all
        audit_analysis_trajectory_metadata
        verify_inventory
        ;;
    panel)
        ensure_dirs
        load_taxa
        check_environment
        check_corpus_structure
        verify_inventory
        audit_reference_trajectories
        run_panel_test
        verify_inventory
        ;;
    all)
        preflight
        create_run_configs
        validate_all
        analyze_all
        audit_analysis_trajectory_metadata
        verify_inventory
        run_panel_test
        verify_inventory
        ;;
esac

log_line "============================================================"
log_line "HARP v4.1 reference26 ${MODE}: PASS"
log_line "Frozen scientific inputs remained SHA-256 identical."
log_line "============================================================"
