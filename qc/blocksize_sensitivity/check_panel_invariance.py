#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from harp.stats import panel_leave_one_out_test


BLOCKS = range(1, 8)


def load_signatures(run_root: Path):
    block_data = {}

    for b in BLOCKS:
        block_dir = run_root / f"block_{b:02d}" / "taxa"
        if not block_dir.is_dir():
            raise SystemExit(f"Missing block directory: {block_dir}")

        taxa = {}
        for taxon_dir in sorted(p for p in block_dir.iterdir() if p.is_dir()):
            summary = taxon_dir / "harp_v4_summary.json"
            if not summary.is_file():
                raise SystemExit(f"Missing summary: {summary}")

            obj = json.loads(summary.read_text(encoding="utf-8"))
            sig = np.asarray(
                obj["descriptive"]["phase_signature"],
                dtype=float,
            )

            if sig.shape != (7,) or not np.all(np.isfinite(sig)):
                raise SystemExit(
                    f"{taxon_dir.name}, block {b}: invalid phase_signature"
                )

            taxa[taxon_dir.name] = sig

        block_data[b] = taxa

    return block_data


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Verify that taxon phase signatures and the observed HARP panel "
            "statistic are invariant across block-size sensitivity runs."
        )
    )
    parser.add_argument("run_root", type=Path)
    parser.add_argument(
        "--atol",
        type=float,
        default=1e-12,
        help="Absolute tolerance for floating-point comparisons.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260801,
        help="Seed passed to the existing panel routine.",
    )
    args = parser.parse_args()

    run_root = args.run_root.resolve()
    data = load_signatures(run_root)

    reference_taxa = list(data[1].keys())
    reference_set = set(reference_taxa)

    rows = []
    failures = []

    for b in BLOCKS:
        if set(data[b]) != reference_set:
            failures.append(
                f"block {b}: taxon set differs from block 1"
            )

    if failures:
        for msg in failures:
            print(f"FAIL\t{msg}", file=sys.stderr)
        raise SystemExit(2)

    # Taxon-level phase-signature invariance
    ref = data[1]

    for taxon in reference_taxa:
        for b in BLOCKS:
            diff = float(
                np.max(np.abs(data[b][taxon] - ref[taxon]))
            )

            status = "PASS" if diff <= args.atol else "FAIL"

            rows.append(
                {
                    "taxon": taxon,
                    "block_size": b,
                    "max_abs_signature_difference": diff,
                    "status": status,
                }
            )

            if status == "FAIL":
                failures.append(
                    f"{taxon}, block {b}: "
                    f"phase_signature differs by {diff:.3e}"
                )

    # Observed panel-statistic invariance
    panel_rows = []

    for b in BLOCKS:
        signatures = np.asarray(
            [data[b][taxon] for taxon in reference_taxa],
            dtype=float,
        )

        # We only need the observed statistic. One permutation is sufficient;
        # the randomized null itself is not being tested here.
        result, _, _ = panel_leave_one_out_test(
            signatures,
            n_permutations=2,
            seed=args.seed,
        )

        panel_rows.append(
            {
                "block_size": b,
                "observed_panel_statistic": float(result.observed),
            }
        )

    reference_panel = panel_rows[0]["observed_panel_statistic"]

    for row in panel_rows:
        diff = abs(
            row["observed_panel_statistic"] - reference_panel
        )
        row["abs_difference_from_block1"] = diff
        row["status"] = "PASS" if diff <= args.atol else "FAIL"

        if row["status"] == "FAIL":
            failures.append(
                f"block {row['block_size']}: observed panel statistic "
                f"differs by {diff:.3e}"
            )

    # Write audit products
    import csv

    sig_out = run_root / "panel_invariance_signatures.tsv"
    with sig_out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "taxon",
                "block_size",
                "max_abs_signature_difference",
                "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    panel_out = run_root / "panel_invariance.tsv"
    with panel_out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "block_size",
                "observed_panel_statistic",
                "abs_difference_from_block1",
                "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(panel_rows)

    report = run_root / "panel_invariance_report.md"

    if failures:
        verdict = "FAIL"
    else:
        verdict = "PASS"

    report.write_text(
        "\n".join(
            [
                "# HARP v4.1 panel invariance QC",
                "",
                f"**Panel invariance verdict: {verdict}**",
                "",
                "Block size affects the taxon-level block-shuffle permutation "
                "null but should not affect the observed taxon phase signatures "
                "or the observed panel leave-one-out similarity statistic.",
                "",
                f"- Taxa checked: **{len(reference_taxa)}**",
                "- Block sizes checked: **1, 2, 3, 4, 5, 6, 7**",
                f"- Comparison tolerance: **{args.atol:g}**",
                "",
                "## Interpretation",
                "",
                (
                    "PASS: all taxon phase signatures and the observed panel "
                    "statistic are invariant across block-size runs."
                    if verdict == "PASS"
                    else
                    "FAIL: at least one quantity expected to be invariant "
                    "changed across block-size runs. This is a QC integrity "
                    "failure and should be investigated before interpreting "
                    "the sensitivity results."
                ),
                "",
                "Detailed comparisons are recorded in "
                "`panel_invariance_signatures.tsv` and `panel_invariance.tsv`.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Panel invariance: {verdict}")
    print(f"Wrote: {sig_out}")
    print(f"Wrote: {panel_out}")
    print(f"Wrote: {report}")

    if failures:
        print("", file=sys.stderr)
        for msg in failures:
            print(f"FAIL: {msg}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
