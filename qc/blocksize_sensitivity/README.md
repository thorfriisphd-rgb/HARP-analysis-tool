# HARP v4.1 — Block-size Sensitivity QC

This folder contains a **standalone release-QC harness**, not a new HARP feature.
It tests whether the taxon-level `period7_max_mode` inference depends materially
on the production contiguous-block shuffle size of 4 residues.

## Location

Place the folder at:

```text
<release-root>/qc/blocksize_sensitivity/
```

so that the script is two directory levels below the release root containing
`run_reference_26taxa.sh`, `src/`, and `reference/`.

## Run

Activate the normal HARP v4.1 environment, then:

```bash
cd <release-root>/qc/blocksize_sensitivity
chmod +x qc_blocksize_sensitivity.sh
./qc_blocksize_sensitivity.sh all
```

Defaults:

```text
block sizes       1 2 3 4 5 6 7
permutations      9999 per taxon statistic
analysis workers  1
reporting alpha   0.05
reference block   4
```

To use two analysis workers:

```bash
HARP_QC_ANALYSIS_JOBS=2 ./qc_blocksize_sensitivity.sh all
```

Increase this cautiously because concurrent MDAnalysis jobs can be memory-heavy.

## What it does

1. Calls the frozen 26-taxon reference runner's `configs` gate, which performs
   the normal preflight and SHA-256 checks and generates the current portable
   reference configs.
2. Creates QC-local copies of those configs for each tested block size.
3. Changes only `statistics.block_size`, `statistics.n_permutations`, and the
   QC-local `output_dir`.
4. Validates and analyses all 26 taxa de novo for every block size.
5. Harvests all three existing HARP taxon statistics, with
   `period7_max_mode` treated as the primary sensitivity target.
6. Verifies that the observed `period7_max_mode` statistic itself is invariant
   to block size; only its null distribution should change.
7. Produces a long-format TSV, compact primary-statistic summary, Markdown
   sensitivity report, heatmap, source hashes, metadata, and full logs.

## Output

Each run is self-contained:

```text
runs/
└── HARP_v4.1_blocksize_sensitivity_n26_<timestamp>/
    ├── block_01/
    ├── block_02/
    ├── ...
    ├── block_07/
    ├── sensitivity_results.tsv
    ├── period7_max_mode_summary.tsv
    ├── sensitivity_report.md
    ├── period7_max_mode_sensitivity.png
    ├── source_hashes.tsv
    ├── run_metadata.txt
    └── qc.log
```

## QC verdict

For reporting only, a taxon is `STABLE` when all tested block sizes lie on the
same side of `p = 0.05`; otherwise it is flagged `CROSSES_ALPHA` for review.
The overall report says `ROBUST` only when no taxon crosses that threshold.
This is **not a new HARP significance rule**—it is a transparent QC summary.

The report also gives the full p-value range, so marginal changes are visible
rather than hidden behind the binary flag.

## Why the panel is not rerun seven times

`block_size` is consumed by the taxon-level contiguous-block permutation test.
The taxon `phase_signature` is constructed from the observed phase enrichment,
not from the taxon null distribution, and the panel uses those signatures in
its independent circular-rotation null. Therefore varying the taxon block size
does not alter the panel input or panel statistic. Re-running the panel seven
times would only duplicate the same panel analysis.
