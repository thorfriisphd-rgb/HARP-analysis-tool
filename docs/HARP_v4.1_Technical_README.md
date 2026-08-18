# HARP v4.1 — Heptad Assignment Register Probe

HARP tests whether trajectory-derived IBAM–MyhT contact organization contains non-random structure relative to an explicit DeepCoil2-anchored heptad register and defined null models.

## What HARP v4.1 does

For each taxon, HARP:

1. parses DeepCoil2 output and establishes a canonical a–g register;
2. validates sequence/register/topology agreement without repairing inputs;
3. derives an MG↔MyhT residue-contact trajectory from MD coordinates;
4. summarizes contact occupancy by heptad phase;
5. evaluates per-taxon statistics under block-shuffle null models; and
6. writes machine-readable results, diagnostic figures, and provenance.

Across completed taxa, `harp panel` tests shared phase alignment. Each taxon's centered seven-phase signature is compared with the consensus of the remaining taxa. The panel null independently circularly rotates every taxon's signature, preserving its shape and strength while breaking common phase registration.

The publication-facing panel-null figure reports the observed alignment, null 95th percentile, permutation p-value, null mean/SD, taxon count, and permutation count. Its x-axis is labelled **Mean cross-taxon phase similarity**. All inferential quantities are calculated upstream; the renderer is presentation-only.

## Installation

HARP v4.1 requires Python 3.11 or later. Conda is recommended for a reproducible scientific environment.

From the repository root:

```bash
conda env create -f environment.yml
conda activate harp41
```

Alternatively, install HARP into an existing Python 3.11+ environment:

```bash
python -m pip install -e .
```

Verify the installation:

```bash
harp --help
```

The full 26-taxon reference analysis additionally requires the separately distributed scientific input corpus described below.

## Core commands

```bash
harp register --deepcoil myht.out --out deepcoil_register.tsv
harp validate --config config.yaml --report validation_report.json
harp analyze --config config.yaml
harp panel --manifest panel_manifest.csv --outdir results \
  --n-permutations 9999 --seed 20260801
```

A complete configuration template is provided in `examples/config.yaml` and a minimal panel manifest in `examples/panel_manifest.csv`.

## Validation philosophy

HARP v4.1 validation is deliberately conservative and diagnostic-only. It does not reverse sequences, trim them, relax thresholds, repair registers, or rewrite scientific selections in order to make a dataset pass.

A canonical DeepCoil2 heptad register requires dominant d anchors at a+3 relative to the dominant a frame.

Reference-corpus-specific constraints, such as the frozen 26-taxon trajectory length and sampling interval, are enforced by the reference runner rather than by generic HARP validation. This allows HARP itself to analyse trajectories of different durations and sampling schemes.

## Outputs

A per-taxon analysis writes, among other files:

* `harp_v4_summary.json`
* `validation_report.json`
* `provenance.json`
* `per_residue_contacts.tsv`
* `phase_summary.tsv`
* `trajectory_contacts.npz`
* per-taxon diagnostic PNGs

A complete 26-taxon reference `all` run creates a timestamped results directory such as:

```text
results/HARP_v4.1_panel_n26_20260814T174346+0200/
```

containing:

```text
taxa/
panel/
```

The `taxa/` directory contains the individual taxon outputs. The `panel/` directory contains the canonical panel summary, manifest used, taxon scores, parameter/register audits, null distribution, provenance, and the panel-null figure.

## 26-taxon reference corpus

The large molecular-dynamics reference corpus is distributed separately from the GitHub source repository.

The GitHub repository retains the small canonical metadata required to identify and verify that corpus:

```text
reference/26taxa/
├── taxa.txt
└── HARP_v4.1_26taxa_input_sha256.tsv
```

The scientific payload consists of 26 taxon directories, each containing:

```text
md.tpr
md.xtc
myht.fa
myht.out
```

The complete corpus therefore contains 104 scientific input files.

After downloading the archived reference corpus from Zenodo, extract it so that the external data directory has the form:

```text
/path/to/reference_data/26taxa/
├── <Taxon_1>/
│   ├── md.tpr
│   ├── md.xtc
│   ├── myht.fa
│   └── myht.out
├── <Taxon_2>/
│   └── ...
└── ...
```

Set `HARP_REFERENCE_DATA_ROOT` to the extracted **`26taxa/` directory itself**, not to its parent:

```bash
export HARP_REFERENCE_DATA_ROOT=/path/to/reference_data/26taxa
```

The reference runner then combines the external scientific payload with the canonical taxon roster, SHA-256 inventory, authoritative configurations, code, and tests retained in the GitHub repository.

## Reproducing the frozen 26-taxon reference analysis

From the HARP repository root, with the `harp41` environment active and `HARP_REFERENCE_DATA_ROOT` set:

```bash
./run_reference_26taxa.sh preflight
```

The preflight checks:

* the 26-taxon corpus structure;
* all 104 scientific files against the canonical SHA-256 inventory;
* all trajectories with MDAnalysis;
* the authoritative per-taxon reference configurations; and
* the packaged pytest suite.

For the frozen reference corpus, each trajectory must contain 1001 frames spanning 0–10,000 ps at 10 ps intervals.

To perform the complete reference analysis:

```bash
./run_reference_26taxa.sh all
```

The runner validates and analyses all 26 taxa, constructs the panel manifest from the newly generated summaries, performs the panel permutation analysis, verifies the scientific inputs again, and applies the frozen numerical regression gate.

## Frozen HARP v4.1 benchmark

The clean 26-taxon reference analysis for HARP v4.1 is frozen at:

```text
taxa:          26
observed:      0.645128594598065
p-value:       0.0001
null mean:    -0.000521984561169
null SD:       0.151023083942446
null q95:      0.225049488444680
permutations:  9999
seed:          20260801
```

`run_reference_26taxa.sh` treats these values as a release regression gate. A canonical reference run must reproduce the frozen benchmark within the numerical tolerance defined by the runner.

The benchmark values are verification metadata only. They do not enter the HARP analysis or alter its scientific calculations.

The runner also verifies the SHA-256 inventory after analysis to confirm that the frozen scientific inputs remained unchanged.

## Tests

Run the packaged test suite with:

```bash
python -m pytest
```

The bundled tests cover deterministic panel permutation behaviour, DeepCoil register inference, provenance naming, and generation of the panel-null figure.

The full 26-taxon trajectory corpus is intentionally not bundled with the source repository because of its size. `run_reference_26taxa.sh` therefore provides the separate integration and scientific reproducibility gate against the Zenodo corpus.

## Legacy reference

`reference/legacy_25taxa/` retains the preceding 25-taxon benchmark:

```text
observed = 0.632109037896326
permutations = 9999
seed = 20260801
```

It is retained as a historical regression reference and is not the HARP v4.1 headline corpus.
