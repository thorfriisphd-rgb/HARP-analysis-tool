## Scientific Context

HARP is part of the broader CCMHCG/IBAM computational ecosystem, alongside PRCO and SWING.

In this framework:

- **PRCO** decodes persistent residue-level contacts from molecular dynamics trajectories.
- **HARP** tests whether those contacts are enriched in particular coiled-coil heptad phases.
- **SWING** tests whether projected interface residues show cross-taxon conservation and biochemical-class convergence.

Together, these methods support the working hypothesis that IBAM/C12orf29 recognises a degenerate coiled-coil interaction grammar associated with contractile systems.

HARP should not be viewed only as an endpoint validation assay. In the May 2026 RefSeq-guided revision, HARP/SWING behaviour helped expose and rationalise improved Myh-tail window selection, including the Ovis Myh7T rescue. This supports the use of HARP as an exploratory interaction-grammar discovery framework.

---

## May 2026 RefSeq-Rescued Dataset

The current canonical dataset incorporates the May 2026 input corrections used across the IBAM computational workflow.

Key updates include:

- correction of a major **Magallana angulata** C12 truncation
    - previous: C12 1–231, MyhT 232–304
    - current: C12 1–345, MyhT 346–413
- extension of the **Ovis aries** Myh7T window
    - previous: C12 1–325, MyhT 326–382
    - current: C12 1–325, MyhT 326–400
- regenerated PRCO tables for updated or harmonised taxa

The canonical May 2026 run is provided under:

```
HARP-output/2026-05-10_08-53-04/
```

---

## Repository Structure

```text
├── data/
│   ├── samples.csv
│   └── prco/
│       └── *_prco.csv
│
├── scripts/
│   ├── run_batch_harp_panel.py
│   ├── generate_harp_plots.py
│   └── harp_run_scaffold.py
│
├── HARP-output/
│   └── 2026-05-10_08-53-04/
│       ├── harp_panel_summary.tsv
│       ├── harp_panel_summary_20260510_085653.tsv
│       ├── harp_rankings.tsv
│       ├── harp_run_config.json
│       ├── harp_run.log
│       └── figures/
│
├── requirements.txt
└── README.md
```
---

## Installation

Tested with Python 3.10+.

Install dependencies:

```
pip install -r requirements.txt
```

---

## Running the Pipeline

From the repository root:

```
python scripts/run_batch_harp_panel.py
```

Outputs are written to:

```
HARP-output/<timestamp>/
```

Each run is self-contained and does not overwrite previous results.

---

## Outputs

A standard HARP run writes:

```text
HARP-output/<timestamp>/
├── harp_panel_summary.tsv
├── harp_panel_summary_<timestamp>.tsv
├── harp_rankings.tsv
├── harp_run_config.json
├── harp_run.log
└── figures/
    ├── 01_observed_vs_null_p95.svg
    ├── 01_observed_vs_null_p95.png
    ├── 01_observed_vs_null_p95.pdf
    ├── 02_score_vs_margin.svg
    ├── 02_score_vs_margin.png
    ├── 02_score_vs_margin.pdf
    ├── 03_phase_score_heatmap.svg
    ├── 03_phase_score_heatmap.png
    ├── 03_phase_score_heatmap.pdf
    └── harp_plots_manifest.txt
```

---
### Summary Tables

- **harp_panel_summary.tsv**
    Canonical output used for downstream plotting

- **harp_panel_summary_.tsv**
    Archived run-specific output

- **harp_rankings.tsv**
    Ranked view for rapid inspection

---
### Figures (auto-generated)

1. **Observed vs null (95th percentile)**
    Tests whether observed phase enrichment exceeds null expectation

2. **Best score vs phase margin**
    Separates strong/specific vs multi-register binding regimes

3. **Phase-score landscape heatmap**
    Full heptad phase distribution across taxa

Figures are generated as SVG, PNG, and PDF.

---
## Interpretation

HARP quantifies whether persistent IBAM–MyhT contacts preferentially occupy particular heptad phases.

Core metrics include:

- **best score** — strongest observed phase-aligned contact enrichment
- **margin** — separation between best and second-best phase scores
- **best-minus-flat** — deviation from a flat/non-enriched profile
- **empirical null p-values** — comparison against shuffled occupied-span nulls
- **phase stability** — consistency of phase calls across cutoffs

HARP distinguishes:

- stable phase-biased interfaces
- competing or multi-register phase profiles
- weak or flat interaction profiles

For the IBAM project, competing phases are not automatically interpreted as failure. A degenerate coiled-coil substrate may retain biologically meaningful multi-register compatibility.

---

## Null Model

The current canonical run uses:

```text
shuffle_occupied_span
```

with:

```text
10,000 permutations
seed = 123
```

This tests whether the observed contact distribution exceeds expectations from shuffled occupied positions while preserving the occupied span structure of the interface.

---

## Relationship to SWING

HARP and SWING test related but distinct features of the IBAM–MyhT interaction grammar.

- HARP asks whether MD-derived contact residues are enriched in coiled-coil heptad phases.
- SWING asks whether projected interface positions show evolutionary conservation and biochemical-class convergence across taxa.

The May 2026 RefSeq-guided rescue, especially the Ovis Myh7T window extension and Magallana truncation correction, strengthened both workflows and supports their use as exploratory tools for interaction-window discovery.

---

## Reproducibility

All included outputs and figures can be regenerated directly:

```bash
python scripts/run_batch_harp_panel.py
```

No external datasets are required for the included 26-taxon panel.

The canonical May 2026 run was generated with:

```text
samples: data/samples.csv
prco_dir: data/prco
null_model: shuffle_occupied_span
shuffle_iter: 10000
seed: 123
```

---
## License

MIT License

Copyright (c) Thor Einar Friis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## Citation

If you use this pipeline, please cite:

>Friis TE. C12orf29 encodes IBAM (In Between Actin and Myosin), a conserved actomyosin-associated protein exhibiting deeply conserved interaction grammar across eukaryotic evolution. Manuscript in preparation.

---


## Author

Thor Einar Friis

[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--4132--4912-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0000-0002-4132-4912)

Independent researcher, Bodø, Norway.
PhD in Molecular Biology, Queensland University of Technology (QUT).
