# HARP — Heptad Assignment Register Probe

HARP is a reproducible computational pipeline for detecting phase-aligned contact enrichment in coiled-coil interactions.  
It is developed and applied here to test the hypothesis that **C12/IBAM is a conserved contractile-system-associated protein**, rather than an RNA ligase.

---

##  Overview

This repository provides:

- A complete HARP analysis pipeline
    
- Precomputed PRCO interface datasets across 26 taxa
    
- Fully reproducible generation of all manuscript figures
    

Running:

```bash
python scripts/run_batch_harp_panel.py
```

will:

- Compute HARP scores across taxa
    
- Perform empirical null-model testing
    
- Generate all manuscript figures (SVG, PNG, PDF)
    

---

##  Repository Structure

```
.├── data/  
│ ├── samples.csv  
│ └── prco/  
│ ├── *_prco.csv  
│  
├── scripts/  
│ ├── run_batch_harp_panel.py  
│ ├── generate_harp_plots.py  
│ └── harp_run_scaffold.py  
│  
├── requirements.txt  
└── README.md
```

---

##  Installation

Tested with Python 3.10+

Install dependencies:

```bash
pip install -r requirements.txt
```

---

##  Running the Pipeline

From the repository root:

```bash
python scripts/run_batch_harp_panel.py
```

Outputs are written to:

```
HARP-output/<timestamp>/
```

Each run is self-contained and does not overwrite previous results.

---

##  Outputs

#### Outputs are generated on run and written to
```
HARP-output/<timestamp>/
├── harp_panel_summary.tsv
├── harp_panel_summary_<timestamp>.tsv
├── harp_rankings.tsv
└── figures/
    ├── 01_observed_vs_null_p95.*
    ├── 02_score_vs_margin.*
    └── 03_phase_score_heatmap.*

```

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
    

---

##  Interpretation

HARP quantifies:

- **Best score** → enrichment of a/d heptad contacts
    
- **Margin** → decisiveness of phase assignment
    
- **Null model** → statistical significance via permutation
    

Together, these distinguish:

- Strong, phase-specific binding
    
- Multi-register compatibility
    
- Weak or flat interaction profiles
    

---

##  Scientific Context

This pipeline supports the interpretation that:

> IBAM (C12orf29) encodes a conserved structural interface that recognises the physicochemical grammar of a degenerate coiled-coil substrate, rather than acting as a sequence-specific RNA ligase.

The evidence integrates:

- Structural prediction (AlphaFold3)
    
- Molecular dynamics simulations
    
- PRCO interface decoding
    
- HARP register analysis
    

---

##  Reproducibility

All results and figures can be regenerated directly:

```bash
python scripts/run_batch_harp_panel.py
```

No external datasets are required.

---

##  License

MIT License

Copyright (c)

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

## 🤝 Contact

Author:  Thor Einar Friis, PhD
ORCID:(https://orcid.org/0000-0002-4132-4912)
