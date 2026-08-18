# HARP v4.1 Release Notes

HARP v4.1 is the frozen GitHub release of the **Heptad Assignment Register Probe**.

This release consolidates the reviewed v4.1 code into an installable and reproducible package. The work since the scientific freeze has been confined to packaging, validation, output organization, documentation, QC, and reproducibility safeguards.

## What changed

* Converted the core modules into an installable `harp` package with the `harp` CLI.
* Separated generic HARP analysis from reference-corpus-specific orchestration.
* Added packaged unit and smoke tests.
* Added SHA-256 input verification and reference trajectory auditing.
* Added timestamped, taxon-labelled reference output organization.
* Retained the publication-facing panel-null renderer with clarified axis labels.
* Established the 26-taxon panel as the canonical HARP v4.1 reference benchmark.
* Added a frozen n=26 numerical regression gate to `run_reference_26taxa.sh`.
* Separated the large 26-taxon scientific payload from the GitHub repository while retaining its canonical roster and SHA-256 inventory in the release.
* Retained the preceding 25-taxon benchmark as a legacy regression reference.

## Frozen reference result

The final clean 26-taxon run reproduced:

```text
observed      0.645128594598065
p-value       0.0001
null mean    -0.000521984561169
null SD       0.151023083942446
null q95      0.225049488444680
permutations  9999
seed          20260801
```

The untouched reference corpus passed structure, SHA-256, MDAnalysis trajectory, authoritative configuration, and packaged-test checks before the complete analysis was run.

## Scientific status

No intended scientific or statistical logic change was introduced during release consolidation.

The clean-corpus reproduction of the frozen 26-taxon benchmark is the final integration check for the HARP v4.1 release.

For installation, reference-data setup, commands, outputs, and reproduction instructions, see `README.md`.

