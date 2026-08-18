# HARP v4.1 26-taxon reference corpus

The current reference panel contains the 26 taxa listed in `taxa.txt`.
Large MD inputs (`.tpr`, `.xtc`) are intentionally not embedded in this
streamlined source tree; copy the frozen corpus into a local reference-data
location before reproducing the complete panel.

Naegleria is not treated as a special statistical class. An earlier
*Naegleria gruberi* MyhT candidate was rejected by HARP because its
DeepCoil2-predicted coiled-coil register did not satisfy the canonical
heptad structure required by the analysis. The accepted 26-taxon reference
corpus is analysed under the same HARP pipeline as the other taxa.

layout. Freeze the 26-taxon reference score and null statistics only after a
clean release-candidate rerun under this packaged version.
