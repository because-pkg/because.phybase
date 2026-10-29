# because.phybase 1.0.2 (2026-10-28)

* **Preserve Observation Row Order in Structure Preparation**: `prepare_structure_data.phylo()` and `prepare_structure_data.multiPhylo()` now align phylogeny tip labels to data rows using the exact original observation row order (`as.character(unique(args$row_ids))`) rather than forcing alphabetical sorting (`levels(as.factor())`). This prevents loss or distortion of phylogenetic signal when data rows are not pre-sorted alphabetically.
* **MultiPhylo Data Compatibility**: Added support for nested data frames when matching tip labels across multiPhylo tree lists.

# because.phybase 1.0.1 (2025-12-22)

* Bump version from 1.0.0 to 1.0.1.
* Add K_tree selection and priors for multiPhylo.
* Print current package version on package load.
