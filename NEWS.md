# Changelog

All notable changes to `because.phybase` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/2.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

# because.phybase 1.1.0 (2026-11-07)

* **Deterministic Phylogenetic Ensemble (`tree_method = "ensemble"`, new default for `multiPhylo`)**: Passing a `multiPhylo` object to `because(..., structure = trees)` now defaults to `tree_method = "ensemble"` with `chains_per_tree = 2L` (Bayesian Model Averaging with guaranteed equal representation across all $M$ trees, following [`brms::brm_multiple`](https://paul-buerkner.github.io/brms/reference/brm_multiple.html) and [*mulTree*](https://doi.org/10.5281/zenodo.12902); Guillerme & Healy 2014). In `engine = "numpyro"`, the pure `NUTS` kernel is compiled only once and executed across all $M \times C$ tree chains in multi-core parallel waves (`n.cores`).
* **Split MCMC vs. Phylogenetic Sensitivity Diagnostics (`check_multiphylo()`)**: Added exported `check_multiphylo()` diagnostic function and updated `summary()` to report both within-tree MCMC convergence (`Rhat` / `Rhat_within`) and between-tree phylogenetic sensitivity (`Rhat_tree`).
* **Spectral Gaussian Marginalization & Eigenvector Sign Alignment for `tree_method = "dcat"`**: When `tree_method = "dcat"` is requested in NumPyro, observed Gaussian equations analytically marginalize out the species random effect in the tree's eigenbasis ($\mathbf{Q}_K^\top (\mathbf{y} - \boldsymbol{\mu})$), and eigenvector signs are aligned across trees relative to Tree 1.

# because.phybase 1.0.2 (2026-10-28)

* **Fix MultiPhylo JAGS Model Definition**: Fixed an uninitialized variable error (`object 'sigma_phylo_param' not found`) in `jags_structure_definition.multiPhylo()` when `use_partitioning` is disabled (such as additive error models or non-Gaussian paths).
* **Preserve Observation Row Order in Structure Preparation**: `prepare_structure_data.phylo()` and `prepare_structure_data.multiPhylo()` now align phylogeny tip labels to data rows using the exact original observation row order (`as.character(unique(args$row_ids))`) rather than forcing alphabetical sorting (`levels(as.factor())`). This prevents loss or distortion of phylogenetic signal when data rows are not pre-sorted alphabetically.
* **MultiPhylo Data Compatibility**: Added support for nested data frames when matching tip labels across multiPhylo tree lists.

# because.phybase 1.0.1 (2025-12-22)

* Bump version from 1.0.0 to 1.0.1.
* Add K_tree selection and priors for multiPhylo.
* Print current package version on package load.
