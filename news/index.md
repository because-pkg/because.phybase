# Changelog

## because.phybase 1.1.0 (2026-11-07)

> All notable changes to `because.phybase` are documented in this file.
> The format is based on [Keep a
> Changelog](https://keepachangelog.com/en/2.0.0/), and this project
> adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

- **Deterministic Phylogenetic Ensemble (`tree_method = "ensemble"`, new
  default for `multiPhylo`)**: Passing a `multiPhylo` object to
  `because(..., structure = trees)` now defaults to
  `tree_method = "ensemble"` with `chains_per_tree = 2L` (Bayesian Model
  Averaging with guaranteed equal representation across all M trees,
  following the multi-tree Bayesian integration framework of
  [*mulTree*](https://doi.org/10.5281/zenodo.12902) \[Guillerme & Healy
  2014\] and [Nakagawa & de Villemereuil
  2019](https://doi.org/10.1093/sysbio/syy089)). In
  `engine = "numpyro"`, the pure `NUTS` kernel is compiled only once and
  executed across all M \times C tree chains in multi-core parallel
  waves (`n.cores`).
- **Split MCMC vs. Phylogenetic Sensitivity Diagnostics
  ([`check_multiphylo()`](https://because-pkg.github.io/because.phybase/reference/check_multiphylo.md))**:
  Added exported
  [`check_multiphylo()`](https://because-pkg.github.io/because.phybase/reference/check_multiphylo.md)
  diagnostic function (providing a domain-specific phylogenetic alias
  for the general covariance engine
  [`because::check_multicov()`](https://rdrr.io/pkg/because/man/check_multicov.html))
  and updated [`summary()`](https://rdrr.io/r/base/summary.html) to
  report both within-tree MCMC convergence (`Rhat` / `Rhat_within`) and
  between-tree phylogenetic sensitivity (`Rhat_tree`).
- **Spectral Gaussian Marginalization & Eigenvector Sign Alignment for
  `tree_method = "dcat"`**: When `tree_method = "dcat"` is requested in
  NumPyro (sampling candidate trees dynamically via a discrete
  categorical prior; [de Villemereuil et
  al. 2012](https://doi.org/10.1186/1471-2148-12-102)), observed
  Gaussian equations analytically marginalize out the species random
  effect in the tree’s eigenbasis (\mathbf{Q}\_K^\top (\mathbf{y} -
  \boldsymbol{\mu})), and eigenvector signs are aligned across trees
  relative to Tree 1.

## because.phybase 1.0.2 (2026-10-28)

- **Fix MultiPhylo JAGS Model Definition**: Fixed an uninitialized
  variable error (`object 'sigma_phylo_param' not found`) in
  [`jags_structure_definition.multiPhylo()`](https://because-pkg.github.io/because.phybase/reference/jags_structure_definition.multiPhylo.md)
  when `use_partitioning` is disabled (such as additive error models or
  non-Gaussian paths).
- **Preserve Observation Row Order in Structure Preparation**:
  [`prepare_structure_data.phylo()`](https://because-pkg.github.io/because.phybase/reference/prepare_structure_data.phylo.md)
  and
  [`prepare_structure_data.multiPhylo()`](https://because-pkg.github.io/because.phybase/reference/prepare_structure_data.multiPhylo.md)
  now align phylogeny tip labels to data rows using the exact original
  observation row order (`as.character(unique(args$row_ids))`) rather
  than forcing alphabetical sorting (`levels(as.factor())`). This
  prevents loss or distortion of phylogenetic signal when data rows are
  not pre-sorted alphabetically.
- **MultiPhylo Data Compatibility**: Added support for nested data
  frames when matching tip labels across multiPhylo tree lists.

## because.phybase 1.0.1 (2025-12-22)

- Bump version from 1.0.0 to 1.0.1.
- Add K_tree selection and priors for multiPhylo.
- Print current package version on package load.
