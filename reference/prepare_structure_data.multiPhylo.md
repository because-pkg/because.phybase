# Prepare Structure Data for Multiple Trees

Implements the `prepare_structure_data` S3 method for `multiPhylo`
objects.

## Usage

``` r
# S3 method for class 'multiPhylo'
prepare_structure_data(structure, data, optimize = TRUE, quiet = FALSE, ...)
```

## Arguments

- structure:

  A `multiPhylo` object from the `ape` package.

- data:

  The model data frame.

- optimize:

  Logical. If TRUE (default), returns precision matrices /
  eigendecompositions.

- quiet:

  Logical. If TRUE, suppresses messages.

- ...:

  Additional arguments passed from `because()`, including:

  - `tree_method`: `"ensemble"` (default; runs `n.chains` MCMC chains
    per tree and combines posteriors with equal weight) or `"dcat"`
    (samples a categorical tree index `K_tree` inside MCMC, using
    eigenvector sign alignment and spectral Gaussian marginalization for
    pure Gaussian models in NumPyro).

  - `n.chains`: Integer number of MCMC chains to run per tree when
    `tree_method = "ensemble"` (default `2L`).

## Value

A list with `structure_object` (standardized trees with `tree_method`
and `n.chains` attributes) and `data_list` (MCMC data).
