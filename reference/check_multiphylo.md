# Check Multi-Tree Convergence and Phylogenetic Sensitivity Diagnostics

Evaluates tree sampling, within-tree MCMC convergence, and between-tree
phylogenetic sensitivity. This is a phylogenetic alias for
[`because::check_multicov()`](https://rdrr.io/pkg/because/man/check_multicov.html).

## Usage

``` r
check_multiphylo(model, plot = TRUE, ...)
```

## Arguments

- model:

  A `because` model object fitted with multiple trees or covariance
  structures.

- plot:

  Logical. If `TRUE` (default), displays a diagnostic plot.

- ...:

  Additional arguments passed to
  [`because::check_multicov()`](https://rdrr.io/pkg/because/man/check_multicov.html).

## Value

An object of class `check_multicov`.
