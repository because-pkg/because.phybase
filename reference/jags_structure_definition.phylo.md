# JAGS Structure Definition for Phylogenetic Trees

Implements the `jags_structure_definition` S3 method for `phylo`
objects. This method is provided by the `because.phybase` package to
extend the base `because` engine with phylogenetic covariance support.

## Usage

``` r
# S3 method for class 'phylo'
jags_structure_definition(
  structure,
  variable_name = "err",
  optimize = TRUE,
  precision_parameter = "lambda",
  ...
)
```

## Arguments

- structure:

  A `phylo` object from the `ape` package.

- variable_name:

  Name of the error variable (default "err").

- optimize:

  Logical. If TRUE (default), uses pre-calculated precision matrix.

- precision_parameter:

  Name of the precision parameter (default "lambda").

- ...:

  Additional arguments (ignored).

## Value

A list with `setup_code` and `error_prior` JAGS code strings.
