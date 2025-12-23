# JAGS Structure Definition for Multiple Trees (Phylogenetic Uncertainty)

Implements the `jags_structure_definition` S3 method for `multiPhylo`
objects.

## Usage

``` r
# S3 method for class 'multiPhylo'
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

  A `multiPhylo` object from the `ape` package.

- variable_name:

  Name of the error variable (default "err").

- optimize:

  Logical. If TRUE (default), uses pre-calculated precision matrices.

- precision_parameter:

  Name of the precision parameter (default "lambda").

- ...:

  Additional arguments (ignored).

## Value

A list with `setup_code` and `error_prior` JAGS code strings.
