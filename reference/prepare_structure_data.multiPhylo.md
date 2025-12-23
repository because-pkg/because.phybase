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

  Logical. If TRUE (default), returns precision matrices.

- quiet:

  Logical. If TRUE, suppresses messages.

- ...:

  Additional arguments (ignored).

## Value

A list with `structure_object` (standardized trees) and `data_list`
(JAGS data).
