# Prepare Structure Data for Phylogenetic Trees

Implements the `prepare_structure_data` S3 method for `phylo` objects.

## Usage

``` r
# S3 method for class 'phylo'
prepare_structure_data(structure, data, optimize = TRUE, quiet = FALSE, ...)
```

## Arguments

- structure:

  A `phylo` object from the `ape` package.

- data:

  The model data frame.

- optimize:

  Logical. If TRUE (default), returns precision matrix.

- quiet:

  Logical. If TRUE, suppresses messages.

- ...:

  Additional arguments (ignored).

## Value

A list with `structure_object` (standardized tree) and `data_list` (JAGS
data).
