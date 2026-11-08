# NumPyro Structure Definition for Phylogenetic Trees

Implements the `numpyro_structure_definition` S3 method for `phylo`
objects. Returns the JAX/NumPyro Python code to compute the Cholesky
factor of the phylogenetic variance-covariance matrix dynamically scaled
by Pagel's lambda.

## Usage

``` r
# S3 method for class 'phylo'
numpyro_structure_definition(structure, variable_name = "err", ...)
```

## Arguments

- structure:

  A `phylo` object.

- variable_name:

  Name of the variable.

- ...:

  Additional arguments.

## Value

A character string of Python code.
