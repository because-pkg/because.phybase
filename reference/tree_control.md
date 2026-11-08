# Configure Multi-Tree (`multiPhylo`) Sampling and Ensemble Options

Attaches multi-tree sampling specifications to a `multiPhylo` object. By
default in `because.phybase` (\>= 1.1.0), `multiPhylo` objects use
`"ensemble"` mode (Bayesian Model Averaging across trees with
deterministic chain-to-tree assignment, defaulting to `n.chains = 2`
chains per tree so within-tree convergence can be separated from
between-tree phylogenetic sensitivity). Users can optionally select
`method = "dcat"` to use the categorical index variable
(`K_tree ~ dcat(p[])`).

## Usage

``` r
tree_control(
  tree,
  method = c("ensemble", "dcat"),
  n.chains = 2L,
  backbone = NULL
)
```

## Arguments

- tree:

  A `multiPhylo` object from the `ape` package, or a named list of
  `multiPhylo` / `phylo` objects representing multiple phylogenetic
  backbones (e.g., `list(Jetz = jetz_trees, McTavish = mctav_trees)`).

- method:

  Character string specifying the multi-tree sampling strategy:

  - `"ensemble"` (default): Runs replicate MCMC chains assigned
    deterministically to individual trees to produce an exact Bayesian
    Model Average (BMA).

  - `"dcat"`: Samples a discrete tree index dynamically via Gibbs
    sampling.

- n.chains:

  Integer. Number of MCMC chains to run per tree in `"ensemble"` mode
  (default 2).

- backbone:

  Optional character vector specifying the backbone name for each tree.
  Automatically inferred if `tree` is a named list.

## Value

The `multiPhylo` object with `"tree_method"` and `"n.chains"` attributes
attached.

The `multiPhylo` object with `"tree_method"`, `"n.chains"`, and optional
`"backbone"` attributes attached.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default ensemble mode with 2 chains per tree
trees <- tree_control(my_multiphylo, method = "ensemble", n.chains = 2)

# Legacy categorical tree-index mode
trees_dcat <- tree_control(my_multiphylo, method = "dcat")
} # }
if (FALSE) { # \dontrun{
# Standard ensemble mode
trees <- tree_control(my_multiphylo, method = "ensemble", n.chains = 2)

# Multi-backbone ensemble mode
trees_bb <- tree_control(list(Jetz = jetz_trees, McTavish = mctav_trees), n.chains = 2)
} # }
```
