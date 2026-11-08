# Attach Evolutionary Model Specifications to Phylogenetic Structures

This helper attaches an `"evo_model"` attribute to a `phylo` or
`multiPhylo` object. By default, `because` assumes a Brownian Motion
(BM) model. You can use this function to specify node-specific
evolutionary models (e.g., Ornstein-Uhlenbeck) for different response
variables within the same causal graph.

## Usage

``` r
evo_model(tree, model = "BM", n_alpha = 50)
```

## Arguments

- tree:

  A `phylo` or `multiPhylo` object.

- model:

  A named character vector mapping response variables to evolutionary
  models. Valid options are `"BM"` (Brownian Motion) and `"OU"`
  (Ornstein-Uhlenbeck). If unnamed and length 1, the model is applied
  globally to all phylogenetically structured responses in the SEM.

- n_alpha:

  Number of points for the OU `alpha` parameter grid (default 50).

## Value

The original tree object with an `"evo_model"` attribute attached.

## Examples

``` r
if (FALSE) { # \dontrun{
# Apply an OU model globally with default grid (50 points)
my_tree <- evo_model(my_tree, "OU")

# Increase grid resolution for better precision
my_tree <- evo_model(my_tree, "OU", n_alpha = 100)

# Apply BM to TraitA and OU to TraitB
my_tree <- evo_model(my_tree, c(TraitA = "BM", TraitB = "OU"))
} # }
```
