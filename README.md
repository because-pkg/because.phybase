# because.phybase

### Phylogenetic Structure Module for 'because'

**because.phybase** extends the [because](https://github.com/because-pkg/because) package with support for phylogenetic covariance structures. It implements S3 methods for `phylo` and `multiPhylo` objects from the `ape` package, enabling Bayesian Structural Equation Models (SEMs) with phylogenetically-structured residuals.

## Features

-   **Phylogenetic Covariance**: Automatically incorporates phylogenetic relatedness into `because` models.
-   **Phylogenetic Uncertainty**: Supports `multiPhylo` objects to account for uncertainty in tree topology or branch lengths.
-   **Seamless Integration**: Designed to work transparently with `because` via S3 method dispatch.

## Installation

To install **because.phybase**, run:

```r
remotes::install_github("because-pkg/because.phybase")
```

Note: This package requires the main **because** package to be installed.

```r
remotes::install_github("because-pkg/because")
```

## Usage

Load the package alongside `because` and pass a `phylo` object to `because_model`:

```r
library(because)
library(because.phybase)
library(ape)

# 1. Load your data and tree
data(my_data)
tree <- read.tree("my_tree.tre")

# 2. Define your structural equations
eqs <- list(
  formula(Trait1 ~ Predictor),
  formula(Trait2 ~ Trait1)
)

# 3. Run model with phylogenetic structure
mod <- because_model(
  equations = eqs,
  data = my_data,
  structures = list(phylo = tree),
  structure_names = "phylo"
)
```

## License

MIT
