# because.phybase

### Phylogenetic Bayesian Structural Equation Models (PhyBaSE) extension for 'because'

`because.phybase` extends the [because](https://because-pkg.github.io) package with support for phylogenetic covariance structures. 
It implements S3 methods for `phylo` and `multiPhylo` objects from the `ape` package, enabling Bayesian Structural Equation Models with phylogeneticaly-structured residuals (PhyBaSE, von Hardenberg & Gonzalez-Voyer, 2025).

## Features

-   **Phylogenetic Covariance**: Incorporates phylogenetic relatedness into `because` models.
-   **Phylogenetic Uncertainty**: Supports `multiPhylo` objects to account for uncertainty in tree topology or branch lengths.
-   **Measurement error**: Accounting for measurement error in traits both providing repeated measures as well as specifying known measurement error variances.
-   **Phylogenetic missing data imputation**: informs imputation of missing data using phylogenetic relationships.
-   **Seamless Integration**: Designed to work transparently with `because` via S3 method dispatch.

## Installation

To install `because.phybase`, run:

``` r
remotes::install_github("because-pkg/because.phybase")
```

`because.phybase` depends on `because` which will be installed automatically as a dependency. 
For information on the functions available in the main `because` package and 
for tutorials on how to fully use all its functionalities please refer to the `because` [documentation](https://because-pkg.github.io)




