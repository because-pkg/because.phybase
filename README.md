# because.phybase

### Phylogenetic Bayesian Structural Equation Models (PhyBaSE) extension for 'because'

`because.phybase` extends the [because](https://because-pkg.github.io) package with support for phylogenetic covariance structures. 
It implements S3 methods for `phylo` and `multiPhylo` objects from the `ape` package, enabling Bayesian Structural Equation Models with phylogenetically-structured residuals (PhyBaSE, von Hardenberg & Gonzalez-Voyer, 2025).

## Features

-   **Phylogenetic Covariance**: Incorporates phylogenetic relatedness into `because` models.
-   **Phylogenetic Uncertainty**: Supports `multiPhylo` objects to account for uncertainty in tree topology or branch lengths.
-   **Measurement error**: Accounting for measurement error in traits both providing repeated measures as well as specifying known measurement error variances.
-   **Phylogenetic missing data imputation**: informs imputation of missing data using phylogenetic relationships.
-   **Seamless Integration**: Designed to work transparently with `because` via S3 method dispatch.

## Installation

To install the **stable release** (`v1.0.2`), run:

``` r
remotes::install_github("because-pkg/because.phybase@v1.0.2")
```

To install the **latest development version**, run:

``` r
remotes::install_github("because-pkg/because.phybase")
```

`because.phybase` depends on `because` which will be installed automatically as a dependency. 
For information on the functions available in the main `because` package and 
for tutorials on how to fully use all its functionalities please refer to the `because` [documentation](https://because-pkg.github.io).

## Citation

To cite `because.phybase` in publications:

* **Journal article**:
  > von Hardenberg, A., & Gonzalez-Voyer, A. (2025). PhyBaSE: A Bayesian structural equation model approach to causal inference in phylogenetic comparative analyses. *Methods in Ecology and Evolution*, 16(6), 1136–1148. <https://doi.org/10.1111/2041-210X.70044>

* **R package**:
  > von Hardenberg, A., & Gonzalez-Voyer, A. (2026). *because.phybase: Phylogenetic Structural Equation Models (PhyBaSE) in R*. R package version 1.0.2, <https://because-pkg.github.io/because.phybase/>.




