#' JAGS Structure Definition for Phylogenetic Trees
#'
#' Implements the `jags_structure_definition` S3 method for `phylo` objects.
#' This method is provided by the `because.phybase` package to extend the base
#' `because` engine with phylogenetic covariance support.
#'
#' @param structure A `phylo` object from the `ape` package.
#' @param variable_name Name of the error variable (default "err").
#' @param optimize Logical. If TRUE (default), uses pre-calculated precision matrix.
#' @param precision_parameter Name of the precision parameter (default "lambda").
#' @param ... Additional arguments (ignored).
#'
#' @return A list with `setup_code` and `error_prior` JAGS code strings.
#'
#' @importFrom because jags_structure_definition
#' @export
jags_structure_definition.phylo <- function(
    structure,
    variable_name = "err",
    optimize = TRUE,
    precision_parameter = "lambda",
    ...
) {
    if (optimize) {
        # Optimized: Use pre-calculated precision matrix (Prec_phylo)
        setup_code <- c(
            "    # Phylogenetic precision matrix (pre-calculated)"
        )

        error_prior <- paste0(
            "    ",
            precision_parameter,
            " ~ dgamma(0.001, 0.001)\n",
            "    ",
            variable_name,
            "[1:N] ~ dmnorm(zeros[1:N], ",
            precision_parameter,
            " * Prec_phylo[1:N, 1:N])"
        )

        return(list(setup_code = setup_code, error_prior = error_prior))
    } else {
        # Unoptimized: Invert VCV in JAGS
        setup_code <- c(
            "    # Phylogenetic VCV setup (Marginal / Unoptimized)",
            "    # Invert VCV to get base precision matrix",
            "    Sigma_phylo[1:N, 1:N] <- inverse(VCV[1:N, 1:N])"
        )

        error_prior <- paste0(
            "    ",
            precision_parameter,
            " ~ dgamma(0.001, 0.001)\n",
            "    ",
            variable_name,
            "[1:N] ~ dmnorm(zeros[1:N], ",
            precision_parameter,
            " * Sigma_phylo[1:N, 1:N])"
        )

        return(list(setup_code = setup_code, error_prior = error_prior))
    }
}

#' Prepare Structure Data for Phylogenetic Trees
#'
#' Implements the `prepare_structure_data` S3 method for `phylo` objects.
#'
#' @param structure A `phylo` object from the `ape` package.
#' @param data The model data frame.
#' @param optimize Logical. If TRUE (default), returns precision matrix.
#' @param quiet Logical. If TRUE, suppresses messages.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with `structure_object` (standardized tree) and `data_list` (JAGS data).
#'
#' @importFrom ape vcv branching.times
#' @importFrom because prepare_structure_data
#' @export
prepare_structure_data.phylo <- function(
    structure,
    data,
    optimize = TRUE,
    quiet = FALSE,
    ...
) {
    if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required for phylogenetic models.")
    }

    # 1. Standardize Tree (scaling branch lengths)
    phylo_tree <- structure
    max_bt <- max(ape::branching.times(phylo_tree))
    if (abs(max_bt - 1.0) > 0.01) {
        if (!quiet) {
            message(sprintf("Standardizing tree (max_bt: %.2f -> 1.0)", max_bt))
        }
        phylo_tree$edge.length <- phylo_tree$edge.length / max_bt
    }

    # 2. Calculate VCV
    vcv_mat <- ape::vcv(phylo_tree)

    # 3. Check dimensions against data
    N <- nrow(data)

    data_list <- list()

    if (optimize) {
        # Return Precision Matrix P = inv(V)
        if (!quiet) {
            message(
                "Calculating phylogenetic precision matrix (optimize=TRUE)..."
            )
        }
        P <- solve(vcv_mat)
        data_list[["Prec_phylo"]] <- P
    } else {
        # Return VCV
        data_list[["VCV"]] <- vcv_mat
    }

    return(list(
        structure_object = phylo_tree, # Return standardised tree
        data_list = data_list
    ))
}

#' JAGS Structure Definition for Multiple Trees (Phylogenetic Uncertainty)
#'
#' Implements the `jags_structure_definition` S3 method for `multiPhylo` objects.
#'
#' @param structure A `multiPhylo` object from the `ape` package.
#' @param variable_name Name of the error variable (default "err").
#' @param optimize Logical. If TRUE (default), uses pre-calculated precision matrices.
#' @param precision_parameter Name of the precision parameter (default "lambda").
#' @param ... Additional arguments (ignored).
#'
#' @return A list with `setup_code` and `error_prior` JAGS code strings.
#'
#' @importFrom because jags_structure_definition
#' @export
jags_structure_definition.multiPhylo <- function(
    structure,
    variable_name = "err",
    optimize = TRUE,
    precision_parameter = "lambda",
    ...
) {
    if (optimize) {
        setup_code <- c(
            "    # Multi-tree phylogenetic Precision setup",
            "    # Prec_multiPhylo[,,k] is passed as data",
            "    # We select the k-th precision matrix: Prec_multiPhylo[1:N, 1:N, K]"
        )
        error_prior <- paste0(
            "    ",
            precision_parameter,
            " ~ dgamma(0.001, 0.001)\n",
            "    ",
            variable_name,
            "[1:N] ~ dmnorm(zeros[1:N], ",
            precision_parameter,
            " * Prec_multiPhylo[1:N, 1:N, K])"
        )
    } else {
        setup_code <- c(
            "    # Multi-tree phylogenetic VCV setup",
            "    # K is the selected tree index (categorical)",
            "    for(k in 1:Ntree) {",
            "       Sigma_phylo[1:N, 1:N, k] <- inverse(multiVCV[1:N, 1:N, k])",
            "    }"
        )
        error_prior <- paste0(
            "    ",
            precision_parameter,
            " ~ dgamma(0.001, 0.001)\n",
            "    ",
            variable_name,
            "[1:N] ~ dmnorm(zeros[1:N], ",
            precision_parameter,
            " * Sigma_phylo[1:N, 1:N, K])"
        )
    }

    return(list(setup_code = setup_code, error_prior = error_prior))
}

#' Prepare Structure Data for Multiple Trees
#'
#' Implements the `prepare_structure_data` S3 method for `multiPhylo` objects.
#'
#' @param structure A `multiPhylo` object from the `ape` package.
#' @param data The model data frame.
#' @param optimize Logical. If TRUE (default), returns precision matrices.
#' @param quiet Logical. If TRUE, suppresses messages.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with `structure_object` (standardized trees) and `data_list` (JAGS data).
#'
#' @importFrom ape vcv branching.times
#' @importFrom because prepare_structure_data
#' @export
prepare_structure_data.multiPhylo <- function(
    structure,
    data,
    optimize = TRUE,
    quiet = FALSE,
    ...
) {
    if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required for phylogenetic models.")
    }

    n_trees <- length(structure)

    # 1. Standardize all trees
    structure <- lapply(structure, function(phylo_tree) {
        max_bt <- max(ape::branching.times(phylo_tree))
        if (abs(max_bt - 1.0) > 0.01) {
            phylo_tree$edge.length <- phylo_tree$edge.length / max_bt
        }
        return(phylo_tree)
    })

    # 2. Compute matrices
    # Use first tree to get N
    first_vcv <- ape::vcv(structure[[1]])
    N <- nrow(first_vcv)

    data_list <- list()
    data_list[["Ntree"]] <- n_trees

    if (optimize) {
        # Array of Precision Matrices
        if (!quiet) {
            message(
                "Calculating array of precision matrices (optimize=TRUE)..."
            )
        }
        Prec_multi <- array(NA, dim = c(N, N, n_trees))
        for (i in 1:n_trees) {
            vcv <- ape::vcv(structure[[i]])
            Prec_multi[,, i] <- solve(vcv)
        }
        data_list[["Prec_multiPhylo"]] <- Prec_multi
    } else {
        # Array of VCV Matrices
        multiVCV <- array(NA, dim = c(N, N, n_trees))
        for (i in 1:n_trees) {
            multiVCV[,, i] <- ape::vcv(structure[[i]])
        }
        data_list[["multiVCV"]] <- multiVCV
    }

    return(list(
        structure_object = structure,
        data_list = data_list
    ))
}
