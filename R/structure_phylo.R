#' Attach Evolutionary Model Specifications to Phylogenetic Structures
#'
#' This helper attaches an `"evo_model"` attribute to a `phylo` or `multiPhylo` object.
#' By default, `because` assumes a Brownian Motion (BM) model. You can use this
#' function to specify node-specific evolutionary models (e.g., Ornstein-Uhlenbeck)
#' for different response variables within the same causal graph.
#'
#' @param tree A `phylo` or `multiPhylo` object.
#' @param model A named character vector mapping response variables to evolutionary models.
#'   Valid options are `"BM"` (Brownian Motion) and `"OU"` (Ornstein-Uhlenbeck).
#'   If unnamed and length 1, the model is applied globally to all phylogenetically
#'   structured responses in the SEM.
#' @param n_alpha Number of points for the OU `alpha` parameter grid (default 50).
#'
#' @return The original tree object with an `"evo_model"` attribute attached.
#' @export
#'
#' @examples
#' \dontrun{
#' # Apply an OU model globally with default grid (50 points)
#' my_tree <- evo_model(my_tree, "OU")
#'
#' # Increase grid resolution for better precision
#' my_tree <- evo_model(my_tree, "OU", n_alpha = 100)
#'
#' # Apply BM to TraitA and OU to TraitB
#' my_tree <- evo_model(my_tree, c(TraitA = "BM", TraitB = "OU"))
#' }
evo_model <- function(tree, model = "BM", n_alpha = 50) {
    if (!inherits(tree, "phylo") && !inherits(tree, "multiPhylo")) {
        stop("'tree' must be a phylo or multiPhylo object.")
    }

    # Validate models
    valid_models <- c("BM", "OU")
    if (any(!model %in% valid_models)) {
        stop(sprintf(
            "Invalid evolutionary model(s). Valid options: %s",
            paste(valid_models, collapse = ", ")
        ))
    }

    attr(tree, "evo_model") <- model
    attr(tree, "n_alpha") <- n_alpha
    return(tree)
}

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
    ...
) {
    args <- list(...)
    loop_bound <- args$loop_bound
    if (is.null(loop_bound)) loop_bound <- "N"
    zeros_name <- args$zeros_name
    if (is.null(zeros_name)) zeros_name <- "zeros"
    i_index <- args$i_index
    if (is.null(i_index)) i_index <- "i"
    s_name <- args$s_name
    if (is.null(s_name)) s_name <- "phylo"

    # Unique random effect node name to avoid collision with response variable
    err_var <- paste0("err_", s_name, "_", variable_name)

    # Unique parameter naming to avoid JAGS collisions
    prec_param <- args$precision_parameter
    if (is.null(prec_param)) {
        prec_param <- paste0("lambda_", s_name, "_", variable_name)
    }
    sig_param <- sub("lambda_", "sigma_", prec_param)
    raw_var <- paste0("err_raw_", s_name, "_", variable_name)

    # Unique loop variable for this structure and response
    j_idx <- paste0("j_", s_name, "_", variable_name)

    evo_model <- attr(structure, "evo_model")
    has_ou <- FALSE
    if (!is.null(evo_model)) {
        has_ou <- "OU" %in% evo_model
    }

    if (optimize) {
        # Optimized: Use pre-calculated precision matrix
        setup_code <- c(
            "    # Phylogenetic precision matrix (pre-calculated)"
        )

        if (variable_name == "err" && has_ou) {
            # Global setup block
            setup_code <- c(
                setup_code,
                "    for(k in 1:N_alpha) {",
                "        prior_alpha_probs[k] <- 1 / N_alpha",
                "    }"
            )
        }

        # Check if called for a specific variable
        is_ou <- FALSE
        if (variable_name != "err" && !is.null(evo_model)) {
            # Extract response name, e.g. "u_std_Y_OU_phylo" -> "Y_OU"
            var_base <- sub("^u_std_(.*)_phylo$", "\\1", variable_name)
            var_base <- sub("^u_(.*)_phylo$", "\\1", var_base)

            if (length(evo_model) == 1 && is.null(names(evo_model))) {
                is_ou <- (evo_model == "OU")
            } else if (var_base %in% names(evo_model)) {
                is_ou <- (evo_model[var_base] == "OU")
            }
        }

        if (is_ou) {
            var_base <- sub("^u_std_(.*)_phylo$", "\\1", variable_name)
            var_base <- sub("^u_(.*)_phylo$", "\\1", var_base)
            setup_code <- c(
                setup_code,
                sprintf("    # Ornstein-Uhlenbeck parameter for %s", var_base),
                sprintf("    idx_alpha_%s ~ dcat(prior_alpha_probs)", var_base),
                sprintf(
                    "    ou_alpha_%s <- alpha_vals[idx_alpha_%s]",
                    var_base,
                    var_base
                )
            )
            model_lines <- paste0(
                "    ", prec_param, " ~ dgamma(0.001, 0.001)\n",
                "    ", sig_param, " <- 1/sqrt(", prec_param, ")\n",
                "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
                "Prec_phylo_OU[1:", loop_bound, ", 1:", loop_bound, ", idx_alpha_", var_base, "])\n",
                "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * ", sig_param, "\n",
                "    }"
            )
            prec_index <- paste0(
                "Prec_phylo_OU[1:", loop_bound, ", 1:", loop_bound, ", idx_alpha_",
                var_base,
                "]"
            )
        } else {
            model_lines <- paste0(
                "    ", prec_param, " ~ dgamma(0.001, 0.001)\n",
                "    ", sig_param, " <- 1/sqrt(", prec_param, ")\n",
                "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
                "Prec_phylo[1:", loop_bound, ", 1:", loop_bound, "])\n",
                "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * ", sig_param, "\n",
                "    }"
            )
            prec_index <- NULL
        }

        term_str <- paste0(err_var, "[", i_index, "]")

        return(list(
            setup_code = setup_code,
            model_lines = model_lines,
            term = term_str,
            prec_index = prec_index
        ))
    } else {
        # Unoptimized: Invert VCV in JAGS
        setup_code <- c(
            "    # Phylogenetic VCV setup (Marginal / Unoptimized)",
            "    # Invert VCV to get base precision matrix",
            "    Sigma_phylo[1:", loop_bound, ", 1:", loop_bound, "] <- inverse(VCV[1:", loop_bound, ", 1:", loop_bound, "])"
        )

        model_lines <- paste0(
            "    ", prec_param, " ~ dgamma(0.001, 0.001)\n",
            "    ", sig_param, " <- 1/sqrt(", prec_param, ")\n",
            "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
            "Sigma_phylo[1:", loop_bound, ", 1:", loop_bound, "])\n",
            "    for(", j_idx, " in 1:", loop_bound, ") {\n",
            "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * ", sig_param, "\n",
            "    }"
        )

        term_str <- paste0(err_var, "[", i_index, "]")

        return(list(setup_code = setup_code, model_lines = model_lines, term = term_str))
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
    N <- nrow(vcv_mat)

    data_list <- list()

    evo_model <- attr(structure, "evo_model")
    has_ou <- "OU" %in% evo_model

    if (optimize) {
        if (!quiet) {
            message(
                "Calculating phylogenetic precision matrix (optimize=TRUE)..."
            )
        }
        P <- solve(vcv_mat)
        data_list[["Prec_phylo"]] <- P

        if (has_ou) {
            if (!quiet) {
                message("Calculating phylogenetic OU precision grid...")
            }
            D <- ape::cophenetic.phylo(phylo_tree)
            N_alpha <- attr(structure, "n_alpha")
            if (is.null(N_alpha)) {
                N_alpha <- 50
            }
            # From large alpha (fast evolution / little correlation) to small alpha (BM)
            # Log-spaced grid
            alpha_vals <- exp(seq(log(0.001), log(10), length.out = N_alpha))

            Prec_phylo_OU <- array(NA, dim = c(N, N, N_alpha))
            for (k in 1:N_alpha) {
                alpha <- alpha_vals[k]
                if (alpha < 1e-5) {
                    V <- vcv_mat
                } else {
                    V <- exp(-alpha * D) *
                        (1 - exp(-2 * alpha * vcv_mat)) /
                        (2 * alpha)
                }

                # Small ridge for numerical stability
                diag(V) <- diag(V) + 1e-6
                Prec_phylo_OU[,, k] <- solve(V)
            }

            data_list[["Prec_phylo_OU"]] <- Prec_phylo_OU
            data_list[["N_alpha"]] <- N_alpha
            data_list[["alpha_vals"]] <- alpha_vals
        }
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
