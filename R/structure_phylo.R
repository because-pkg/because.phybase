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
    precision_parameter = "lambda",
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
    err_var <- paste0("err_", variable_name, "_", s_name)

    # Unique parameter naming to avoid JAGS collisions - UNIFIED naming: tau_u_Structure_Variable
    prec_param <- if (precision_parameter == "lambda") {
        # [UNIFICATION] Modern naming: tau_u_phylo_Variable
        paste0("tau_u_", s_name, "_", variable_name)
    } else {
        precision_parameter
    }
    sig_param <- sub("tau_u_", "sigma_", prec_param)
    raw_var <- paste0("err_raw_", variable_name, "_", s_name)

    # Unique loop variable for this structure and response
    j_idx <- paste0("j_", variable_name, "_", s_name)

    evo_model <- attr(structure, "evo_model")
    has_ou <- FALSE
    if (!is.null(evo_model)) {
        has_ou <- "OU" %in% evo_model
    }

    use_partitioning <- isTRUE(args$use_partitioning)

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

        engine <- args$engine %||% "jags"

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
            lambda_param      <- paste0("lambda_", variable_name)
            sigma_total_param <- paste0("sigma_total_", variable_name)
            sigma_phylo_param <- paste0("sigma_", s_name, "_", variable_name)
            sigma_res_param   <- paste0("sigma_", variable_name, "_res")
            tau_res_param     <- paste0("tau_res_", variable_name)

            if (use_partitioning && engine == "nimble") {
                # [NIMBLE NON-CENTERED] z ~ N(0,I), err_raw = L_phylo %*% z
                # L_phylo = t(chol(VCV)) pre-computed in R, passed as data.
                # ESS on z ~ N(0,I) is maximally efficient (prior = proposal).
                # This avoids the dmnorm funnel that causes poor mixing for large N.
                z_var <- paste0("z_", variable_name, "_", s_name)
                iz_idx <- paste0("i_z_", variable_name, "_", s_name)
                model_lines <- paste0(
                    "    ", lambda_param, " ~ dunif(0, 1) # Pagel lambda for ", variable_name, "\n",
                    "    ", sigma_total_param, " ~ dunif(0, 10) # total phylogenetic SD\n",
                    "    ", prec_param, " <- 1 / (", lambda_param, " * ", sigma_total_param,
                        " * ", sigma_total_param, ") # phylogenetic precision\n",
                    "    ", tau_res_param, " <- 1 / ((1 - ", lambda_param, ") * ",
                        sigma_total_param, " * ", sigma_total_param, ") # residual precision\n",
                    "    ", sigma_phylo_param, " <- sqrt(", lambda_param, ") * ", sigma_total_param, " # phylogenetic SD\n",
                    "    ", sigma_res_param,   " <- sqrt(1 - ", lambda_param, ") * ", sigma_total_param, " # residual SD\n",
                    "    for(", iz_idx, " in 1:", loop_bound, ") {\n",
                    "        ", z_var, "[", iz_idx, "] ~ dnorm(0, 1)\n",
                    "    }\n",
                    "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                    "        ", raw_var, "[", j_idx, "] <- inprod(L_phylo[", j_idx, ", 1:", loop_bound, "], ",
                            z_var, "[1:", loop_bound, "])\n",
                    "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * (1/sqrt(", prec_param, "))\n",
                    "    }"
                )
            } else if (use_partitioning) {
                # [JAGS / NIMBLE-OU] MEE-paper style: lambda + sigma_total + dmnorm(Prec_phylo)
                # NOTE: Do NOT start model_lines with a comment.
                model_lines <- paste0(
                    "    ", lambda_param, " ~ dunif(0, 1) # Pagel lambda for ", variable_name, "\n",
                    "    ", sigma_total_param, " ~ dunif(0, 10) # total phylogenetic SD\n",
                    "    ", prec_param, " <- 1 / (", lambda_param, " * ", sigma_total_param,
                        " * ", sigma_total_param, ") # phylogenetic precision\n",
                    "    ", tau_res_param, " <- 1 / ((1 - ", lambda_param, ") * ",
                        sigma_total_param, " * ", sigma_total_param, ") # residual precision\n",
                    "    ", sigma_phylo_param, " <- sqrt(", lambda_param, ") * ", sigma_total_param, " # phylogenetic SD\n",
                    "    ", sigma_res_param,   " <- sqrt(1 - ", lambda_param, ") * ", sigma_total_param, " # residual SD\n",
                    "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
                    "Prec_phylo[1:", loop_bound, ", 1:", loop_bound, "])\n",
                    "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                    "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * (1/sqrt(", prec_param, "))\n",
                    "    }"
                )
            } else {
                # [ADDITIVE] Standard random effect prior (fallback when multiple structures)
                model_lines <- paste0(
                    "    ", prec_param, " ~ dgamma(1, 1)\n",
                    "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
                    "Prec_phylo[1:", loop_bound, ", 1:", loop_bound, "])\n",
                    "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                    "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * (1/sqrt(", prec_param, "))\n",
                    "    }"
                )
            }
            prec_index <- NULL
        }

        term_str <- paste0(err_var, "[", i_index, "]")

        return(list(
            setup_code = setup_code,
            model_lines = model_lines,
            term = term_str,
            prec_index = prec_index,
            partition_handled = use_partitioning,
            variable_name = variable_name
        ))

    } else {
        # Unoptimized: Invert VCV in JAGS
        setup_code <- c(
            "    # Phylogenetic VCV setup (Marginal / Unoptimized)",
            "    # Invert VCV to get base precision matrix",
            "    Sigma_phylo[1:", loop_bound, ", 1:", loop_bound, "] <- inverse(VCV[1:", loop_bound, ", 1:", loop_bound, "])"
        )

        model_lines <- paste0(
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
    args <- list(...)
    engine <- args$engine %||% "jags"
    if (engine == "numpyro") {
        optimize <- FALSE
    }
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

    # 2b. Align VCV matrix with data factor levels
    # Try to find target_order either from args$row_ids or from data columns
    target_order <- NULL
    
    if (!is.null(args$row_ids) && setequal(phylo_tree$tip.label, as.character(args$row_ids))) {
        target_order <- levels(as.factor(args$row_ids))
        if (!quiet) message("Aligning phylogeny to provided row_ids (id_col).")
    } else if (!is.null(data)) {
        matching_col <- NULL
        for (col_name in names(data)) {
            if (is.character(data[[col_name]]) || is.factor(data[[col_name]])) {
                # Check if unique elements exactly match tip labels
                if (setequal(phylo_tree$tip.label, as.character(data[[col_name]]))) {
                    matching_col <- col_name
                    break
                }
            }
        }
        if (!is.null(matching_col)) {
            target_order <- levels(as.factor(data[[matching_col]]))
            if (!quiet) message(sprintf("Aligning phylogeny to data column '%s'", matching_col))
        }
    }
    
    if (!is.null(target_order)) {
        vcv_mat <- vcv_mat[target_order, target_order, drop = FALSE]
    } else {
        if (!quiet) {
            warning("Could not find a matching data column to align phylogeny tip labels. Phylogenetic signal may be lost if data is not ordered alphabetically.")
        }
    }

    # 3. Check dimensions against data
    N <- nrow(vcv_mat)

    data_list <- list()

    evo_model <- attr(structure, "evo_model")
    has_ou <- "OU" %in% evo_model

    if (optimize) {
        if (engine == "nimble" && !has_ou) {
            # NIMBLE: pre-compute lower Cholesky factor for non-centered parameterization.
            # z ~ N(0,I) is sampled, err_raw = L_phylo %*% z gives err_raw ~ N(0, VCV).
            # ESS (Elliptical Slice Sampler) on z ~ N(0,I) is far more efficient
            # than ESS on err_raw ~ dmnorm(0, Prec_phylo) for large N.
            if (!quiet) message("Calculating lower Cholesky factor for NIMBLE (optimize=TRUE)...")
            data_list[["L_phylo"]] <- t(chol(vcv_mat))
        } else {
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
    args <- list(...)
    loop_bound <- args$loop_bound
    if (is.null(loop_bound)) loop_bound <- "N"
    zeros_name <- args$zeros_name
    if (is.null(zeros_name)) zeros_name <- "zeros"
    i_index <- args$i_index
    if (is.null(i_index)) i_index <- "i"
    s_name <- args$s_name
    if (is.null(s_name)) s_name <- "multiPhylo"
    engine <- args$engine
    if (is.null(engine)) engine <- "jags"

    # Unique random effect node name to avoid collision with response variable
    err_var <- paste0("err_", variable_name, "_", s_name)

    # Unique parameter naming to avoid JAGS collisions
    prec_param <- if (precision_parameter == "lambda") {
        paste0("tau_u_", s_name, "_", variable_name)
    } else {
        precision_parameter
    }
    sig_param <- sub("tau_u_", "sigma_", prec_param)
    raw_var <- paste0("err_raw_", variable_name, "_", s_name)

    # Unique loop variable for this structure and response
    j_idx <- paste0("j_", variable_name, "_", s_name)

    use_partitioning <- isTRUE(args$use_partitioning)

    if (optimize) {
        if (engine == "nimble") {
            setup_code <- c(
                "    # Multi-tree phylogenetic Cholesky setup (NIMBLE optimized)",
                "    # L_multiPhylo[,,k] is passed as data",
                paste0("    # We select the k-th lower Cholesky factor: L_multiPhylo[1:", loop_bound, ", 1:", loop_bound, ", K]")
            )

            z_var <- paste0("z_", variable_name, "_", s_name)

            model_lines_common <- paste0(
                "    for (i_z_", variable_name, " in 1:", loop_bound, ") {\n",
                "        ", z_var, "[i_z_", variable_name, "] ~ dnorm(0, 1)\n",
                "    }\n",
                "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                "        ", raw_var, "[", j_idx, "] <- inprod(L_multiPhylo[", j_idx, ", 1:", loop_bound, ", K], ", z_var, "[1:", loop_bound, "])\n",
                "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * (1/sqrt(", prec_param, "))\n",
                "    }"
            )

            if (use_partitioning) {
                model_lines <- paste0(
                    "    # Pagel's Lambda Partitioning (Scaling only) for ", variable_name, " and multiPhylo (NIMBLE)\n",
                    model_lines_common
                )
            } else {
                model_lines <- paste0(
                    "    ", prec_param, " ~ dgamma(1, 1)\n",
                    model_lines_common
                )
            }
            prec_index <- paste0("L_multiPhylo[1:", loop_bound, ", 1:", loop_bound, ", K]")
        } else {
            setup_code <- c(
                "    # Multi-tree phylogenetic Precision setup (JAGS Memory Safe)",
                "    # Prec_multiPhylo[,,k] is passed as data",
                paste0("    # We select the k-th precision matrix: Prec_multiPhylo[1:", loop_bound, ", 1:", loop_bound, ", K]")
            )

            if (use_partitioning) {
                # MEE-paper style: direct lambda + sigma_total, pre-computed Prec_multiPhylo[,,K]
                # NOTE: First line must NOT be a comment (safe_add_lines would skip declared_nodes registration).
                lambda_param      <- paste0("lambda_", variable_name)
                sigma_total_param <- paste0("sigma_total_", variable_name)
                sigma_phylo_param <- paste0("sigma_", s_name, "_", variable_name)
                sigma_res_param   <- paste0("sigma_", variable_name, "_res")
                tau_res_param     <- paste0("tau_res_", variable_name)
                model_lines <- paste0(
                    "    ", lambda_param, " ~ dunif(0, 1) # Pagel lambda for ", variable_name, "\n",
                    "    ", sigma_total_param, " ~ dunif(0, 10) # total phylogenetic SD\n",
                    "    ", prec_param, " <- 1 / (", lambda_param, " * ", sigma_total_param,
                        " * ", sigma_total_param, ") # phylogenetic precision\n",
                    "    ", tau_res_param, " <- 1 / ((1 - ", lambda_param, ") * ",
                        sigma_total_param, " * ", sigma_total_param, ") # residual precision\n",
                    "    ", sigma_phylo_param, " <- sqrt(", lambda_param, ") * ", sigma_total_param, " # phylogenetic SD\n",
                    "    ", sigma_res_param, " <- sqrt(1 - ", lambda_param, ") * ", sigma_total_param, " # residual SD\n",
                    "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
                    "Prec_multiPhylo[1:", loop_bound, ", 1:", loop_bound, ", K])\n",
                    "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                    "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * (1/sqrt(", prec_param, "))\n",
                    "    }"
                )
            } else {
                model_lines <- paste0(
                    "    ", prec_param, " ~ dgamma(1, 1)\n",
                    "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
                    "Prec_multiPhylo[1:", loop_bound, ", 1:", loop_bound, ", K])\n",
                    "    for(", j_idx, " in 1:", loop_bound, ") {\n",
                    "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * (1/sqrt(", prec_param, "))\n",
                    "    }"
                )
            }
            prec_index <- paste0("Prec_multiPhylo[1:", loop_bound, ", 1:", loop_bound, ", K]")
        }

        term_str <- paste0(err_var, "[", i_index, "]")

        return(list(
            setup_code = setup_code,
            model_lines = model_lines,
            term = term_str,
            prec_index = prec_index,
            partition_handled = use_partitioning,
            variable_name = variable_name
        ))

    } else {
        setup_code <- c(
            "    # Multi-tree phylogenetic VCV setup",
            "    # K is the selected tree index (categorical)",
            "    for(k in 1:Ntree) {",
            paste0("       Sigma_phylo[1:", loop_bound, ", 1:", loop_bound, ", k] <- inverse(multiVCV[1:", loop_bound, ", 1:", loop_bound, ", k])"),
            "    }"
        )
        model_lines <- paste0(
            "    ", prec_param, " ~ dgamma(10, 10)\n",
            "    ", raw_var, "[1:", loop_bound, "] ~ dmnorm(", zeros_name, "[1:", loop_bound, "], ",
            "Sigma_phylo[1:", loop_bound, ", 1:", loop_bound, ", K])\n",
            "    for(", j_idx, " in 1:", loop_bound, ") {\n",
            "        ", err_var, "[", j_idx, "] <- ", raw_var, "[", j_idx, "] * (1/sqrt(", prec_param, "))\n",
            "    }"
        )

        term_str <- paste0(err_var, "[", i_index, "]")

        return(list(setup_code = setup_code, model_lines = model_lines, term = term_str))
    }
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
    args <- list(...)
    engine <- args$engine %||% "jags"
    if (engine == "numpyro") {
        if (!quiet) message("Using DiscreteHMCGibbs for multiPhylo marginalization in NumPyro...")
    }
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
    class(structure) <- "multiPhylo"

    # 2. Compute matrices
    # Use first tree to get N
    first_vcv <- ape::vcv(structure[[1]])
    
    # Align tree tips with data factor levels
    target_order <- NULL
    
    if (!is.null(args$row_ids) && setequal(structure[[1]]$tip.label, as.character(args$row_ids))) {
        target_order <- levels(as.factor(args$row_ids))
        if (!quiet) message("Aligning multiPhylo to provided row_ids (id_col).")
    } else if (!is.null(data)) {
        matching_col <- NULL
        for (col_name in names(data)) {
            if (is.character(data[[col_name]]) || is.factor(data[[col_name]])) {
                if (setequal(structure[[1]]$tip.label, as.character(data[[col_name]]))) {
                    matching_col <- col_name
                    break
                }
            }
        }
        if (!is.null(matching_col)) {
            target_order <- levels(as.factor(data[[matching_col]]))
            if (!quiet) message(sprintf("Aligning multiPhylo to data column '%s'", matching_col))
        }
    }
    
    if (!is.null(target_order)) {
        first_vcv <- first_vcv[target_order, target_order, drop = FALSE]
    } else {
        if (!quiet) {
            warning("Could not find a matching data column to align multiPhylo tip labels. Phylogenetic signal may be lost if data is not ordered alphabetically.")
        }
    }

    N <- nrow(first_vcv)

    data_list <- list()
    data_list[["Ntree"]] <- n_trees

    if (engine == "numpyro") {
        if (!quiet) message("Calculating array of eigen decompositions for NumPyro...")
        eigvals_multi <- array(NA, dim = c(N, n_trees))
        eigvecs_multi <- array(NA, dim = c(N, N, n_trees))
        for (i in 1:n_trees) {
            vcv <- ape::vcv(structure[[i]])
            if (!is.null(target_order)) vcv <- vcv[target_order, target_order, drop=FALSE]
            eig <- eigen(vcv, symmetric = TRUE)
            eigvals_multi[, i] <- eig$values
            eigvecs_multi[,, i] <- eig$vectors
        }
        data_list[["multiPhylo"]] <- list(
            eigvals = eigvals_multi,
            eigvecs = eigvecs_multi
        )
    } else if (optimize) {
        if (engine == "nimble") {
            # Array of Cholesky Factors for NIMBLE
            if (!quiet) message("Calculating array of lower Cholesky factors for NIMBLE (optimize=TRUE)...")
            L_multi <- array(NA, dim = c(N, N, n_trees))
            for (i in 1:n_trees) {
                vcv <- ape::vcv(structure[[i]])
                if (!is.null(target_order)) vcv <- vcv[target_order, target_order, drop=FALSE]
                L_multi[,, i] <- t(chol(vcv))
            }
            data_list[["L_multiPhylo"]] <- L_multi
        } else {
            # Array of Precision Matrices for JAGS
            if (!quiet) message("Calculating array of precision matrices for JAGS (optimize=TRUE)...")
            Prec_multi <- array(NA, dim = c(N, N, n_trees))
            for (i in 1:n_trees) {
                vcv <- ape::vcv(structure[[i]])
                if (!is.null(target_order)) vcv <- vcv[target_order, target_order, drop=FALSE]
                Prec_multi[,, i] <- solve(vcv)
            }
            data_list[["Prec_multiPhylo"]] <- Prec_multi
        }
    } else {
        # Array of VCV Matrices
        multiVCV <- array(NA, dim = c(N, N, n_trees))
        for (i in 1:n_trees) {
            vcv <- ape::vcv(structure[[i]])
            if (!is.null(target_order)) vcv <- vcv[target_order, target_order, drop=FALSE]
            multiVCV[,, i] <- vcv
        }
        data_list[["VCV_multiPhylo"]] <- multiVCV
    }

    return(list(
        structure_object = structure,
        data_list = data_list
    ))
}

#' NumPyro Structure Definition for Phylogenetic Trees
#'
#' @description
#' Implements the `numpyro_structure_definition` S3 method for `phylo` objects.
#' Returns the JAX/NumPyro Python code to compute the Cholesky factor of the
#' phylogenetic variance-covariance matrix dynamically scaled by Pagel's lambda.
#'
#' @param structure A `phylo` object.
#' @param variable_name Name of the variable.
#' @param ... Additional arguments.
#'
#' @return A character string of Python code.
#' @exportS3Method because::numpyro_structure_definition phylo
numpyro_structure_definition.phylo <- function(structure, variable_name = "err", ...) {
    py_code <- "
def phylo_transform(numpyro, jnp, jax, dist, var, group_name, num_groups, matrix, z_raw, sigma, shared_state):
    import numpy as np
    
    # Eigendecomposition is computed EXACTLY ONCE during model tracing
    # This avoids doing an O(N^3) Cholesky decomposition at every leapfrog step!
    eigvals, eigvecs = np.linalg.eigh(np.array(matrix))
    eigvals = jnp.array(eigvals)
    eigvecs = jnp.array(eigvecs)
    
    # Prevent negative eigenvalues due to numerical errors
    eigvals_safe = jnp.maximum(eigvals, 1e-8)
    
    # Scale z_raw by the structural standard deviations (sqrt of eigenvalues) and rotate back
    z_scaled = z_raw * jnp.sqrt(eigvals_safe)
    z_group = jnp.dot(eigvecs, z_scaled)
    
    # Return z_group and the UNCHANGED sigma.
    # because_py will do: u_group = z_group * sigma
    # This results in covariance = sigma^2 * Phylo_Cov
    return z_group, sigma
"
    return(py_code)
}

#' Implements the `numpyro_structure_definition` S3 method for `multiPhylo` objects.
#'
#' @param structure The `multiPhylo` object.
#' @param variable_name The response variable name.
#' @param ... Additional arguments.
#'
#' @exportS3Method because::numpyro_structure_definition multiPhylo
numpyro_structure_definition.multiPhylo <- function(structure, variable_name = "err", ...) {
    args <- list(...)
    group_name <- args$s_name
    if (is.null(group_name)) group_name <- "multiPhylo"
    
    py_code <- "
def multiPhylo_transform(numpyro, jnp, jax, dist, var, group_name, num_groups, matrix_dict, z_raw, sigma, shared_state):
    eigvals_all = jnp.array(matrix_dict['eigvals'])
    eigvecs_all = jnp.array(matrix_dict['eigvecs'])
    n_trees = eigvals_all.shape[-1]
    
    k_name = f'K_tree_{group_name}'
    if k_name not in shared_state:
        shared_state[k_name] = numpyro.sample(k_name, dist.Categorical(probs=jnp.ones(n_trees)/n_trees))
    K = shared_state[k_name]
    
    # Safe indexing for arbitrary batch dimensions from enumeration
    # eigvals_all shape: (N, n_trees) -> (N, ...)
    eigvals = jnp.take(eigvals_all, K, axis=-1)
    # Move batch dimensions to the front if they exist: shape becomes (..., N)
    if eigvals.ndim > 1:
        eigvals = jnp.moveaxis(eigvals, 0, -1)
    
    # eigvecs_all shape: (N, N, n_trees) -> (N, N, ...)
    eigvecs = jnp.take(eigvecs_all, K, axis=-1)
    # Move batch dimensions to the front: shape becomes (..., N, N)
    if eigvecs.ndim > 2:
        eigvecs = jnp.moveaxis(eigvecs, [0, 1], [-2, -1])
    
    eigvals_safe = jnp.maximum(eigvals, 1e-8)
    
    # z_raw has shape (..., N). eigvals_safe has shape (..., N).
    z_scaled = z_raw * jnp.sqrt(eigvals_safe)
    
    # Batched matrix-vector multiplication
    z_group = jnp.einsum('...ij,...j->...i', eigvecs, z_scaled)
    
    # Return z_group and the UNCHANGED sigma.
    return z_group, sigma
"
    return(py_code)
}

#' Get Structure Name Hook for multiPhylo
#'
#' @param structure The multiPhylo structure.
#' @param ... Additional arguments.
#' @exportS3Method because::get_structure_name_hook multiPhylo
get_structure_name_hook.multiPhylo <- function(structure, ...) {
    return("multiPhylo")
}
