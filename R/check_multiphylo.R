#' Check Multi-Tree Convergence and Phylogenetic Sensitivity Diagnostics
#'
#' Evaluates tree sampling, within-tree MCMC convergence, and between-tree
#' phylogenetic sensitivity. This is a phylogenetic alias for [because::check_multicov()].
#'
#' @param model A `because` model object fitted with multiple trees or covariance structures.
#' @param plot Logical. If `TRUE` (default), displays a diagnostic plot.
#' @param ... Additional arguments passed to [because::check_multicov()].
#'
#' @return An object of class `check_multicov`.
#' @export
check_multiphylo <- function(model, plot = TRUE, ...) {
  because::check_multicov(model = model, plot = plot, ...)
}
