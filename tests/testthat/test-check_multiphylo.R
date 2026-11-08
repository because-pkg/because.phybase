library(because)
library(because.phybase)
library(ape)

testthat::test_that("check_multiphylo runs on multi-tree and multi-backbone ensembles", {
  testthat::skip_if_not_installed("rjags")
  
  set.seed(123)
  N <- 8
  t1 <- rtree(N)
  t2 <- rtree(N)
  t3 <- rtree(N)
  t4 <- rtree(N)
  t1$tip.label <- t2$tip.label <- t3$tip.label <- t4$tip.label <- paste0("sp", 1:N)

  dat <- data.frame(sp = paste0("sp", 1:N), x = rnorm(N), y = rnorm(N))

  fit <- because(
    equations = list(y ~ x),
    data = dat,
    structure = list(
      Jetz = list(t1, t2),
      McTavish = list(t3, t4)
    ),
    id_col = "sp",
    engine = "jags",
    n.chains = 2,
    n.iter = 300,
    n.burnin = 100,
    quiet = TRUE
  )

  diag <- check_multiphylo(fit, plot = FALSE)
  testthat::expect_s3_class(diag, "check_multicov")
  testthat::expect_true(!is.null(diag$variance_decomposition))
  testthat::expect_true(all(c("Parameter", "SD", "Pct_MCMC", "Pct_Struct", "Pct_Group") %in% names(diag$variance_decomposition)))

  pct_sum <- diag$variance_decomposition$Pct_MCMC + 
             diag$variance_decomposition$Pct_Struct + 
             diag$variance_decomposition$Pct_Group
  testthat::expect_true(all(abs(pct_sum - 100.0) <= 0.2))
})
