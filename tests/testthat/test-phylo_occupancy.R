library(because)
# ------------------------------------------------------------
#  Multi‑species phylogenetic occupancy model (50 species)
# ------------------------------------------------------------
library(because) # main modelling package
library(ape) # for phylogenetic trees
set.seed(123)
# 1. Simulate a phylogenetic tree for 50 species
N_species <- 50
tree <- rtree(N_species) # random ultrametric tree
species_names <- tree$tip.label
# 2. Site‑level covariate (e.g., habitat quality)
N_sites <- 200
J_visits <- 3
Habitat <- rnorm(N_sites, mean = 0, sd = 1) # centred & scaled
# 3. Species‑level trait (e.g., body size)
Trait <- rnorm(N_species, mean = 0, sd = 1) # centred & scaled
# 4. Phylogenetic covariance matrix (Brownian motion)
C_phylo <- vcv(tree, corr = TRUE) # N_species × N_species
draw_phylo_coefs <- function(sd = 1) {
  MASS::mvrnorm(1, mu = rep(0, N_species), Sigma = sd^2 * C_phylo)
}
# 5. True phylogenetically correlated parameters
alpha_psi <- draw_phylo_coefs(sd = 1.0) # occupancy intercepts
beta_psi <- draw_phylo_coefs(sd = 0.5) # occupancy slope for Habitat
alpha_p <- draw_phylo_coefs(sd = 1.0) # detection intercepts
beta_p <- draw_phylo_coefs(sd = 0.5) # detection slope for Trait
gamma_p <- 0.8 # interaction Trait × Habitat
# 6. Simulate latent occupancy (z) for each species/site
psi_mat <- matrix(NA, nrow = N_species, ncol = N_sites)
z_mat <- matrix(NA, nrow = N_species, ncol = N_sites)
for (s in 1:N_species) {
  logit_psi <- alpha_psi[s] + beta_psi[s] * Habitat
  psi_mat[s, ] <- plogis(logit_psi)
  z_mat[s, ] <- rbinom(N_sites, 1, psi_mat[s, ])
}
# 7. Simulate detection probabilities (p) for each species/site
p_mat <- matrix(NA, nrow = N_species, ncol = N_sites)
for (s in 1:N_species) {
  logit_p <- alpha_p[s] + beta_p[s] * Trait[s] + gamma_p * Trait[s] * Habitat
  p_mat[s, ] <- plogis(logit_p)
}
# 8. Generate detection histories (Y) for each species
Y_obs <- vector("list", N_species)
for (s in 1:N_species) {
  Y_mat <- matrix(0, nrow = N_sites, ncol = J_visits)
  for (i in 1:N_sites) {
    for (j in 1:J_visits) {
      Y_mat[i, j] <- rbinom(1, 1, z_mat[s, i] * p_mat[s, i])
    }
  }
  Y_obs[[s]] <- Y_mat
}
names(Y_obs) <- species_names

# 10. Model equations
eqs <- list(
  Y ~ Habitat, # occupancy part (state)
  p_Y ~ Trait + Trait:Habitat # detection part (observation)
)

# 12. Fit the model
res <- because(
  equations = eqs,
  data = list(
    Y = Y_obs, # list of detection matrices (species dimension)
    Habitat = Habitat,
    Trait = Trait
  ),
  structure = tree, # phylogenetic covariance
  family = c(Y = "occupancy"),
  n.iter = 6000,
  n.burnin = 2000,
  n.thin = 5,
  quiet = FALSE
)
# 13. Summarise results
summ <- summary(res)
alpha_psi_est <- summ$results[
  grep("^alpha_Y", rownames(summ$results)),
  ,
  drop = FALSE
]
alpha_det_est <- summ$results[grep("^alpha_p_Y", rownames(summ$results)), ]
# Quick sanity check: plot true vs. estimated occupancy intercepts
plot(
  alpha_psi,
  alpha_psi_est$Mean,
  xlab = "True occupancy intercept (α_ψ)",
  ylab = "Posterior mean estimate",
  main = "Phylogenetic occupancy: true vs. estimated"
)
abline(0, 1, col = "red", lwd = 2)
