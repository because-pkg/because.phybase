library(because)
# library(because)
devtools::load_all(".")
library(ape)

set.seed(42)

# 1. Simulate Tree
N_species <- 50 # Increase to 50 as per user scenario
tree <- rtree(N_species)
species_names <- tree$tip.label
tree$node.label <- NULL

# 2. Simulate Data
N_sites <- 20
J_visits <- 3

Habitat <- rnorm(N_sites)
Trait <- rnorm(N_species)
names(Trait) <- species_names

# True Parameters (Logit scale, reasonably small)
alpha_psi_mean <- 0 # Occupancy intercept (~0.5 prob)
alpha_p_mean <- 0 # Detection intercept (~0.5 prob)
beta_psi_hab <- 1.0 # Habitat effect
beta_p_trait <- -0.5 # Trait effect

# Phylo signal on Occupancy
alpha_psi <- MASS::mvrnorm(1, rep(alpha_psi_mean, N_species), 1.0 * vcv(tree))
alpha_p <- rep(alpha_p_mean, N_species) # No phylo on detection for simplicity

# Generate flattened data
Y_flat <- matrix(NA, nrow = N_species * N_sites, ncol = J_visits)
row_idx <- 0

# Store true z for checking
z_true <- numeric(N_species * N_sites)

for (s in 1:N_species) {
    p_sp <- plogis(alpha_p[s] + beta_p_trait * Trait[s])
    # Species-specific occupancy probability (varying by site)

    for (i in 1:N_sites) {
        psi_sp <- plogis(alpha_psi[s] + beta_psi_hab * Habitat[i])

        z <- rbinom(1, 1, psi_sp)

        row_idx <- row_idx + 1
        z_true[row_idx] <- z
        Y_flat[row_idx, ] <- rbinom(J_visits, 1, z * p_sp)
    }
}

df <- expand.grid(Site = 1:N_sites, Species = species_names)
df$Habitat <- Habitat[df$Site]
df$Trait <- Trait[df$Species]
df$Y <- I(Y_flat)

print("Fitting Model with Strong Priors...")

# Define custom priors to prevent explosion
# Key: Use N(0, 1) or N(0, 1.5) for logit-scale parameters
my_priors <- list(
    alpha_Y = "dnorm(0, 1)",
    beta_Y_Habitat = "dnorm(0, 1)",
    alpha_p_Y = "dnorm(0, 1)",
    beta_p_Y_Trait = "dnorm(0, 1)"
)

tryCatch(
    {
        res <- because(
            equations = list(
                Y ~ Habitat,
                p_Y ~ Trait
            ),
            data = df,
            family = c(Y = "occupancy"),
            structure = tree,
            id_col = "Species",
            priors = my_priors, # APPLYING FIX HERE
            n.iter = 500,
            quiet = TRUE, # Should be silent now
            verbose = FALSE # Should be silent
        )

        print("Success! Checking summary stats...")
        summ <- summary(res)$results
        print(summ[, c("Mean", "SD", "2.5%", "97.5%", "Rhat")])

        if (abs(summ["alpha_Y", "Mean"]) > 10) {
            stop("Parameter EXPLOSION detected despite priors!")
        } else {
            print("Parameters are stable.")
        }
    },
    error = function(e) {
        print(paste("Error:", e$message))
    }
)
