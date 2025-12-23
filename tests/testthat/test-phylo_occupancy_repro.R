library(because)
library(because)
library(ape)

set.seed(42)

# 1. Simulate Tree
N_species <- 10
tree <- rtree(N_species)
species_names <- tree$tip.label
tree$node.label <- NULL

# 2. Simulate Data (Flattened)
N_sites <- 20
J_visits <- 3

Habitat <- rnorm(N_sites)
Trait <- rnorm(N_species)
names(Trait) <- species_names

alpha_psi <- rnorm(N_species)
alpha_p <- rnorm(N_species)

# Generate flattened data
Y_flat <- matrix(NA, nrow = N_species * N_sites, ncol = J_visits)
row_idx <- 0
for (s in 1:N_species) {
    p_sp <- plogis(alpha_p[s] + 0.5 * Trait[s])
    psi_sp <- plogis(alpha_psi[s] + 0.5 * Habitat)
    z <- rbinom(N_sites, 1, psi_sp)

    for (i in 1:N_sites) {
        row_idx <- row_idx + 1
        Y_flat[row_idx, ] <- rbinom(J_visits, 1, z[i] * p_sp)
    }
}

df <- expand.grid(Site = 1:N_sites, Species = species_names)
df$Habitat <- Habitat[df$Site]
df$Trait <- Trait[df$Species]
df$Species_ID <- df$Species

# Calculate manual Workaround items
C <- vcv(tree, corr = TRUE)
Prec_Species <- solve(C)
zeros_Species <- rep(0, nrow(C))

data_list <- list(
    Y = Y_flat,
    Habitat = df$Habitat,
    Trait = df$Trait,
    Species = df$Species,
    N = nrow(df),
    N_Species = length(unique(df$Species)), # Workaround 1
    zeros_Species = zeros_Species, # Workaround 2
    Prec_Species = Prec_Species # Workaround 3
)

print("Fitting model with workarounds...")

tryCatch(
    {
        res <- because(
            equations = list(
                # Occupancy:
                # uses id_col="Species" + structure=tree for Phylo effect
                Y ~ Habitat,

                # Detection:
                # Random effect (1|Species). Since 'Species' matches tips and 'structure' is present,
                # because generates phylogenetic RE code. We supply the data manually.
                p_Y ~ Trait + (1 | Species)
            ),
            data = data_list,
            family = c(Y = "occupancy"),
            structure = tree,
            id_col = "Species",
            n.iter = 100,
            n.chains = 1,
            quiet = FALSE
        )
        print("Success!")
        print(summary(res))
    },
    error = function(e) {
        print(paste("Error:", e$message))
    }
)
