test_that("package loads correctly", {
    expect_true(requireNamespace("because.phybase", quietly = TRUE))
    expect_true(requireNamespace("because", quietly = TRUE))
})

test_that("S3 methods are registered", {
    # These methods should exist if ape and because are correctly linked
    expect_true(exists("jags_structure_definition.phylo"))
    expect_true(exists("prepare_structure_data.phylo"))
})
