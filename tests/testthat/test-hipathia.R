
## Marta R. Hidalgo

library(hipathia)
context("Hipathia function")

data("exp_data")
mgi <- load.pathways("hsa", pathways.list = c("hsa03320", "hsa04012"))
results <- hipathia(exp_data[,1, drop = FALSE], mgi)

test_that("Classes are correct", {
    expect_is(results, "list")
    expect_is(results$all$path.vals, "matrix")
    expect_is(results$all$nodes.vals, "matrix")
    expect_is(results$by.path, "list")
    expect_equal(length(results$by.path), 2)
})

test_that("Genes.vals & metaginfo parameters are needed", {
    expect_error(hipathia())
    expect_error(hipathia(exp_data))
    expect_error(hipathia(metaginfo = mgi))
})

test_that("Results are complete", {
    expect_equal(2, length(results))
    expect_equal(2, length(results$all))
    expect_equal(length(mgi$pathigraphs),
                 length(results$by.path))
    expect_equal(ncol(results$all$path.vals), 1)
    expect_equal(ncol(results$all$nodes.vals), 1)
})

test_that("Values are in rank", {
    pv <- results$all$path.vals
    expect_true(all(pv <= 1 & pv >= 0))
    nv <- results$all$nodes.vals
    expect_true(all(nv <= 1 & nv >= 0))
})

test_that("Extreme values are preserved", {
    genes <- mgi$all.genes
    m0 <- matrix(0, ncol = 1, nrow = length(genes), dimnames = list(genes, "0"))
    m1 <- matrix(1, ncol = 1, nrow = length(genes), dimnames = list(genes, "1"))
    m <- cbind(m0, m1)
    results <- hipathia(m, mgi)
    expect_true(all(results$all$path.vals %in% c(0,1)))
})

