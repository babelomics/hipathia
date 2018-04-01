
## Marta R. Hidalgo

library(hipathia)
context("Hipathia function")

data("exp_data")
mgi <- load_pathways("hsa", pathways_list = c("hsa03320", "hsa04012"))
results <- hipathia(exp_data[,1, drop = FALSE], mgi)

test_that("Classes are correct", {
    expect_is(results, "MultiAssayExperiment")
    expect_is(results[["paths"]], "SummarizedExperiment")
    expect_is(results[["nodes"]], "SummarizedExperiment")
    expect_is(colData(results), "DataFrame")
    expect_is(assay(results, "paths"), "matrix")
    expect_is(assay(results, "nodes"), "matrix")
})

test_that("Genes.vals & metaginfo parameters are needed", {
    expect_error(hipathia())
    expect_error(hipathia(exp_data))
    expect_error(hipathia(metaginfo = mgi))
})

test_that("Results are complete", {
    expect_equal(length(assays(results)), 2)
    expect_equal(1, ncol(results[["paths"]]))
    expect_equal(1, ncol(results[["nodes"]]))
    expect_equal(sum(sapply(mgi$pathigraphs, function(p){
                    length(p$effector.subgraphs)})),
                 nrow(results[["paths"]]))
    expect_equal(nrow(mgi$all.labelids), nrow(results[["nodes"]]))
})

test_that("Values are in rank", {
    pv <- assay(results[["paths"]])
    expect_true(all(pv <= 1 & pv >= 0))
    nv <- assay(results[["nodes"]])
    expect_true(all(nv <= 1 & nv >= 0))
})

test_that("Extreme values are preserved", {
    genes <- mgi$all.genes
    m0 <- matrix(0, ncol = 1, nrow = length(genes), dimnames = list(genes, "0"))
    m1 <- matrix(1, ncol = 1, nrow = length(genes), dimnames = list(genes, "1"))
    m <- cbind(m0, m1)
    results <- hipathia(m, mgi)
    expect_true(all(assay(results[["paths"]]) %in% c(0,1)))
})

