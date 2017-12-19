
## Marta R. Hidalgo

library(hipathia)
context("Quantify terms")

data("results")
pathways <- load.pathways("hsa", pathways.list = c("hsa03320", "hsa04012"))
uniprot_vals <- quantify.terms(results, pathways, "uniprot")
go_vals <- quantify.terms(results, pathways, "GO")

test_that("Resulting object is a matrix", {
    expect_is(uniprot_vals, "matrix")
    expect_is(go_vals, "matrix")
})

test_that("Colnames are preserved", {
    expect_equal(colnames(results$all$path.vals), colnames(uniprot_vals))
    expect_equal(colnames(results$all$path.vals), colnames(go_vals))
})

test_that("Values are in rank", {
    expect_true(all(uniprot_vals <= 1 & uniprot_vals >= 0))
    expect_true(all(go_vals <= 1 & go_vals >= 0))
})


