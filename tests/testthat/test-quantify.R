
## Marta R. Hidalgo

library(hipathia)
context("Quantify terms")

data("results")
pathways <- load_pathways("hsa", pathways_list = c("hsa03320", "hsa04012"))
uniprot_vals <- quantify_terms(results, pathways, "uniprot")
go_vals <- quantify_terms(results, pathways, "GO")
uniprot_vals_mat <- quantify_terms(results, pathways, "uniprot", out_matrix = TRUE)
go_vals_mat <- quantify_terms(results, pathways, "GO", out_matrix = TRUE)

test_that("Resulting object is a matrix", {
    expect_is(uniprot_vals, "SummarizedExperiment")
    expect_is(go_vals, "SummarizedExperiment")
    expect_is(uniprot_vals_mat, "matrix")
    expect_is(go_vals_mat, "matrix")
})

test_that("Colnames are preserved", {
    expect_equal(colnames(results)[["paths"]], colnames(uniprot_vals))
    expect_equal(colnames(results)[["paths"]], colnames(go_vals))
})

test_that("Values are in rank", {
    expect_true(all(uniprot_vals_mat <= 1 & uniprot_vals_mat >= 0))
    expect_true(all(go_vals_mat <= 1 & go_vals_mat >= 0))
    expect_true(all(assay(uniprot_vals) <= 1 & assay(uniprot_vals) >= 0))
    expect_true(all(assay(go_vals) <= 1 & assay(go_vals) >= 0))
})


