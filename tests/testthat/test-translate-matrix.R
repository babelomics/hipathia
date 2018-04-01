
## Marta R. Hidalgo

library(hipathia)
context("Label translation function")

data("brca_data")
trans_data <- translate_matrix(brca_data, "hsa")
mgi <- load_pathways("hsa")

test_that("Resulting object is a matrix", {
    expect_is(trans_data, "matrix")
})

test_that("Colnames are preserved", {
    expect_equal(colnames(brca_data), colnames(trans_data))
})

test_that("No NA returned", {
    expect_equal(0, sum(is.na(trans_data)))
})

test_that("All rownames are Entrez genes", {
    expect_true(sum(rownames(trans_data) %in% mgi$all.genes) ==
                    nrow(trans_data))
})

