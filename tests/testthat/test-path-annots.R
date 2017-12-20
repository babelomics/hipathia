
## Marta R. Hidalgo

library(hipathia)
context("Get pathway annotations")

data("comp")
mgi <- load.pathways("hsa", pathways.list = c("hsa03320", "hsa04012"))
names <- get.pathways.annotations(rownames(comp), mgi, "uniprot")
names_col <- get.pathways.annotations(rownames(comp), mgi, "uniprot", 
                                      collapse = T)

test_that("Classes are correct", {
    expect_is(names, "data.frame")
    expect_is(names_col, "data.frame")
})

test_that("Dims are correct", {
    expect_equal(ncol(names), 2)
    expect_equal(ncol(names_col), 2)
    expect_equal(length(unique(names$paths)), nrow(comp))
    expect_equal(nrow(names_col), nrow(comp))
})

test_that("Annotated functions are correct", {
    expect_true(all(names$paths %in% rownames(comp)))
    expect_true(all(names_col$paths %in% rownames(comp)))
})
