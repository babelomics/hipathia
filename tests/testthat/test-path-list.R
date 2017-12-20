
## Marta R. Hidalgo

library(hipathia)
context("Get pathway list")

data("comp")
mgi <- load.pathways("hsa", pathways.list = c("hsa03320", "hsa04012"))
paths <- get.pathways.list(mgi)

test_that("Class is correct", {
    expect_is(paths, "character")
})

test_that("Length is correct", {
    expect_equal(length(paths), length(mgi$pathigraphs))
})

test_that("Pathways are correct", {
    expect_true(all(paths %in% names(mgi$pathigraphs)))
})
