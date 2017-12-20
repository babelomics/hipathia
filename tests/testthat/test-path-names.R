
## Marta R. Hidalgo

library(hipathia)
context("Get pathway names")

data("comp")
mgi <- load.pathways("hsa", pathways.list = c("hsa03320", "hsa04012"))
names <- get.path.names(mgi, rownames(comp))

test_that("Class is correct", {
    expect_is(names, "character")
})

test_that("Length is correct", {
    expect_equal(length(names), nrow(comp))
})

test_that("Values are not empty", {
    expect_true(!all(names == ":"))
})
