
## Marta R. Hidalgo

library(hipathia)
context("Get pathway names")

data("comp")
mgi <- load_pathways("hsa", pathways_list = c("hsa03320", "hsa04012"))
names <- get_path_names(mgi, rownames(comp))

test_that("Class is correct", {
    expect_is(names, "character")
})

test_that("Length is correct", {
    expect_equal(length(names), nrow(comp))
})

test_that("Values are not empty", {
    expect_true(!all(names == ":"))
})
