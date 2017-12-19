
## Marta R. Hidalgo

library(hipathia)
context("Pathways summary")

data("comp")
pathways_list <- c("hsa03320", "hsa04012")
mgi <- load.pathways("hsa", pathways.list = pathways_list)
pathsum <- get.pathways.summary(comp, pathways)

test_that("Resulting object is a data.frame", {
    expect_is(pathsum, "data.frame")
})

test_that("Pathways summary is correct", {
    # Pathways are only the necessary
    expect_equal(character(0), setdiff(pathsum[,1], pathways_list))
    expect_equal(nrow(pathsum), length(pathways_list))
    expect_equal(ncol(pathsum), 8)
    expect_equal(pathsum$percent.up.paths + pathsum$percent.down.paths,
                 pathsum$percent.significant.paths)
    expect_equal(pathsum$num.up.paths + pathsum$num.down.paths,
                 pathsum$num.significant.paths)
})

