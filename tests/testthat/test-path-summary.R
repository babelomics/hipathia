
## Marta R. Hidalgo

library(hipathia)
context("Pathways summary")

data("comp")
pathways_list <- c("hsa03320", "hsa04012")
mgi <- load_pathways("hsa", pathways_list = pathways_list)
pathsum <- get_pathways_summary(comp, mgi)

test_that("Resulting object is a data.frame", {
    expect_is(pathsum, "data.frame")
})

test_that("Pathways summary is correct", {
    # Pathways are only the necessary
    expect_equal(character(0), setdiff(pathsum[,1], pathways_list))
    expect_equal(nrow(pathsum), length(pathways_list))
    expect_equal(ncol(pathsum), 8)
    expect_equal(pathsum$percent_up_paths + pathsum$percent_down_paths,
                 pathsum$percent_significant_paths)
    expect_equal(pathsum$num_up_paths + pathsum$num_down_paths,
                 pathsum$num_significant_paths)
})

