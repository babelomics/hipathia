
## Marta R. Hidalgo

library(hipathia)
context("Load pathways")

pathways_list <- c("hsa03320", "hsa04012")
mgi <- load_pathways("hsa", pathways_list = pathways_list)

test_that("Species parameter is needed", {
    expect_error(load_pathways())
})

test_that("Pathways object is a list", {
    expect_is(mgi, "list")
})

test_that("Pathways object is correct", {
    expect_equal(length(mgi), 6)
    expect_equal(setdiff(names(mgi),
                         c("all.labelids", "pathigraphs", "all.genes",
                           "path.norm", "eff.norm", "species")),
                 character(0))
    expect_true(is_accepted_species(mgi$species))
    expect_true(all(mgi$all.labelids[,"path.id"] %in% pathways_list))
    expect_true(all(mgi$eff.norm >= 0 & mgi$eff.norm <= 1))
    expect_true(all(mgi$path.norm >= 0 & mgi$path.norm <= 1))
    expect_true(all(sapply(mgi$pathigraphs,
                           function(pg) class(pg$graph) == "igraph")))
    all.effs <- unlist(sapply(mgi$pathigraphs,
                              function(pg) names(pg$effector.subgraphs)))
    all.paths <- unlist(sapply(mgi$pathigraphs,
                              function(pg) names(pg$subgraphs)))
    expect_true(all(all.effs == names(mgi$eff.norm)))
    expect_true(all(all.paths == names(mgi$path.norm)))
    expect_true(all(sapply(mgi$pathigraphs, function(pg) length(pg)) == 11))
})

