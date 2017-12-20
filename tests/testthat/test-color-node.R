
## Marta R. Hidalgo

library(hipathia)
context("Color node per DE")

data("results")
data("brca_design")
sample_group <- brca_design[colnames(results$all$path.vals),"group"]
mgi <- load.pathways("hsa", pathways.list = c("hsa03320", "hsa04012"))
colors_de <- node.color.per.differential.expression(results,
                                                    mgi,
                                                    sample_group,
                                                    "Tumor",
                                                    "Normal")

test_that("Classes are correct", {
    expect_is(colors_de, "list")
    expect_true(all(sapply(colors_de, class) == "character"))
    expect_equal(sum(sapply(unlist(colors_de), is.null)), 0)
})

test_that("Names are correct", {
    expect_equal(length(setdiff(names(colors_de), names(mgi$pathigraphs))), 0)
    expect_equal(length(setdiff(names(mgi$pathigraphs), names(colors_de))), 0)
})
