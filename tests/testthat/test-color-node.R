
## Marta R. Hidalgo

library(hipathia)
context("Color node per DE")

data("results")
data("brca")
sample_group <- colData(brca)[,1]
mgi <- load_pathways("hsa", pathways_list = c("hsa03320", "hsa04012"))
colors_de <- node_color_per_de(results, mgi, sample_group, "Tumor", "Normal")

test_that("Classes are correct", {
    expect_is(colors_de, "list")
    expect_true(all(sapply(colors_de$colors, function(x) is(x, "character"))))
    expect_equal(sum(sapply(unlist(colors_de), is.null)), 0)
})

test_that("Names are correct", {
    expect_equivalent(length(setdiff(names(colors_de$colors), names(mgi$pathigraphs))), 0)
})
