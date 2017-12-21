#
# ## Marta R. Hidalgo
#
# library(hipathia)
# context("Create report")
#
# data("exp_data")
# data("brca_design")
# data("comp")
# extdata_folder <- system.file("extdata", package = "hipathia")
# output_folder <- paste0(extdata_folder, "/save_noColors")
# mgi <- load.pathways("hsa", pathways.list = c("hsa03320", "hsa04012"))
# results <- hipathia(exp_data, mgi)
# sample_group <- brca_design[colnames(path_vals),"group"]
# colors_de <- node.color.per.differential.expression(results,
#                                                     mgi,
#                                                     sample_group,
#                                                     "Tumor",
#                                                     "Normal")
# create.report(results, comp, mgi, output_folder)
#
# test_that("Folders are created", {
#     expect_true(file.exists(paste0(output.folder, "/pathway-viewer")))
#     expect_true(file.exists(paste0(output.folder, "/pathway-viewer/conf")))
#     expect_true(file.exists(paste0(output.folder, "/pathway-viewer/css")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/fontawesome")))
#     expect_true(file.exists(paste0(output.folder, "/pathway-viewer/fonts")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/normalize-css")))
#     expect_true(file.exists(paste0(output.folder, "/pathway-viewer/pathways")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/webcomponentsjs")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/index.html")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/pathway-viewer.html")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/report_legend.png")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/report_logo.png")))
#     expect_true(file.exists(paste0(output.folder,
#                                    "/pathway-viewer/report_title.png")))
# })
#
# test_that("All pathway files have been created", {
#     expect_equal(length(list.files(paste0(output.folder,
#                                           "/pathway-viewer/pathways"))),
#                  3 * length(mgi$pathigraphs) + 1)
# })
#
