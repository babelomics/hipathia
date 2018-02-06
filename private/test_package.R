##
## test_package.R
## Script for testing all main functionalittis of the package
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##

library(hipathia)


# Preprocess
#---------------------------------

# load and filter graphs
pathways <- load.pathways(species = "hsa")
pathways_only2 <- load.pathways(species = "hsa",
                                pathways.list = c("hsa03320", "hsa04012"))

# prepare data
data("brca_data")
data("brca_design")
trans_data <- translate.matrix(brca_data, "hsa")
exp_data <- normalize.data(trans_data)


##############################
# Hipathia
##############################

results <- hipathia(exp_data, pathways, verbose = TRUE)


# Descriptive plots
#---------------------------------
path_vals <- get.paths.matrix(results)

# Define groups to compare
sample_group <- brca_design[colnames(path_vals),"group"]

# Plot heatmap for all paths
heatmap.plot(path_vals, sample.type = sample_group)
heatmap.plot(path_vals, sample.type = sample_group, colors = "hipathia",
             variable.clust = TRUE)

# Select sample colors
sample_colors <- c("#e6612c", "#b6ebe7")
names(sample_colors) <- unique(sample_group)
heatmap.plot(path_vals, sample_group, colors = "hipathia",
             sample.colors = sample_colors)


# Perform comparison
#---------------------------
comp <- do.wilcoxon(path_vals, sample_group, g1 = "Tumor", g2 = "Normal")
path_names <- get.path.names(pathways, rownames(comp))
comp <- cbind(path_names, comp)

pathways_summary <- get.pathways.summary(comp, pathways)
head(pathways_summary, n = 15)

ranked_path_vals <- path_vals[order(comp$p.value, decreasing = FALSE),]
ranked_comp <- comp[order(comp$p.value, decreasing = FALSE), ]

table(comp$FDRp.value<0.05)

# Plot heatmap for most relevant paths
heatmap.plot(ranked_path_vals[comp$FDRp.value < 0.05,],
             sample.type = sample_group, variable.clust = TRUE)

# Plot heatmap for most relevant paths
heatmap.plot(ranked_path_vals[1:15,], sample.type = sample_group,
             variable.clust = TRUE)
heatmap.plot(ranked_path_vals[1:15,], sample_group, colors = "hipathia",
             variable.clust = TRUE)


# PCA
#--------------
# Perform PCA model
pca_model <- do.pca(ranked_path_vals[1:ncol(ranked_path_vals),])

# plot PCA
pca.plot(pca_model, sample_group)
multiple.pca.plot(pca_model, sample_group, cex = 3, plot.variance = TRUE)

# Select sample colors
pca.plot(pca_model, sample_group, sample_colors)


##############################
# Pathway visualization
##############################

# Define node colors
colors_de <- node.color.per.de(results, pathways, sample_group, "Tumor",
                               "Normal")
colors_de_hipathia <- node.color.per.de(results, pathways, sample_group,
                                        "Tumor", "Normal", colors = "hipathia")
# Node colors with Uniprot grouping
colors_uni <- node.color.per.de(results, pathways, sample_group, "Tumor",
                               "Normal", group.by = "uniprot")
# Node colors with genes grouping
colors_gen <- node.color.per.de(results, pathways, sample_group, "Tumor",
                                "Normal", group.by = "genes")

# Visualize comparison as static image
#-------------------------------------------
pathway.comparison.plot(comp, metaginfo = pathways, pathway = "hsa03320",
                        node.colors = colors_de)
pathway.comparison.plot(comp, metaginfo = pathways, pathway = "hsa04014",
                        node.colors = colors_de_hipathia, colors = "hipathia")

# Save & Visualize comparison in server
#-------------------------------------------
create.report(comp, pathways, "save_noColors")
visualize.report("save_noColors")

create.report(comp, pathways, "save_colors", node.colors = colors_de)
visualize.report("save_colors", port = 4001)

# Visualize with grouping by Uniprot functions
create.report(comp, pathways, group.by = "uniprot",
              "save_colors_uniprot", node.colors = colors_uni)
visualize.report("save_colors_uniprot", port = 4002)

# Visualize with grouping by genes
create.report(comp, pathways, group.by = "genes",
              "save_colors_genes", node.colors = colors_gen)
visualize.report("save_colors_genes", port = 4003)

# Stop server
servr::daemon_stop()







# Decomposed pathways
#---------------------------

# Hipathia (for 2 pathways)
results_decomposed <- hipathia(exp_data, pathways_only2, verbose=TRUE,
                               decompose = TRUE)
dec_vals <- get.paths.matrix(results_decomposed)

# Comparison
comp_dec <- do.wilcoxon(dec_vals, sample_group, g1 = "Tumor", g2 = "Normal")
dec.names <- get.path.names(pathways_only2, rownames(comp_dec))
comp_dec <- cbind(dec.names, comp_dec)

# Save & visualize report
create.report(comp_dec, pathways_only2, "save_decomposed")
visualize.report("save_decomposed")
servr::daemon_stop()






##############################
# Function analysis
##############################

uniprot.annots <- get.pathways.annotations(rownames(comp), pathways, "uniprot")

# Uniprot Keywords
#------------------
uniprot_vals <- quantify.terms(results, pathways, "uniprot")

# Heatmap
heatmap.plot(uniprot_vals, sample.type = sample_group, variable.clust = TRUE)
# Comparison
comp_uni <- do.wilcoxon(uniprot_vals, sample_group, g1 = "Tumor", g2 = "Normal")


# GO terms
#-----------
go_vals <- quantify.terms(results, pathways, "GO")

# GO terms heatmap
heatmap.plot(go_vals, sample.type = sample_group, colors = "hipathia")
heatmap.plot(go_vals, sample.type = sample_group, colors = "hipathia",
             variable.clust = TRUE)

# Select most relevant GO terms
comp.go <- do.wilcoxon(go_vals, sample_group, g1 = "Tumor", g2 = "Normal")
ranked_go_vals <- go_vals[order(comp.go$p.value, decreasing = FALSE),]

# Number of significant GO terms
table(comp.go$FDRp.value < 0.05)

# 15 most significant GO terms heatmap
heatmap.plot(ranked_go_vals[1:15,], sample.type = sample_group,
             variable.clust = TRUE)
heatmap.plot(ranked_go_vals[1:15,], sample_group, colors="hipathia",
             variable.clust = TRUE)

# PCA
pca_model_go <- do.pca(ranked_go_vals[1:ncol(ranked_go_vals),])
# Plot PCA
pca.plot(pca_model_go, sample_group)
multiple.pca.plot(pca_model_go, sample_group)

# Create table of results with GO ancestors
go_table <- paths.to.go.ancestor(pathways, comp.paths = comp, comp.go = comp.go)
