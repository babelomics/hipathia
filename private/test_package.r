
library(hipathia)

data("brca_data")
data("brca_design")

# load and filter graphs
pathways <- load.pathways(species = "hsa")

# prepare data
trans.data <- translate.matrix(brca_data, "hsa")
exp.data <- normalize.data(trans.data)


# Hipathia
results <- hipathia(exp.data, pathways, verbose=TRUE)


# Descriptive statistics
#---------------------------------
path.vals <- get.paths.matrix(results)

sample.group <- brca_design[colnames(path.vals),"group"]

# Plot heatmap for all paths
heatmap.plot(path.vals, sample.type=sample.group)
heatmap.plot(path.vals, sample.type=sample.group, colors="hipathia", variable.clust = TRUE)

# select most relevant
comp <- do.wilcoxon(path.vals, sample.group, g1 = "Tumor", g2 = "Normal")
path.names <- get.path.names(pathways, rownames(comp))
comp <- cbind(path.names, comp)

pathways.summary <- get.pathways.summary(comp, pathways)
head(pathways.summary, n = 15)

ranked.path.vals <- path.vals[order(comp$p.value, decreasing = FALSE),]
ranked.comp <- comp[order(comp$p.value, decreasing = FALSE), ]

table(comp$FDRp.value<0.05)

# Plot heatmap for most relevant paths
heatmap.plot(ranked.path.vals[comp$FDRp.value<0.05,], sample.type=sample.group, variable.clust = TRUE)

# Plot heatmap for most relevant paths
heatmap.plot(ranked.path.vals[1:15,], sample.type=sample.group, variable.clust = TRUE)
heatmap.plot(ranked.path.vals[1:15,], sample.group, colors="hipathia", variable.clust = TRUE)


# Visualization
#------------------------
# Visualize comparison
colors.de <- node.color.per.differential.expression(results, pathways, sample.group, "Tumor", "Normal")
pathway.comparison.plot(comp, metaginfo = pathways, pathway = "hsa03320", node.colors = colors.de)
colors.de.hipathia <- node.color.per.differential.expression(results, pathways, sample.group, "Tumor", "Normal", colors = "hipathia")
pathway.comparison.plot(comp, metaginfo = pathways, pathway = "hsa04014", node.colors = colors.de.hipathia, colors = "hipathia")

# Visualize comparison in server
create.report(results, comp, pathways, "save_noColors")
create.report(results, comp, pathways, "save_colors", colors = colors.de)

visualize.report("save_colors")
visualize.report("save_noColors")


 # PCA
pca.model <- do.pca(ranked.path.vals[1:ncol(ranked.path.vals),])
# plot PCA
pca.plot(pca.model, sample.group)
multiple.pca.plot(pca.model, sample.group, cex=3, plot.variance = TRUE)

# # select sample colors
# sample.colors <- c("#e6612c", "#b6ebe7")
# names(sample.colors) <- unique(sample.group)
# pca.plot(pca.model, sample.group, sample.colors)
# heatmap.plot(path.vals, sample.group, colors="hipathia", sample.colors=sample.colors)



# Decomposed pathways
pathways.decomposed <- load.pathways(species = "hsa", pathways.list = c("hsa03320", "hsa04024"))
results.decomposed <- hipathia(exp.data, pathways.decomposed, verbose=TRUE, decompose = TRUE)
dec.vals <- get.paths.matrix(results.decomposed)

comp.dec <- do.wilcoxon(dec.vals, sample.group, g1 = "Tumor", g2 = "Normal")
dec.names <- get.path.names(pathways.decomposed, rownames(comp.dec))
comp.dec <- cbind(dec.names, comp.dec)

create.report(results.decomposed, comp.dec, pathways.decomposed, "save_decomposed")








# Function analysis
#---------------------------------
# Uniprot Keywords
uniprot.vals <- quantify.terms(results, pathways, "uniprot")
# plot
heatmap.plot(uniprot.vals, sample.type=sample.group, variable.clust = TRUE)
comp.uni <- do.wilcoxon(uniprot.vals, sample.group, g1 = "Tumor", g2 = "Normal")


# GO terms
go.vals <- quantify.terms(results, pathways, "GO")
# plot heatmap
heatmap.plot(go.vals, sample.type=sample.group, colors = "hipathia")

# Select most relevant GO terms
comp.go <- do.wilcoxon(go.vals, sample.group, g1 = "Tumor", g2 = "Normal")
ranked.go.vals <- go.vals[order(comp.go$p.value, decreasing = FALSE),]

table(comp.go$FDRp.value<0.05)

heatmap.plot(ranked.go.vals[1:15,], sample.type=sample.group, variable.clust = TRUE)
heatmap.plot(ranked.go.vals[1:15,], sample.group, colors="hipathia", variable.clust = TRUE)

# PCA
pca.model.go <- do.pca(ranked.go.vals[1:ncol(ranked.go.vals),])
# plot PCA
pca.plot(pca.model.go, sample.group)
multiple.pca.plot(pca.model.go, sample.group)



# Create table of results with GO ancestors 
go.table <- paths.to.go.ancestor(pathways, comp.paths = comp, comp.go = comp.go)
