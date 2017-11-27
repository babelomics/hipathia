
load("../private/data/brca_design.RData")
load("../private/data/brca_data__all_genes.RData")
load("../private/extdata/meta_graph_info_hsa.RData")

brca_data <- brca_data__all_genes
comm <- intersect(rownames(brca_data), metaginfo$all.genes)
brca_data <- brca_data[comm,]

brca <- NULL
brca$data <- brca_data
brca$design <- brca_design

# save(brca, file="../private/data/brca.RData")
save(brca_data, file="../private/data/brca_data.RData")
save(brca_design, file="../private/data/brca_design.RData")

