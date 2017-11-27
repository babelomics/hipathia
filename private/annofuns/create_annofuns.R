
library(igraph)
source("../R/functions.R")

all.species <- c("hsa","mmu","rno")

for(species in all.species){
  
  load(paste0("extdata/meta_graph_info_", species, ".RData"))
  entrez2hgnc <- read.table(paste0("extdata/entrez_hgnc_", species, ".annot"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  
  # Uniprot
  
  uniannot <- load.annot.file(paste0("extdata/uniprot_keywords_", species, "__biological_process.annot"))
  annofuns <- annotate.paths(metaginfo$pathigraphs, uniannot, entrez2hgnc)
  save(annofuns, file=paste0("annofuns/annofuns_uniprot_", species, ".RData"))
  save(annofuns, file=paste0("extdata/annofuns_uniprot_", species, ".RData"))
  
  # GO
  
  goannot <- load.annot.file(paste0("extdata/go_bp_", species, ".annot"))
  annofuns <- annotate.paths(metaginfo$pathigraphs, goannot, entrez2hgnc)
  save(annofuns, file=paste0("annofuns/annofuns_GO_", species, ".RData"))
  save(annofuns, file=paste0("extdata/annofuns_GO_", species, ".RData"))
  
}

