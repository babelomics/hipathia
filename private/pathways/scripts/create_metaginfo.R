
library(KEGGgraph)
library(igraph)
library(graph)

hipath <- getwd()
source(paste0(hipath, "/R/main.R"))
source(paste0(hipath, "/R/functions.R"))
source(paste0(hipath, "/R/utils.R"))
source(paste0(hipath, "/private/pathways/scripts/graphs.R"))
source(paste0(hipath, "/private/pathways/scripts/KEGG_net.R"))
source(paste0(hipath, "/private/pathways/scripts/layout.R"))

# Parameters
species <- "hsa" # "rno" # "mmu" #

# set folders 
kgml.folder <- paste0(hipath, "/private/pathways/", species, "/kgml/")
sif.folder <- paste0(hipath, "/private/pathways/", species, "/sif/")
tmp.folder <- paste0(hipath, "/private/pathways/", species, "/temp/")
ammend.file <- paste0(hipath, "/private/pathways/sif_amendments.txt")
pathway.names <- unique(gsub(".xml", "", list.files(kgml.folder, pattern="xml")))

# Load annotations
dbannot <- load.annot.file(paste0(hipath, "/private/annotations/", species, "/uniprot_keywords_", species, "__biological_process.annot"))
entrez2hgnc <- utils::read.table(paste0(hipath, "/private/annotations/", species, "/entrez_hgnc_", species, ".annot"),header=F,sep="\t",stringsAsFactors=F)


# Process KGML files
#-------------------------------------------------
# create name_pathways.txt file
create.pathway.names.file(pathway.names, species, kgml.folder, sif.folder)

# Transform KGML to SIF files
transform.XML.to.SIF(pathway.names, kgml.folder, sif.folder)

# Load pathways from created SIF files
pgs <- load.graphs(sif.folder, species)
save(pgs, file=paste0(tmp.folder, "/pgs.RData"))

# Ammend pathways
apgs <- amend.kegg.pathways(ammend.file, pgs, species)
save(apgs, file=paste0(tmp.folder, "/apgs.RData"))
     
# Add final functions to the pathways
fpgs <- add.functions.to.pathigraphs(apgs, entrez2hgnc, dbannot, maxiter = 1000)
save(fpgs, file=paste0(tmp.folder, "/fpgs.RData"))

# Compute Path Normalization Values
metaginfo <- create.metaginfo.object(fpgs)
save(metaginfo, file=paste0(tmp.folder, "/meta_graph_info_", species, ".RData"))




