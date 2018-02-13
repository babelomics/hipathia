##
## load.R
## Load from hpAnnot package functions
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##


#' Loads annotations object
#'
#' @param db Database to be used. Either "GO" or "uniprot".
#' @param species Species of the samples.
#'
#' #@examples
#' #load_annofuns("GO", "hsa")
#' #load_annofuns("uniprot", "hsa")
#'
#' @return Annotations object
#'
load_annofuns <- function(db, species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_database(db))
        stop("Database not accepted")
    file <- paste0("annofuns_", db, "_", species, ".rda")
    af <- load(system.file("extdata", file, package = "hpAnnot"), 
               envir = environment())
    annofuns <- get(af)
    return(annofuns)
}


#' Loads object with graph information
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load_mgi("hsa")
#'
#' @return Graph information object
#'
load_mgi <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    file <- paste0("meta_graph_info_", species, ".rda")
    mgi <- load(system.file("extdata", file, package = "hpAnnot"), 
                envir = environment())
    meta_graph_info <- get(mgi)
    return(meta_graph_info)
}


#' Loads object with pseudo graph information
#'
#' @param species Species of the samples.
#' @param group_by How to group the subpathways to be visualized. By default
#' they are grouped by the pathway to which they belong. Available groupings
#' include "uniprot", to group subpathways by their annotated Uniprot functions,
#' "GO", to group subpathways by their annotated GO terms, and "genes", to group
#' subpathways by the genes they include.
#'
#' #@examples
#' #load_pseudo_mgi("hsa", "uniprot")
#'
#' @return Pseudo graph information object
#'
load_pseudo_mgi <- function(species, group_by){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_grouping(group_by))
        stop("Grouping not accepted")
    file <- paste0("pmgi_", species, "_", group_by, ".rda")
    pseudo <- load(system.file("extdata", file, package = "hpAnnot"), 
         envir = environment())
    pmgi <- get(pseudo)
    return(pmgi)
}


#' Loads table of references
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load_xref("hsa")
#'
#' @return Table of references
#'
load_xref <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    file <- paste0("xref_", species, ".rda")
    xr <- load(system.file("extdata", file, package = "hpAnnot"), 
         envir = environment())
    xref <- get(xr)
    return(xref)
}


#' Loads table of translation from HGNC to Entrez
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load_entrez_hgnc("hsa")
#'
#' @return Table of translation from HGNC to Entrez
#'
load_entrez_hgnc <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    file <- paste0("entrez_hgnc_", species, ".rda")
    eh <- load(system.file("extdata", file, package = "hpAnnot"), 
         envir = environment())
    entrez_hgnc <- get(eh)
    return(entrez_hgnc)
}


#' Loads functional annotations to genes
#'
#' Loads functional annotations from HGNC to the selected database.
#'
#' @param db Database to be used. Either "GO" or "uniprot".
#' @param species Species of the samples.
#'
#' #@examples
#' #load_annots("GO", "hsa")
#'
#' @return Functional annotations from HGNC to the selected database.
#'
load_annots <- function(db, species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_database(db))
        stop("Database not accepted")
    file <- paste0("annot_", db, "_", species, ".rda")
    ann <- load(system.file("extdata", file, package = "hpAnnot"), 
                  envir = environment())
    annot <- get(ann)
    return(annot)
}


#' Loads GO graph information
#'
#' #@examples
#' #load_gobp_frame()
#'
#' @return GO graph information
#'
load_gobp_frame <- function(){
    gbf <- load(system.file("extdata", "go_bp_frame.rda", package = "hpAnnot"), 
                envir = environment())
    go_bp_frame <- get(gbf)
    return(go_bp_frame)
}


#' Loads GO graph
#'
#' #@examples
#' #load_gobp_net()
#'
#' @return GO graph
#'
load_gobp_net <- function(){
    gbn <- load(system.file("extdata", "go_bp_net.rda", package = "hpAnnot"), 
                envir = environment())
    go_bp_net <- get(gbn)
    return(go_bp_net)
}

