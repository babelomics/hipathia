##
## load.R
## Load objects from hpAnnot package through AnnotationHub
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##


get_hpannot_version <- function(){
    return("v2")
}


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
#' @import AnnotationHub
#'
load_annofuns <- function(db, species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_database(db))
        stop("Database not accepted")
    v <- get_hpannot_version()
    file <- paste0("annofuns_", db, "_", species, "_", v, ".rda")
    hp <- hub()
    annofuns <- suppressMessages(hp[[names(hp)[hp$title == file]]])
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
#' @import AnnotationHub
#'
load_mgi <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    v <- get_hpannot_version()
    file <- paste0("meta_graph_info_", species, "_", v, ".rda")
    hp <- hub()
    mgi <- suppressMessages(hp[[names(hp)[hp$title == file]]])
    return(mgi)
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
#' @import AnnotationHub
#'
load_pseudo_mgi <- function(species, group_by){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_grouping(group_by))
        stop("Grouping not accepted")
    v <- get_hpannot_version()
    file <- paste0("pmgi_", species, "_", group_by, "_", v, ".rda")
    hp <- hub()
    pmgi <- suppressMessages(hp[[names(hp)[hp$title == file]]])
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
#' @import AnnotationHub
#'
load_xref <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    v <- get_hpannot_version()
    file <- paste0("xref_", species, "_", v, ".rda")
    hp <- hub()
    xref <- suppressMessages(hp[[names(hp)[hp$title == file]]])
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
#' @import AnnotationHub
#'
load_entrez_hgnc <- function(species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    v <- get_hpannot_version()
    file <- paste0("entrez_hgnc_", species, "_", v, ".rda")
    hp <- hub()
    entrez_hgnc <- suppressMessages(hp[[names(hp)[hp$title == file]]])
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
#' @import AnnotationHub
#'
load_annots <- function(db, species){
    if(!is_accepted_species(species))
        stop("Species not accepted")
    if(!is_accepted_database(db))
        stop("Database not accepted")
    v <- get_hpannot_version()
    file <- paste0("annot_", db, "_", species, "_", v, ".rda")
    hp <- hub()
    annot <- suppressMessages(hp[[names(hp)[hp$title == file]]])
    return(annot)
}


#' Loads GO graph information
#'
#' #@examples
#' #load_gobp_frame()
#'
#' @return GO graph information
#' @import AnnotationHub
#'
load_gobp_frame <- function(){
    hp <- hub()
    v <- get_hpannot_version()
    file <- paste0("go_bp_frame_", v, ".rda")
    gbf <- suppressMessages(hp[[names(hp)[hp$title == file]]])
    return(gbf)
}


#' Loads GO graph
#'
#' #@examples
#' #load_gobp_net()
#'
#' @return GO graph
#' @import AnnotationHub
#'
load_gobp_net <- function(){
    hp <- hub()
    v <- get_hpannot_version()
    file <- paste0("go_bp_net_", v, ".rda")
    gbn <- suppressMessages(hp[[names(hp)[hp$title == file]]])
    return(gbn)
}


# Package-global cache of the Annotation-hub() object
hub = local({
    hp = NULL
    function(){
        if(is.null(hp))
            hp <- query(AnnotationHub(), "hpAnnot")
        hp
    }
})
