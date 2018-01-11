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
#' #load.annofuns("GO", "hsa")
#' #load.annofuns("uniprot", "hsa")
#'
#' @return Annotations object
#'
load.annofuns <- function(db, species){
    if(!is.accepted(species))
        stop("Species not accepted")
    if(db == "GO"){
        ag <- utils::data("annofuns_GO", envir = environment())
        annofuns_GO <- get(ag)
        annofuns <- annofuns_GO[[species]]
    }else if(db == "uniprot"){
        au <- utils::data("annofuns_uniprot", envir = environment())
        annofuns_uniprot <- get(au)
        annofuns <- annofuns_uniprot[[species]]
    }else{
        stop("Database not accepted")
    }
    return(annofuns)
}


#' Loads object with graph information
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load.mgi("hsa")
#'
#' @return Graph information object
#'
load.mgi <- function(species){
    if(!is.accepted(species))
        stop("Species not accepted")
    if(species == "hsa"){
        mgi <- utils::data("meta_graph_info_hsa", envir = environment())
        meta_graph_info_hsa <- get(mgi)
        return(meta_graph_info_hsa)
    }else if(species == "rno"){
        mgi <- utils::data("meta_graph_info_rno", envir = environment())
        meta_graph_info_rno <- get(mgi)
        return(meta_graph_info_rno)
    }else if(species == "mmu"){
        mgi <- utils::data("meta_graph_info_mmu", envir = environment())
        meta_graph_info_mmu <- get(mgi)
        return(meta_graph_info_mmu)
    }else{
        stop("Species not accepted")
    }
}


#' Loads object with pseudo graph information
#'
#' @param species Species of the samples.
#' @param group.by How to group the subpathways to be visualized. By default
#' they are grouped by the pathway to which they belong. Available groupings
#' include "uniprot", to group subpathways by their annotated Uniprot functions,
#' "GO", to group subpathways by their annotated GO terms, and "genes", to group
#' subpathways by the genes they include.
#'
#' #@examples
#' #load.pseudo.mgi("hsa", "uniprot")
#'
#' @return Pseudo graph information object
#'
load.pseudo.mgi <- function(species, group.by){
    if(!is.accepted(species))
        stop("Species not accepted")
    if(!is.accepted.grouping(group.by))
        stop("Grouping not accepted")
    if(species == "hsa"){
        if(group.by == "uniprot"){
            pmgi <- utils::data("pmgi_hsa_uniprot", envir = environment())
            pmgi_hsa_uniprot <- get(pmgi)
            return(pmgi_hsa_uniprot)
        }else if(group.by == "GO"){
            pmgi <- utils::data("pmgi_hsa_GO", envir = environment())
            pmgi_hsa_GO <- get(pmgi)
            return(pmgi_hsa_GO)
        }else if(group.by == "genes"){
            pmgi <- utils::data("pmgi_hsa_genes", envir = environment())
            pmgi_hsa_genes <- get(pmgi)
            return(pmgi_hsa_genes)
        }
    }else if(species == "rno"){
        if(group.by == "uniprot"){
            pmgi <- utils::data("pmgi_rno_uniprot", envir = environment())
            pmgi_rno_uniprot <- get(pmgi)
            return(pmgi_rno_uniprot)
        }else if(group.by == "GO"){
            pmgi <- utils::data("pmgi_rno_GO", envir = environment())
            pmgi_rno_GO <- get(pmgi)
            return(pmgi_rno_GO)
        }else if(group.by == "genes"){
            pmgi <- utils::data("pmgi_rno_genes", envir = environment())
            pmgi_rno_genes <- get(pmgi)
            return(pmgi_rno_genes)
        }
    }else if(species == "mmu"){
        if(group.by == "uniprot"){
            pmgi <- utils::data("pmgi_mmu_uniprot", envir = environment())
            pmgi_mmu_uniprot <- get(pmgi)
            return(pmgi_mmu_uniprot)
        }else if(group.by == "GO"){
            pmgi <- utils::data("pmgi_mmu_GO", envir = environment())
            pmgi_mmu_GO <- get(pmgi)
            return(pmgi_mmu_GO)
        }else if(group.by == "genes"){
            pmgi <- utils::data("pmgi_mmu_genes", envir = environment())
            pmgi_mmu_genes <- get(pmgi)
            return(pmgi_mmu_genes)
        }
    }else{
        stop("Species not accepted")
    }
}


#' Loads table of references
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load.xref("hsa")
#'
#' @return Table of references
#'
load.xref <- function(species){
    if(!is.accepted(species))
        stop("Species not accepted")
    xr <- utils::data("xref", envir = environment())
    xref_spe <- get(xr)[[species]]
    return(xref_spe)
}


#' Loads table of translation from HGNC to Entrez
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #load.entrez.hgnc("hsa")
#'
#' @return Table of translation from HGNC to Entrez
#'
load.entrez.hgnc <- function(species){
    if(!is.accepted(species))
        stop("Species not accepted")
    eh <- utils::data("entrez_hgnc", envir = environment())
    entrez_hgnc <- get(eh)[[species]]
    return(entrez_hgnc)
}


#' Loads GO Annotations
#'
#' @param species Species of the samples.
#'
#' @return GO Annotations
#'
load.gobp <- function(species){
    if(!is.accepted(species))
        stop("Species not accepted")
    gba <- utils::data("go_bp_annots", envir = environment())
    go_bp_annot <- get(gba)[[species]]
    return(go_bp_annot)
}


#' Loads Uniprot annotations
#'
#' @param species Species of the samples.
#'
#' @return Uniprot annotations
#'
load.unibp <- function(species){
    if(!is.accepted(species))
        stop("Species not accepted")
    uba <- utils::data("uni_bp_annots", envir = environment())
    uni_bp_annot <- get(uba)[[species]]
    return(uni_bp_annot)
}


#' Loads functional annotations
#'
#' Loads functional annotations from HGNC to the selected database.
#'
#' @param db Database to be used. Either "GO" or "uniprot".
#' @param species Species of the samples.
#'
#' #@examples
#' #load.annots("GO", "hsa")
#'
#' @return Functional annotations from HGNC to the selected database.
#'
load.annots <- function(db, species){
    if(!is.accepted(species))
        stop("Species not accepted")
    if(db == "GO"){
        annofuns <- load.gobp(species)
    }else if(db == "uniprot"){
        annofuns <- load.unibp(species)
    }else{
        stop("Database not accepted")
    }
    return(annofuns)
}


#' Loads GO graph information
#'
#' #@examples
#' #load.gobp.frame()
#'
#' @return GO graph information
#'
load.gobp.frame <- function(){
    gbf <- utils::data("go_bp_frame", envir = environment())
    go_bp_frame <- get(gbf)
    return(go_bp_frame)
}


#' Loads GO graph
#'
#' #@examples
#' #load.gobp.net()
#'
#' @return GO graph
#'
load.gobp.net <- function(){
    gbn <- utils::data("go_bp_net", envir = environment())
    go_bp_net <- get(gbn)
    return(go_bp_net)
}

