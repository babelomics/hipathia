##
## utils.R
## Utility functions for package Hipathia
##
## Written by Marta R. Hidalgo, Jose Carbonell-Caballero
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##



translate_ids <- function(ids, xref){

    # get translation
    ltids <- xref[ids]
    tids <- sapply(ltids, function(x){
        if(is.null(x)){
            return(NA)
        } else {
            return(x[1])
        }
    })
    names(tids) <- ids

    # check problems
    is_na <- is.na(tids)
    tids_count <- table(tids)
    is_dup <- tids_count[as.character(tids)] > 1 & is_na == FALSE

    # separate ids
    translated_ids <- unique(ids[!is_na])
    untranslated_ids <- unique(ids[is_na])
    duplicated_ids <- unique(ids[is_dup])

    # compute counts and ratios
    translated_ids_count <- length(translated_ids)
    untranslated_ids_count <- length(untranslated_ids)
    duplicated_ids_count <- length(duplicated_ids)
    n <- length(ids)
    translated_ids_ratio <- translated_ids_count/n
    untranslated_ids_ratio <- untranslated_ids_count/n
    duplicated_ids_ratio <- duplicated_ids_count/n

    return(list(
        translation = tids,
        is_na = is_na,
        is_dup = is_dup,
        translated_ids = translated_ids,
        untranslated_ids = untranslated_ids,
        duplicated_ids = duplicated_ids,
        translated_ids_count = translated_ids_count,
        untranslated_ids_count = untranslated_ids_count,
        duplicated_ids_count = duplicated_ids_count,
        translated_ids_ratio = translated_ids_ratio,
        untranslated_ids_ratio = untranslated_ids_ratio,
        duplicated_ids_ratio = duplicated_ids_ratio
    )
    )

}




#' @title
#' Translation of the rownames IDs of a SummarizedExperiment to Entrez IDs.
#'
#' @description
#' Translates the IDs in the rownames of a SummarizedExperiment to Entrez IDs.
#' For accepted IDs to be transformed see the DOCUMENTATION.
#'
#' @param data Either a SummarizedExperiment object or a matrix of gene 
#' expression.
#' @param species Species of the samples.
#' @param sel_assay Character or integer, indicating the assay to be translated 
#' in the SummarizedExperiment. Default is 1.
#' @param verbose Boolean, whether to show details about the results of the
#' execution.
#'
#' @examples data("brca_data")
#' trans_data <- translate_data(brca_data, "hsa")
#'
#' @return Either a SummarizedExperiment or a matrix (depending on the input 
#' type) of gene expression with Entrez IDs as rownames.
#'
#' @export
#' @import SummarizedExperiment
#' @import AnnotationHub
#' @importFrom methods is
#'
translate_data <- function(data, species, sel_assay = 1, verbose=TRUE){
    
    if(is(data, "SummarizedExperiment")){
        se_flag <- TRUE
        mat <- assay(data, sel_assay)
    }else if(is(data, "matrix")){
        se_flag <- FALSE
        mat <- data
    }else{
        stop("Only SummarizedExperiment or matrix classes accepted as data")
    }
    trans_mat <- translate_matrix(mat, species)
    if(se_flag == TRUE)
        trans_mat <- SummarizedExperiment(list(trans = trans_mat), 
                                     colData = colData(data))
    return(trans_mat)
}


#' @title
#' Translation of the rownames IDs of a matrix to Entrez IDs.
#'
#' @description
#' Translates the IDs in the rownames of a matrix to Entrez IDs.
#' For accepted IDs to be transformed see the DOCUMENTATION.
#'
#' @param exp Matrix of gene expression.
#' @param species Species of the samples.
#' @param verbose Boolean, whether to show details about the results of the
#' execution.
#'
#' @return Matrix of gene expression with Entrez IDs as rownames.
#'
#' @import AnnotationHub
#'
translate_matrix <- function(exp, species, verbose = TRUE){

    xref <- load_xref(species)

    new_ids <- gsub("\\.[0123456789]$", "", rownames(exp))

    # translate ids
    tt <- translate_ids(new_ids, xref)

    # filter untranslated ids
    exp2 <- exp[!tt$is_na,,drop=FALSE]
    valid_translation <- tt$translation[!tt$is_na]

    # average duplicated ids
    raw_exp3 <- by(exp2, valid_translation, colMeans, na.rm = TRUE)
    if(ncol(exp2) > 1){
        exp3 <- do.call("rbind", raw_exp3)
    } else {
        exp3 <- matrix(raw_exp3, ncol = 1)
        rownames(exp3) <- names(raw_exp3)
        colnames(exp3) <- colnames(exp2)
    }

    if(verbose==TRUE){
        cat("translated ids = ", tt$translated_ids_count, " (",
            format(digits = 2, tt$translated_ids_ratio), ") \n", sep = "")
        cat("untranslated ids = ", tt$untranslated_ids_count, " (",
            format(digits = 2, tt$untranslated_ids_ratio), ") \n", sep = "")
        cat("multihit ids = ", sum(tt$duplicated_ids_count), " (",
            format(digits = 2, tt$duplicated_ids_ratio), ") \n", sep = "")
    }

    attr(exp3,"translation") <- tt
    return(exp3)

}



#' Tranlates path IDs to path names
#'
#' @description
#' Translates the subpathway IDs to readable and comprensible names.
#'
#' For effector subpathways, the names of the subpathways are encoded
#' as "pathway: effector_protein", where "pathway" is the pathway to
#' which the subpathway belongs and "effector_protein" is the name of
#' the last node in the subpathway.
#'
#' For decomposed subpathways, the names of the subpathways are encoded
#' as "pathway: receptor_protein - effector_protein", where "pathway" is
#' the pathway to which the subpathway belongs, "receptor_protein" is the
#' name of the initial node of the subpathway and "effector_protein" is
#' the name of the last node in the subpathway.
#'
#' @param metaginfo Pathways object
#' @param names Character vector with the subpathway IDs to be translated
#' @param maxchar Integer, describes the number of maximum characters to
#' be shown. By default no filter is applied.
#'
#' @examples
#' data(path_vals)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' translated_names <- get_path_names(pathways, rownames(path_vals))
#'
#' @return A character vector including the readable names of the
#' subpathways IDs, in the same order as provided.
#' @export
#'
get_path_names <- function(metaginfo, names, maxchar=NULL){

    labels <- metaginfo$all.labelids

    prettynames <- unlist(lapply(names, function(name){
        strname <- unlist(strsplit(name, "\\-"))

        if(length(strname) == 4){
            ininode <- paste("N", strname[2], strname[3], sep = "-")
            effnode <- paste("N", strname[2], strname[4], sep = "-")
            inilabel <- as.matrix(labels)[ininode, "label"]
            efflabel <- as.matrix(labels)[effnode, "label"]
            pathlabel <- as.matrix(labels)[ininode, "path.name"]
            label <- paste0(pathlabel, ": ", inilabel, " -> ", efflabel )
            label
        }else if(length(strname) == 3){
            effnode <- paste("N", strname[2], strname[3], sep = "-")
            label <- as.matrix(labels)[effnode, c("path.name", "label")]
            label <- paste(label, collapse = ": ")
            label
        }else{
            stop("Not recognized name")
        }
    }))

    if(!is.null(maxchar))
        prettynames <- clip_names(prettynames, maxchar = maxchar)

    return(prettynames)

}



#' @title
#' Upgrade igraphs to current version
#'
#' @description
#' Upgrades the \code{igraph} objects in metaginfo object to the corresponding
#' version of the \code{igraph} package.
#'
#' @import igraph
#'
#' @param metaginfo Pathways object
#'
#' @return The pathways object with the upgraded igraph objects
#'
igraphs_upgrade <- function(metaginfo){

    for(pw in names(metaginfo$pathigraphs)){

        ug <- upgrade_graph(metaginfo$pathigraphs[[pw]]$graph)
        metaginfo$pathigraphs[[pw]]$graph <- ug

        for(sp in names(metaginfo$pathigraphs[[pw]]$effector.subgraphs)){
            ug <- upgrade_graph(metaginfo$pathigraphs[[pw]]$
                                    effector.subgraphs[[sp]])
            metaginfo$pathigraphs[[pw]]$effector.subgraphs[[sp]] <- ug
        }

        for(sp in names(metaginfo$pathigraphs[[pw]]$subgraphs)){
            ug <- upgrade_graph(metaginfo$pathigraphs[[pw]]$subgraphs[[sp]])
            metaginfo$pathigraphs[[pw]]$subgraphs[[sp]] <- ug
        }
    }
    return(metaginfo)
}


#' @title
#' Loads the pathways object.
#'
#' @description
#' Loads the pathways object, which includes information about the pathways
#' to be analyzed.
#'
#' @import igraph
#'
#' @details
#' The object of pathways includes information about the pathways and the
#' subpathways which will be analyzed. This object must be provided to some
#' of the functions (like \code{hipathia} or \code{quantify_terms}) in the
#' package. These functions will analyze all the pathways included in this
#' object. By default, all available pathways are load. In order to restrict
#' the analysis to a predefined set of pathways, specify the set of pathways
#' to load with the parameter \code{pathways_list}.
#'
#' @param species Species of the samples.
#' @param pathways_list Vector of the IDs of the pathways to load. By default
#' all available pathways are load.
#'
#' @examples
#' pathways <- load_pathways("hsa")   # Loads all pathways for human
#' pathways <- load_pathways("mmu", c("mmu03320", "mmu04024", "mmu05200"))
#'    # Loads pathways 03320, 04024 and 05200 for mouse
#'
#' @return An pathways object including
#' * \code{species} Species to which the pathways are related.
#' * \code{pathigraphs} List of Pathigraph objects. Each Pathigraph contains
#' the necessary information of a pathway for it to be analyzed
#' with \code{Hipathia}.
#' * \code{all_genes} List of all the genes included in the selection of
#' pathways stored in \code{pathigraphs}.
#' * \code{eff_norm} Vector of normalization values for effector subpathways.
#' * \code{path_norm} Vector of normalization values for decomposed
#' subpathways.
#'
#' @export
#' @import AnnotationHub
#'
load_pathways <- function(species, pathways_list = NULL){
    metaginfo <- load_mgi(species)
    metaginfo <- filter_pathways(metaginfo, pathways_list = pathways_list)
    metaginfo <- igraphs_upgrade(metaginfo)
    message("Loaded ", length(metaginfo$pathigraphs), " pathways")
    return(metaginfo)
}




#' Lists the IDs of the pathways in a pathways object
#'
#' Lists the IDs of the pathways included in the pathways object
#' \code{metaginfo}
#'
#' @param metaginfo Pathways object
#' @return List of the pathway IDs included in the pathways object
#'
#' @examples
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' pathways_list <- get_pathways_list(pathways)
#'
#' @export
#'
get_pathways_list <- function(metaginfo){
    return(names(metaginfo$pathigraphs))
}


filter_pathways <- function(metaginfo, pathways_list = NULL){
    if(!is.null(pathways_list)){
        metaginfo$pathigraphs <- metaginfo$pathigraphs[pathways_list]
        metaginfo$all.genes <- all_needed_genes(metaginfo$pathigraphs)
        pn_paths <- sapply(strsplit(names(metaginfo$path.norm), "-"), "[[", 2)
        pn_idx <- pn_paths %in% pathways_list
        metaginfo$path.norm <- metaginfo$path.norm[pn_idx]
        en_paths <- sapply(strsplit(names(metaginfo$eff.norm), "-"), "[[", 2)
        en_idx <- en_paths %in% pathways_list
        metaginfo$eff.norm <- metaginfo$eff.norm[en_idx]
        metaginfo$all.labelids <- metaginfo$all.labelids[
            metaginfo$all.labelids[,"path.id"] %in% pathways_list,]
    }
    return(metaginfo)
}


all_needed_genes <- function(pathigraphs){
    genes <- unique(unlist(sapply(pathigraphs, function(x){
        unique(unlist(V(x$graph)$genesList))
    })))
    return(genes[!is.na(genes) & genes!="/"])
}



#' Gets the object of subpathway activation values
#'
#' @description
#' This function returns the object with the levels of activation of each
#' subpathway for each sample. Rows represent the subpathways and columns
#' represent the samples. Each cell is the value of activation of a subpathway
#' in a sample.
#'
#' Rownames are the IDs of the subpathways. In order to transform IDs into
#' readable names, use \code{get_path_names}.
#'
#' Effector subpathways are subgraphs of a pathway including all the paths
#' leading to an effector protein. Effector proteins are defined as final
#' nodes in the graph. Each effector protein (final node) in a pathway
#' defines its own effector subpathway as the nodes and edges in a path leading
#' to it.
#'
#' Decomposed subpathways are subgraphs of a pathway including all the paths
#' starting in a receptor protein and ending in an effector protein. Receptor
#' proteins are defined as initial nodes and effector proteins are defined
#' as final nodes in the graph. Each effector subpathway can be decomposed
#' in as many decomposed subpathways as initial nodes it includes.
#'
#' @param results Results object as returned by \code{hipathia}.
#' @param matrix Boolean, if TRUE the function returns a matrix object, if 
#' FALSE (as default) returns a SummarizedExperiment object.
#'
#' @examples
#' data(results)
#' path_vals <- get_paths_data(results)
#'
#' @return Object, either a SummarizedExperiment or a matrix, with the levels 
#' of activation of each decomposed subpathway for each sample.
#' 
#' @export
#' @import SummarizedExperiment
#'
get_paths_data <- function(results, matrix = FALSE){
    if(matrix == TRUE){
        return(assay(results[["paths"]]))
    }else{
        return(results[["paths"]])
    }
}




clip_names <- function(snames, maxchar = 30){
    n <- nchar(snames)
    idx <- n > maxchar
    n[idx] <- maxchar - 3
    paste0(substr(snames, 1, n), ifelse(idx, "...", ""))
}



#' @importFrom stats median
add_missing_genes <- function(exp_data, genes, default = NULL){
    if(is.null(default))
        default <- stats::median(exp_data)
    missing_genes <- setdiff(genes, rownames(exp_data))
    if(length(missing_genes > 0)){
        fakemat <- default + matrix(0, nrow = length(missing_genes),
                                    ncol = ncol(exp_data))
        rownames(fakemat) <- missing_genes
        colnames(fakemat) <- colnames(exp_data)
        exp_data <- rbind(exp_data, fakemat)
        # message("----------------------------------------------------")
        message("Added missing genes: ",
                length(missing_genes),
                " (",
                round(length(missing_genes)/nrow(exp_data) * 100, digits = 2),
                "%)")
        # message("----------------------------------------------------")
    }
    return(exp_data)
}




#' Head function for SummarizedExperiment, data.frames and matrix objects
#'
#' Shows the first \code{n} rows and the first \code{n} columns of a matrix,
#' in case the matrix has more than \code{n+5} rows or columns.
#' Otherwise, it shows all the rows or columns, respectively.
#'
#' @param mat Object to be shown
#' @param n Number of rows and columns
#' @param sel_assay Character or integer, indicating the assay to be translated 
#' in the SummarizedExperiment. Default is 1.
#'
#' @examples mat <- matrix(rnorm(100), ncol = 10)
#' hhead(mat)
#' hhead(mat, 3)
#' hhead(mat, 7)
#'
#' @return Matrix with as much as \code{n} rows and \code{n} columns.
#'
#' @importFrom utils head
#' @importFrom methods is
#' @export
hhead <- function(mat, n = 5, sel_assay = 1){
    if(is(mat, "SummarizedExperiment"))
        mat <- assay(mat, sel_assay)
    if(!is.null(ncol(mat))){
        if(ncol(mat) >= n+5){
            if(nrow(mat) >= n+5){
                mat[seq_len(n),seq_len(n)]
            }else{
                mat[,seq_len(n)]
            }
        }else{
            if(nrow(mat) >= n+5){
                mat[seq_len(n),]
            }else{
                mat
            }
        }
    }else{
        utils::head(mat, n)
    }
}



get_effnode_id <- function(path_name){
    path_split <- unlist(strsplit(path_name, split="-"))
    if(!path_split[1] == "P") return(NA)
    if(length(path_split) == 4){
        paste("N", path_split[2], path_split[4], sep="-")
    }else if(length(path_split) == 3){
        paste("N", path_split[2], path_split[3], sep="-")
    }else{
        return(NA)
    }
}


get_ininode_id <- function(path_name){
    path_split <- unlist(strsplit(path_name, split="-"))
    if(length(path_split) > 3){
        paste("N", path_split[2], path_split[3], sep="-")
    }else{
        return(NA)
    }
}


get_effpath_id <- function(node_name){
    node_split <- unlist(strsplit(node_name, split="-"))
    if(!node_split[1] == "N") return(NA)
    gsub("N", "P", node_name)
}


#' Get Pathways functional annotations
#'
#' Get functional annotation of the pathways, either for a particular
#' annotation or a stored one.
#'
#' @param pathway_names Character vector of the names of the pathways
#' @param metaginfo Pathways object
#' @param dbannot Either a string indicating which precomputed annotation
#' to use ("uniprot" for Uniprot Keywords or "GO" for Gene Ontology terms),
#' or a dataframe with the annotation of the genes to the functions. First
#' column are gene symbols, second column the functions.
#' @param collapse Boolean, whether to collapse all functions of the same
#' path in a single character string.
#'
#' @return 2-columns matrix with the annotations of each pathway ID in the
#' annotation \code{dbannot}.
#'
#' @examples
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' pathway_names <- c("P-hsa03320-37", "P-hsa03320-61", "P-hsa03320-46",
#' "P-hsa03320-57", "P-hsa03320-64", "P-hsa03320-47", "P-hsa03320-65")
#' get_pathways_annotations(pathway_names, pathways, "GO")
#' get_pathways_annotations(pathway_names, pathways, "uniprot")
#'
#' @export
#' @import AnnotationHub
#'
get_pathways_annotations <- function(pathway_names, metaginfo, dbannot,
                                     collapse = FALSE){

    if(is.character(dbannot)){
        annofuns <- load_annofuns(dbannot, metaginfo$species)
    }else{
        annofuns <- annotate_paths(metaginfo, dbannot)
    }

    annofuns$funs[is.na(annofuns$funs)] <- ""
    decomposed <- is_decomposed(pathway_names)
    if(decomposed == TRUE)
        pathway_names <- sapply(pathway_names, function(n)
            paste(unlist(strsplit(n, split="-"))[c(1,2,4)], collapse = "-"))

    if(collapse == TRUE){
        miniaf <- do.call("rbind", lapply(pathway_names, function(path){
            path_id <- annofuns$paths == path
            data.frame(effector.nodes = annofuns$effector.nodes[path_id][1],
                       paths = path,
                       funs = paste(annofuns$funs[path_id], collapse = ", "),
                       stringsAsFactors = FALSE)
        }))
    }else{
        miniaf <- do.call("rbind", lapply(pathway_names, function(path){
            annofuns[annofuns$paths == path,]
        }))
    }
    rownames(miniaf) <- NULL

    return(miniaf[,c(2,3)])
}



#' Get highest common GO ancestor of GO annotations
#'
#' @param go_terms GO terms for which the highest common ancestors are
#' to be looked for.
#' @param go_comp Wilcoxon comparison of the matrix of GO values as returned
#' by \code{do_wilcoxon}.
#' @param metaginfo Pathways object
#' @param unique Boolean, whether to return only one highest significant GO
#' ancestor or all of them. By default, TRUE.
#' @param pval P-value cut-off. Default values is set to 0.05.
#'
#' @return highest common ancestors
#'
#' #@export
#' @import AnnotationHub
#'
get_highest_sig_ancestor <- function(go_terms, go_comp, metaginfo,
                                     unique = TRUE, pval = 0.05){

    go_bp_frame <- load_gobp_frame()
    go_bp_net <- load_gobp_net()

    # Relacionar GO term con su etiqueta
    go_labels <- sapply(go_terms, function(term) {
        go_bp_frame[go_bp_frame$name == term, "id"]
    })
    go_comp$labels <- sapply(rownames(go_comp), function(term){
        go_bp_frame[go_bp_frame$name == term, "id"]
    })
    go_comp$terms <- rownames(go_comp)
    rownames(go_comp) <- go_comp$labels
    sig_go_labels <- go_comp[go_comp$FDRp.value < pval, "labels"]

    # Encontrar GO superior
    sup <- do.call("rbind", lapply(go_labels, function(label){
        short <- shortest.paths(go_bp_net, label, mode = "in")[1,]
        ancestors <- names(short)[!short == "Inf"]
        sig_ancs <- intersect(ancestors, sig_go_labels)
        if(length(sig_ancs) > 0 ){
            if(length(sig_ancs) > 1)
                sig_ancs <- setdiff(sig_ancs, label)
            sig_levels <- go_bp_frame[sig_ancs, "level"]
            min_sl_idx <- sig_levels == min(sig_levels)
            highest_ancestors <- go_bp_frame[sig_ancs, "id"][min_sl_idx]
            if(unique == TRUE)
                highest_ancestors <- highest_ancestors[1]
            df <- data.frame(GO_term = label,
                             GO_name = go_bp_frame[label, "name"],
                             GO_adj_pval = go_comp[label, "FDRp.value"],
                             Highest_Significant_Ancestor = highest_ancestors,
                             HSA_name = go_bp_frame[highest_ancestors, "name"],
                             HSA_adj_pval = go_comp[highest_ancestors,
                                                    "FDRp.value"],
                             stringsAsFactors = FALSE)
        }else{
            df <- data.frame(GO_term = label,
                             GO_name = go_bp_frame[label, "name"],
                             GO_adj_pval = go_comp[label, "FDRp.value"],
                             Highest_Significant_Ancestor = "",
                             HSA_name = "",
                             HSA_adj_pval = "",
                             stringsAsFactors = FALSE)
        }
    }))

    return(sup)
}


#' Create path results table with highest significant GO ancestors
#'
#' Create table of results with the comparison of the paths together with
#' the GO functional annotation and the highest significant GO ancestor
#' (HSGOA).
#'
#' The table returns in each row: the name of a pathway and its Wilcoxon
#' comparison information (direction, adjusted p-value), the GO term to which
#' the path is related (not necessarily unique), the Wilcoxon comparison
#' informationfor this GO (direction, adjusted p-value), the HSGOA of this
#' GO and its Wilcoxon comparison information (direction, adjusted p-value).
#'
#' The HSGOA is computed as the GO term with minimum level from all the
#' significant (with respect to value \code{pval}) ancestors of a GO.
#' The level of a GO term is computed as the number of nodes in the shortest
#' path from this GO term to the term "GO:0008150". The ancestors of a node
#' are defined as all the nodes from which a path can be defined from the
#' ancestor to the node.
#'
#' @param pathways Pathways object
#' @param comp_paths Wilcoxon comparison of the matrix of pathways values
#' as returned by \code{do_wilcoxon}.
#' @param comp_go Wilcoxon comparison of the matrix of GO values as
#' returned by \code{do_wilcoxon}.
#' @param pval P-value cut-off. Default values is set to 0.05.
#'
#' @return Table of comparisons with Highest common ancestors
#'
#' @examples
#' data(comp)
#' data(go_vals)
#' data(brca_design)
#' data(path_vals)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' comp_go <- do_wilcoxon(go_vals, sample_group, g1 = "Tumor", g2 = "Normal")
#' \dontrun{pathways <- load_pathways(species = "hsa", pathways_list =
#' c("hsa03320", "hsa04012"))
#' table <- paths_to_go_ancestor(pathways, comp, comp_go)}
#'
#' @export
#'
paths_to_go_ancestor <- function(pathways, comp_paths, comp_go, pval = 0.05){
    path_names <- get_path_names(pathways, rownames(comp_paths))
    names(path_names) <- rownames(comp_paths)
    path_annot <- get_pathways_annotations(pathway_names = rownames(comp_paths),
                                           pathways, "GO", collapse = FALSE)
    big_table <- do.call("rbind", lapply(rownames(comp_paths), function(path){
        gos <- path_annot[path_annot$paths == path, "funs"]
        if(length(gos) > 1 || !gos == ""){
            gos_ancst <- get_highest_sig_ancestor(go_terms = gos,
                                                  comp_go,
                                                  metaginfo = pathways)
        }else{
            gos_ancst <- data.frame(GO_term = "",
                                    GO_name = "",
                                    GO_adj_pval = "",
                                    Highest_Significant_Ancestor = "",
                                    HSA_name = "",
                                    HSA_adj_pval = "",
                                    stringsAsFactors = FALSE)
        }
        df <- cbind(path.id = path,
                    path.name = path_names[path],
                    comp_paths[path,c("UP/DOWN", "FDRp.value")],
                    gos_ancst)
        return(df)
    }))
    rownames(big_table) <- NULL

    return(big_table)
}


#' Normalize the pathway matrix by rows
#'
#' Due to the nature of the Hipathia method, the length of a pathway may
#' influence its signal rank. In order to compare signal values among
#' subpathways, we strongly recommend to normalize the matrix with this
#' normalization.
#'
#' This function removes the bias caused by the length of the subpathways
#' by dividing by the value obtained from running the method with a basal
#' value of 0.5 at each node.
#'
#' @param path_vals Matrix of the pathway values
#' @param metaginfo Pathways object
#'
#' @return Matrix of normalized pathway values
#'
#' @examples
#' data(path_vals)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' path_normalized <- normalize_paths(path_vals, pathways)
#'
#' @export
#'
normalize_paths <- function(path_vals, metaginfo){
    decomposed <- is_decomposed_matrix(path_vals)
    if(decomposed == TRUE){
        norm_factors <- metaginfo$path.norm[rownames(path_vals)]
        path_norm <- normalize_data(path_vals/(norm_factors*0.99+0.01),
                                    by_quantiles = FALSE,
                                    by_gene = FALSE,
                                    percentil = FALSE)
    }else{
        norm_factors <- metaginfo$eff.norm[rownames(path_vals)]
        path_norm <- normalize_data(path_vals/(norm_factors*0.99+0.01),
                                    by_quantiles = FALSE,
                                    by_gene = FALSE,
                                    percentil = FALSE)
    }
    return(path_norm)
}


is_decomposed_matrix <- function(mat){
    decomposed <- is_decomposed(rownames(mat))
    return(decomposed)
}

is_decomposed <- function(ids){
    lens <- lengths(sapply(ids, strsplit, "-"))
    if(length(unique(lens)) > 1)
        stop("Not unique type of labels")
    decomposed <- !(lens[1] == 3)
    return(decomposed)
}

