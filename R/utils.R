# utils.R
# Written by Marta R. Hidalgo


translate.ids <- function(ids, xref){

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
#' Translation of the rownames IDs to Entrez IDs.
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
#' @examples data("brca_data")
#' trans.data <- translate.matrix(brca_data, "hsa")
#'
#' @return Matrix of gene expression with Entrez IDs as rownames.
#'
#' @export
#' @import hpAnnot
#'
translate.matrix <- function(exp, species, verbose=TRUE){

    xref <- load.xref(species)

    new_ids <- gsub("\\.[0123456789]$", "", rownames(exp))

    # translate ids
    tt <- translate.ids(new_ids, xref)

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
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' translated.names <- get.path.names(pathways, rownames(path_vals))
#'
#' @return A character vector including the readable names of the
#' subpathways IDs, in the same order as provided.
#' @export
#'
get.path.names <- function(metaginfo, names, maxchar=NULL){

    pathigraphs <- metaginfo$pathigraphs

    prettynames <- unlist(lapply(names, function(name){
        strname <- unlist(strsplit(name, "\\-"))
        pathway <- strname[2]
        labelid <- pathigraphs[[pathway]]$label.id
        if(length(strname) > 3){
            newini <- labelid[which(labelid[,1] ==
                                        paste("N", strname[2], strname[3],
                                              sep = "-")),2]
            newend <- labelid[which(labelid[,1] ==
                                        paste("N", strname[2], strname[4],
                                              sep = "-")),2]
            if(grepl("_func", strname[4])){
                otroend <- labelid[which(labelid[,1] ==
                                             paste("N", strname[2], strname[3],
                                                   sep = "-")),2]
                name <- paste0(pathigraphs[[pathway]]$path.name,
                               ": ", newini, " -> ", otroend, " -> ", newend )
            }else{
                name <- paste0(pathigraphs[[pathway]]$path.name,
                               ": ", newini, " -> ", newend )
            }
        }else{
            path.name <- pathigraphs[[pathway]]$path.name
            node <- labelid[which(labelid[,1] ==
                                      paste("N", strname[2], strname[3],
                                            sep = "-")),2]
            name <- paste0(path.name, ": ", node)
        }
    }))

    if(!is.null(maxchar))
        prettynames <- clip.names(prettynames, maxchar = maxchar)

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
igraphs.upgrade <- function(metaginfo){

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
#' of the functions (like \code{hipathia} or \code{quantify.terms}) in the
#' package. These functions will analyze all the pathways included in this
#' object. By default, all available pathways are load. In order to restrict
#' the analysis to a predefined set of pathways, specify the set of pathways
#' to load with the parameter \code{pathways.list}.
#'
#' @param species Species of the samples.
#' @param pathways.list Vector of the IDs of the pathways to load. By default
#' all available pathways are load.
#'
#' @examples
#' pathways <- load.pathways("hsa")   # Loads all pathways for human
#' pathways <- load.pathways("mmu", c("mmu03320", "mmu04024", "mmu05200"))
#'    # Loads pathways 03320, 04024 and 05200 for mouse
#'
#' @return An pathways object including
#' * \code{species} Species to which the pathways are related.
#' * \code{pathigraphs} List of Pathigraph objects. Each Pathigraph contains
#' the necessary information of a pathway for it to be analyzed
#' with \code{Hipathia}.
#' * \code{all.genes} List of all the genes included in the selection of
#' pathways stored in \code{pathigraphs}.
#' * \code{eff.norm} Vector of normalization values for effector subpathways.
#' * \code{path.norm} Vector of normalization values for decomposed
#' subpathways.
#'
#' @export
#' @import hpAnnot
#'
load.pathways <- function(species, pathways.list = NULL){
    metaginfo <- load.mgi(species)
    metaginfo <- filter.pathways(metaginfo, pathways.list = pathways.list)
    metaginfo <- igraphs.upgrade(metaginfo)
    cat(paste0("Loaded ", length(metaginfo$pathigraphs), " pathways\n"))
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
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' pathways.list <- get.pathways.list(pathways)
#'
#' @export
#'
get.pathways.list <- function(metaginfo){
    return(names(metaginfo$pathigraphs))
}


filter.pathways <- function(metaginfo, pathways.list = NULL){
    if(!is.null(pathways.list)){
        metaginfo$pathigraphs <- metaginfo$pathigraphs[pathways.list]
        metaginfo$all.genes <- all.needed.genes(metaginfo$pathigraphs)
        metaginfo$path.norm <- metaginfo$path.norm[
            sapply(names(metaginfo$path.norm), function(x){
                unlist(strsplit(x, split = "-"))[2]}) %in% pathways.list]
        metaginfo$eff.norm <- metaginfo$eff.norm[
            sapply(names(metaginfo$eff.norm), function(x){
                unlist(strsplit(x, split = "-"))[2]}) %in% pathways.list]
    }
    return(metaginfo)
}


all.needed.genes <- function(pathigraphs){
    genes <- unique(unlist(sapply(pathigraphs, function(x){
        unique(unlist(V(x$graph)$genesList))
    })))
    return(genes[!is.na(genes) & genes!="/"])
}



#' Gets the matrix of subpathway activation values
#'
#' @description
#' This function returns the matrix with the levels of activation of each
#' subpathway for each sample. Rows represent the subpathways and columns
#' represent the samples. Each cell is the value of activation of a subpathway
#' in a sample.
#'
#' Rownames are the IDs of the subpathways. In order to transform IDs into
#' readable names, use \code{get.path.names}.
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
#'
#' @examples
#' data(results)
#' path.vals <- get.paths.matrix(results)
#'
#' @return Matrix with the levels of activation of each decomposed subpathway
#' for each sample.
#' @export
#'
get.paths.matrix <- function(results){
    return(results$all$path.vals)
}





clip.names <- function(snames, maxchar = 30){
    sapply(snames, function(x) {
        if(nchar(x) > maxchar){
            return(paste0(substr(x, 1 , maxchar - 3), "..."))
        } else {
            return(x)
        }
    })
}



#' @importFrom stats median
add.missing.genes <- function(exp.data, genes, default=NULL){
    if(is.null(default))
        default <- stats::median(exp.data)
    missing_genes <- setdiff(genes, rownames(exp.data))
    if(length(missing_genes > 0)){
        if(ncol(exp.data) == 1){
            fakemat <- default +
                as.matrix(mat.or.vec(nr = length(missing_genes),
                                     nc = ncol(exp.data)))
        } else {
            fakemat <- default + mat.or.vec(nr = length(missing_genes),
                                            nc = ncol(exp.data))
        }
        rownames(fakemat) <- missing_genes
        colnames(fakemat) <- colnames(exp.data)
        exp.data <- rbind(exp.data,fakemat)
        #message("----------------------------------------------------")
        message("Added missing genes: ",
                length(missing_genes),
                " (",
                round(length(missing_genes)/nrow(exp.data)*100, digits = 2),
                "%)")
        message("----------------------------------------------------")
    }
    return(exp.data)
}




#' Head function for matrices
#'
#' Shows the first \code{n} rows and the first \code{n} columns of a matrix,
#' in case the matrix has more than \code{n+5} rows or columns.
#' Otherwise, it shows all the rows or columns, respectively.
#'
#' @param mat Matrix to be shown
#' @param n Number of rows and columns
#'
#' @examples mat <- matrix(rnorm(100), ncol = 10)
#' hhead(mat)
#' hhead(mat, 3)
#' hhead(mat, 7)
#'
#' @return Matrix with as much as \code{n} rows and columns.
#'
#' @export
hhead <- function(mat, n = 5){
    if(ncol(mat) >= n+5){
        if(nrow(mat) >= n+5){
            mat[1:n,1:n]
        }else{
            mat[,1:n]
        }
    }else{
        if(nrow(mat) >= n+5){
            mat[1:n,]
        }else{
            mat
        }
    }
}



get.effnode.id <- function(path.name){
    path.split <- unlist(strsplit(path.name, split="-"))
    if(!path.split[1] == "P") return(NA)
    if(length(path.split) > 3){
        paste("N", path.split[2], path.split[4], sep="-")
    }else{
        paste("N", path.split[2], path.split[3], sep="-")
    }
}


get.ininode.id <- function(path.name){
    path.split <- unlist(strsplit(path.name, split="-"))
    if(length(path.split) > 3){
        paste("N", path.split[2], path.split[3], sep="-")
    }else{
        return(NA)
    }
}


get.effpath.id <- function(node.name){
    node.split <- unlist(strsplit(node.name, split="-"))
    if(!node.split[1] == "N") return(NA)
    gsub("N", "P", node.name)
}


#' Get Pathways functional annotations
#'
#' Get functional annotation of the pathways, either for a particular
#' annotation or a stored one.
#'
#' @param pathway.names Character vector of the names of the pathways
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
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' pathway.names <- c("P-hsa03320-37", "P-hsa03320-61", "P-hsa03320-46",
#' "P-hsa03320-57", "P-hsa03320-64", "P-hsa03320-47", "P-hsa03320-65")
#' get.pathways.annotations(pathway.names, pathways, "GO")
#' get.pathways.annotations(pathway.names, pathways, "uniprot")
#'
#' @export
#' @import hpAnnot
#'
get.pathways.annotations <- function(pathway.names, metaginfo, dbannot,
                                     collapse = TRUE){

    if(is.character(dbannot)){
        annofuns <- load.annofuns(dbannot, metaginfo$species)
    }else{
        annofuns <- annotate.paths(metaginfo, dbannot)
    }

    annofuns$funs[is.na(annofuns$funs)] <- ""
    decomposed <- length(unlist(strsplit(pathway.names[1], split = "-"))) == 4
    if(decomposed == TRUE)
        pathway.names <- sapply(pathway.names, function(n)
            paste(unlist(strsplit(n, split="-"))[c(1,2,4)], collapse = "-"))

    if(collapse == TRUE){
        miniaf <- do.call("rbind", lapply(pathway.names, function(path){
            af <- annofuns[annofuns$paths == path,]
            data.frame(effector.nodes =
                           annofuns$effector.nodes[annofuns$paths == path][1],
                       paths = path,
                       funs = paste(annofuns$funs[annofuns$paths == path],
                                    collapse = ", "),
                       stringsAsFactors = FALSE)
        }))
    }else{
        miniaf <- do.call("rbind", lapply(pathway.names, function(path){
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
#' by \code{do.wilcoxon}.
#' @param metaginfo Pathways object
#' @param unique Boolean, whether to return only one highest significant GO
#' ancestor or all of them. By default, TRUE.
#' @param pval P-value cut-off. Default values is set to 0.05.
#'
#' @return highest common ancestors
#'
#' #@export
#' @import hpAnnot
#'
get.highest.sig.ancestor <- function(go_terms, go_comp, metaginfo,
                                     unique = TRUE, pval = 0.05){

    go_bp_frame <- load.gobp.frame()
    go_bp_net <- load.gobp.net()

    # Relacionar GO term con su etiqueta
    go_labels <- sapply(go_terms,
                        function(term) go_bp_frame[go_bp_frame$name ==
                                                       term, "id"])
    go_comp$labels <- sapply(rownames(go_comp),
                             function(term) go_bp_frame[go_bp_frame$name ==
                                                            term, "id"])
    go_comp$terms <- rownames(go_comp)
    rownames(go_comp) <- go_comp$labels
    sig_go_labels <- go_comp[go_comp$FDRp.value < pval, "labels"]

    # Encontrar GO superior
    sup <- do.call("rbind", lapply(go_labels, function(label){
        # print(label)
        short <- shortest.paths(go_bp_net, label, mode = "in")[1,]
        ancestors <- names(short)[!short == "Inf"]
        sig.ancestors <- intersect(ancestors, sig_go_labels)
        if(length(sig.ancestors) > 0 ){
            if(length(sig.ancestors) > 1)
                sig.ancestors <- setdiff(sig.ancestors, label)
            sig.levels <- go_bp_frame[sig.ancestors, "level"]
            if(unique == TRUE){
                highest.ancestors <- go_bp_frame[sig.ancestors,
                                                 "id"][sig.levels ==
                                                           min(sig.levels)][1]
            }else{
                highest.ancestors <- go_bp_frame[sig.ancestors,
                                                 "id"][sig.levels ==
                                                           min(sig.levels)]
            }
            df <- data.frame(GO_term = label,
                             GO_name = go_bp_frame[label, "name"],
                             GO_adj_pval = go_comp[label, "FDRp.value"],
                             Highest_Significant_Ancestor = highest.ancestors,
                             HSA_name = go_bp_frame[highest.ancestors, "name"],
                             HSA_adj_pval = go_comp[highest.ancestors,
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
#' @param comp.paths Wilcoxon comparison of the matrix of pathways values
#' as returned by \code{do.wilcoxon}.
#' @param comp.go Wilcoxon comparison of the matrix of GO values as
#' returned by \code{do.wilcoxon}.
#' @param pval P-value cut-off. Default values is set to 0.05.
#'
#' @return Table of comparisons with Highest common ancestors
#'
#' @examples
#' data(comp)
#' data(go_vals)
#' data(brca_design)
#' data(path_vals)
#' sample.group <- brca_design[colnames(path_vals),"group"]
#' comp.go <- do.wilcoxon(go_vals, sample.group, g1 = "Tumor", g2 = "Normal")
#' \dontrun{pathways <- load.pathways(species = "hsa", pathways.list =
#' c("hsa03320", "hsa04012"))
#' table <- paths.to.go.ancestor(pathways, comp, comp.go)}
#'
#' @export
#'
paths.to.go.ancestor <- function(pathways, comp.paths, comp.go, pval = 0.05){
    path.names <- get.path.names(pathways, rownames(comp.paths))
    names(path.names) <- rownames(comp.paths)
    path.annot <- get.pathways.annotations(pathway.names = rownames(comp.paths),
                                           pathways, "GO", collapse = FALSE)
    big.table <- do.call("rbind", lapply(rownames(comp.paths), function(path){
        gos <- path.annot[path.annot$paths == path, "funs"]
        if(length(gos) > 1 || !gos == ""){
            gos.ancst <- get.highest.sig.ancestor(go_terms = gos,
                                                  comp.go,
                                                  metaginfo = pathways)
        }else{
            gos.ancst <- data.frame(GO_term = "",
                                    GO_name = "",
                                    GO_adj_pval = "",
                                    Highest_Significant_Ancestor = "",
                                    HSA_name = "",
                                    HSA_adj_pval = "",
                                    stringsAsFactors = FALSE)
        }
        df <- cbind(path.id = path,
                    path.name = path.names[path],
                    comp.paths[path,c("UP/DOWN", "FDRp.value")],
                    gos.ancst)
        return(df)
    }))
    rownames(big.table) <- NULL

    return(big.table)
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
#' @param path.vals Matrix of the pathway values
#' @param metaginfo Pathways object
#'
#' @return Matrix of normalized pathway values
#'
#' @examples
#' data(path_vals)
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' path.normalized <- normalize.paths(path_vals, pathways)
#'
#' @export
#'
normalize.paths <- function(path.vals, metaginfo){
    decomposed <- !all(sapply(sapply(rownames(path.vals),
                                     strsplit, "-"), length) == 3)
    if(decomposed == TRUE){
        norm.factors <- metaginfo$path.norm[rownames(path.vals)]
        path.norm <- normalize.data(path.vals/(norm.factors*0.99+0.01),
                                    by.quantiles = FALSE,
                                    by.gene = FALSE,
                                    percentil = FALSE)
    }else{
        norm.factors <- metaginfo$eff.norm[rownames(path.vals)]
        path.norm <- normalize.data(path.vals/(norm.factors*0.99+0.01),
                                    by.quantiles = FALSE,
                                    by.gene = FALSE,
                                    percentil = FALSE)
    }
    return(path.norm)
}

