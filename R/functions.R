##
## functions.R
## Enrichment functions
##
## Written by Marta R. Hidalgo, Jose Carbonell-Caballero
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##


#' Annotates functions to pathways
#'
#' Annotates functions from a database to each pathway
#'
#' @param metaginfo Pathways object
#' @param dbannot Either a string indicating which precomputed annotation to
#' use ("uniprot" for Uniprot Keywords or "GO" for Gene Ontology terms), or
#' a dataframe with the annotation of the genes to the functions. First
#' column are gene symbols, second column the functions.
#'
#' @return Object of annotations from pathways to functions
#'
#' #@examples
#' #pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' #"hsa04012"))
#' #annotate_paths(pathways, "GO")
#'
#' #@export
#'
annotate_paths <- function(metaginfo, dbannot){

    pathigraphs <- metaginfo$pathigraphs
    entrez2hgnc <- load_entrez_hgnc(metaginfo$species)
    if(is.character(dbannot) & length(dbannot) == 1)
        dbannot <- load_annots(dbannot, metaginfo$species)

    annofuns <- do.call("rbind", lapply(pathigraphs,function(pathigraph){

        new_pathigraph <- pathigraph
        vs <- V(new_pathigraph$graph)$name
        new_pathigraph$graph <- induced.subgraph(new_pathigraph$graph,
                                                 vs[!grepl("_func", vs)])
        funs <- get_pathway_functions(new_pathigraph,
                                      dbannot,
                                      entrez2hgnc,
                                      use_last_nodes = TRUE,
                                      unique = FALSE)
        paths <- lapply(names(funs), function(path){
            rep(path, times = length(funs[[path]]))
        })
        df <- data.frame(effector_nodes = unlist(paths),
                         paths = gsub("N", "P", unlist(paths)),
                         funs = unlist(funs),
                         stringsAsFactors = FALSE)

    }))

    return(annofuns)

}


#' Computes the level of activation of the functions related to the
#' previously computed subpathways
#'
#' @param results List of results as returned by the \code{hipathia} function
#' @param metaginfo Pathways object
#' @param dbannot Either a string indicating which precomputed annotation to
#' use ("uniprot" for Uniprot Keywords or "GO" for Gene Ontology terms), or
#' a dataframe with the annotation of the genes to the functions. First
#' column are gene symbols, second column the functions.
#' @param out_matrix Boolean, whther the output object should be a matrix 
#' object. Default is FALSE, returning a SummarizedExperiment object.
#' @param normalize Boolean, whether to normalize the matrix of pathway
#' values with \code{normalize_paths} before quantifying the signal. Due to
#' the nature of the Hipathia method, in which the length of each pathway may
#' alter its signal rank, we strongly recommend to perform this normalization.
#' This normalization removes the bias. Default is set to TRUE.
#'
#' @return Matrix with the level of activation of the functions in
#' \code{dbannot}
#'
#' @examples
#' data(results)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' go_values <- quantify_terms(results, pathways, "GO")
#' uniprot_values <- quantify_terms(results, pathways, "uniprot")
#'
#' @export
#' @importFrom matrixStats colMeans2
#' @importFrom matrixStats colProds
#'
quantify_terms <- function(results, metaginfo, dbannot, out_matrix = FALSE, 
                           normalize = TRUE){

    method = "mean"
    species <- metaginfo$species
    path_vals <- assay(results[["paths"]])
    first <- rownames(path_vals)[1]
    decomposed <- is_decomposed_matrix(path_vals)
    if(decomposed == TRUE)
        stop("Function terms can only be computed from NOT decomposed subpaths")
    
    if(is.character(dbannot) & length(dbannot) == 1){
        annofuns <- load_annofuns(dbannot, species)
    }else{
        annofuns <- annotate_paths(metaginfo, dbannot)
    }
    annofuns <- annofuns[!is.na(annofuns$funs),]
    annofuns$pathway <- sapply(strsplit(annofuns$paths, split = "-"), "[[", 2)
    annofuns <- annofuns[annofuns$pathway %in% names(metaginfo$pathigraphs),]
    fun_names <- unique(annofuns$funs)

    if(normalize == TRUE)
        path_vals <- normalize_paths(path_vals, metaginfo)

    fun_vals <- matrix(0,
                       ncol = ncol(path_vals),
                       nrow = length(fun_names),
                       dimnames = list(fun_names, colnames(path_vals)))
    message("Quantified functions: ", nrow(fun_vals))
    for(fun in fun_names){
        paths <- annofuns$paths[annofuns$funs == fun]
        if(method == "mean"){
            fun_vals[fun,] <- colMeans2(path_vals[paths,,drop = FALSE], 
                                        na.rm = TRUE)
        }else if(method == "signal"){
            minimat <- 1 - path_vals[paths,,drop = FALSE]
            fun_vals[fun,] <- 1 - colProds(minimat, na.rm = TRUE)
        }
    }
    if(out_matrix == FALSE){
        cd <- colData(results[["paths"]])
        if(is.character(dbannot) & length(dbannot) == 1 & dbannot == "GO"){
            rd <- DataFrame(feat.name = get_go_names(rownames(fun_vals), 
                                                     metaginfo$species))
            fun_vals <- SummarizedExperiment(list(terms = fun_vals), 
                                             colData = cd, 
                                             rowData = rd)
        }else{
            fun_vals <- SummarizedExperiment(list(terms = fun_vals), 
                                             colData = cd)
        }
    }
    return(fun_vals)
    
}

get_entrez_function <- function(entrezs, entrez2hgnc, dbannot){
    hgncs <- entrez2hgnc[entrez2hgnc[,1] %in% entrezs,2]
    path_functions <- dbannot[dbannot$gene %in% hgncs,]

    return(path_functions)
}


#' @importFrom stats fisher.test
#' @importFrom stats p.adjust
enrichment <- function(path_functions, dbannot, na_rm = TRUE){
    if(na_rm){
        path_functions <- path_functions[!is.na(path_functions[,1]),]
    }
    hgncs <- unique(path_functions[,1])
    compdb <- dbannot[-match(hgncs,dbannot[,1]) ,]

    term_counts <- table(path_functions[,2])
    terms <- names(term_counts)
    n1 <- length(hgncs)

    miniannot <- dbannot[dbannot[,2] %in% terms,]
    all_term_counts <- table(miniannot[,2])[terms]
    n2 <- length(unique(compdb[,1]))

    term_pvalues <- sapply(terms,function(x){
        mat <- matrix(byrow=TRUE,
                      c(term_counts[x], n1-term_counts[x],
                        all_term_counts[x], n2-all_term_counts[x]),
                      nrow = 2,
                      dimnames = list(c("path","rest"), c("yes","no")))
        test <- stats::fisher.test(mat, alternative = "greater")
        return(test$p.value)
    })

    pa <- stats::p.adjust(term_pvalues, "fdr")
    out <- data.frame(term = names(term_pvalues),
                      pvalue = term_pvalues,
                      adj.pvalue = pa)
    out <- out[order(out$pvalue),]

    return(out)

}

#' Returns functions related to a pathway
#'
#' @param pathigraph Pathway object
#' @param dbannot Dataframe with the annotation of the genes to the functions.
#' First column are gene symbols, second column the functions.
#' @param entrez2hgnc Relation between Entrez and HGNC genes.
#' @param use_last_nodes Boolean, whether to annotate functions to the last
#' nodes of the pathways or not. If FALSE, functions will refer to all the nodes
#' of the pathway.
#' @param unique Boolean, whether to return the first function for each path.
#'
#' @return List of annotations from pathways to functions
#'
get_pathway_functions <- function(pathigraph, dbannot, entrez2hgnc,
                                  use_last_nodes = TRUE, unique = TRUE){

    g <- pathigraph$graph
    last_nodes <- gsub("_func", "", get_last_node(g))
    if(use_last_nodes == TRUE){
        ebn <- V(g)$genesList[which(V(g)$name %in% last_nodes)]
    } else {
        ebn <- get.vertex.attribute(g, "genesList")
    }
    entrezs <- unique(unlist(ebn))
    entrezs <- entrezs[!is.na(entrezs)]
    gpf <- unique(get_entrez_function(entrezs, entrez2hgnc, dbannot))
    if(nrow(gpf) == 0){
        last_node_functions <- rep(NA, length(last_nodes))
        names(last_node_functions) <- last_nodes
        last_node_functions <- as.list(last_node_functions)
    } else {
        fe <- enrichment(gpf,dbannot)
        #enriched <- as.character(fe$term[fe$pvalue < 0.05])

        filtdbannot <- dbannot#[ dbannot[,2] %in% enriched,]
        filtdbannot$rank <- fe[match(filtdbannot[,2], fe$term), "pvalue"]
        filtdbannot$rank[is.na(filtdbannot$rank)] <- 1
        filtdbannot <- filtdbannot[order(filtdbannot$rank),]

        last_node_functions <- lapply(last_nodes,
                                      get_best_node_functions,
                                      pathigraph,
                                      entrez2hgnc,
                                      filtdbannot,
                                      unique = unique)
        names(last_node_functions) <- last_nodes
        idx <- lengths(last_node_functions) == 0
        last_node_functions[idx] <- NA
    }
    return(last_node_functions)

}



get_best_node_functions <- function(node, pathigraph, entrez2hgnc, filtdbannot,
                                    unique = TRUE){
    g <- pathigraph$graph
    entrezs <- V(g)$genesList[[which(V(g)$name == node)]]
    node_funcs <- unique(unlist(lapply(entrezs,
                                       get_best_gene_functions,
                                       entrez2hgnc = entrez2hgnc,
                                       filtdbannot = filtdbannot,
                                       unique = unique)))
    node_funcs <- node_funcs[!is.na(node_funcs)]
    return(node_funcs)
}

get_best_gene_functions <- function(entrez, entrez2hgnc, filtdbannot,
                                    unique = TRUE){
    hgncs <- entrez2hgnc[entrez2hgnc[,1] %in% entrez,2]
    funcs <- filtdbannot[filtdbannot$gene %in% hgncs,]
    funcs <- unique(funcs)
    if(unique == TRUE){
        return(funcs[1,2])
    }else{
        return(funcs[,2])
    }
}

get_last_node <- function(subgraph){
    V(subgraph)$name[!(V(subgraph)$name %in% get.edgelist(subgraph)[,1])]
}



