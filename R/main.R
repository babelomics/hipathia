##
## main.R
## Core functions of package Hipathia
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##

#' Computes the level of activation of the subpathways for each 
#' of the samples
#'
#' #@importFrom igraph
#'
#' @param genes_vals A SummarizedExperiment or matrix with the normalized 
#' expression values of the genes. Rows represent genes and columns represent 
#' samples. Rownames() must be accepted gene IDs.
#' @param metaginfo Pathways object
#' @param sel_assay Character or integer, indicating the assay to be processed 
#' in the SummarizedExperiment. Only applied if \code{genes_vals} is a 
#' \code{SummarizedExperiment}.Default is 1.
#' @param decompose Boolean, whether to compute the values for the decomposed
#' subpathways. By default, effector subpathways are computed.
#' @param maxnum Number of maximum iterations when iterating the signal
#' through the loops into the pathways
#' @param verbose Boolean, whether to show details about the results of
#' the execution of hipathia
#' @param tol Tolerance for the difference between two iterations when
#' iterating the signal through the loops into the pathways
#' @param test Boolean, whether to test the input objects. Default is TRUE.
#'
#' @return A MultiAssayExperiment object with the level of activation of the 
#' subpathways from
#' the pathways in \code{pathigraphs} for the experiment
#' with expression values in \code{genes_vals}.
#'
#' @examples
#' data(exp_data)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' results <- hipathia(exp_data, pathways, verbose = TRUE)
#' \dontrun{results <- hipathia(exp_data, pathways, decompose = TRUE, 
#' verbose = FALSE)}
#'
#' @export
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#'
hipathia <- function(genes_vals, metaginfo, sel_assay = 1, decompose = FALSE, 
                     maxnum = 100, verbose = TRUE, tol = 0.000001, test = TRUE){

    if(is(genes_vals, "SummarizedExperiment")){
        coldata <- colData(genes_vals)
        genes_vals <- assay(genes_vals, sel_assay)
    }else{
        cols <- colnames(genes_vals)
        coldata <- data.frame(cols = cols, stringsAsFactors = FALSE)
    }
    if(test == TRUE){
        if(is.null(genes_vals))
            stop("Missing input matrix")
        if(is.null(metaginfo))
            stop("Missing pathways object")
        test_matrix(genes_vals)
        test_pathways_object(metaginfo)
        test_tolerance(tol)
    }
    pathigraphs <- metaginfo$pathigraphs
    genes_vals <- add_missing_genes(genes_vals, genes = metaginfo$all.genes)
    results <- list()

    if(verbose == TRUE)
        cat("HiPathia processing...\n")

    results$by.path <- lapply(pathigraphs, function(pathigraph){

        if(verbose == TRUE)
            cat(pathigraph$path.id, "-", pathigraph$path.name, "\n")

        res <- list()
        res$nodes.vals <- nodes_values_from_genes(genes_vals, pathigraph$graph)

        if(decompose == FALSE){
            respaths <- all_path_values( res$nodes.vals,
                                         pathigraph$effector.subgraphs,
                                         maxnum = maxnum,
                                         tol = tol )
        }else{
            respaths <- all_path_values( res$nodes.vals,
                                         pathigraph$subgraphs,
                                         maxnum = maxnum,
                                         tol = tol )
        }
        res$path.vals <- respaths[[1]]
        res$convergence <- respaths[[2]]

        return(res)
    })

    paths <- do.call("rbind", lapply(results$by.path, function(x) x$path.vals))
    nodes <- do.call("rbind", lapply(results$by.path, function(x) x$nodes.vals))
    
    paths_rd <- DataFrame(feat.ID = rownames(paths), 
                          feat.name = get_path_names(metaginfo, 
                                                        rownames(paths)),
                          decomposed = decompose)
    nodes_rd <- as.data.frame(metaginfo$all.labelids[rownames(nodes),], 
                              stringsAsFactors = FALSE)
        
    paths_se <- SummarizedExperiment(list(paths = paths), rowData = paths_rd, 
                                     colData = coldata)
    nodes_se <- SummarizedExperiment(list(nodes = nodes), rowData = nodes_rd, 
                                     colData = coldata)
    resmae <- MultiAssayExperiment(list(paths = paths_se, nodes = nodes_se))
    
    return(resmae)
}


#' @importFrom DelayedArray colMins
nodes_values_from_genes <- function(genes_vals, ig, summ = "per90"){
    genes_list <- V(ig)$genesList
    names(genes_list) <- V(ig)$name
    genes_list <- genes_list[!grepl("_func", names(genes_list))]
    nodes_vals <- matrix(NA,
                         nrow = length(names(genes_list)),
                         ncol = ncol(genes_vals),
                         dimnames = list(names(genes_list),
                                         colnames(genes_vals)))
    for (node_name in names(genes_list)){
        genes <- genes_list[[node_name]]
        if( "/" %in% genes ){ #Then the node is a protein complex
            lists <- get_genes_lists( genes )
            probabilities_mat <- matrix(NA, nrow = 0, ncol = ncol(genes_vals))
            for( list1 in lists ){
                if( length(list1) > 1 ){
                    gv_l1 <- genes_vals[list1,,drop = FALSE]
                    prob <- summarize_probabilities(gv_l1, summ)
                }else{
                    prob <- genes_vals[list1,,drop = FALSE]
                }
                probabilities_mat <- rbind(probabilities_mat, prob)
            }
            nodes_vals[node_name,] <- colMins(probabilities_mat, na.rm = TRUE)
        }else{
            glist <- genes_list[[node_name]]
            if(length(glist) > 1){
                gv <- genes_vals[glist,,drop = FALSE]
                nodes_vals[node_name,] <- summarize_probabilities(gv, summ)
            }else if (length(glist) == 1 && !is.na(glist)){
                dm <- data.matrix(genes_vals[glist,,drop = FALSE])
                nodes_vals[node_name,] <- dm
            }else{
                nodes_vals[node_name,] <- rep(1, ncol(nodes_vals))
            }
        }
    }
    return(nodes_vals)
}


#' @importFrom matrixStats colMedians
#' @importFrom matrixStats colMeans2
#' @importFrom stats quantile
#' @importFrom DelayedArray colMaxs
#' @importFrom DelayedArray colMins
summarize_probabilities <- function(probabilities, summ = "per90"){
    if (summ == "mean"){
        prob <- colMeans2(probabilities, na.rm = TRUE)
    }else if(summ == "median"){
        prob <- colMedians(probabilities, na.rm = TRUE)
    }else if (summ == "max"){
        prob <- colMaxs(probabilities, na.rm = TRUE)
    }else if (summ == "min"){
        prob <- colMins(probabilities, na.rm = TRUE)
    }else if (summ == "per90"){
        prob <- apply(probabilities, 2, stats::quantile, 0.9, na.rm = TRUE)
    }else if (summ == "per95"){
        prob <- apply(probabilities, 2, stats::quantile, 0.95, na.rm = TRUE)
    }else if (summ == "per99"){
        prob <- apply(probabilities, 2, stats::quantile, 0.99, na.rm = TRUE)
    }else{
        stop("Summarizing probabilities option", summ, "is not valid")
    }
    return(prob)
}


get_genes_lists <- function( genes_list ){
    g_list <- NULL
    g_list_list <- list()
    while( length(genes_list) > 0 ){
        if( genes_list[[1]] != "/" ){
            g_list <- c( g_list, genes_list[[1]])
        }
        else{
            g_list_list[[length(g_list_list) + 1]] <- g_list
            g_list <- NULL
        }
        genes_list <- genes_list[-1]
    }
    g_list_list[[length(g_list_list) + 1]] <- g_list
    return(g_list_list)
}



all_path_values <- function( nodes_vals, subgraph_list, method = "maxmin",
                             maxnum = 100, tol = 0.000001, divide = FALSE,
                             response_tol = 0 ){
    path_vals <- matrix(0,
                        ncol = ncol(nodes_vals),
                        nrow = length(subgraph_list),
                        dimnames = list(names(subgraph_list),
                                        colnames(nodes_vals)))
    signal_dif <- list()
    for( path in names(subgraph_list)){
        dec_name <- unlist(strsplit(path, "\\-"))
        if(length(dec_name) == 4){
            ininodes <- paste("N", dec_name[2], dec_name[3], sep = "-")
            endnode <- paste("N", dec_name[2], dec_name[4], sep = "-")
        }else if(length(dec_name) == 3){
            endnode <- paste("N", dec_name[2], dec_name[3], sep = "-")
            sl <- subgraph_list[[path]]
            ininodes <- V(sl)$name[!V(sl)$name %in% get.edgelist(sl)[,2]]
        }else{
            stop("Error: Unknown path ID")
        }
        res <- path_value(nodes_vals,
                          subgraph_list[[path]],
                          ininodes,
                          endnode,
                          method,
                          maxnum = maxnum,
                          tol = tol,
                          divide = divide,
                          response_tol = response_tol)
        path_vals[path,] <- res[[1]]
        signal_dif[[path]] <- res[[2]]
    }
    return(list(path_vals, signal_dif))
}



path_value <- function( nodes_vals, subgraph, ininodes, endnode,
                        method = "maxmin", maxnum = 100, tol = 0.000001,
                        divide = FALSE, response_tol = 0 ){

    # Initialize lists
    ready <- ininodes
    # Initialize node values
    node_signal <- matrix(NA,
                          ncol = ncol(nodes_vals),
                          nrow = length(V(subgraph)),
                          dimnames = list(V(subgraph)$name,
                                          colnames(nodes_vals)))
    endnode_signal_dif <- 10

    num <- 0
    reached_last <- FALSE
    while( length(ready) > 0 && num <= maxnum){
        num <- num + 1
        actnode <- ready[[1]]
        old_signal <- node_signal[actnode,]

        # Compute node signal
        if(divide && actnode != endnode){
            nfol <- length(incident(subgraph, actnode, mode = "out"))
        }
        else{
            nfol <- 1
        }
        node_signal[actnode,] <- compute_node_signal(actnode,
                                                     nodes_vals[actnode,],
                                                     node_signal,
                                                     subgraph,
                                                     method,
                                                     response_tol) / nfol

        # Transmit signal
        nextnodes <- get.edgelist(subgraph)[incident(subgraph,
                                                     actnode, mode = "out"),2]
        dif <- old_signal - node_signal[actnode,]

        if(actnode == endnode){
            reached_last <- TRUE
            if(!all(is.na(dif)))
                endnode_signal_dif <- c(endnode_signal_dif, sqrt(sum(dif^2)))
            #num <- num+1
        }
        if(all(is.na(old_signal)) ||
           endnode_signal_dif[length(endnode_signal_dif)] > tol )
            ready <- unique(c(ready, nextnodes))
        ready <- ready[-1]
    }
    if(reached_last == FALSE){
        endnode_signal_dif <- NA
    }
    return(list(node_signal[endnode,], endnode_signal_dif))
}


#' @importFrom matrixStats colProds
#' @importFrom DelayedArray colMins
#' @importFrom DelayedArray colMaxs
compute_node_signal <- function(actnode, node_val, node_signal, subgraph,
                                method="maxmin", response_tol = 0){

    incis <- incident(subgraph, actnode, mode = "in")

    if(length(incis)==0){
        signal <- rep(1, length(node_val))

    } else {

        # get activators and inhibitors signal
        prevs <- get.edgelist(subgraph)[incis,1]
        input_signals <- node_signal[prevs,,drop = FALSE]
        nas <- is.na(input_signals[,1])
        prevs <- prevs[!nas]
        incis <- incis[!nas]
        input_signals <- input_signals[!nas,,drop = FALSE]
        typeincis <- E(subgraph)$relation[incis]
        activators <- typeincis == 1
        nactivators <- sum(activators)
        inhibitors <- typeincis == -1
        ninhibitors <- sum(inhibitors)
        activator_signals <- input_signals[activators,,drop = FALSE]
        inhibitor_signals <- input_signals[inhibitors,,drop = FALSE]

        if( method == "sum"){
            s1 <- prettyifelse(nactivators > 0,
                               colSums(activator_signals),
                               rep(1,length(node_val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               colSums(inhibitor_signals),
                               rep(0,length(node_val)))
            signal <- s1-s2
        }
        else if( method == "maxmin"){
            s1 <- prettyifelse(nactivators > 0,
                               colProds(1- activator_signals, na.rm = TRUE),
                               rep(0, length(node_val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               # colProds(rowMaxs(inhibitor_signals) +
                               #          rowMaxs(inhibitor_signals) -
                               #          inhibitor_signals),
                               apply(apply(inhibitor_signals, 1, max) +
                                         apply(inhibitor_signals, 1, min) -
                                         inhibitor_signals, 2, prod),
                               rep(1, length(node_val)))
            signal <- (1-s1)*s2
        }
        else if( method == "pond"){
            s1 <- prettyifelse(nactivators > 0,
                               colProds(1 - activator_signals),
                               rep(0, length(node_val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               colProds(1 - inhibitor_signals),
                               rep(1, length(node_val)))
            signal <- (1-s1)*s2
        }
        else if( method == "min"){
            s1 <- prettyifelse(nactivators > 0,
                               colMins(activator_signals),
                               rep(1, length(node_val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               1 - colMaxs(inhibitor_signals),
                               rep(1, length(node_val)))
            signal <- s1*s2
        }
        else {
            stop("Unknown propagation rule")
        }

        # If signal too low, signal do not propagate
        if(sum(nas) == 0 && signal < response_tol)
            signal <- rep(0, length(node_val))

    }

    signal[signal>1] <- 1
    signal[signal<0] <- 0
    signal <- signal * node_val

    return(signal)
}


prettyifelse <- function(test, una, olaotra){
    if(test){
        return(una)
    } else {
        return(olaotra)
    }
}



