##
## main.R
## Core functions of package Hipathia
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##

#' Computes the level of activation of the subpathways for each of the samples
#'
#' #@importFrom igraph
#'
#' @param genes.vals A matrix with the normalized expression values of the
#' genes. Rows represent genes and columns represent samples.
#' Rownames() must be accepted gene IDs.
#' @param metaginfo Pathways object
#' @param decompose Boolean, whether to compute the values for the decomposed
#' subpathways. By default, effector subpathways are computed.
#' @param maxnum Number of maximum iterations when iterating the signal
#' through the loops into the pathways
#' @param verbose Boolean, whether to show details about the results of
#' the execution of hipathia
#' @param tol Tolerance for the difference between two iterations when
#' iterating the signal through the loops into the pathways
#'
#' @return An object with the level of activation of the subpathways from
#' the pathways in \code{pathigraphs} for the experiment
#' with expression values in \code{genes.vals}.
#'
#' @examples
#' data(exp_data)
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' results <- hipathia(exp_data, pathways, verbose = TRUE)
#' results <- hipathia(exp_data, pathways, decompose = TRUE, verbose = FALSE)
#'
#' @export
#'
hipathia <- function(genes.vals, metaginfo, decompose = FALSE, maxnum = 100,
                     verbose = TRUE, tol = 0.000001, test = TRUE){

    if(test == TRUE){
        if(is.null(genes.vals))
            stop("Missing input matrix")
        if(is.null(metaginfo))
            stop("Missing pathways object")
        test.matrix(genes.vals)
        test.pathways.object(metaginfo)
        test.tolerance(tol)
    }
    pathigraphs <- metaginfo$pathigraphs
    genes.vals <- add.missing.genes(genes.vals, genes = metaginfo$all.genes)
    results <- list()

    if(verbose == TRUE)
        cat("HiPathia processing...\n")

    results$by.path <- lapply(pathigraphs, function(pathigraph){

        if(verbose == TRUE)
            cat(pathigraph$path.id, "-", pathigraph$path.name, "\n")

        res <- list()
        res$nodes.vals <- nodes.values.from.genes(genes.vals, pathigraph$graph)

        if(decompose == FALSE){
            respaths <- all.path.values( res$nodes.vals,
                                         pathigraph$effector.subgraphs,
                                         maxnum = maxnum,
                                         tol = tol )
        }else{
            respaths <- all.path.values( res$nodes.vals,
                                         pathigraph$subgraphs,
                                         maxnum = maxnum,
                                         tol = tol )
        }
        res$path.vals <- respaths[[1]]
        res$convergence <- respaths[[2]]

        return(res)
    })

    results$all$path.vals <- do.call("rbind", lapply(results$by.path, 
                                                     function(x) x$path.vals))
    results$all$nodes.vals <- do.call("rbind", lapply(results$by.path, 
                                                      function(x) x$nodes.vals))
    return(results)
}


nodes.values.from.genes <- function(genes.vals, ig, summ = "per90"){
    genes.list <- V(ig)$genesList
    names(genes.list) <- V(ig)$name
    genes.list <- genes.list[!grepl("_func", names(genes.list))]
    nodes.vals <- matrix(NA,
                         nrow = length(names(genes.list)),
                         ncol = ncol(genes.vals),
                         dimnames = list(names(genes.list),
                                         colnames(genes.vals)))
    for (node.name in names(genes.list)){
        genes <- genes.list[[node.name]]
        if( "/" %in% genes ){ #Then the node is a protein complex
            lists <- get.genes.lists( genes )
            probabilities.mat <- matrix(NA, nrow = 0, ncol = ncol(genes.vals))
            for( list1 in lists ){
                if( length(list1) > 1 ){
                    prob <- summarize.probabilities(genes.vals[list1,,
                                                               drop = FALSE],
                                                    summ)
                }else{
                    prob <- genes.vals[list1,,drop = FALSE]
                }
                probabilities.mat <- rbind( probabilities.mat, prob )
            }
            nodes.vals[node.name,] <- apply(probabilities.mat, 2, min)
        }else{
            glist <- genes.list[[node.name]]
            if (length(glist) > 1){
                nodes.vals[node.name,] <-
                    summarize.probabilities(genes.vals[glist,,drop = FALSE],
                                            summ)
            }else if (length(glist) == 1 && !is.na(glist)){
                nodes.vals[node.name,] <- data.matrix(genes.vals[glist,,
                                                                 drop = FALSE])
            }else{
                nodes.vals[node.name,] <- rep(1, ncol(nodes.vals))
            }
        }
    }
    return(nodes.vals)
}


#' @importFrom stats median
#' @importFrom stats quantile
summarize.probabilities <- function(probabilities, summ = "per90"){
    if (summ == "mean"){
        prob <- apply(probabilities, 2, mean, na.rm = TRUE)
    }else if(summ == "median"){
        prob <- apply(probabilities, 2, stats::median, na.rm = TRUE)
    }else if (summ == "max"){
        prob <- apply(probabilities, 2, max, na.rm = TRUE)
    }else if (summ == "min"){
        prob <- apply(probabilities, 2, min, na.rm = TRUE)
    }else if (summ == "per90"){
        prob <- apply(probabilities, 2, stats::quantile, 0.9, na.rm = TRUE)
    }else if (summ == "per95"){
        prob <- apply(probabilities, 2, stats::quantile, 0.95, na.rm = TRUE)
    }else if (summ == "per99"){
        prob <- apply(probabilities, 2, stats::quantile, 0.99, na.rm = TRUE)
    }else{
        print (paste("The option for summarizing the probabilities",
                     summ, "is not valid"))
    }
    return(prob)
}


get.genes.lists <- function( genes.list ){
    g.list <- NULL
    g.list.list <- list()
    while( length(genes.list) > 0 ){
        if( genes.list[[1]] != "/" ){
            g.list <- c( g.list, genes.list[[1]])
        }
        else{
            g.list.list[[length(g.list.list) + 1]] <- g.list
            g.list <- NULL
        }
        genes.list <- genes.list[-1]
    }
    g.list.list[[length(g.list.list) + 1]] <- g.list
    return(g.list.list)
}



all.path.values <- function( nodes.vals, subgraph.list, method = "maxmin",
                             maxnum = 100, tol = 0.000001, divide = FALSE,
                             response.tol = 0 ){
    path.vals <- matrix(0,
                        ncol = ncol(nodes.vals),
                        nrow = length(subgraph.list),
                        dimnames = list(names(subgraph.list),
                                        colnames(nodes.vals)))
    signal.dif <- list()
    for( path in names(subgraph.list)){
        dec.name <- unlist(strsplit(path, "\\-"))
        if(length(dec.name) == 4){
            ininodes <- paste("N", dec.name[2], dec.name[3], sep = "-")
            endnode <- paste("N", dec.name[2], dec.name[4], sep = "-")
        }else if(length(dec.name) == 3){
            endnode <- paste("N", dec.name[2], dec.name[3], sep = "-")
            sl <- subgraph.list[[path]]
            ininodes <- V(sl)$name[!V(sl)$name %in% get.edgelist(sl)[,2]]
        }else{
            stop("Error: Unknown path ID")
        }
        res <- path.value(nodes.vals,
                          subgraph.list[[path]],
                          ininodes,
                          endnode,
                          method,
                          maxnum = maxnum,
                          tol = tol,
                          divide = divide,
                          response.tol = response.tol)
        path.vals[path,] <- res[[1]]
        signal.dif[[path]] <- res[[2]]
    }
    return(list(path.vals, signal.dif))
}



path.value <- function( nodes.vals, subgraph, ininodes, endnode,
                        method = "maxmin", maxnum = 100, tol = 0.000001,
                        divide = FALSE, response.tol = 0 ){

    # Initialize lists
    ready <- ininodes
    # Initialize node values
    node.signal <- matrix(NA,
                          ncol = ncol(nodes.vals),
                          nrow = length(V(subgraph)),
                          dimnames = list(V(subgraph)$name,
                                          colnames(nodes.vals)))
    endnode.signal.dif <- 10

    num <- 0
    reached_last <- FALSE
    while( length(ready) > 0 && num <= maxnum){
        num <- num + 1
        actnode <- ready[[1]]
        old.signal <- node.signal[actnode,]

        # Compute node signal
        if(divide && actnode != endnode){
            nfol <- length(incident(subgraph, actnode, mode = "out"))
        }
        else{
            nfol <- 1
        }
        node.signal[actnode,] <- compute.node.signal(actnode,
                                                     nodes.vals[actnode,],
                                                     node.signal,
                                                     subgraph,
                                                     method,
                                                     response.tol) / nfol

        # Transmit signal
        nextnodes <- get.edgelist(subgraph)[incident(subgraph,
                                                     actnode, mode = "out"),2]
        dif <- old.signal - node.signal[actnode,]

        if(actnode == endnode){
            reached_last <- TRUE
            if(!all(is.na(dif)))
                endnode.signal.dif <- c(endnode.signal.dif, sqrt(sum(dif^2)))
            #num <- num+1
        }
        if(all(is.na(old.signal)) ||
           endnode.signal.dif[length(endnode.signal.dif)] > tol )
            ready <- unique(c(ready, nextnodes))
        ready <- ready[-1]
    }
    if(reached_last == FALSE){
        endnode.signal.dif <- NA
    }
    return(list(node.signal[endnode,], endnode.signal.dif))
}


compute.node.signal <- function(actnode, node.val, node.signal, subgraph,
                                method="maxmin", response.tol = 0){

    incis <- incident(subgraph, actnode, mode = "in")

    if(length(incis)==0){
        signal <- rep(1, length(node.val))

    } else {

        # get activators and inhibitors signal
        prevs <- get.edgelist(subgraph)[incis,1]
        input_signals <- node.signal[prevs,,drop = FALSE]
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
                               rep(1,length(node.val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               colSums(inhibitor_signals),
                               rep(0,length(node.val)))
            signal <- s1-s2
        }
        else if( method == "maxmin"){
            s1 <- prettyifelse(nactivators > 0,
                               apply(1- activator_signals, 2, prod),
                               rep(0, length(node.val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               apply(apply(inhibitor_signals, 1, max) +
                                         apply(inhibitor_signals, 1, min) -
                                         inhibitor_signals, 2, prod),
                               rep(1, length(node.val)))
            signal <- (1-s1)*s2
        }
        else if( method == "pond"){
            s1 <- prettyifelse(nactivators > 0,
                               apply(1 - activator_signals, 2, prod),
                               rep(0, length(node.val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               apply(1 - inhibitor_signals, 2, prod),
                               rep(1, length(node.val)))
            signal <- (1-s1)*s2
        }
        else if( method == "min"){
            s1 <- prettyifelse(nactivators > 0,
                               apply(activator_signals,2,min),
                               rep(1, length(node.val)))
            s2 <- prettyifelse(ninhibitors > 0,
                               1 - apply(inhibitor_signals,2,max),
                               rep(1, length(node.val)))
            signal <- s1*s2
        }
        else {
            stop("Unknown propagation rule")
        }

        # If signal too low, signal do not propagate
        if(sum(nas) == 0 && signal < response.tol)
            signal <- rep(0, length(node.val))

    }

    signal[signal>1] <- 1
    signal[signal<0] <- 0
    signal <- signal * node.val

    return(signal)
}


prettyifelse <- function(test, una, olaotra){
    if(test){
        return(una)
    } else {
        return(olaotra)
    }
}



