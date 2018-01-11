##
## testing.R
## Testing functions for package Hipathia
##
## Written by Marta R. Hidalgo
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##


test.matrix <- function(mat){
    is01 <- all(mat >= 0 & mat <= 1)
    nona <- sum(is.na(mat)) == 0
    ismatrix <- class(mat) == "matrix"
    if(!is01 == TRUE)
        stop("Matrix rank exceeds from [0,1]")
    if(!nona == TRUE)
        stop("NAs found in matrix")
    if(!ismatrix == TRUE)
        stop("Input values must be of class 'matrix'")
}

test.pathways.object <- function(pathways){
    hasall <- length(pathways) == 6
    spec <- is.accepted(pathways$species)
    isigraph <- class(pathways$pathigraphs[[1]]$graph) == "igraph"
    if(!hasall == TRUE | !spec == TRUE | !isigraph == TRUE)
        stop("Pathways object not allowed")
}

test.tolerance <- function(tol){
    isless1 <- tol < 1
    ispositive <- tol > 0
    if(!ispositive == TRUE)
        stop("Tolerance parameter must be positive")
    if(!isless1 == TRUE)
        stop("Tolerance parameter must be lower than 1")
}


#' Checks whether a species is accepted
#'
#' @param species Species of the samples.
#'
#' #@examples
#' #is.accepted("hsa")
#' #is.accepted("fca")
#'
#' @return Boolean, whether \code{species} is accepted or not.
#'
is.accepted <- function(species){
    isacc <- species %in% c("hsa", "mmu", "rno")
    return(isacc)
}

is.accepted.grouping <- function(group){
    isacc <- group %in% c("uniprot", "GO", "genes")
    return(isacc)
}

