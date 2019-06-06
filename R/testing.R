##
## testing.R
## Testing functions for package Hipathia
##
## Written by Marta R. Hidalgo
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##


test_matrix <- function(mat){
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

test_pathways_object <- function(pathways){
    hasall <- length(pathways) == 7 | length(pathways) == 8
    if(length(pathways) == 7){
        byuser <- FALSE
    }else if(length(pathways) == 8){
        byuser <- pathways$by.user
    }
    spec <- is_accepted_species(pathways$species) || byuser == TRUE
    isigraph <- class(pathways$pathigraphs[[1]]$graph) == "igraph"
    if(!hasall == TRUE | !spec == TRUE | !isigraph == TRUE)
        stop("Pathways object not allowed")
}

test_tolerance <- function(tol){
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
#' #is_accepted_species("hsa")
#' #is_accepted_species("fca")
#'
#' @return Boolean, whether \code{species} is accepted or not.
#'
is_accepted_species <- function(species){
    isacc <- species %in% c("hsa", "mmu", "rno")
    return(isacc)
}

is_accepted_grouping <- function(group){
    isacc <- group %in% c("uniprot", "GO", "genes")
    return(isacc)
}

is_accepted_database <- function(db){
    isacc <- db %in% c("uniprot", "GO")
    return(isacc)
}
