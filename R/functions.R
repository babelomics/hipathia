
load.annot.file <- function(file, header=FALSE){
  data <- utils::read.table(file,
                            header = header,
                            stringsAsFactors = FALSE,
                            sep="\t",
                            quote="", 
                            comment.char="")
  data <- data[,1:2]
  colnames(data) <- c("gene","function")
  return(data)
}

#' @export
#'
annotate.paths <- function(pathigraphs, dbannot, entrez2hgnc){
  
  annofuns <- do.call("rbind", lapply(pathigraphs,function(pathigraph){
    
    new.pathigraph <- pathigraph
    vs <- V(new.pathigraph$graph)$name
    new.pathigraph$graph <- induced.subgraph(new.pathigraph$graph, 
                                             vs[!grepl("_func", vs)])
    funs <- get.pathway.functions(new.pathigraph, 
                                  dbannot, 
                                  entrez2hgnc, 
                                  use.last.nodes=TRUE, 
                                  unique=FALSE)
    paths <- lapply(names(funs), function(path){
      rep(path, times=length(funs[[path]]))
      })
    df <- data.frame(effector.nodes = unlist(paths), 
                     paths = gsub("N", "P", unlist(paths)), 
                     funs=unlist(funs), 
                     stringsAsFactors=FALSE)
    
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
#' @param normalize Boolean, whether to normalize the matrix of pathway
#' values with \code{normalize.paths} before quantifying the signal. Due to 
#' the nature of the Hipathia method, in which the length of each pathway may
#' alter its signal rank, we strongly recommend to perform this normalization.
#' This normalization removes the bias. Default is set to TRUE.
#' 
#' @return Matrix with the level of activation of the functions in 
#' \code{dbannot}
#' 
#' @examples 
#' data(brca_data)
#' trans.data <- translate.matrix(brca_data, "hsa")
#' exp.data <- normalize.data(trans.data)
#' pathways <- load.pathways(species = "hsa")
#' results <- hipathia(exp.data, pathways, verbose=FALSE)
#' go.values <- quantify.terms(results, pathways, "GO")
#' uniprot.values <- quantify.terms(results, pathways, "uniprot")
#' # annot.values <- quantify.terms(results, pathways, annotation)
#' 
#' @export
#'
quantify.terms <- function(results, metaginfo, dbannot, normalize = TRUE){
  
  method="mean"

  afu.file <- paste0("annofuns_uniprot_", metaginfo$species, ".RData")
  afg.file <- paste0("annofuns_GO_", metaginfo$species, ".RData")
  eh.file <- paste0("entrez_hgnc_", metaginfo$species, ".annot")
    
  if(!file.exists(system.file("extdata", afu.file, package="hipathia")) | 
     !file.exists(system.file("extdata", afg.file, package="hipathia")) |
     !file.exists(system.file("extdata", eh.file, package="hipathia")))
    retrieve.pathways.info(metaginfo$species)
  if(is.character(dbannot)){
    if(dbannot == "uniprot"){
      if(file.exists(system.file("extdata", afu.file, package="hipathia")) == TRUE){
        load(system.file("extdata", afu.file, package="hipathia"))
      }else{
        stop(paste0("File ", system.file("extdata",package="hipathia"), "/" , 
                    afu.file, " does not exist"))
      }
    }else if(dbannot == "GO"){
      if(file.exists(system.file("extdata", afg.file, package="hipathia")) == TRUE){
        load(system.file("extdata", afg.file, package="hipathia"))
      }else{
        stop(paste0("File ", system.file("extdata",package="hipathia"), "/", 
                    afg.file, " does not exist"))
      }
    }else{
      stop("dbannot parameter not accepted")
    }
  }else{
    entrez2hgnc <- utils::read.table(
      system.file("extdata", eh.file, package="hipathia"),
      header=FALSE,
      sep="\t",
      stringsAsFactors=FALSE)
    annofuns <- annotate.paths(metaginfo$pathigraphs, dbannot, entrez2hgnc)
  }

  annofuns <- annofuns[!is.na(annofuns$funs),]
  annofuns$pathway <- sapply(strsplit(annofuns$paths, split="-"), "[[", 2)
  annofuns <- annofuns[annofuns$pathway %in% names(metaginfo$pathigraphs),]
  fun.names <- unique(annofuns$funs)
  
  path.vals <- results$all$path.vals
  if(normalize == TRUE)
    path.vals <- normalize.paths(path.vals, metaginfo)
  
  fun.vals <- matrix(0, ncol=ncol(path.vals), 
                     nrow = length(fun.names), 
                     dimnames = list(fun.names, 
                                     colnames(path.vals)))
  cat(paste("Quantified functions:", nrow(fun.vals)))
  for( fun in fun.names){
    paths <- annofuns$paths[annofuns$funs==fun]
    if(method == "mean"){
      fun.vals[fun,] <- apply(path.vals[paths,,drop=FALSE],2,mean)
    }else if(method == "signal"){
      minimat <- 1 - path.vals[paths,,drop=FALSE]
      fun.vals[fun,] <- 1 - apply(minimat, 2, prod)
    }
  }
  return(fun.vals)
  
}

get.entrez.function <- function(entrezs,entrez2hgnc,dbannot){
  hgncs <- entrez2hgnc[entrez2hgnc[,1] %in% entrezs,2]
  path_functions <- dbannot[dbannot$gene%in% hgncs,]
  
  return(path_functions)
}


#' @importFrom stats fisher.test
#' @importFrom stats p.adjust
enrichment <- function(path_functions,dbannot,na.rm=TRUE){
  if(na.rm){
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
                  c(term_counts[x],n1-term_counts[x],
                    all_term_counts[x],n2-all_term_counts[x]),
                  nrow=2,
                  dimnames=list(c("path","rest"),c("yes","no")))
    test <- stats::fisher.test(mat,alternative = "greater")
    return(test$p.value)
  })
  
  pa <- stats::p.adjust(term_pvalues,"fdr")
  out <- data.frame(term=names(term_pvalues), 
                    pvalue=term_pvalues,
                    adj.pvalue=pa)
  out <- out[order(out$pvalue),]
  
  return(out)
  
}

#' @export
#'
get.pathway.functions <- function(pathigraph, dbannot, entrez2hgnc, use.last.nodes=TRUE, unique=TRUE){
  
  g <- pathigraph$graph
  last_nodes <- gsub("_func","",get.last.node(g))
  if(use.last.nodes==TRUE){
    ebn <- V(g)$genesList[which(V(g)$name %in% last_nodes)]
  } else {
    ebn <- get.vertex.attribute(g,"genesList")
  }
  entrezs <- unique(unlist(ebn))
  entrezs <- entrezs[!is.na(entrezs)]
  gpf <- unique(get.entrez.function(entrezs,entrez2hgnc,dbannot))
  if(nrow(gpf)==0){
    #warning(paste0("No functions found for any node in ",pathigraph$path.id))
    last_node_functions <- rep(NA,length(last_nodes))
    names(last_node_functions) <- last_nodes
    last_node_functions <- as.list(last_node_functions)
  } else {
    fe <- enrichment(gpf,dbannot)
    enriched <- as.character(fe$term[fe$pvalue<0.05])
    
    filtdbannot <- dbannot#[ dbannot[,2] %in% enriched,]
    filtdbannot$rank <- fe[match(filtdbannot[,2],fe$term),"pvalue"]
    filtdbannot$rank[is.na(filtdbannot$rank)] <- 1
    filtdbannot <- filtdbannot[order(filtdbannot$rank),]
    
    last_node_functions <- lapply(last_nodes,
                                  get.best.node.functions,
                                  pathigraph,
                                  entrez2hgnc,
                                  filtdbannot, 
                                  unique=unique)
    names(last_node_functions) <- last_nodes
    last_node_functions[sapply(last_node_functions,function(x) length(x))==0] <- NA
  }  
  return(last_node_functions)
  
}



get.best.node.functions <- function(node,pathigraph,entrez2hgnc,filtdbannot, unique=TRUE){
  g <- pathigraph$graph
  entrezs <- V(g)$genesList[[which(V(g)$name==node)]]
  node_funcs <- unique(unlist(lapply(entrezs,
                                     get.best.gene.functions,
                                     entrez2hgnc=entrez2hgnc,
                                     filtdbannot=filtdbannot, 
                                     unique=unique)))
  node_funcs <- node_funcs[!is.na(node_funcs)]
  return(node_funcs)
}

get.best.gene.functions <- function(entrez,entrez2hgnc,filtdbannot, unique=TRUE){
  hgncs <- entrez2hgnc[entrez2hgnc[,1] %in% entrez,2]
  funcs <- filtdbannot[filtdbannot$gene%in% hgncs,]
  funcs <- unique(funcs)
  if(unique == TRUE){
    return(funcs[1,2])
  }else{
    return(funcs[,2])
  }
}

get.last.node <- function(subgraph){
  V(subgraph)$name[!(V(subgraph)$name %in% get.edgelist(subgraph)[,1])]
}



