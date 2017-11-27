
#_____________________________________________________
#
#              SIFS
#_____________________________________________________
#

load.graphs <- function(input.folder,species,pathway.names=NULL){
  
  nam <- read.delim(file=paste(input.folder, "/name.pathways_", species, ".txt", sep=""), comment.char="", sep="\t", header=F, stringsAsFactors=F, row.names=1, colClasses="character")
  
  if(is.null(pathway.names)){
    #pathway.names <- unique(gsub(".att", "", gsub(".sif","",list.files(input.folder,pattern=species))))
    pathway.names <- paste(species,rownames(nam),sep="")
  } 
  
  pathigraphs <- list()
  for(pathway in pathway.names){
    print(pathway)
    pathigraphs[[pathway]] <- list()
    pathigraphs[[pathway]]$graph <- sif2patIGraph(paste(input.folder, "/", pathway, sep=""))
    subs <- create.subgraphs(pathigraphs[[pathway]]$graph)
    pathigraphs[[pathway]]$subgraphs <- subs[[1]]
    pathigraphs[[pathway]]$subgraphs.mean.length <- subs[[2]]
    pathigraphs[[pathway]]$effector.subgraphs <- create.effector.subgraphs(pathigraphs[[pathway]])
    pathigraphs[[pathway]]$path.name <- nam[gsub(species, "", pathway),1]
    pathigraphs[[pathway]]$path.id <- pathway
    pathigraphs[[pathway]]$label.id <- cbind(V(pathigraphs[[pathway]]$graph)$name, V(pathigraphs[[pathway]]$graph)$label)
    colnames(pathigraphs[[pathway]]$label.id) <- c("name", "label")
  }
  
  return(pathigraphs)
  
}


sif2patIGraph <- function(pref.file){  
  # Read files
  sif <- utils::read.table(paste0(pref.file, ".sif"), sep="\t", fill=T, stringsAsFactors=F, header=F)
  att <- utils::read.table(paste0(pref.file, ".att"), sep="\t", fill=T, stringsAsFactors=F, header=T)
  rownames(att) <- att[,1]
  
  # Create igraph
  ig <- graph.data.frame(sif[,c(1,3)], directed=T)
  
  # Add attributes
  E(ig)$relation <- unlist(sapply(sif[,2], function(x){ifelse(x=="activation", 1, -1)}))
  E(ig)$curved <- F
  V(ig)$nodeX <- att[V(ig)$name,"X"]
  V(ig)$nodeY <- att[V(ig)$name,"Y"] # max(att$Y) - att[V(ig)$name,"Y"] # 
  V(ig)$shape <- att[V(ig)$name,"shape"]
  V(ig)$label.cex <- att[V(ig)$name,"label.cex"]
  V(ig)$label.color <- att[V(ig)$name,"label.color"]
  V(ig)$width <- att[V(ig)$name,"width"]
  V(ig)$height <- att[V(ig)$name,"height"]
  glist <- sapply(as.character(att[,"genesList"]), function(x){unlist(strsplit(x, split=","))})
  names(glist) <- att[,"ID"]
  V(ig)$genesList <- glist[V(ig)$name]
  if(!is.null(att$label))
    V(ig)$label <- att[V(ig)$name,"label"]
  
  return(ig)
}


# Write SIF files from KGML files
transform.XML.to.SIF <- function(pathway.names, kgml.folder, sif.folder, metabolite = T){
  g <- sapply(pathway.names, function(name.pathway){
    print(name.pathway)
    ig <- parse.kegg.to.igraph( name.pathway, kgml.folder, genesOnly= !metabolite )
    write.sif.files(ig, rep("white", length(V(ig))), paste0(sif.folder, name.pathway))
  })
}


write.sif.files <- function(ig, color, exit.file){
  write.sif.file(ig, paste(exit.file, ".sif", sep=""))
  write.nodes.file(ig, color, paste(exit.file, ".att", sep=""))
}

# Creates a .SIF file from the graph "rendergraph" with "NODE1 - REL- NODE2" format
#      rgraph: pathigraph with render info (output from plotPathigraph)
#      exit.file: root where the file will be stored
write.sif.file <- function(ig, exit.file){
  sif <- get.edgelist(ig)
  sif <- cbind(sif, sif[,2])
  sif[,2] <- unlist(lapply(E(ig)$relation, function(x){ifelse(x==1, "activation", "inhibition")}))
  write.table(sif, file=exit.file, row.names=F, col.names=F, quote=F, sep="\t")
}


# Creates a .TXT file with the attributes of the nodes. Included info: X and Y position, color, letter color and shape.
#      rgraph: pathigraph with render info (output from plotPathigraph)
#      exit.file: root where the file will be stored
write.nodes.file <- function(ig, color, exit.file){
  
  attrs <- cbind(V(ig)$name, gsub("\n"," / ",V(ig)$label), V(ig)$nodeX, V(ig)$nodeY, color, V(ig)$shape, V(ig)$type, V(ig)$label.cex, V(ig)$label.color, V(ig)$width, V(ig)$height, sapply(V(ig)$genesList, paste, collapse=","))
  write.table(t(c("ID","label","X", "Y", "color", "shape", "type", "label.cex", "label.color", "width", "height", "genesList")), file=exit.file, row.names=F, col.names=F, quote=F, sep="\t")
  write.table(attrs, append=TRUE, file=exit.file, row.names=F, col.names=F, quote=F, sep="\t")
}





#_____________________________________________________
#
#              PATHIGRAPHS
#_____________________________________________________


create.subgraphs <- function(ig){
  # Compute subpaths
  in.nodes <- V(ig)$name[!V(ig)$name%in%get.edgelist(ig)[,2]]
  out.nodes <- V(ig)$name[!V(ig)$name%in%get.edgelist(ig)[,1]]
  possible.paths <- find.possible.paths(ig, in.nodes, out.nodes)
  subgraphs <- NULL
  sublen <- NULL
  for( ininode in names(possible.paths)){
    subs <- lapply(possible.paths[[ininode]], function(x){induced.subgraph(ig, unique(unlist(x)))})
    pathway.name <- unlist(strsplit(ininode, split="\\-"))[2]
    ininode.id <- unlist(strsplit(ininode, split="\\-"))[3]
    endnodes.id <- sapply(names(subs), function(name) unlist(strsplit(name, split="\\-"))[3])
    names(subs) <- paste("P", pathway.name, ininode.id, endnodes.id, sep="-")
    subgraphs <- c(subgraphs, subs)
    lens <- lapply(possible.paths[[ininode]], function(x){mean(unlist(lapply(x, length)))})
    names(lens) <- names(subs)
    sublen <- c(sublen, lens)
  }
  return(list(subgraphs, sublen))
}



create.effector.subgraphs <- function(pathigraph, funs=F){
  if(funs == T){
    pathisubs <- pathigraph$subgraphs_funs
  }else{
    pathisubs <- pathigraph$subgraphs
  }
  pathisubs.effectors <- sapply(names(pathisubs), function(name) paste(unlist(strsplit(name, split="\\-"))[c(1, 2, 4)], collapse = "-"))
  effectors <- unique(pathisubs.effectors)
  effector.subs <- list()
  for(effector in effectors){
    subs <- pathisubs[effector == pathisubs.effectors]
    effector.subs[[effector]] <- induced.subgraph(pathigraph$graph, V(graph.union(subs))$name)
  }
  return(effector.subs)
}


find.possible.paths <- function( ig, in.nodes, out.nodes ){
  paths <- list()
  for( in.node in in.nodes ){
    all.pathsi <- find.all.paths.from( in.node, ig, NULL )
    endnodes <- sapply(all.pathsi, function(path){path[length(path)]})
    cycles <- all.pathsi[!endnodes %in% out.nodes]
    for(cycle in cycles){
      cycle.end <- cycle[length(cycle)]
      contain <- which(sapply(all.pathsi, function(p) cycle.end %in% p))
      for(i in contain)
        all.pathsi[[i]] <- unique(unlist(c(all.pathsi[[i]], cycle)))
    }
    uniend <- unique(endnodes)
    uniend <- uniend[uniend %in% out.nodes]
    if(length(uniend)>0){
      paths[[in.node]] <- lapply(uniend, function(end){all.pathsi[endnodes==end]})
      names(paths[[in.node]]) <- uniend
    }
  }
  return(paths)
}



# Returns a list containing all paths from the given "node" to any final node, in list format
find.all.paths.from <- function( node, ig, visited ){
  paths <- list()
  edges <- incident(ig, node, mode="out")
  edges <- setdiff(edges, visited)
  for( ed in edges ){
    visited <- c(visited, ed)
    # Find next node (opposite)
    op <- get.edgelist(ig)[ed,2]
    # IF op is a final node
    if( length(incident(ig, op, mode="out")) == 0 ){
      path <- c(node, op)
      paths[[length(paths)+1]] <- path
    }else{
      new.paths <- find.all.paths.from(op, ig, visited)
      for( j in 1:length( new.paths)){
        new.paths[[j]] <- c( node, new.paths[[j]] )
        paths[[length(paths)+1]] <- new.paths[[j]]
      }
    }
  }
  if( length(paths)== 0)
    paths[[1]] <- node
  return( paths )
}



#_____________________________________________________
#
#              OTHERS
#_____________________________________________________


create.pathway.names.file <- function(pathway.names, species, kgml.folder, sif.folder){
  titles <- sapply( pathway.names, function(name){
    kg <- parseKGML(paste0(kgml.folder, "/", name, ".xml"))
    kg@pathwayInfo@title
  })
  titles.df <- cbind(gsub(species, "", names(titles)), titles)
  write.table(titles.df, file=paste0(sif.folder, "/name.pathways_", species, ".txt"), sep="\t", quote=F, col.names=F, row.names=F)
}



#_____________________________________________________
#
#              AMMENDMENTS
#_____________________________________________________


amend.kegg.pathways <- function(amend.file, pathigraphs, species, verbose=T){
  if(!file.exists(amend.file)){
    warning("amend.file does not exist")
  }else{
    amend.df <- utils::read.table(amend.file, header=T, sep="\t", stringsAsFactors=F, colClasses = "character")
    amend.df$path <- paste0(species, amend.df$path)
    amend.df$node1 <- apply(amend.df, 1, function(v)paste("N", v[1], v[3], sep="-"))
    amend.df$node2 <- apply(amend.df, 1, function(v)paste("N", v[1], v[5], sep="-"))
    tochange <- intersect(names(pathigraphs), unique(amend.df$path))
    if(verbose) cat("Amending pathway...\n")
    for(path in tochange){
      if(verbose)
        cat(paste(path, "-", pathigraphs[[path]]$path.name, "\n"))
      mini.df <- amend.df[amend.df$path==path,]
      for(i in 1:nrow(mini.df)){
        if(mini.df$action[i] == "+"){
          pathigraphs[[path]] <- amend.add.edge(pathigraphs[[path]], mini.df[i,])
        }else if(mini.df$action[i] == "-"){
          pathigraphs[[path]] <- amend.delete.edge(pathigraphs[[path]], mini.df[i,])
        }
      }
      todelete <- V(pathigraphs[[path]]$graph)$name[!V(pathigraphs[[path]]$graph)$name%in% unique(as.vector(get.edgelist(pathigraphs[[path]]$graph)))]
      for(node in todelete)
        pathigraphs[[path]] <- amend.delete.node(pathigraphs[[path]], node)
      subs <- create.subgraphs(pathigraphs[[path]]$graph)
      pathigraphs[[path]]$subgraphs <- subs[[1]]
      pathigraphs[[path]]$subgraphs.mean.length <- subs[[2]]
      pathigraphs[[path]]$effector.subgraphs <- create.effector.subgraphs(pathigraphs[[path]])
      pathigraphs[[path]]$label.id <- cbind(V(pathigraphs[[path]]$graph)$name, V(pathigraphs[[path]]$graph)$label)
      colnames(pathigraphs[[path]]$label.id) <- c("name", "label")    
    }
  }
  return(pathigraphs)
}


amend.add.edge <- function(pathigraph, amend){
  if(amend$node1 %in% V(pathigraph$graph)$name & amend$node2 %in% V(pathigraph$graph)$name & sum(duplicated(rbind(as.character(amend[,c(3,5)]), get.edgelist(pathigraph$graph))))<1)
    pathigraph$graph <- pathigraph$graph + edge(amend$node1, amend$node2, curved=F, relation=ifelse(amend$relation=="activation", 1, -1))
  return(pathigraph)
}

amend.delete.edge <- function(pathigraph, amend){
  if(amend$node1 %in% V(pathigraph$graph)$name & amend$node2 %in% V(pathigraph$graph)$name & sum(duplicated(rbind(as.character(amend[,c(3,5)]), get.edgelist(pathigraph$graph))))==1)
    pathigraph$graph <- pathigraph$graph - edge(paste0(amend$node1, "|", amend$node2))
  return(pathigraph)
}

amend.delete.node <- function(pathigraph, node.name){
  pathigraph$graph <- pathigraph$graph - node.name
  return(pathigraph)
}



#_____________________________________________________
#
#              CREATE METAGINFO
#_____________________________________________________

create.metaginfo.object <- function(fpgs, basal.value = 0.5){
  
  pathigraph.genes <- all.needed.genes(fpgs)
  
  # Todos los genes a valor 0.5
  genes.vals.05 <- matrix(basal.value, ncol=2, nrow=length(pathigraph.genes), dimnames=list(pathigraph.genes))
  meta.05 <- NULL
  meta.05$pathigraphs <- fpgs
  results.05 <- hipathia(genes.vals.05, meta.05)
  results.dec.05 <- hipathia(genes.vals.05, meta.05, decompose = TRUE)

  # Create metaginfo object
  metaginfo <- NULL
  metaginfo[["all.genes"]] <- pathigraph.genes
  metaginfo[["path.norm"]] <- results.dec.05$all$path.vals[,1]
  metaginfo[["eff.norm"]] <- results.05$all$path.vals[,1]
  metaginfo[["pathigraphs"]] <- fpgs
  
  return(metaginfo)
}
