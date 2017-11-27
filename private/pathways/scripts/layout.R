
#_____________________________________________________
#
#              LAYOUTS
#_____________________________________________________


refine.layout <- function(gg,coords,fixed,e=10e-7,maxiter=500,w=0.5,align.final.at.right=T,p0q=0.01,p0mult=5,use.edges=T){
  
  cs <- coords
  
  V(gg)$x <- cs[,1]
  V(gg)$y <- cs[,2]
  
  minx <- min(cs[fixed==T,1])
  maxx <- max(cs[fixed==T,1])
  miny <- min(cs[fixed==T,2])
  maxy <- max(cs[fixed==T,2])
  
  #cs[fixed==F,1] <- sample(seq(minx,maxx,length.out=1000),1) #minx
  #cs[fixed==F,2] <- sample(seq(miny,maxy,length.out=1000),1) #miny
  
  n <- length(V(gg))
  #maxdis <- sqrt((maxx-minx)^2+(maxy-miny)^2)
  maxdis <- max(abs(maxx-minx),4*abs(maxy-miny))
  
  # get neighbour
  global_dis <- shortest.paths(gg)
  global_dis[is.infinite(global_dis)] <- max(global_dis[!is.infinite(global_dis)]) + 1
  nglobal_dis <- (max(global_dis)-global_dis)*(1-diag(nrow(global_dis)))
  nglobal_dis <- nglobal_dis/max(nglobal_dis)  
  adjm <- (global_dis==1)+0
  dists <- c()
  for(i in 1:nrow(adjm)){
    anodes <- which(adjm[i,]==1)
    for(j in anodes){
      dists <- c(dists,sqrt((cs[i,1]-cs[j,1])^2 + (cs[i,2]-cs[j,2])^2))
    }
  }  
  p0 <- stats::quantile(dists[dists>0],p0q)
  p00 <- p0*p0mult
  
  if(align.final.at.right==T){
    nonfixed <- intersect(which(fixed==F),which(!(V(gg)$name %in% get.edgelist(gg)[,1])))    
    for(nf in nonfixed){
      friends <- which(global_dis[nf,]==1)
      x <- max(V(gg)$x[friends])+p0*3
      y <- mean(V(gg)$y[friends])
      cs[nf,] <- c(x,y)
    }
  }
  
  # Move also nodes with duplicated coordinates
  fixed[duplicated(cs)] <- F  
  
  iter <- 1
  medif <- e+1
  locations <- list()
  
  while(medif>e & iter<maxiter){
    
    cs2 <- cs 
    
    for(i in 1:n){
      
      if(fixed[i]==F){
        
        # node forces
        node_forces <- mat.or.vec(nr=n,nc=2)
        
        for(j in 1:n){
          
          if(i!=j){
            
            is_only_repeller <- !(nglobal_dis[i,j]==1)
            pp <- c(p00,p0)[is_only_repeller+1]
            node_forces[j,] <- get.node.vectorial.force(cs[i,1],cs[i,2],cs[j,1],cs[j,2],pp,maxdis,repeller=is_only_repeller)
            if(global_dis[i,j]==1){
              node_forces[j,] <- node_forces[j,]*10
            }
            
          }
        }
        
        #print(".")
        
        # edge forces
        if(use.edges==T){
          edges <- get.edgelist(gg)
          ne <- nrow(edges)
          edge_forces <- mat.or.vec(nr=ne,nc=2)
          for(j in 1:nrow(edges)){
            if(V(gg)$name[i]!=edges[j,1] & V(gg)$name[i]!=edges[j,2]){
              c0 <- cs[i,]
              c1 <- cs[which(V(gg)$name==edges[j,1]),]
              c2 <- cs[which(V(gg)$name==edges[j,2]),]
              edge_forces[j,] <- get.edge.vectorial.force(c0,c1,c2,p0,maxdis)
            }          
          }
          
          forces <- rbind(node_forces,edge_forces)
        } else {
          forces <- node_forces
        }
        cs2[i,] <- cs2[i,] + colMeans(forces)*w
        
        #print("*")
        
        #if(i==61) stop("pepe")
        
      }
      
    }
    
    
    medif <- max(abs(cs-cs2))
    
    cs <- cs2
    locations[[iter]] <- cs
    
    iter <- iter + 1
    # cat(iter,medif,"\n")
  }
  
  return(list(    
    layout=cs,
    locations=locations,
    iter=iter,
    maxiter=maxiter,
    w=w,
    e=e,
    laste=medif,
    p0=p0,
    p00=p00,
    maxdis=maxdis,
    dists=dists,
    fixed=fixed
  ))
  
}

plot.node.trip <- function(locations,fixed){
  nodes <- which(fixed==F)
  for(node in nodes){
    lines(t(sapply(locations,function(x) c(x[node,1],x[node,2]))),col="red",lty=2)
  }
}

get.node.vectorial.force <- function(x1,y1,x2,y2,p0,maxdis,repeller=F){
  
  dis <- sqrt((x1-x2)^2+(y1-y2)^2)
  dis <- max(abs(x1-x2),4*abs(y1-y2))
  
  if(repeller==T){
    if(dis>p0){
      force <- 0
    } else {
      ndis <- (p0-dis)/p0
      force <- -ndis
    }
  } else {
    if(dis>p0){
      ndis <- (dis-p0)/(maxdis-p0)      
      force <- ndis
    } else {    
      ndis <- (p0-dis)/p0      
      force <- -ndis
    }
  }
  
  v <- c(x2-x1,y2-y1)
  if(v[1]==0 & v[2]==0){
    v <- rnorm(2)
  }
  v <- v/sqrt(v[1]^2+v[2]^2)
  v <- v*force
  
  return(v)
  
}

get.edge.vectorial.force <- function(c0,c1,c2,p0,maxdis){
  
  x0 <- c0[1]
  y0 <- c0[2]
  x1 <- c1[1]
  y1 <- c1[2]
  x2 <- c2[1]
  y2 <- c2[2]
  
  mod <- sqrt( (x2-x1)^2 + (y2-y1)^2 )
  s <- ( (x0-x1)*(x2-x1) + (y0-y1)*(y2-y1)) / mod
  sx <- x1 + (s/mod)*(x2-x1)
  sy <- y1 + (s/mod)*(y2-y1)
  
  if(x1<x2) {
    xmin <- x1
    xmax <- x2
  } else {
    xmin <- x2
    xmax <- x1
  }
  if(y1<y2) {
    ymin <- y1
    ymax <- y2
  } else {
    ymin <- y2
    ymax <- y1
  }
  if(sx<xmin | sx>xmax | sy<ymin | sy>ymax){
    force <- 0
  } else {
    dis <- sqrt((x0-sx)^2+(y0-sy)^2)
    if(dis>p0){
      force <- 0
    } else {
      ndis <- (p0-dis)/p0
      force <- -ndis
    }
  }  
  
  v <- c(sx-x0,sy-y0)
  if(v[1]==0 & v[2]==0){
    v <- rnorm(2)
  }
  v <- v/sqrt(v[1]^2+v[2]^2)
  v <- v*force
  
  return(v)
  
}




#_____________________________________________________
#
#          INTEGRATE LAYOUT AND FUNCTIONS
#_____________________________________________________


add.param <- function(param,graph,newgraph,common,init=T){
  if(init==T){
    current <- get.vertex.attribute(graph,param)
    if(is.numeric(current)){
      newgraph <- set.vertex.attribute(newgraph,param,value=mean(current))
    } else {
      newgraph <- set.vertex.attribute(newgraph,param,value=names(sort(table(current),decreasing=T))[1])
    }  
  }
  newgraph <- set.vertex.attribute(newgraph,param,index = common, get.vertex.attribute(graph,param,index = common))  
}


add.functions.to.pathigraph <- function(pathigraph,entrez2hgnc,dbannot,maxiter=500,verbose=T,w=10,e=0.001,use.last.nodes=T,use.edges=F,p0q=0.01,p0mult=5){
  
  newpathigraph <- pathigraph
  
  last_node_funcs_list <- get.pathway.functions(pathigraph,dbannot,entrez2hgnc,use.last.nodes=use.last.nodes)
  if(sum(is.na(last_node_funcs_list))>0) last_node_funcs_list <- last_node_funcs_list[!is.na(last_node_funcs_list)]
  last_node_funcs <- sapply(last_node_funcs_list,function(x) paste(x,collapse="\n"))  
  nff <- data.frame(node=names(last_node_funcs),functions=last_node_funcs, stringsAsFactors=F)
  
  V(pathigraph$graph)$x <- V(pathigraph$graph)$nodeX
  V(pathigraph$graph)$y <- V(pathigraph$graph)$nodeY
  
  if(nrow(nff)==0){  
    return(pathigraph)
  } else {
    # create new graph
    if(verbose==T) cat("Creating new graph...\n")
    newgraph <- graph.data.frame(rbind(get.edgelist(pathigraph$graph),cbind(nff$node,paste0(nff$node,"_func"))))
    common <- V(pathigraph$graph)$name[V(pathigraph$graph)$name %in% V(newgraph)$name]      
    newgraph <- add.param("x",pathigraph$graph,newgraph,common)
    newgraph <- add.param("y",pathigraph$graph,newgraph,common)
    newgraph <- add.param("shape",pathigraph$graph,newgraph,common)
    newgraph <- add.param("size",pathigraph$graph,newgraph,common)
    newgraph <- add.param("size2",pathigraph$graph,newgraph,common)
    newgraph <- add.param("label",pathigraph$graph,newgraph,common,init=F)
    E(newgraph)$arrow.size <- .2
    V(newgraph)$label.cex <- .62
    newgraph <- set.vertex.attribute(newgraph,"label",index = paste0(nff$node,"_func"),value = nff$functions)
    newgraph <- set.vertex.attribute(newgraph,"color",value="cyan")
    
    # recalculate layout
    if(verbose==T) cat("Recomputing graph layout...\n")
    fixed <- rep(T,length(V(newgraph)))
    fixed[grep("func",V(newgraph)$name)] <- F  
    lay <- cbind(V(newgraph)$x,V(newgraph)$y)
    rl <- refine.layout(newgraph,lay,fixed,w=w,maxiter = maxiter,e = e, use.edges=use.edges,p0q=p0q,p0mult=p0mult)
    rl$iter
    
    V(newgraph)$nodeX <- rl$layout[,1]
    V(newgraph)$x <- rl$layout[,1]
    V(newgraph)$nodeY <- rl$layout[,2]
    V(newgraph)$y <- rl$layout[,2]
    
    gl <- as.list(V(pathigraph$graph)$genesList)
    names(gl) <- V(pathigraph$graph)$name    
    V(newgraph)$genesList <- gl[V(newgraph)$name]
    
    negl <- length(E(newgraph))-length(E(pathigraph$graph))    
    E(newgraph)$relation <- c(E(pathigraph$graph)$relation,rep(1,negl))
    
    newpathigraph$graph <- newgraph
    newpathigraph$subgraphs_funs <- create.subgraphs(newgraph)[[1]]
    newpathigraph$effector.subgraphs_funs <- create.effector.subgraphs(newpathigraph, funs=T)
    newpathigraph$subgraphs <- pathigraph$subgraphs
    newpathigraph$rl <- rl
    newpathigraph$fixed <- rl$fixed
    
    return(newpathigraph)
  }
  
}


add.functions.to.pathigraphs <- function(apgs, entrez2hgnc, dbannot, maxiter= 100, p0mult= 4){
  cat("Adding functions to pathways...\n")
  fpgs <- lapply(apgs,function(x) {
    cat(x$path.id)
    cat(" - ")
    cat(x$path.name) 
    cat("\n")
    add.functions.to.pathigraph(x,entrez2hgnc,dbannot, maxiter=maxiter, p0mult=p0mult, verbose = FALSE)
  })
  return(fpgs)
}
