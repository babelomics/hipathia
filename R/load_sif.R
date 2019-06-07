

#' Create a Pathways object from SIF files
#'
#' Creates a Pathways object from the information of a pathway stored in a SIF
#' file with some attributes. This pathways object can be used by function
#' \code{hipathia} to analyze data.
#'
#' @param sif.folder Path to the folder in which SIF and ATT files are stored.
#' @param spe Species
#' @param entrez_symbol Relation between Entrez (NCBI) genes and gene symbols.
#' Data.frame with 2 columns: First column is the EntrezGene ID, second column
#' is the gene Symbol. The genes in the nodes of the pathways should be defined
#' by Entrez IDs in the SIF and ATT files of the pathways. In order to be more
#' readable, gene names are used when plotting the pathways.
#' @param dbannot Functional annotation of the genes in the pathways to create
#' function nodes.
#'
#' @return A pathways object with the same structure of that returned by
#' function \code{load_pathways}.
#'
#' @export
#'
mgi_from_sif <- function(sif.folder, spe, entrez_symbol = NULL, dbannot = NULL){

    ## NO USAR ESTA FUNCIÓN PARA HACER LOS PATHWAYS DE KEGG: NO HACE AMMENDMENTS

    message("Loading graphs...")
    pgs <- load_graphs(sif.folder, spe)
    if(!is.null(dbannot) & !is.null(entrez_symbol)){
        message("Adding functions to pathways...")
        pgs <- add_funs_to_pgs(pgs, entrez_symbol, dbannot, maxiter = 1000)
    }
    message("Creating MGI...")
    metaginfo <- create_metaginfo_object(pgs, spe, by.user = TRUE)

    message("Created MGI with ", length(metaginfo$pathigraphs), " pathway(s)")

    return(metaginfo)
}



load_graphs <- function(input.folder, species, pathway.names = NULL,
                        verbose = FALSE){

    file <- paste(input.folder, "/name.pathways_", species, ".txt", sep="")
    nam <- utils::read.delim(file, comment.char = "", sep = "\t",
                             header = FALSE, stringsAsFactors = FALSE,
                             row.names = 1, colClasses = "character")

    if(is.null(pathway.names))
        pathway.names <- paste0(species, rownames(nam))

    pathigraphs <- list()
    for(pathway in pathway.names){
        if(verbose == TRUE)
            print(pathway)
        pathigraphs[[pathway]] <- list()
        s2p <- sif_to_pathigraph(paste(input.folder, "/", pathway, sep=""))
        pathigraphs[[pathway]]$graph <- s2p
        subs <- create_subgraphs(pathigraphs[[pathway]]$graph)
        pathigraphs[[pathway]]$subgraphs <- subs[[1]]
        pathigraphs[[pathway]]$subgraphs.mean.length <- subs[[2]]
        ces <- create_effector_subgraphs(pathigraphs[[pathway]])
        pathigraphs[[pathway]]$effector.subgraphs <- ces
        pathigraphs[[pathway]]$path.name <- nam[gsub(species, "", pathway),1]
        pathigraphs[[pathway]]$path.id <- pathway
        labs <- cbind(V(pathigraphs[[pathway]]$graph)$name,
                      V(pathigraphs[[pathway]]$graph)$label)
        pathigraphs[[pathway]]$label.id <- labs
        colnames(pathigraphs[[pathway]]$label.id) <- c("name", "label")
    }

    return(pathigraphs)

}



sif_to_pathigraph <- function(pref.file){
    # Read files
    sif <- utils::read.table(paste0(pref.file, ".sif"), sep = "\t", fill = TRUE,
                             stringsAsFactors = FALSE, header = FALSE)
    att <- utils::read.table(paste0(pref.file, ".att"), sep = "\t", fill = TRUE,
                             stringsAsFactors = FALSE, header = TRUE)
    rownames(att) <- att[,1]

    # Create igraph
    ig <- graph.data.frame(sif[,c(1,3)], directed = TRUE)

    # Add attributes
    E(ig)$relation <- unlist(sapply(sif[,2], function(x){
        ifelse(x == "activation", 1, -1)
    }))
    E(ig)$curved <- FALSE
    V(ig)$nodeX <- att[V(ig)$name, "X"]
    V(ig)$nodeY <- att[V(ig)$name, "Y"] # max(att$Y) - att[V(ig)$name,"Y"] #
    V(ig)$shape <- att[V(ig)$name, "shape"]
    V(ig)$label.cex <- att[V(ig)$name, "label.cex"]
    V(ig)$label.color <- att[V(ig)$name, "label.color"]
    V(ig)$width <- att[V(ig)$name, "width"]
    V(ig)$height <- att[V(ig)$name, "height"]
    V(ig)$tooltip <- att[V(ig)$name, "tooltip"]
    glist <- sapply(as.character(att[, "genesList"]), function(x){
        unlist(strsplit(x, split = ","))
    })
    names(glist) <- att[,"ID"]
    V(ig)$genesList <- glist[V(ig)$name]
    if(!is.null(att$label))
        V(ig)$label <- att[V(ig)$name,"label"]

    return(ig)
}




#_____________________________________________________
#
#              PATHIGRAPHS
#_____________________________________________________


create_subgraphs <- function(ig){
    # Compute subpaths
    in.nodes <- V(ig)$name[!V(ig)$name%in%get.edgelist(ig)[,2]]
    out.nodes <- V(ig)$name[!V(ig)$name%in%get.edgelist(ig)[,1]]
    possible.paths <- find_possible_paths(ig, in.nodes, out.nodes)
    subgraphs <- NULL
    sublen <- NULL
    for(ininode in names(possible.paths)){
        funinduced <- function(x){induced.subgraph(ig, unique(unlist(x)))}
        subs <- lapply(possible.paths[[ininode]], funinduced)
        pathname <- unlist(strsplit(ininode, split="\\-"))[2]
        ininode.id <- unlist(strsplit(ininode, split="\\-"))[3]
        endnodes.id <- sapply(strsplit(names(subs), split="\\-"), "[[", 3)
        names(subs) <- paste("P", pathname, ininode.id, endnodes.id, sep="-")
        subgraphs <- c(subgraphs, subs)
        funmeanlength <- function(x){mean(unlist(lapply(x, length)))}
        lens <- lapply(possible.paths[[ininode]], funmeanlength)
        names(lens) <- names(subs)
        sublen <- c(sublen, lens)
    }
    return(list(subgraphs, sublen))
}



create_effector_subgraphs <- function(pathigraph, funs=FALSE){
    if(funs == TRUE){
        pathisubs <- pathigraph$subgraphs_funs
    }else{
        pathisubs <- pathigraph$subgraphs
    }
    pathisubs.effectors <- sapply(names(pathisubs),
                                  function(name) paste(unlist(
                                      strsplit(name, split="\\-"))[c(1, 2, 4)],
                                      collapse = "-"))
    effectors <- unique(pathisubs.effectors)
    effector.subs <- list()
    for(effector in effectors){
        subs <- pathisubs[effector == pathisubs.effectors]
        effector.subs[[effector]] <- induced.subgraph(pathigraph$graph,
                                                      V(graph.union(subs))$name)
    }
    return(effector.subs)
}


find_possible_paths <- function( ig, in.nodes, out.nodes ){
    paths <- list()
    for( in.node in in.nodes ){
        all.pathsi <- find_all_paths_from( in.node, ig, NULL )
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
        if(length(uniend) > 0){
            fununiend <- function(end){all.pathsi[endnodes == end]}
            paths[[in.node]] <- lapply(uniend, fununiend)
            names(paths[[in.node]]) <- uniend
        }
    }
    return(paths)
}



# Returns a list containing all paths from the given "node" to any final node,
# in list format
find_all_paths_from <- function(node, ig, visited){
    paths <- list()
    edges <- incident(ig, node, mode="out")
    edges <- setdiff(edges, visited)
    for(ed in edges){
        visited <- c(visited, ed)
        # Find next node (opposite)
        op <- get.edgelist(ig)[ed,2]
        # IF op is a final node
        if(length(incident(ig, op, mode="out")) == 0){
            path <- c(node, op)
            paths[[length(paths)+1]] <- path
        }else{
            new.paths <- find_all_paths_from(op, ig, visited)
            for(j in seq_along(new.paths)){
                new.paths[[j]] <- c(node, new.paths[[j]])
                paths[[length(paths)+1]] <- new.paths[[j]]
            }
        }
    }
    if(length(paths) == 0)
        paths[[1]] <- node
    return( paths )
}





#_____________________________________________________
#
#          INTEGRATE LAYOUT AND FUNCTIONS
#_____________________________________________________


add_param <- function(param, graph, newgraph, common, init = TRUE){
    if(init == TRUE){
        current <- get.vertex.attribute(graph, param)
        if(is.numeric(current)){
            newgraph <- set.vertex.attribute(newgraph, param,
                                             value = mean(current))
        } else {
            newgraph <- set.vertex.attribute(newgraph, param,
                                             value = names(sort(
                                                 table(current),
                                                 decreasing = TRUE))[1])
        }
    }
    newgraph <- set.vertex.attribute(newgraph, param, index = common,
                                     get.vertex.attribute(graph, param,
                                                          index = common))
}


add_funs_to_pg <- function(pathigraph, entrez_symbol, dbannot, maxiter = 500,
                           verbose = TRUE, w = 10, e = 0.001,
                           use.last.nodes = TRUE, use.edges = FALSE, p0q = 0.01,
                           p0mult = 5){

    newpathigraph <- pathigraph

    V(pathigraph$graph)$x <- V(pathigraph$graph)$nodeX
    V(pathigraph$graph)$y <- V(pathigraph$graph)$nodeY

    lastNodeFuns_list <- get_pathway_functions(pathigraph,
                                               dbannot,
                                               entrez_symbol,
                                               use_last_nodes = use.last.nodes)
    if(sum(is.na(lastNodeFuns_list)) > 0)
        lastNodeFuns_list <- lastNodeFuns_list[!is.na(lastNodeFuns_list)]
    if(length(lastNodeFuns_list) == 0){
        return(pathigraph)
    } else {
        lastNodeFuns <- sapply(lastNodeFuns_list,
                                  function(x) paste(x, collapse = "\n"))
        nff <- data.frame(node = names(lastNodeFuns),
                          node_funs = paste0(names(lastNodeFuns), "_func"),
                          functions = lastNodeFuns,
                          stringsAsFactors = FALSE)

        # create new graph
        if(verbose == TRUE)
            cat("Creating new graph...\n")
        newgraph <- graph.data.frame(rbind(get.edgelist(pathigraph$graph),
                                           cbind(nff$node,
                                                 paste0(nff$node, "_func"))))
        common <- V(pathigraph$graph)$name[V(pathigraph$graph)$name %in%
                                               V(newgraph)$name]
        V(pathigraph$graph)$size <- 15
        V(pathigraph$graph)$size2 <- 5
        newgraph <- add_param("x", pathigraph$graph, newgraph, common)
        newgraph <- add_param("y", pathigraph$graph, newgraph, common)
        newgraph <- add_param("shape", pathigraph$graph, newgraph, common)
        newgraph <- add_param("size", pathigraph$graph, newgraph, common)
        newgraph <- add_param("size2", pathigraph$graph, newgraph, common)
        # newgraph <- add_param("label", pathigraph$graph, newgraph, common,
        #                       init = FALSE)
        E(newgraph)$arrow.size <- .2
        V(newgraph)$label.cex <- .62
        # Labels
        labels <- sapply(V(newgraph)$name, function(name){
            if(grepl("_func", name)){
                nff[gsub("_func", "", name), "functions"]
            }else{
                V(pathigraph$graph)$label[which(V(pathigraph$graph)$name ==
                                                    name)]
            }
        })
        newgraph <- set.vertex.attribute(newgraph, "label", value = labels)
        newgraph <- set.vertex.attribute(newgraph, "color", value = "cyan")
        tooltips <- sapply(V(newgraph)$name, function(name){
            if(grepl("_func", name)){
                gsub("\n", "<br>", nff[gsub("_func", "", name), "functions"])
            }else{
                V(pathigraph$graph)$tooltip[which(
                    V(pathigraph$graph)$name == name)]
            }
        })
        newgraph <- set.vertex.attribute(newgraph, "tooltip", value = tooltips)

        # recalculate layout
        if(verbose == TRUE)
            cat("Recomputing graph layout...\n")
        fixed <- rep(TRUE, length(V(newgraph)))
        fixed[grep("func", V(newgraph)$name)] <- FALSE
        lay <- cbind(V(newgraph)$x, V(newgraph)$y)
        rl <- refine_layout(newgraph, lay, fixed, w = w, maxiter = maxiter,
                            e = e, use.edges = use.edges, p0q = p0q,
                            p0mult = p0mult)
        if(verbose == TRUE)
            rl$iter

        V(newgraph)$nodeX <- rl$layout[,1]
        V(newgraph)$x <- rl$layout[,1]
        V(newgraph)$nodeY <- rl$layout[,2]
        V(newgraph)$y <- rl$layout[,2]

        gl <- as.list(V(pathigraph$graph)$genesList)
        names(gl) <- V(pathigraph$graph)$name
        V(newgraph)$genesList <- gl[V(newgraph)$name]

        negl <- length(E(newgraph)) - length(E(pathigraph$graph))
        E(newgraph)$relation <- c(E(pathigraph$graph)$relation, rep(1, negl))

        newpathigraph$graph <- newgraph
        newpathigraph$subgraphs_funs <- create_subgraphs(newgraph)[[1]]
        newpathigraph$effector.subgraphs_funs <-
            create_effector_subgraphs(newpathigraph, funs = TRUE)
        newpathigraph$subgraphs <- pathigraph$subgraphs
        newpathigraph$rl <- rl
        newpathigraph$fixed <- rl$fixed

        return(newpathigraph)
    }

}


add_funs_to_pgs <- function(apgs, entrez_symbol, dbannot, maxiter = 100,
                            p0mult = 4, verbose = FALSE){
    if(verbose == TRUE)
        cat("Adding functions to pathways...\n")
    fpgs <- lapply(apgs,function(x) {
        if(verbose == TRUE)
            cat(x$path.id, " - ", x$path.name, "\n")
        add_funs_to_pg(x, entrez_symbol, dbannot, maxiter = maxiter,
                       p0mult = p0mult, verbose = FALSE)
    })
    return(fpgs)
}




#_____________________________________________________
#
#              LAYOUTS
#_____________________________________________________


refine_layout <- function(gg, coords, fixed, e = 10e-7, maxiter = 500, w = 0.5,
                          align.final.at.right = TRUE, p0q = 0.01, p0mult = 5,
                          use.edges = TRUE){

    cs <- coords

    V(gg)$x <- cs[,1]
    V(gg)$y <- cs[,2]

    minx <- min(cs[fixed == TRUE, 1])
    maxx <- max(cs[fixed == TRUE, 1])
    miny <- min(cs[fixed == TRUE, 2])
    maxy <- max(cs[fixed == TRUE, 2])

    n <- length(V(gg))
    #maxdis <- sqrt((maxx-minx)^2+(maxy-miny)^2)
    maxdis <- max(abs(maxx - minx), 4*abs(maxy - miny))

    # get neighbour
    global_dis <- shortest.paths(gg)
    global_dis[is.infinite(global_dis)] <-
        max(global_dis[!is.infinite(global_dis)]) + 1
    nglobal_dis <- (max(global_dis) - global_dis)*(1 - diag(nrow(global_dis)))
    nglobal_dis <- nglobal_dis/max(nglobal_dis)
    adjm <- (global_dis == 1) + 0
    dists <- c()
    for(i in seq_len(nrow(adjm))){
        anodes <- which(adjm[i,] == 1)
        for(j in anodes){
            dists <- c(dists, sqrt((cs[i,1] - cs[j,1])^2 +
                                       (cs[i,2] - cs[j,2])^2))
        }
    }
    p0 <- stats::quantile(dists[dists > 0], p0q)
    p00 <- p0 * p0mult

    if(align.final.at.right == TRUE){
        nonfixed <- intersect(which(fixed == FALSE),
                              which(!(V(gg)$name %in% get.edgelist(gg)[,1])))
        for(nf in nonfixed){
            friends <- which(global_dis[nf,] == 1)
            x <- max(V(gg)$x[friends]) + p0 * 3
            y <- mean(V(gg)$y[friends])
            cs[nf,] <- c(x,y)
        }
    }

    # Move also nodes with duplicated coordinates
    fixed[duplicated(cs)] <- FALSE

    iter <- 1
    medif <- e + 1
    locations <- list()

    while(medif > e & iter < maxiter){
        cs2 <- cs
        for(i in seq_len(n)){
            if(fixed[i] == FALSE){
                # node forces
                node_forces <- mat.or.vec(nr = n, nc = 2)
                for(j in seq_len(n)){
                    if(i != j){
                        is_uniq_rep <- !(nglobal_dis[i,j] == 1)
                        pp <- c(p00,p0)[is_uniq_rep + 1]
                        node_forces[j,] <- get_node_vecForce(cs[i,1], cs[i,2],
                                                             cs[j,1], cs[j,2],
                                                             pp, maxdis,
                                                             repeller =
                                                                 is_uniq_rep)
                        if(global_dis[i,j] == 1){
                            node_forces[j,] <- node_forces[j,] * 10
                        }
                    }
                }
                # edge forces
                if(use.edges == TRUE){
                    edges <- get.edgelist(gg)
                    ne <- nrow(edges)
                    edge_forces <- mat.or.vec(nr = ne, nc = 2)
                    for(j in seq_len(nrow(edges))){
                        if(V(gg)$name[i] != edges[j,1] &
                           V(gg)$name[i] != edges[j,2]){
                            c0 <- cs[i,]
                            c1 <- cs[which(V(gg)$name == edges[j,1]),]
                            c2 <- cs[which(V(gg)$name == edges[j,2]),]
                            edge_forces[j,] <- get_edge_vectorial_force(
                                c0, c1, c2, p0, maxdis)
                        }
                    }

                    forces <- rbind(node_forces, edge_forces)
                } else {
                    forces <- node_forces
                }
                cs2[i,] <- cs2[i,] + colMeans(forces) * w

                #print("*")

                #if(i==61) stop("pepe")

            }

        }


        medif <- max(abs(cs - cs2))

        cs <- cs2
        locations[[iter]] <- cs

        iter <- iter + 1
        # cat(iter,medif,"\n")
    }

    return(list(
        layout = cs,
        locations = locations,
        iter = iter,
        maxiter = maxiter,
        w = w,
        e = e,
        laste = medif,
        p0 = p0,
        p00 = p00,
        maxdis = maxdis,
        dists = dists,
        fixed = fixed
    ))

}

# #
# plot.node.trip <- function(locations, fixed){
#     nodes <- which(fixed == FALSE)
#     for(node in nodes){
#         lines(t(sapply(locations, function(x) c(x[node,1], x[node,2]))),
#               col = "red",
#               lty = 2)
#     }
# }


# get.node.vectorial.force
get_node_vecForce <- function(x1, y1, x2, y2, p0, maxdis, repeller = FALSE){

    dis <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
    dis <- max(abs(x1 - x2), 4 * abs(y1 - y2))

    if(repeller == TRUE){
        if(dis > p0){
            force <- 0
        } else {
            ndis <- (p0 - dis)/p0
            force <- -ndis
        }
    } else {
        if(dis > p0){
            ndis <- (dis - p0)/(maxdis - p0)
            force <- ndis
        } else {
            ndis <- (p0 - dis)/p0
            force <- -ndis
        }
    }

    v <- c(x2 - x1, y2 - y1)
    if(v[1] == 0 & v[2] == 0){
        v <- stats::rnorm(2)
    }
    v <- v/sqrt(v[1]^2 + v[2]^2)
    v <- v * force

    return(v)

}

get_edge_vectorial_force <- function(c0, c1, c2, p0, maxdis){

    x0 <- c0[1]
    y0 <- c0[2]
    x1 <- c1[1]
    y1 <- c1[2]
    x2 <- c2[1]
    y2 <- c2[2]

    mod <- sqrt( (x2 - x1)^2 + (y2 - y1)^2 )
    s <- ( (x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1)) / mod
    sx <- x1 + (s/mod)*(x2 - x1)
    sy <- y1 + (s/mod)*(y2 - y1)

    if(x1 < x2) {
        xmin <- x1
        xmax <- x2
    } else {
        xmin <- x2
        xmax <- x1
    }
    if(y1 < y2) {
        ymin <- y1
        ymax <- y2
    } else {
        ymin <- y2
        ymax <- y1
    }
    if(sx < xmin | sx > xmax | sy < ymin | sy > ymax){
        force <- 0
    } else {
        dis <- sqrt((x0 - sx)^2 + (y0 - sy)^2)
        if(dis > p0){
            force <- 0
        } else {
            ndis <- (p0 - dis)/p0
            force <- -ndis
        }
    }

    v <- c(sx - x0, sy - y0)
    if(v[1] == 0 & v[2] == 0){
        v <- stats::rnorm(2)
    }
    v <- v/sqrt(v[1]^2 + v[2]^2)
    v <- v * force

    return(v)

}



#_____________________________________________________
#
#              CREATE METAGINFO
#_____________________________________________________

create_metaginfo_object <- function(fpgs, species, by.user = FALSE,
                                    basal.value = 0.5){

    pathigraph.genes <- all_needed_genes(fpgs)

    # Create all.labelids table
    labelids <- lapply(fpgs, function(pg) cbind(pg$label.id,
                                                path.id = pg$path.id,
                                                path.name = pg$path.name))
    labelids <- do.call("rbind", labelids)
    rownames(labelids) <- labelids[,1]
    # labelids <- as.data.frame(labelids, strignsAsFactors = FALSE)

    # Todos los genes a valor 0.5
    genes.vals.05 <- matrix(basal.value, ncol=2, nrow=length(pathigraph.genes),
                            dimnames=list(pathigraph.genes, c("1", "2")))
    meta.05 <- NULL
    meta.05$pathigraphs <- fpgs
    meta.05$all.labelids <- labelids
    results.05 <- hipathia(genes.vals.05, meta.05, test = FALSE,
                           verbose = FALSE)
    results.dec.05 <- hipathia(genes.vals.05, meta.05, decompose = TRUE,
                               test = FALSE, verbose = FALSE)

    # Create metaginfo object
    metaginfo <- NULL
    metaginfo$species <- species
    metaginfo$all.genes <- pathigraph.genes
    metaginfo$path.norm <- assay(results.dec.05, "paths")[,1]
    metaginfo$eff.norm <- assay(results.05, "paths")[,1]
    metaginfo$pathigraphs <- fpgs
    metaginfo$all.labelids <- labelids
    metaginfo$group.by <- "pathways"
    metaginfo$by.user <- by.user

    return(metaginfo)
}



#_____________________________________________________
#
#              OTHERS
#_____________________________________________________


all_needed_genes <- function(pathigraphs){
    genes <- unique(unlist(sapply(pathigraphs, function(x){
        unique(unlist(V(x$graph)$genesList))
    })))
    return(genes[!is.na(genes) & genes!="/"])
}

