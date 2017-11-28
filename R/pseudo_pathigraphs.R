######## ######## ######## ######## ######## ######## ######## ########
######## PSEUDO PATHIGRAPHS
######## ######## ######## ######## ######## ######## ######## ########

get.pseudo.pathigraphs <- function(subgraphs, annots){

  categories <- unique(annots[,2])
  pseudo_pathigraphs <- list()
  for(i in 1:length(categories)){
    cat("Processing ", categories[i], " (", i, " of ", length(categories),
        ")\n", sep = "")
    selnames <- annots[which(annots[,2] == categories[i]),1]
    selsub <- subgraphs[selnames]
    cat("    found ", length(selsub), " subgraphs...\n", sep = "")
    path.id <- paste0("term_", i)
    cpfs <- create.pathigraph.from.subgraphs(selsub, categories[i], path.id)
    pseudo_pathigraphs[[path.id]] <- cpfs
  }

  return(pseudo_pathigraphs)

}


create.pathigraph.from.subgraphs <- function(selsub, path.name, path.id){

  # # relabel nodes by pathway
  # for(i in 1:length(selsub)){
  #   pathway <- strsplit(names(selsub)[i],"__")[[1]][1]
  #   V(selsub[[i]])$name <- paste0(pathway,"__",V(selsub[[i]])$name)
  # }

  # overlapping between paths
  imat <- data.matrix(mat.or.vec(nr = length(selsub), nc = length(selsub)))
  rownames(imat) <- names(selsub)
  colnames(imat) <- names(selsub)
  for(i in 1:length(selsub)){
    for(j in 1:length(selsub)){
      imat[i,j] <- length(intersect(V(selsub[[i]])$name, V(selsub[[j]])$name))
    }
  }
  ga <- graph.adjacency(imat > 0)
  cga <- clusters(ga)
  cga$groups <- unique(cga$membership)

  # recalculate Y coordinates
  path_coords <- do.call("rbind", lapply(selsub,function(x) range(V(x)$nodeY)))
  init_path_coords <- cbind(path_coords, 0, 0, 0, 0, 0)
  group_coords <- do.call("rbind",by(path_coords, cga$membership, range))

  ymargin <- 50
  lasty <- 0
  for(i in 1:cga$no){
    m <- cga$groups[i]
    indexes <- which(cga$membership == m)
    for(j in indexes){
      V(selsub[[j]])$nodeY <- V(selsub[[j]])$nodeY - group_coords[i,1] + lasty
      init_path_coords[j,3] <- group_coords[i,1]
      init_path_coords[j,4] <- lasty
      init_path_coords[j,5] <- m
      init_path_coords[j,6] <- path_coords[j,1] - group_coords[i,1] + lasty
      init_path_coords[j,7] <- path_coords[j,2] - group_coords[i,1] + lasty
    }
    lasty <- lasty + (group_coords[i,2] - group_coords[i,1]) + ymargin
  }

  # create general graphs

  node_atts <- unique(do.call("rbind",lapply(selsub,function(x) {
    data.frame(name = V(x)$name,
               label = V(x)$label,
               shape = V(x)$shape,
               x = V(x)$nodeX,
               y = V(x)$nodeY,
               width = V(x)$width,
               height = V(x)$height,
               label.color = V(x)$label.color,
               label.cex = V(x)$label.cex,
               stringsAsFactors = F)
  })))
  rownames(node_atts) <- node_atts$name

  edges <- unique(do.call("rbind", lapply(selsub, function(x) {
      el <- get.edgelist(x)
      data.frame(source = el[,1], relation = E(x)$relation, target = el[,2])
      })))

  supergraph <- graph.data.frame(edges[,c(1,3)], directed = TRUE)
  E(supergraph)$relation <- edges[,2]
  for(k in 1:ncol(node_atts)){
    field <- colnames(node_atts)[k]
    supergraph <- set.vertex.attribute(supergraph,
                                       name = field,
                                       value = node_atts[V(supergraph)$name,
                                                         field])
  }
  V(supergraph)$nodeX <- V(supergraph)$x
  V(supergraph)$nodeY <- V(supergraph)$y
  V(supergraph)$genesList <- rep(NA, length(V(supergraph)$name))

  pathigraph <- list()
  pathigraph$graph <- supergraph
  pathigraph$path.name <- path.name
  pathigraph$path.id <- path.id
  pathigraph$effector.subgraphs <- selsub
  pathigraph$path_coords <- path_coords
  pathigraph$init_path_coords <- init_path_coords
  pathigraph$group_coords <- group_coords
  pathigraph$ga <- ga
  pathigraph$cga <- cga

  return(pathigraph)

}

create.html.report3 <- function(pseudo_pathigraphs, metaginfo,comp, home,
                                output.folder, effector = FALSE, conf = 0.05,
                                extra_javascript = "", after_html = "",
                                before_html = "", clean_out_folder = TRUE,
                                template_name = "index_template.html",
                                output_name = "index.html"){


  # if(clean_out_folder==T & file.exists(paste0(output.folder,"/report"))){
  #   unlink(paste0(output.folder,"/report"),recursive=T)
  # }
  # if(!file.exists(paste0(output.folder,"/report")))
    # dir.create(paste0(output.folder,"/report"))
  # cp_command <- paste0("cp -r ", home, "/network-viewer/* ",
  #                      output.folder, "/report/")
  # system(cp_command)
  #
  # index <- scan(paste0(output.folder, "/report/", template_name),
  #               comment.char = "", sep = "\n", what = "character")


  if(clean_out_folder == TRUE &
     file.exists(paste0(output.folder, "/pathway-viewer"))){
    unlink(paste0(output.folder, "/pathway-viewer"), recursive = TRUE)
    unlink(paste0(output.folder, "/index.html"), recursive = TRUE)
  }
  if(!file.exists(paste0(output.folder,"/pathway-viewer")))
      dir.create(paste0(output.folder,"/pathway-viewer"))
  cp_command <- paste0("cp -r ", home, "/pathway-viewer/* '", output.folder,
                       "/pathway-viewer/'")
  system(cp_command)
  mv_legend_command <- paste0("cp  ", home, "/report-files/pretty_legend.png '",
                              output.folder, "'")
  system(mv_legend_command)

  index <- scan(paste0(home, '/report-files/', template_name),
                comment.char = "", sep = "\n", what = "character")


  global_div <- c()
  global_define <- c()
  global_init <- c()

  global_div <- c(global_div, paste0("<div id='pathway_viewer'>"))
  global_div <- c(global_div, paste0("  <div id='mini_pathway_title'>",
                                     "Pathways list</div>"))
  global_div <- c(global_div, paste0("  <div id='pathway_title'></div>"))
  global_div <- c(global_div, paste0("  <div id='network_viewer_scroll'>",
                                     "<div id='network_viewer'></div></div>"))
  global_div <- c(global_div, paste0("  <div id='pathway_selector'>"))
  for(pathway in names(pseudo_pathigraphs)){
    global_div <- c(global_div,
                    paste0("  <div class='pathway_selector_item aaa' data-",
                           "pathway='", pathway, "' data-pathwayname='",
                           pseudo_pathigraphs[[pathway]]$path.name,
                           "' onclick='handlePathwayClick(\"", pathway, "\",\"",
                           pseudo_pathigraphs[[pathway]]$path.name,
                           "\")'>",
                           pseudo_pathigraphs[[pathway]]$path.name,
                           "</div>"))
  }
  global_div <- c(global_div, paste0("  </div>"))
  global_div <- c(global_div, paste0(" <div id='mini_path_title'>",
                                     "Available subpaths</div>"))
  global_div <- c(global_div, paste0("  <input id='pathwayFilterInput' ",
                                     "type='text' placeholder='Search...'/>"))
  global_div <- c(global_div, paste0("  <div id='path_selector'>"))

  global_div <- c(global_div, paste0("  </div>"))
  global_div <- c(global_div, paste0("</div>"))

  #
  for(pathway in names(pseudo_pathigraphs)){

    # define
    raw_sif <- read.table(paste0(output.folder, "/pathways/", pathway, ".sif"),
                          sep = "\n", header = FALSE, stringsAsFactors = FALSE)
    sif <- paste0("sif: '", paste(raw_sif$V1, collapse = "\\\\n"), "'")

    # node att
    raw_att <- read.table(paste0(output.folder, "/pathways/", pathway, ".natt"),
                          sep = "\n", header = FALSE, stringsAsFactors = FALSE,
                          comment.char = "")
    raw_att$V1 <- gsub("\\\\n", " & ", raw_att$V1)
    att <- paste0("natt: '#", paste(raw_att$V1, collapse = "\\\\n"), "'")

    # edge att
    raw_eatt <- read.table(paste0(output.folder, "/pathways/",
                                  pathway, ".eatt"),
                           sep = "\n", header = FALSE,
                           stringsAsFactors = F, comment.char = "")
    raw_eatt$V1 <- gsub("\\\\n", " & ", raw_eatt$V1)
    eatt <- paste0("eatt: '#", paste(raw_eatt$V1, collapse = "\\\\n"), "'")

    # path list
    if(effector == TRUE){
      s <- pseudo_pathigraphs[[pathway]]$effector.subgraphs
    }else{
      s <- pseudo_pathigraphs[[pathway]]$subgraphs
    }

    pathlist <- paste0("plist: ",
                       paste("[", paste("\"", names(s), "\"", collapse = ",",
                                        sep = ""), "]"))
    prettynames <- gsub("\\*", "", sapply(strsplit(get.path.names(
        metaginfo, names(s)), ": "), "[[", 2))
    prettynames <- gsub("\n", " - ", prettynames)
    prettypathlist <- paste0("pretty_plist: ",
                             paste("[", paste("\"", prettynames, "\"",
                                              collapse = ",", sep = ""), "]"))

    pcolors <- c("UP"="darkred","DOWN"="darkblue")[comp[names(s),"UP/DOWN"]]
    pcolors[comp[names(s),"FDRp.value"]>=conf] <- "darkgrey"

    prettycolors <- paste0("pretty_color_plist: ", paste(
        "[", paste("\"", pcolors, "\"", collapse = ",", sep = ""), "]"))


    pathstr <- paste0(pathway, ":{\n",
                      paste(sif, att, eatt, pathlist, prettypathlist,
                            prettycolors, sep = ",\n"), "\n}")

    global_define <- c(global_define,pathstr)

  }

  global_init <- c(global_init, paste0("this.loadPathway(\"",
                                       names(pseudo_pathigraphs)[1], "\",\"",
                                       pseudo_pathigraphs[[1]]$path.name,
                                       "\",el);"))

  new_index <- gsub("PUT_HERE_YOUR_ELEMENTS", paste(global_div,collapse="\n"),
                    gsub("PUT_HERE_YOUR_NETWORKS",
                         paste(global_define, collapse = ",\n"),
                         gsub("INIT_HERE_YOUR_NETWORKS",
                              paste(global_init, collapse = "\n"),index)))

  new_index <- gsub("PUT_HERE_EXTRA_JAVASCRIPT", extra_javascript, new_index)
  new_index <- gsub("PUT_HERE_BEFORE_HTML", before_html, new_index)
  new_index <- gsub("PUT_HERE_AFTER_HTML", after_html, new_index)



  write(paste(new_index, collapse = "\n"),
        file = paste0(output.folder, "/report/", output_name))

}



create.pseudo.path.info <- function(all_comp, metaginfo.pseudo, metaginfo){
  fpgs <- metaginfo.pseudo$pathigraphs
  path_info<-lapply(fpgs,function(fpg) all_comp[names(fpg$effector.subgraphs),])
  path_json_list <- lapply(names(path_info), function(x){
    out <- paste0("{\n\t\"id\":\"", x, "\",\n")
    out <- paste0(out, "\t\"name\":\"", fpgs[[x]]$path.name, "\",\n")
    anysig <- FALSE
    anyup <- FALSE
    anydown <- FALSE
    anysigup <- FALSE
    anysigdown <- FALSE
    anychanged <- FALSE
    for(i in 1:nrow(path_info[[x]])){
      if(path_info[[x]]$has_change[i] == TRUE)
          anychanged <- TRUE
      if(path_info[[x]]$FDRp.value[i] <= 0.05) {
        anysig <- TRUE
        if(path_info[[x]]$status[i] == "UP")
            anysigup <- TRUE
        if(path_info[[x]]$status[i] == "DOWN")
            anysigdown <- TRUE
      }
      if(path_info[[x]]$status[i] == "UP")
          anyup <- TRUE
      if(path_info[[x]]$status[i] == "DOWN")
          anydown <- TRUE
    }
    out <- paste0(out, "\t\"haschanged\":", tolower(anychanged), ",\n")
    out <- paste0(out, "\t\"sig\":", tolower(anysig), ",\n")
    out <- paste0(out, "\t\"up\":", tolower(anyup), ",\n")
    out <- paste0(out, "\t\"down\":", tolower(anydown), ",\n")
    out <- paste0(out, "\t\"upsig\":", tolower(anysigup), ",\n")
    out <- paste0(out, "\t\"downsig\":", tolower(anysigdown), ",\n")
    out <- paste0(out, "\t\"paths\":[\n")
    for(i in 1:nrow(path_info[[x]])){
      out <- paste0(out, "\t\t{")
      out <- paste0(out, "\"id\":\"", rownames(path_info[[x]])[i], "\", ")
      out <- paste0(out, "\"name\":\"",
                    get.path.names(metaginfo,
                                   rownames(path_info[[x]])[i]),"\", ")
      # out <- paste0(out, "\"shortname\":\"",
      #               gsub("\\*", "", strsplit(get.path.names(
      #                   metaginfo, rownames(path_info[[x]])[i]),
      #                   ": ")[[1]][2]),
      #               "\", ")
      out <- paste0(out, "\"shortname\":\"",
                    gsub("\\*", "",
                         get.path.names(metaginfo,
                                        rownames(path_info[[x]])[i])), "\", ")
      out <- paste0(out, "\"pvalue\":", path_info[[x]]$FDRp.value[i], ", ")
      out <- paste0(out, "\"status\":\"", path_info[[x]]$status[i], "\", ")
      out <- paste0(out, "\"sig\":\"",
                    tolower(path_info[[x]]$FDRp.value[i] < 0.05), "\", ")
      out <- paste0(out, "\"haschanged\":",
                    tolower(path_info[[x]]$has_change[i]), ", ")
      out <- paste0(out, "\"up\":",
                    tolower(path_info[[x]]$status[i] == "UP"), ", ")
      out <- paste0(out, "\"down\":",
                    tolower(path_info[[x]]$status[i] == "DOWN"), ", ")
      out <- paste0(out, "\"upsig\":",
                    tolower(path_info[[x]]$status[i] == "UP" &
                                path_info[[x]]$FDRp.value[i] < 0.05),", ")
      out <- paste0(out, "\"downsig\":",
                    tolower(path_info[[x]]$status[i] == "DOWN" &
                                path_info[[x]]$FDRp.value[i] < 0.05),", ")
      out <- paste0(out, "\"color\":\"", path_info[[x]]$color[i], "\"")
      out <- paste0(out, "}")
      if(i == nrow(path_info[[x]])){
        out <- paste0(out, "\n")
      } else {
        out <- paste0(out, ",\n")
      }
    }
    out <- paste0(out, "\t]\n")
    out <- paste0(out, "}")
    out
  })
  path_json <- paste0("[\n", paste(path_json_list, collapse = ","),"\n]")
  return(path_json)
}



create.report3 <- function(results, comp, pseudo, metaginfo, output.folder,
                           home, node.colors = NULL, conf = 0.05,
                           verbose = FALSE){

    if(!is.null(node.colors)){
        moreatts <- summarize.atts(list(node.colors), c("color"))
    }else{
        moreatts <- NULL
    }

    # save.results(results, comp, pseudo, output.folder)

    if(!file.exists(paste0(output.folder, "/pathways/")))
        dir.create(paste0(output.folder, "/pathways/"))

    for(pathway in names(pseudo$pathigraphs)){
        if(verbose == TRUE) cat(pathway)
        write.attributes(comp, pathway, pseudo,
                         paste0(output.folder, "/pathways/", pathway),
                         moreatts_pathway = moreatts[[pathway]],
                         conf = conf)
        # write.sif.files( pathigraphs[[pathway]]$graph,
        #                  nodeColor,
        #                  paste0(output.folder, "/pathways/", pathway))
    }

    pv.path <- paste0(system.file("extdata", package = "hipathia"))
    create.html.report3(pseudo$pathigraphs,
                        metaginfo = metaginfo,
                        comp = comp,
                        home = home,
                        output.folder = output.folder,
                        effector = T )
    # create.html.report3(metaginfo,
    #                     comp,
    #                     pv.path,
    #                     output.folder,
    #                     template_name = "index_template.html",
    #                     output_name = "index.html",
    #                     clean_out_folder = FALSE)

    # pathways.folder <- paste0(output.folder, "/pathway-viewer/pathways/")
    # if(!file.exists(pathways.folder))
    #     dir.create(pathways.folder)
    # for(pathway in names(metaginfo$pathigraphs)){
    #     if(verbose == TRUE) cat(pathway)
    #     write.attributes(comp,
    #                      pathway,
    #                      metaginfo,
    #                      paste0(pathways.folder, pathway),
    #                      moreatts_pathway=moreatts[[pathway]],
    #                      conf=conf)
    # }

    comp$status <- comp$"UP/DOWN"
    comp$has_changed <- TRUE
    path_json <- create.pseudo.path.info(comp, pseudo.meta, pathways )
    write(path_json, file = paste0(output.folder, "/pathways/path_info.json"))
}


create.html.report3 <- function(pseudo_pathigraphs, metaginfo,comp, home,
                                output.folder, effector = FALSE, conf = 0.05,
                                extra_javascript = "", after_html = "",
                                before_html = "", clean_out_folder = TRUE,
                                template_name = "index_template.html",
                                output_name = "index.html"){


    # if(clean_out_folder==T & file.exists(paste0(output.folder,"/report"))){
    #   unlink(paste0(output.folder,"/report"),recursive=T)
    # }
    # if(!file.exists(paste0(output.folder,"/report")))
    #     dir.create(paste0(output.folder,"/report"))
    # cp_command <- paste0("cp -r ",home,"/network-viewer/* ",
    #                      output.folder,"/report/")
    # system(cp_command)
    #
    # index <- scan(paste0(output.folder,"/report/",template_name),
    #               comment.char = "",sep="\n",what="character")


    if(clean_out_folder == TRUE &
       file.exists(paste0(output.folder, "/pathway-viewer"))){
        unlink(paste0(output.folder,"/pathway-viewer"), recursive = TRUE)
        unlink(paste0(output.folder,"/index.html"), recursive = TRUE)
    }
    if(!file.exists(paste0(output.folder, "/pathway-viewer")))
        dir.create(paste0(output.folder, "/pathway-viewer"))
    cp_command <- paste0("cp -r ", home, "/pathway-viewer/* '", output.folder,
                         "/pathway-viewer/'")
    system(cp_command)
    mv_legend_command <- paste0("cp  ", home,
                                "/report-files/pretty_legend.png '",
                                output.folder, "'")
    system(mv_legend_command)

    index <- scan(paste0(home, '/report-files/', template_name),
                  comment.char = "", sep = "\n", what = "character")


    global_div <- c()
    global_define <- c()
    global_init <- c()

    global_div <- c(global_div, paste0("<div id='pathway_viewer'>"))
    global_div <- c(global_div,
                    paste0("  <div id='mini_pathway_title'>",
                           "Pathways list</div>"))
    global_div <- c(global_div, paste0("  <div id='pathway_title'></div>"))
    global_div <- c(global_div,
                    paste0("  <div id='network_viewer_scroll'>",
                           "<div id='network_viewer'></div></div>"))
    global_div <- c(global_div, paste0("  <div id='pathway_selector'>"))
    for(pathway in names(pseudo_pathigraphs)){
        global_div <- c(global_div,
                        paste0("  <div class='pathway_selector_item aaa' ",
                               "data-pathway='", pathway,
                               "' data-pathwayname='",
                               pseudo_pathigraphs[[pathway]]$path.name,
                               "' onclick='handlePathwayClick(\"", pathway,
                               "\",\"", pseudo_pathigraphs[[pathway]]$path.name,
                               "\")'>", pseudo_pathigraphs[[pathway]]$path.name,
                               "</div>"))
    }
    global_div <- c(global_div, paste0("  </div>"))
    global_div <- c(global_div, paste0(" <div id='mini_path_title'>",
                                       "Available subpaths</div>"))
    global_div <- c(global_div, paste0("  <input id='pathwayFilterInput' ",
                                       "type='text' placeholder='Search...'/>"))
    global_div <- c(global_div, paste0("  <div id='path_selector'>"))

    global_div <- c(global_div, paste0("  </div>"))
    global_div <- c(global_div, paste0("</div>"))

    #
    for(pathway in names(pseudo_pathigraphs)){

        # define
        raw_sif <- read.table(paste0(output.folder, "/pathways/",
                                     pathway, ".sif"),
                              sep = "\n", header = FALSE,
                              stringsAsFactors = FALSE)
        sif <- paste0("sif: '",paste(raw_sif$V1,collapse="\\\\n"),"'")

        # node att
        raw_att <- read.table(paste0(output.folder, "/pathways/",
                                     pathway, ".natt"),sep="\n",
                              header = FALSE,
                              stringsAsFactors = FALSE, comment.char = "")
        raw_att$V1 <- gsub("\\\\n"," & ",raw_att$V1)
        att <- paste0("natt: '#", paste(raw_att$V1, collapse = "\\\\n"), "'")

        # edge att
        raw_eatt <- read.table(paste0(output.folder, "/pathways/",
                                      pathway, ".eatt"),
                               sep = "\n", header = FALSE,
                               stringsAsFactors = FALSE, comment.char = "")
        raw_eatt$V1 <- gsub("\\\\n"," & ",raw_eatt$V1)
        eatt <- paste0("eatt: '#",paste(raw_eatt$V1, collapse = "\\\\n"), "'")

        # path list
        if(effector==T){
            s <- pseudo_pathigraphs[[pathway]]$effector.subgraphs
        }else{
            s <- pseudo_pathigraphs[[pathway]]$subgraphs
        }

        pathlist <- paste0("plist: ",
                           paste("[",
                                 paste("\"", names(s), "\"",
                                       collapse = ",", sep = ""), "]"))
        prettynames <- gsub("\\*", "",
                            sapply(strsplit(get.path.names(metaginfo,names(s)),
                                            ": "),"[[",2))
        prettynames <- gsub("\n"," - ",prettynames)
        prettypathlist <- paste0("pretty_plist: ",
                                 paste("[", paste("\"", prettynames,
                                                  "\"", collapse = ",",
                                                  sep = ""),
                                       "]"))

        pcolors <- c("UP"="darkred","DOWN"="darkblue")[comp[names(s),"UP/DOWN"]]
        pcolors[comp[names(s),"FDRp.value"]>=conf] <- "darkgrey"

        prettycolors <- paste0("pretty_color_plist: ",
                               paste("[", paste("\"", pcolors,
                                                "\"", collapse = ",", sep = ""),
                                     "]"))


        pathstr <- paste0(pathway, ":{\n",
                          paste(sif, att, eatt, pathlist, prettypathlist,
                                prettycolors, sep = ",\n"),
                          "\n}")

        global_define <- c(global_define,pathstr)

    }

    global_init <- c(global_init,
                     paste0("this.loadPathway(\"", names(pseudo_pathigraphs)[1],
                            "\",\"", pseudo_pathigraphs[[1]]$path.name,
                            "\",el);"))

    new_index <- gsub("PUT_HERE_YOUR_ELEMENTS",
                      paste(global_div, collapse = "\n"),
                      gsub("PUT_HERE_YOUR_NETWORKS",
                           paste(global_define, collapse = ",\n"),
                           gsub("INIT_HERE_YOUR_NETWORKS",
                                paste(global_init, collapse = "\n"), index)))

    new_index <- gsub("PUT_HERE_EXTRA_JAVASCRIPT", extra_javascript, new_index)
    new_index <- gsub("PUT_HERE_BEFORE_HTML", before_html, new_index)
    new_index <- gsub("PUT_HERE_AFTER_HTML", after_html, new_index)



    write(paste(new_index, collapse = "\n"),
          file = paste0(output.folder, "/", output_name))

}



create.node.and.edge.attributes <- function(comp, pathway, metaginfo,
                                            moreatts_pathway = NULL,
                                            conf = 0.05, effector = TRUE,
                                            reverse_xref = NULL, exp = NULL){
    print(pathway)
    pathigraphs <- metaginfo$pathigraphs
    pcomp <- comp[grep(paste0(pathway,"__"), rownames(comp)),]
    if(effector){
        s <- pathigraphs[[pathway]]$effector.subgraphs
    }else{
        s <- pathigraphs[[pathway]]$subgraphs
    }
    ig <- pathigraphs[[pathway]]$graph

    if(is.null(V(ig)$type))
        V(ig)$type <- "node"
    if(is.null(V(ig)$width))
        V(ig)$width <- 15
    if(is.null(V(ig)$height))
        V(ig)$height <- 5
    if(is.null(V(ig)$label.color))
        V(ig)$label.color <- "black"
    if(is.null(V(ig)$label.cex))
        V(ig)$label.cex <- 0.7
    V(ig)$stroke.color <- "lightgrey"
    #V(ig)$stroke.color <- find.node.colors(pcomp, s, ig, conf)[V(ig)$name]

    #V(ig)$color <- "white"
    V(ig)$width[V(ig)$shape=="circle"] <- 15
    V(ig)$width[V(ig)$shape!="circle"] <- 22
    V(ig)$shape[V(ig)$shape=="rectangle" &
                    !grepl("func", V(ig)$name)] <- "ellipse"
    V(ig)$shape[V(ig)$shape=="rectangle" &
                    grepl("func", V(ig)$name)] <- "rectangle"

    natt <- cbind(V(ig)$name, V(ig)$label, 10, V(ig)$nodeX, V(ig)$nodeY,
                  V(ig)$stroke.color, V(ig)$shape, V(ig)$type, V(ig)$label.cex,
                  V(ig)$label.color, V(ig)$width, V(ig)$height,
                  sapply(V(ig)$genesList, paste, collapse=","))
    colnames(natt) <- c("ID","label","labelSize","X", "Y", "color", "shape",
                        "type", "label.cex", "label.color", "width", "height",
                        "genesList")
    rownames(natt) <- natt[,1]
    natt[,"label"] <- gsub("\n", " / ", natt[,"label"])
    if(!is.null(moreatts_pathway)){
        natt <- cbind(natt, moreatts_pathway)
    }
    node_path_assoc <- matrix(0, nrow = nrow(natt), ncol = length(s))
    colnames(node_path_assoc) <- names(s)
    natt <- cbind(natt,node_path_assoc)

    sif <- c()
    eatt <- c()
    epath_assoc <- c()

    for(i in 1:length(s)){

        # get subgraph
        subgraph <- s[[i]]
        name <- names(s)[i]
        #cname <- paste(pathway,name,sep="__")
        cname <- name
        pname <- get.path.names(metaginfo, name)

        # sif
        raw_edges <- get.edgelist(subgraph)
        type <- c("activation","inhibition")[(E(subgraph)$relation == -1) + 1]
        edges <- cbind(raw_edges[,1], type, raw_edges[,2])
        sif <- rbind(sif,edges)

        # edge attributes
        eids <- apply(edges, 1, function(x) paste0(x, collapse = "_"))
        status <- comp[cname, "UP/DOWN"]
        if("color" %in% colnames(comp)){
            color <- comp[cname, "color"]
        } else {
            if( comp[cname, "FDRp.value"] < conf){
                color <- c("#1f78b4", "#e31a1c")[(status == "UP") + 1]
            } else {
                color <- "darkgrey"
            }
        }
        path_assoc <- matrix(0, nrow = nrow(edges), ncol = length(s))
        colnames(path_assoc) <- names(s)
        path_assoc[,name] <- 1
        edges_atts <- cbind(id = eids,
                            status = status,
                            color = color,
                            name = name,
                            cname = cname,
                            pname = pname,
                            pvalue = comp[cname,"p.value"],
                            adj.pvalue = comp[cname,"FDRp.value"])
        eatt <- rbind(eatt,edges_atts)

        epath_assoc <- rbind(epath_assoc,path_assoc)

        # node attributes
        natt[get.vertex.attribute(subgraph, "name"), name] <- 1

    }

    # melt multi path interactions
    unique_edges <- unique(eatt[,1])
    def_eatt <- c()
    def_sif <- c()
    def_epath_assoc <- c()

    for(ue in unique_edges){

        indexes <- which(eatt[,1] == ue)
        subeatt <- eatt[indexes,,drop = FALSE]
        subepath_assoc <- epath_assoc[indexes,,drop = FALSE]
        subsif <- sif[indexes,,drop = FALSE]

        anysig <- F
        # up regulated
        upsig <- which(subeatt[,"status"] == "UP" &
                           as.numeric(subeatt[,"adj.pvalue"]) < conf)
        if(length(upsig ) > 0){

            selected_subsif <- subsif[1,]
            selected_subsif[2] <- paste0(selected_subsif[2], ".up")
            def_sif <- rbind(def_sif, selected_subsif)

            mini_subeatt <- subeatt[upsig,, drop = FALSE]
            selected_subeatt <- mini_subeatt[1, c("id", "status", "color",
                                                  "pvalue", "adj.pvalue")]
            selected_subeatt["id"] <- paste(selected_subsif, collapse = "_")
            def_eatt <- rbind(def_eatt, selected_subeatt)

            selected_subepath_assoc <- subepath_assoc[upsig,,drop = FALSE]
            def_epath_assoc <- rbind(def_epath_assoc,
                                     colSums(selected_subepath_assoc) > 0)

            anysig <- T
        }

        # down regulated
        downsig <- which(subeatt[,"status"] == "DOWN" &
                             as.numeric(subeatt[,"adj.pvalue"]) < conf)
        if(length(downsig) > 0){

            selected_subsif <- subsif[1,]
            selected_subsif[2] <- paste0(selected_subsif[2], ".down")
            def_sif <- rbind(def_sif, selected_subsif)

            mini_subeatt <- subeatt[downsig,, drop = FALSE]
            selected_subeatt <- mini_subeatt[1, c("id", "status", "color",
                                                  "pvalue", "adj.pvalue")]
            selected_subeatt["id"] <- paste(selected_subsif, collapse = "_")
            def_eatt <- rbind(def_eatt, selected_subeatt)

            selected_subepath_assoc <- subepath_assoc[downsig,,drop = FALSE]
            def_epath_assoc <- rbind(def_epath_assoc,
                                     colSums(selected_subepath_assoc) > 0)

            anysig <- TRUE
        }

        # no sigs
        nosigs <- which(as.numeric(subeatt[,"adj.pvalue"]) >= conf)
        if(length(nosigs)>0){

            selected_subsif <- subsif[1,]
            def_sif <- rbind(def_sif, selected_subsif)

            mini_subeatt <- subeatt[nosigs,, drop = FALSE]
            selected_subeatt <- mini_subeatt[1,c("id", "status", "color")]
            def_eatt <- rbind(def_eatt, selected_subeatt)

            selected_subepath_assoc <- subepath_assoc[nosigs,, drop = FALSE]
            def_epath_assoc <- rbind(def_epath_assoc,
                                     colSums(selected_subepath_assoc) > 0)

        }




        #     if(length(indexes)==1){
        #       def_eatt <- rbind(def_eatt,subeatt[,c("id","status","color")])
        #       def_sif <- rbind(def_sif,subsif)
        #       def_epath_assoc <- rbind(def_epath_assoc,subepath_assoc)
        #     } else {
        #
        #       down_ind <- which(subeatt[,"status"]=="DOWN")
        #       if(length(down_ind)>0){
        #         down_subeatt <- subeatt[down_ind,,drop=F]
        #
        #         selected_subsif <- subsif[1,]
        #         selected_subsif[2] <- paste0(selected_subsif[2],".down")
        #
        #         selected_down_subeatt <- down_subeatt[1,c("id","status","color")]
        #         selected_down_subeatt["id"] <- paste0(selected_subsif,collapse="_")
        #
        #         def_eatt <- rbind(def_eatt,selected_down_subeatt)
        #         def_sif <- rbind(def_sif,selected_subsif)
        #
        #         down_subepath_assoc <- subepath_assoc[down_ind,,drop=F]
        #         def_epath_assoc <- rbind(def_epath_assoc,colSums(down_subepath_assoc)>0)
        #       }
        #
        #       up_ind <- which(subeatt[,"status"]=="UP")
        #       if(length(up_ind)>0){
        #         up_subeatt <- subeatt[up_ind,,drop=F]
        #
        #         selected_subsif <- subsif[1,]
        #         selected_subsif[2] <- paste0(selected_subsif[2],".up")
        #
        #         selected_up_subeatt <- up_subeatt[1,c("id","status","color")]
        #         selected_up_subeatt["id"] <- paste0(selected_subsif,collapse="_")
        #
        #         def_eatt <- rbind(def_eatt,selected_up_subeatt)
        #         def_sif <- rbind(def_sif,selected_subsif)
        #
        #         up_subepath_assoc <- subepath_assoc[up_ind,,drop=F]
        #         def_epath_assoc <- rbind(def_epath_assoc,colSums(up_subepath_assoc)>0)
        #       }
        #
        #     }

    }

    rownames(def_eatt) <- NULL
    def_eatt <- as.data.frame(def_eatt, stringsAsFactors = FALSE)
    def_epath_assoc <- as.data.frame(def_epath_assoc, stringsAsFactors = FALSE)
    rownames(def_sif) <- NULL
    def_sif <- as.data.frame(def_sif, stringsAsFactors = FALSE)

    def_eatt$head <- c("inhibited", "directed")[grepl("activation",
                                                      def_sif[,2]) + 1]

    def_eatt <- cbind(def_eatt, (def_epath_assoc == TRUE) + 0)

    natt[,"label"] <- gsub("\\*", "", natt[,"label"])

    # Add functions
    #---------------------
    left <- which(grepl("func", get.edgelist(ig)[,2]))
    if(length(left) > 0 ){
        if(length(left) == 1){
            ids <- paste(get.edgelist(ig)[left,1], "activation",
                         get.edgelist(ig)[left,2], sep = "_")
        }else{
            ids <- apply(get.edgelist(ig)[left,], 1, function(x){
                paste(x[1], "activation", x[2], sep="_")
                })
        }
        funejes <- as.data.frame(matrix(0,
                                        nrow = length(ids),
                                        ncol = ncol(def_eatt)),
                                 stringsAsFactors = FALSE)
        colnames(funejes) <- colnames(def_eatt)
        rownames(funejes) <- ids
        funejes$id <- ids
        funejes$status <- "DOWN"
        funejes$color <- "darkgrey"
        if("pvalue" %in% colnames(funejes))
            funejes$pvalue <- ids
        if("adj.pvalue" %in% colnames(funejes))
            funejes$adj.pvalue <- "DOWN"
        funejes$head <- "directed"
        nods <- get.edgelist(ig)[left,1]
        names(nods) <- ids
        names(ids) <- nods
        if(effector == T){
            funs <- t(apply(funejes, 1, function(x){
                if(any(colnames(funejes) == nods[as.character(x[1])])){
                    x[which(colnames(funejes) == nods[as.character(x[1])])] <- 1
                    x
                }else{
                    x
                }
            }))
        }else{
            funs <- t(apply(funejes, 1, function(x){
                lastnodes <- sapply(colnames(funejes), function(col){
                    unlist(strsplit(col, split=" - "))[2]
                    })
                if(any(lastnodes==nods[x[1]])){
                    x[which(lastnodes == nods[x[1]])] <- 1
                    x
                }else{
                    x
                }
            }))
        }
        funs <- as.data.frame(funs, stringsAsFactors = FALSE)
        sif_funs <- data.frame(V1 = get.edgelist(ig)[left,1],
                               type = rep("activation", times=length(left)),
                               V3 = get.edgelist(ig)[left,2],
                               stringsAsFactors = FALSE)

        def_sif <- rbind(def_sif, sif_funs)
        def_eatt <- rbind(def_eatt, funs)
    }

    fun_indexes <- grep("_func", rownames(natt))
    fun_names <- rownames(natt)[fun_indexes]
    if(length(fun_indexes) > 0){
        for(i in 1:length(fun_names)){
            pp <- gsub("_func", "", fun_names[i])
            if(effector == TRUE){
                natt[fun_names[i],pp] <- 1
            } else {
                natt[fun_names[i], grep(paste0("- ", pp), colnames(natt))] <- 1
            }
        }
    }

    if(!is.null(reverse_xref)){
        sids <- strsplit(as.character(natt[,"genesList"]), split = ",")
        translate_ids <- function(ids){
            if(length(ids) > 0){
                ids <- setdiff(ids, "/")
                tids <- sapply(reverse_xref[ids], function(x){
                    if(is.null(x)){
                        return("?")
                    }else{
                        return(x)
                    }
                })
                return(paste(tids, collapse = ","))
            } else {
                return("?")
            }
        }
        natt <- cbind(natt, tids = sapply(sids, translate_ids))
    }
    if(!is.null(exp)){
        sids <- strsplit(as.character(natt[,"genesList"]), split = ",")
        ids_list <- as.list(1:nrow(exp))
        names(ids_list) <- rownames(exp)
        get_expr_ids <- function(ids){
            if(length(ids) > 0){
                ids <- setdiff(ids, "/")
                exp_values <- sapply(ids_list[ids],function(x){
                    if(is.null(x)){
                        return("?")
                    }else{
                        return(exp[x,])
                    }
                })
                return(paste(exp_values, collapse = ","))
            } else {
                return("?")
            }
        }
        natt <- cbind(natt, exp_values = sapply(sids, get_expr_ids))
    }

    return(list(sif = def_sif, edge_att = def_eatt, node_att = natt))
}
