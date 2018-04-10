##
## save.R
## Saving functions for package Hipathia
##
## Written by Marta R. Hidalgo, Jose Carbonell-Caballero
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##

##################################
# Save results
##################################

#' Save results to folder
#'
#' Saves results to a folder. In particular, it saves the matrix of subpathway
#' values, a table with the results of the provided comparison,
#' the accuracy of the results and the .SIF and attributes of the pathways.
#'
#' @param results Results object as returned by the \code{hipathia} function.
#' @param comp Comparison as returned by the \code{do_wilcoxon} function.
#' @param metaginfo Pathways object
#' @param output_folder Name of the folder in which the results will be stored.
#' @param path Absolute path to the parent directory in which `output_folder` 
#' will be saved. If it is not provided, it will be created in a temp folder.
#'
#' @return Creates a folder in disk in which all the information to browse the
#' pathway results is stored.
#'
#' @examples
#' data(results)
#' data(comp)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' save_results(results, comp, pathways, "output_results")
#'
#' @export
#'
save_results <- function(results, comp, metaginfo, output_folder = NULL, 
                         path = NULL){

    if(is.null(path))
        path <- tempdir()
    if(is.null(output_folder)){
        n <- length(list.files(path, pattern = "hipathia_results")) + 1
        output_folder <- paste0("hipathia_results_", n)
    }
    output_folder <- paste0(path, "/", output_folder)
    if(!file.exists(output_folder))
        dir.create(output_folder)

        # Write files
    utils::write.table(results$all$path.vals,
                       file = paste0(output_folder,"/all_path_vals.txt"),
                       col.names = TRUE,
                       row.names = TRUE,
                       quote = FALSE,
                       sep="\t")
    comp$path.name <- get_path_names(metaginfo, rownames(comp))
    utils::write.table(comp,
                       file = paste0(output_folder,"/all_path_stats.txt"),
                       col.names = TRUE,
                       row.names = TRUE,
                       quote = FALSE,
                       sep = "\t")

    if(!is.null( results$all$accuracy )){
        accu <- c(results$all$accuracy$total, results$all$accuracy$percent,
                  results$all$accuracy$by.path)
        names(accu) <- c("Accuracy", "Percent", names(accu)[3:length(accu)])
        utils::write.table(accu,
                           file = paste0(output_folder,"/accuracy.txt"),
                           col.names = TRUE,
                           row.names = TRUE,
                           quote = FALSE,
                           sep="\t")
    }
    return(output_folder)
}


write_attributes <- function(this_comp, pathway, metaginfo, prefix,
                             moreatts_pathway = NULL, conf = 0.05,
                             reverse_xref = NULL, exp = NULL){
    atts <- create_node_and_edge_attributes(this_comp, pathway, metaginfo,
                                            moreatts_pathway = moreatts_pathway,
                                            conf = conf,
                                            reverse_xref = reverse_xref,
                                            exp = exp)
    utils::write.table(atts$sif, file = paste0(prefix, ".sif"),
                       row.names = FALSE,
                       col.names = FALSE,
                       quote = FALSE,
                       sep = "\t")
    utils::write.table(atts$node_att, file = paste0(prefix, ".natt"),
                       row.names = FALSE,
                       col.names = TRUE,
                       quote = FALSE,
                       sep = "\t")
    utils::write.table(atts$edge_att, file = paste0(prefix, ".eatt"),
                       row.names = FALSE,
                       col.names = TRUE,
                       quote = FALSE,
                       sep = "\t")
}


create_node_and_edge_attributes <- function(comp, pathway, metaginfo,
                                            moreatts_pathway = NULL, conf=0.05,
                                            reverse_xref = NULL, exp = NULL){


    pathigraphs <- metaginfo$pathigraphs
    effector <- length(unlist(strsplit(rownames(comp)[1], split="-"))) == 3
    if(effector == TRUE){
        s <- pathigraphs[[pathway]]$effector.subgraphs
    }else{
        s <- pathigraphs[[pathway]]$subgraphs
    }
    ig <- pathigraphs[[pathway]]$graph

    if(is.null(V(ig)$type)) V(ig)$type <- "node"
    if(is.null(V(ig)$width)) V(ig)$width <- 15
    if(is.null(V(ig)$height)) V(ig)$height <- 5
    if(is.null(V(ig)$label.color)) V(ig)$label.color <- "black"
    if(is.null(V(ig)$label.cex)) V(ig)$label.cex <- 0.7
    # V(ig)$stroke.color <- find_node_colors(comp, s, ig, conf)[V(ig)$name]
    V(ig)$stroke.color <- "lightgrey"
    #V(ig)$stroke.color[grepl("func",V(ig)$name)] <- "white"
    V(ig)$stroke.size <- 2
    #V(ig)$stroke.size[grepl("func",V(ig)$name)] <- 0

    V(ig)$color <- "white"
    V(ig)$width[V(ig)$shape=="circle"] <- 15
    V(ig)$width[V(ig)$shape!="circle"] <- 22
    V(ig)$shape[V(ig)$shape=="rectangle" & !grepl("func", V(ig)$name)] <-
        "ellipse"
    V(ig)$shape[V(ig)$shape=="rectangle" & grepl("func", V(ig)$name)] <-
        "rectangle"
    V(ig)$width[grepl("func",V(ig)$name)] <- -1

    V(ig)$tooltip <- sapply(seq_along(V(ig)), function(i){
        if(V(ig)$shape[i] == "ellipse"){
            paste(sapply(V(ig)$genesList[[i]], function(gen){
                if(!gen == "/"){
                    paste0("<a target='_blank' ",
                           "href='http://www.genome.jp/dbget-bin/www_bget?",
                           metaginfo$species, ":", gen, "'>", gen, "</a>")
                }else{
                    ""
                }
            }), collapse="<br>")
        }else if(V(ig)$shape[i] == "rectangle"){
            if(grepl("\n", V(ig)$label[i])){
                gsub("\n", "<br>", V(ig)$label[i])
            }else{""}
        }else if(V(ig)$shape[i] == "circle"){
            paste0("<a target='_blank' ",
                   "href='http://www.genome.jp/dbget-bin/www_bget?",
                   V(ig)$label[i], "'>", V(ig)$label[i], "</a>")
        }
    })

    natt <- cbind(V(ig)$name,
                  V(ig)$label,
                  10,
                  V(ig)$nodeX,
                  V(ig)$nodeY,
                  V(ig)$color,
                  V(ig)$stroke.color,
                  V(ig)$stroke.size,
                  V(ig)$shape,
                  V(ig)$type,
                  V(ig)$label.cex,
                  V(ig)$label.color,
                  V(ig)$width,
                  sapply(V(ig)$genesList, paste, collapse=","),
                  V(ig)$tooltip)
    colnames(natt) <- c("ID",
                        "label",
                        "labelSize",
                        "X",
                        "Y",
                        "color",
                        "strokeColor",
                        "strokeSize",
                        "shape",
                        "type",
                        "labelCex",
                        "labelColor",
                        "size",
                        "genesList",
                        "tooltip")
    rownames(natt) <- natt[,1]
    natt[,"label"] <- sapply(natt[,"label"], function(x){
        ul <- unlist(strsplit(x, split="\n"))
        if(length(ul) > 1){
            paste0(ul[1], ", ...")
        }else{
            ul[1]
        }
    })
    # Add
    if(!is.null(moreatts_pathway)){
        common_col_idx <- colnames(moreatts_pathway) %in% colnames(natt)
        common_col <- colnames(moreatts_pathway)[common_col_idx]
        not_common_col_idx <- !colnames(moreatts_pathway) %in% colnames(natt)
        not_common_col <- colnames(moreatts_pathway)[not_common_col_idx]
        for(col in common_col)
            natt[,col] <- moreatts_pathway[,col]
        if(!"strokeColor" %in% common_col){
            natt[,"strokeColor"] <- natt[,"color"]
            natt[natt[,"strokeColor"] == "white","strokeColor"] <- "lightgrey"
        }
        natt <- cbind(natt, moreatts_pathway[,not_common_col])
    }
    node_path_assoc <- matrix(0, nrow = nrow(natt), ncol = length(s))
    colnames(node_path_assoc) <- names(s)
    natt <- cbind(natt, node_path_assoc)

    sif <- c()
    eatt <- c()
    epath_assoc <- c()

    for(i in seq_along(s)){

        # get subgraph
        subgraph <- s[[i]]
        name <- names(s)[i]
        pname <- get_path_names(metaginfo, name)

        # sif
        raw_edges <- get.edgelist(subgraph)
        type <- c("activation","inhibition")[(E(subgraph)$relation == -1) + 1]
        edges <- cbind(raw_edges[,1], type, raw_edges[,2])
        sif <- rbind(sif, edges)

        # edge attributes
        eids <- apply(edges, 1, function(x) paste0(x, collapse = "_"))
        status <- comp[name,"UP/DOWN"]
        if("color" %in% colnames(comp)){
            color <- comp[name,"color"]
        } else {
            if( comp[name,"FDRp.value"] < conf){
                color <- c("#1f78b4","#e31a1c")[(status == "UP") + 1]
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
                            pname = pname,
                            pvalue = comp[name,"p.value"],
                            adj.pvalue = comp[name,"FDRp.value"])
        eatt <- rbind(eatt, edges_atts)

        epath_assoc <- rbind(epath_assoc, path_assoc)

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

        # up regulated
        upsig <- which(subeatt[,"status"] == "UP" &
                           as.numeric(subeatt[,"adj.pvalue"]) < conf)
        if(length(upsig) > 0){

            selected_subsif <- subsif[1,]
            selected_subsif[2] <- paste0(selected_subsif[2], ".up")
            def_sif <- rbind(def_sif, selected_subsif)

            mini_subeatt <- subeatt[upsig,,drop = FALSE]
            selected_subeatt <- mini_subeatt[1, c("id", "status", "color",
                                                  "pvalue", "adj.pvalue")]
            selected_subeatt["id"] <- paste(selected_subsif, collapse = "_")
            def_eatt <- rbind(def_eatt, selected_subeatt)

            selected_subepath_assoc <- subepath_assoc[upsig,,drop = FALSE]
            def_epath_assoc <- rbind(def_epath_assoc,
                                     colSums(selected_subepath_assoc) > 0)
        }

        # down regulated
        downsig <- which(subeatt[,"status"] == "DOWN" &
                             as.numeric(subeatt[,"adj.pvalue"]) < conf)
        if(length(downsig) > 0){

            selected_subsif <- subsif[1,]
            selected_subsif[2] <- paste0(selected_subsif[2], ".down")
            def_sif <- rbind(def_sif, selected_subsif)

            mini_subeatt <- subeatt[downsig,,drop = FALSE]
            selected_subeatt <- mini_subeatt[1,c("id",
                                                 "status",
                                                 "color",
                                                 "pvalue",
                                                 "adj.pvalue")]
            selected_subeatt["id"] <- paste(selected_subsif, collapse = "_")
            def_eatt <- rbind(def_eatt, selected_subeatt)

            selected_subepath_assoc <- subepath_assoc[downsig,,drop = FALSE]
            def_epath_assoc <- rbind(def_epath_assoc,
                                     colSums(selected_subepath_assoc) > 0)
        }

        # no sigs
        nosigs <- which(as.numeric(subeatt[,"adj.pvalue"]) >= conf)
        if(length(nosigs) > 0){

            selected_subsif <- subsif[1,]
            def_sif <- rbind(def_sif, selected_subsif)

            mini_subeatt <- subeatt[nosigs,,drop = FALSE]
            selected_subeatt <- mini_subeatt[1,c("id",
                                                 "status",
                                                 "color",
                                                 "pvalue",
                                                 "adj.pvalue")]
            def_eatt <- rbind(def_eatt, selected_subeatt)

            selected_subepath_assoc <- subepath_assoc[nosigs,,drop = FALSE]
            def_epath_assoc <- rbind(def_epath_assoc,
                                     colSums(selected_subepath_assoc) > 0)
        }
    }

    rownames(def_eatt) <- NULL
    def_eatt <- as.data.frame(def_eatt, stringsAsFactors = FALSE)
    def_epath_assoc <- as.data.frame(def_epath_assoc, stringsAsFactors = FALSE)
    rownames(def_sif) <- NULL
    def_sif <- as.data.frame(def_sif, stringsAsFactors = FALSE)

    def_eatt$shape <- c("inhibited", "directed")[grepl("activation",
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
                paste(x[1], "activation", x[2], sep = "_")
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
        funejes$shape <- "directed"
        nods <- get.edgelist(ig)[left,1]
        names(nods) <- ids
        names(ids) <- nods
        funs <- t(apply(funejes, 1, function(x){
            lastnodes <- sapply(colnames(funejes), get_effnode_id)
            if(any(lastnodes == nods[x[[1]]])){
                x[lastnodes == nods[x[[1]]]] <- 1
                x
            }else{
                x
            }
        }))
        funs <- as.data.frame(funs, stringsAsFactors = FALSE)
        sif_funs <- data.frame(V1 = get.edgelist(ig)[left,1],
                               type = rep("activation", times = length(left)),
                               V3 = get.edgelist(ig)[left,2],
                               stringsAsFactors = FALSE)

        def_sif <- rbind(def_sif, sif_funs)
        def_eatt <- rbind(def_eatt, funs)
    }

    fun_indexes <- grep("_func", rownames(natt))
    fun_names <- rownames(natt)[fun_indexes]
    if(length(fun_indexes) > 0){
        for(i in seq_along(fun_names)){
            pp <- gsub("N", "P", gsub("_func", "", fun_names[i]))
            if(effector == TRUE){
                natt[fun_names[i], pp] <- 1
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
                tids <- sapply(reverse_xref[ids],function(x){
                    if(is.null(x)){
                        return("?")
                    } else {
                        return(x)
                    }})
                return(paste(tids, collapse = ","))
            } else {
                return("?")
            }
        }
        natt <- cbind(natt, tids = sapply(sids, translate_ids))
    }
    if(!is.null(exp)){
        sids <- strsplit(as.character(natt[,"genesList"]), split = ",")
        ids_list <- as.list(seq_len(nrow(exp)))
        names(ids_list) <- rownames(exp)
        get_expr_ids <- function(ids){
            if(length(ids) > 0){
                ids <- setdiff(ids, "/")
                exp_values <- sapply(ids_list[ids],function(x){
                    if(is.null(x)){
                        return("?")
                    }else{
                        return(exp[x,])
                    }})
                return(paste(exp_values, collapse = ","))
            } else {
                return("?")
            }
        }
        natt <- cbind(natt, exp_values = sapply(sids, get_expr_ids))
    }

    return(list(sif = def_sif,
                edge_att = def_eatt,
                node_att = natt))
}



create_path_info <- function(all_comp, metaginfo){
    fpgs <- metaginfo$pathigraphs
    effector <- length(unlist(strsplit(rownames(all_comp)[1], split="-"))) == 3
    path_info <- lapply(fpgs, function(fpg){
        if(effector == TRUE){
            all_comp[names(fpg$effector.subgraphs),]        
        }else{
            all_comp[names(fpg$subgraphs),]        
        }
    })
    
    path_json_list <- lapply(names(path_info),function(x){
        out <- paste0("{\n\t\"id\":\"", x, "\",\n")
        out <- paste0(out, "\t\"name\":\"", fpgs[[x]]$path.name, "\",\n")
        anysig <- FALSE
        anyup <- FALSE
        anydown <- FALSE
        anysigup <- FALSE
        anysigdown <- FALSE
        anychanged <- FALSE
        for(i in seq_len(nrow(path_info[[x]]))){
            if(path_info[[x]]$has_changed[i] == TRUE)
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
        for(i in seq_len(nrow(path_info[[x]]))){
            out <- paste0(out, "\t\t{")
            out <- paste0(out, "\"id\":\"", rownames(path_info[[x]])[i], "\", ")
            out <- paste0(out, "\"name\":\"",
                          get_path_names(metaginfo,
                                         rownames(path_info[[x]])[i]), "\", ")
            if(grepl("term_", metaginfo$pathigraphs[[1]]$path.id) == TRUE){
                out <- paste0(out, "\"shortname\":\"",
                              get_path_names(metaginfo,
                                             rownames(path_info[[x]])[i]),
                              "\", ")
            }else{
                out <- paste0(out, "\"shortname\":\"" ,
                              gsub("\\*", "", strsplit(get_path_names(
                                  metaginfo,
                                  rownames(path_info[[x]])[i]),": ")[[1]][2]),
                              "\", ")
            }
            out <- paste0(out, "\"pvalue\":", path_info[[x]]$FDRp.value[i],
                          ", ")
            out <- paste0(out, "\"status\":\"", path_info[[x]]$status[i],
                          "\", ")
            out <- paste0(out, "\"sig\":\"",
                          tolower(path_info[[x]]$FDRp.value[i] < 0.05), "\", ")
            out <- paste0(out, "\"haschanged\":",
                          tolower(path_info[[x]]$has_changed[i]), ", ")
            out <- paste0(out, "\"up\":",
                          tolower(path_info[[x]]$status[i] == "UP"), ", ")
            out <- paste0(out, "\"down\":",
                          tolower(path_info[[x]]$status[i] == "DOWN"), ", ")
            out <- paste0(out, "\"upsig\":",
                          tolower(path_info[[x]]$status[i] == "UP" &
                                      path_info[[x]]$FDRp.value[i] < 0.05),
                          ", ")
            out <- paste0(out, "\"downsig\":",
                          tolower(path_info[[x]]$status[i] == "DOWN" &
                                      path_info[[x]]$FDRp.value[i] < 0.05),
                          ", ")
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
    path_json <- paste0("[\n", paste(path_json_list, collapse = ","), "\n]")
    return(path_json)
}


create_report_folders <- function(output_folder, home, clean_out_folder = TRUE){

    pv_folder <- paste0(output_folder,"/pathway-viewer")

    if(clean_out_folder == TRUE & file.exists(pv_folder)){
        unlink(pv_folder, recursive = TRUE)
        unlink(paste0(output_folder, "/index.html"), recursive = TRUE)
    }
    file.copy(paste0(home,"/pathway-viewer/"), output_folder, recursive = TRUE)
    report_path <- paste0(home, "/report-files/")
    png_files_copy <- list.files(path = report_path, pattern = ".png")
    png_files_copy <- paste0(home, "/report-files/", png_files_copy)
    file.copy(png_files_copy, pv_folder)

}

create_pathways_folder <- function(output_folder, metaginfo, comp, moreatts,
                                   conf, verbose = FALSE){

    pathways_folder <- paste0(output_folder, "/pathway-viewer/pathways/")
    if(!file.exists(pathways_folder))
        dir.create(pathways_folder)
    for(pathway in names(metaginfo$pathigraphs)){
        if(verbose == TRUE)
            cat(pathway)
        write_attributes(comp,
                         pathway,
                         metaginfo,
                         paste0(pathways_folder, pathway),
                         moreatts_pathway = moreatts[[pathway]],
                         conf = conf)
    }

    comp$status <- comp$"UP/DOWN"
    comp$has_changed <- TRUE
    path_json <- create_path_info(comp, metaginfo)
    write(path_json, file = paste0(output_folder,
                                   "/pathway-viewer/pathways/path_info.json"))

}


create_html_index <- function(home, output_folder,
                                template_name = "index_template.html",
                                output_name = "index.html"){


    index <- scan(paste0(home,'/report-files/',template_name),
                  comment.char = "", sep = "\n", what = "character",
                  quiet = TRUE)

    global_div <- c()

    global_div <- c(global_div, paste0("<pathway-viewer id='pathway-viewer'",
                                       " path-type='url' path='pathways'>",
                                       "</pathway-viewer>"))

    new_index <- gsub("PUT_HERE_YOUR_ELEMENTS",
                      paste(global_div, collapse = "\n"),
                      index)

    write(paste(new_index, collapse = "\n"),
          file = paste0(output_folder,"/pathway-viewer/",output_name))
}



#' Create visualization HTML
#'
#' Saves the results of a Wilcoxon comparison for the Hipathia pathway values
#' into a folder, and creates a HTML from which to visualize the results on
#' top of the pathways. The results are stored into the specified folder.
#' If this folder does not exist, it will be created. The parent folder must
#' exist.
#'
#' @examples
#' data(results)
#' data(comp)
#' data(brca)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' create_report(comp, pathways, "save_results")
#' 
#' \dontrun{
#' sample_group <- colData(brca)[,1]
#' colors_de <- node_color_per_de(results, pathways,
#' sample_group, "Tumor", "Normal")
#' create_report(comp, pathways, "save_results",
#' node_colors = colors_de)
#'}
#'
#' @param comp Comparison object as given by the \code{do_wilcoxon} function
#' @param metaginfo Pathways object as returned by the \code{load_pathways}
#' function
#' @param output_folder Name of the folder in which the report will be stored.
#' @param path Absolute path to the parent directory in which `output_folder` 
#' will be saved. If it is not provided, it will be created in a temp folder.
#' @param node_colors List of colors with which to paint the nodes of the
#' pathways, as returned by the
#' \code{node_color_per_de} function. Default is white.
#' @param group_by How to group the subpathways to be visualized. By default
#' they are grouped by the pathway to which they belong. Available groupings
#' include "uniprot", to group subpathways by their annotated Uniprot functions,
#' "GO", to group subpathways by their annotated GO terms, and "genes", to group
#' subpathways by the genes they include. Default is set to "pathway".
#' @param conf Level of significance. By default 0.05.
#' @param verbose Boolean, whether to show details about the results of the
#' execution
#'
#' @return Saves the results and creates a report to visualize them through
#' a server in the specified \code{output_folder}. Returns the folder where
#' the report has been stored.
#'
#' @export
#'
create_report <- function(comp, metaginfo, output_folder = NULL, path = NULL,
                          node_colors = NULL,
                          group_by = "pathway", conf = 0.05, verbose = FALSE){

    if(group_by != "pathway" &
       length(unlist(strsplit(rownames(comp)[1], split = "-"))) == 4)
        stop("Grouping only available for effector subgraphs")

    if(!is.null(node_colors)){
        if(node_colors$group_by != group_by)
            stop("Grouping in node.colors must agree with group_by")
        moreatts <- summarize_atts(list(node_colors$colors), c("color"))
    }else{
        moreatts <- NULL
    }

    if(group_by != "pathway"){
        message("Creating groupings by ", group_by, "...")
        metaginfo <- get_pseudo_metaginfo(metaginfo, group_by = group_by)
    }

    if(is.null(path))
        path <- tempdir()
    if(is.null(output_folder)){
        n <- length(list.files(path, pattern = "hipathia_report")) + 1
        output_folder <- paste0("hipathia_report_", n)
    }
    output_folder <- paste0(path, "/", output_folder)
    if(!file.exists(output_folder))
        dir.create(output_folder)
    
    pv_path <- paste0(system.file("extdata", package="hipathia"))

    message("Creating report folders...")
    create_report_folders(output_folder, pv_path, clean_out_folder = FALSE)

    message("Creating pathways folder...")
    create_pathways_folder(output_folder, metaginfo, comp, moreatts, conf,
                           verbose)

    message("Creating HTML index...")
    create_html_index(pv_path,
                      output_folder,
                      template_name = "index_template.html",
                      output_name = "index.html")

    return(output_folder)
}


summarize_atts <- function(att_list, att_names){
    df_list <- c()
    for(pathway in names(att_list[[1]])){
        df <- sapply(att_list, function(l){l[[pathway]]})
        colnames(df) <- att_names
        df_list[[pathway]] <- df
    }
    return(df_list)
}

#'
#' Visualize a HiPathia report
#'
#' @param output_folder Folder in which results to visualize are stored
#' @param port Port to use
#'
#' @return The instructions to visualize a HiPathia report in a web browser
#'
#' @examples
#' data(comp)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' report <- create_report(comp, pathways, "save_results")
#' visualize_report(report)
#' 
#' \dontrun{
#' data(results)
#' data(brca)
#' sample_group <- colData(brca)[,1]
#' colors_de <- node_color_per_de(results, pathways,
#' sample_group, "Tumor", "Normal")
#' report <- create_report(comp, pathways, "save_results",
#' node_colors = colors_de)
#' visualize_report(report)
#' visualize_report(report, port = 5000)
#' }
#' \dontshow{servr::daemon_stop()}
#'
#' @import servr
#' @export
#'
visualize_report <- function(output_folder, port = 4000){
    servr::httd(paste0(output_folder, "/pathway-viewer"),
                port = port, browser = FALSE, daemon = TRUE)
    cat("Open a web browser and go to URL http://127.0.0.1:", port, "\n")
}




###########################################

# PSEUDO META_GRAPH_INFORMATION

get_pseudo_metaginfo <- function(pathways, group_by){
    pseudo <- load_pseudo_mgi(pathways$species, group_by)
    rownames(pseudo$all.labelids) <- pseudo$all.labelids[,1]
    pathways_list <- names(pathways$pathigraphs)
    if(!all(unique(pseudo$all.labelids[,"path.id"]) %in% pathways_list))
        pseudo <- filter_pseudo_mgi(pseudo, pathways_list)
    return(pseudo)
}

filter_pseudo_mgi <- function(pseudo_meta, pathways_list){
    num_nodes <- sapply(names(pseudo_meta$pathigraphs), function(term){
        graph <- pseudo_meta$pathigraphs[[term]]$graph
        idx <- unlist(lapply(pathways_list, grep, V(graph)$name))
        vs <- V(graph)[idx]
        length(vs)
    })
    tofilter <- names(pseudo_meta$pathigraphs)[num_nodes >= 1]
    mini_pathigraphs <- lapply(pseudo_meta$pathigraphs[tofilter],
                               function(pg){
        minipg <- NULL
        graph <- pg$graph
        idx <- unlist(lapply(pathways_list, grep, V(graph)$name))
        vs <- V(graph)[idx]
        minipg$graph <- igraph::induced_subgraph(graph, vs)
        minipg$path.name <- pg$path.name
        minipg$path.id <- pg$path.id
        es_ind <- unlist(lapply(pathways_list, grep, pg$effector.subgraphs))
        minipg$effector.subgraphs <- pg$effector.subgraphs[es_ind]
        minipg
                               })
    names(mini_pathigraphs) <- tofilter

    all_labels <- pseudo_meta$all.labelids
    lab_in_pl <- all_labels[,"path.id"] %in% pathways_list
    filter_labelids <- all_labels[lab_in_pl,]

    mini_pseudo <- NULL
    mini_pseudo$pathigraphs <- mini_pathigraphs
    mini_pseudo$species <- pseudo_meta$species
    mini_pseudo$all.labelids <- filter_labelids

    return(mini_pseudo)
}

