##
## chart.R
## Plot functions
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##

#' Plots subpathways heatmap
#'
#' Plots a heatmap with the values of the subpathways.
#'
#' @param data Either a SummarizedExperiment or a matrix with the values to be 
#' plotted. Rows are features and columns are samples.
#' @param sel_assay Character or integer, indicating the assay to be normalized 
#' in the SummarizedExperiment. Default is 1.
#' @param group Either a character indicating the name of the column in 
#' colData 
#' including the classes to plot, or a character vector with the class to 
#' which each sample belongs. Samples must be ordered as in \code{data}. 
#' By default, all samples will be assigned to the same class.
#' @param colors Either a character vector with colors or a key name
#' indicating the color scheme to be used in the heatmap.
#' If a character vector is provided, it is recommended to provide at
#' least 3 colors. Three different predefined color schemes may be
#' selected by providing a key name. Options are:
#' * \code{classic} Blue for lower values, white for medium values,
#' red for higher values.
#' * \code{hipathia} Hipathia predefined color scheme: Green for lower
#' values, white for medium values, orange for higher values.
#' * \code{redgreen} Green for lower values, black for medium values,
#' red for higher values.
#' By default \code{classic} color scheme is applied.
#' @param sample_clust Boolean, whether to cluster samples (columns).
#' By default TRUE.
#' @param variable_clust Boolean, whether to cluster variables (rows).
#' By default FALSE. If TRUE, rows with 0 variance are removed.
#' @param labRow,labCol Character vectors with row and column labels
#' to be used. By default rownames(data) or colnames(data)
#' are used, respectively.
#' @param sample_colors Named character vector of colors. The names of
#' the colors must be the classes in \code{group}. Each sample will
#' be assigned the color corresponding to its class, taken from the
#' \code{group} vector. By default a color will be assigned
#' automatically to each class.
#' @param scale Boolean, whether to scale each row to the interval [0,1].
#' Default is TRUE.
#' @param save_png Path to the file where the image as PNG will be saved.
#' By default, the image is not saved.
#' @param legend Boolean, whether to display a legend.
#' @param legend_xy Position for the legend, in case \code{legend} is TRUE.
#' @param pch Graphical parameter from \code{par()} function.
#' @param main Main title of the image
#'
#' @return Heatmap of the values of the subpathways
#'
#' @examples
#' data(brca_design)
#' data(path_vals)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' heatmap_plot(path_vals, group = sample_group)
#' heatmap_plot(path_vals, group = sample_group, colors = "hipathia",
#' variable_clust = TRUE)
#'
#' @export
#' @import grDevices graphics
#' @importFrom matrixStats rowVars
#' @importFrom DelayedArray rowMins
#' @importFrom DelayedArray rowMaxs
#' @importFrom stats heatmap
#' @importFrom methods is
#'
heatmap_plot <- function(data, group = NULL, sel_assay = 1, 
                         colors = "classic",
                         sample_clust = TRUE, variable_clust = FALSE,
                         labRow = NULL, labCol = NULL, sample_colors = NULL,
                         scale = TRUE, save_png = NULL, legend = TRUE,
                         legend_xy = "topright", pch = 15, main = NULL){
    if(is(data, "SummarizedExperiment")){
        se <- TRUE
        if(is(group, "character") & length(group) == 1)
            if(group %in% colnames(colData(data))){
                group <- colData(data)[[group]]
            }else{
                stop("group variable must be a column in colData(data)")
            }
        vals <- assay(data)
    }else if(is(data, "matrix")){
        vals <- data
    }else{
        stop("Only SummarizedExperiment or matrix classes accepted as data")
    }
    if(length(colors) == 1){
        if(colors == "hipathia"){
            colors <- c("#007462", "white", "#e66430")
        }else if(colors == "classic"){
            colors <- c("blue","gray","red")
        }else if(colors == "redgreen"){
            colors <- c("green","black","red")}
    }
    if(is.null(group))
        group <- rep("A", ncol(vals))
    if(sample_clust==FALSE){
        colv <- NA
    } else {
        colv <- TRUE
    }
    if(variable_clust==FALSE){
        rowv <- NA
    } else {
        vars <- matrixStats::rowVars(vals)
        vals <- vals[!is.na(vars) & vars != 0,]
        rowv <- TRUE
    }
    if(is.null(labRow)){
        if(nrow(vals) < 50){
            labRow <- rownames(vals)
        }else{
            labRow <- FALSE
        }
    }
    if(is.null(labCol)){
        if(ncol(vals) < 50){
            labCol <- colnames(vals)
        }else{
            labCol <- FALSE
        }
    }
    if(is.null(sample_colors)){
        if(length(unique(group)) <= 8){
            sample_colors <- c("#50b7ae", "#b6ebe7", "#e66430",
                               "#305f59", "#ffc868", "#152e2b",
                               "#a0170e", "#f9b493")[seq_along(unique(group))]
        }else{
            sample_colors <- c("#50b7ae", "#b6ebe7", "#e66430", "#305f59",
                               "#ffc868", "#152e2b", "#a0170e", "#f9b493",
                               grDevices::terrain.colors(length(unique(
                                   group)) - 8))
        }
        names(sample_colors) <- unique(group)
    }
    if(scale == TRUE){
        min <- rowMins(vals, na.rm = TRUE)
        max <- rowMaxs(vals, na.rm = TRUE)
        vals <- (vals - min)/(max - min)
    }
    if(!is.null(save_png))
        grDevices::png(filename = save_png)
    if(!is.null(main))
        graphics::par(oma=c(1,0,3,1))
    stats::heatmap(vals,
                   margins = c(10,10),
                   labRow = labRow,
                   labCol = labCol,
                   scale = "none",
                   Rowv = rowv,
                   Colv = colv,
                   ColSideColors = sample_colors[group],
                   col = grDevices::colorRampPalette(colors)(256))
    if(legend == TRUE)
        legend(legend_xy,
               legend = unique(group),
               col = sample_colors[unique(group)],
               pch = pch,
               xpd = TRUE,
               cex = 1,
               border = 0)
    if(!is.null(main)){
        title(main=main, outer = TRUE)
        graphics::par(oma = c(0,0,0,0))
    }
    if(!is.null(save_png))
        grDevices::dev.off()
}


#'
#' Plots two components of a PCA
#'
#' Plots two components of a PCA computed with \code{do_pca}
#'
#' @param fit princomp object as returned by \code{do_pca}
#' @param group Vector with the group to which each sample belongs.
#' The samples must be ordered as in \code{rownames(fit$scores)}.
#' By default, all samples will be assigned to the same class.
#' @param sample_colors Named character vector of colors. The names of
#' the colors must be the classes in \code{group}. Each sample will be
#' assigned the color corresponding to its class, taken from the
#' \code{group} vector. By default a color will be assigned
#' automatically to each class.
#' @param cp1 Integer, number of the component in the X-axis.
#' Default is 1, the first component.
#' @param cp2 Integer, number of the component in the Y-axis.
#' Default is 2, the second component.
#' @param legend Boolean, whether to plot a legend in the plot.
#' Default is TRUE.
#' @param legend_xy Situation of the legend in the plot. Available
#' options are: "bottomright", "bottom", "bottomleft", "left",
#' "topleft", "top", "topright", "right" and "center".
#' @param cex Graphical parameter from \code{par()} function.
#' @param pch Graphical parameter from \code{par()} function.
#' @param mgp Graphical parameter from \code{par()} function.
#' @param main Title of the graphics
#' @param save_png Path to the file where the image as PNG will be saved.
#' By default, the image is not saved.
#'
#' @return Plots two components of a PCA
#'
#' @examples
#' data(path_vals)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' pca_model <- do_pca(path_vals[seq_len(ncol(path_vals)),])
#' pca_plot(pca_model, sample_group)
#'
#' @export
#' @import grDevices graphics
#'
pca_plot <- function(fit, group = NULL, sample_colors = NULL, cp1 = 1,
                     cp2 = 2, legend = TRUE, legend_xy = "bottomleft", cex = 2,
                     pch = 20, mgp = c(3,1,0), main = "PCA plot",
                     save_png = NULL){
    if(is.null(group)) group <- rep("A", fit$n.obs)
    if(is.null(sample_colors)){
        if(length(unique(group)) <= 8){
            sample_colors <- c("#50b7ae", "#b6ebe7", "#e66430",
                               "#305f59", "#ffc868", "#152e2b",
                               "#a0170e", "#f9b493")[seq_along(unique(group))]
        }else{
            sample_colors <- c("#50b7ae", "#b6ebe7", "#e66430", "#305f59",
                               "#ffc868", "#152e2b", "#a0170e", "#f9b493",
                               grDevices::topo.colors(length(unique(
                                   group)) - 8))
        }
        names(sample_colors) <- unique(group)
    }
    cpv1 <- fit$scores[,cp1]
    cpv2 <- fit$scores[,cp2]
    if(!is.null(save_png))
        grDevices::png(filename = save_png)
    graphics::plot(cpv1,
                   cpv2,
                   xlab = paste("PC", cp1),
                   ylab = paste("PC", cp2),
                   col = sample_colors[group],
                   pch = pch,
                   cex = cex,
                   main = main,
                   mgp = mgp)
    if(legend == TRUE){
        legend(legend_xy,
               legend = unique(group),
               col = sample_colors[unique(group)],
               pch = pch,
               xpd = TRUE,
               cex = 1,
               border = 0)
    }
    if(!is.null(save_png))
        grDevices::dev.off()
}


#'
#' Plots multiple components of a PCA
#'
#' Plots multiple components of a PCA analysis computed with \code{do_pca}
#'
#' @param fit princomp object as returned by \code{do_pca}
#' @param group Vector with the group to which each sample belongs.
#' The samples must be ordered as in \code{path_vals}.
#' By default, all samples will be assigned to the same class.
#' @param sample_colors Named character vector of colors. The names of the
#' colors must be the classes in \code{group}. Each sample will be
#' assigned the color corresponding to its class, taken from the
#' \code{group} vector. By default a color will be assigned
#' automatically to each class.
#' @param comps Vector with the components to be plot
#' @param plot_variance Logical, whether to plot the cumulative variance.
#' @param legend Boolean, whether to plot a legend in the plot.
#' Default is TRUE.
#' @param cex Graphical parameter from \code{par()} function.
#' @param pch Graphical parameter from \code{par()} function.
#' @param main Main title of the image
#' @param save_png Path to the file where the image as PNG will be saved.
#' By default, the image is not saved.
#'
#' @return Plots multiple components of a PCA
#'
#' @examples
#' data(path_vals)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' pca_model <- do_pca(path_vals[seq_len(ncol(path_vals)),])
#' multiple_pca_plot(pca_model, sample_group, cex = 3, plot_variance = TRUE)
#'
#' @export
#' @import graphics grDevices
#'
multiple_pca_plot <- function(fit, group = NULL, sample_colors = NULL,
                              comps = seq_len(3), plot_variance = FALSE, 
                              legend = TRUE,
                              cex = 2, pch = 20, main = "Multiple PCA plot",
                              save_png = NULL){
    combs <- utils::combn(comps, 2)
    ncombs <- ncol(combs)
    nn <- ncombs
    if(!is.null(legend))
        nn <- nn + 1
    if(plot_variance==TRUE)
        nn <- nn + 1

    nr <- floor(sqrt(nn))
    nc <- ceiling((nn)/nr)
    oldmfrow <- par("mfrow")
    graphics::par(mfrow=c(nr, nc))
    oldmar <- par("mar")
    graphics::par(mar=c(4, 4, 1, 1))
    oldoma <- par("oma")
    graphics::par(oma=c(0, 0, 2, 0))
    if(!is.null(save_png))
        grDevices::png(filename = save_png)
    for(i in seq_len(ncombs)){
        pca_plot(fit,
                 group = group,
                 sample_colors = sample_colors,
                 cp1 = combs[1,i],
                 cp2 = combs[2,i],
                 cex = cex,
                 pch = pch,
                 main = NULL,
                 legend = FALSE,
                 mgp = c(2.5,1,0))
    }
    if(legend == TRUE){
        if(is.null(sample_colors)){
            ug <- unique(group)
            if(length(unique(group)) <= 8){
                sample_colors <- c("#50b7ae", "#b6ebe7", "#e66430",
                                   "#305f59", "#ffc868", "#152e2b",
                                   "#a0170e", "#f9b493")[seq_along(ug)]
            }else{
                sample_colors <- c("#50b7ae", "#b6ebe7", "#e66430", "#305f59",
                                   "#ffc868", "#152e2b", "#a0170e", "#f9b493",
                                   grDevices::topo.colors(length(ug) - 8))
            }
            names(sample_colors) <- unique(group)
        }
        graphics::plot(1, type="n", axes = FALSE, xlab = "", ylab = "")
        legend("center",
               legend = unique(group),
               col = sample_colors[unique(group)],
               pch = pch,
               lwd = 2,
               xpd = TRUE,
               cex = 1,
               border = NA,
               pt.cex = 1.2)
    }
    if(plot_variance == TRUE){
        plot_pca_variance(fit, acum = TRUE, thresh = 0.1)
    }
    title(main, outer = TRUE)
    graphics::par(mfrow = oldmfrow)
    graphics::par(oma = oldoma)
    graphics::par(mar = oldmar)
    if(!is.null(save_png))
        grDevices::dev.off()
}


plot_pca_variance <- function(fit, thresh = 0, acum = FALSE, minnum = 5){
    if(acum==FALSE){
        comptoplot <- fit$explain_var > thresh
        if(sum(comptoplot) < minnum)
            comptoplot <- seq_len(5)
        graphics::barplot(fit$explain_var[ comptoplot ],
                          ylab = "explain variance",
                          xlab = "",
                          las = 2,
                          cex.names = 0.5,
                          ylim = c(0,1))
    } else {
        comptoplot <- fit$acum_explain_var < (1 - thresh)
        if(sum(comptoplot) < minnum) comptoplot <- seq_len(5)
        graphics::barplot(fit$acum_explain_var[ comptoplot ],
                          ylab = "acum explain variance",
                          xlab = "",
                          las = 2,
                          cex.names = 0.5,
                          ylim = c(0,1))
    }
}



# PLOT RESULTS
##############################

#' @import igraph
plot_pathigraph <- function(g, node_color = NULL, edge_lty = 1, main = "" ){
    V(g)$shape[V(g)$shape == "rectangle"] <- "crectangle"
    V(g)$shape[grepl("_func", V(g)$name)] <- "rectangle"
    V(g)$color[grepl("_func", V(g)$name)] <- "white"
    V(g)$frame.color <- "darkgrey"
    V(g)$frame.color[grepl("_func", V(g)$name)] <- "white"
    V(g)$size <- 10
    V(g)$size[V(g)$shape == "rectangle"] <- 15
    V(g)$size[V(g)$shape == "circle"] <- 5
    V(g)$size[grepl("_func", V(g)$name)] <- 20
    V(g)$size2 <- 5
    if(is.null(node_color)){
        V(g)$color <- "white"}else{V(g)$color <- node_color}
    E(g)$lty <- edge_lty
    plot.igraph(g,
                layout = cbind(V(g)$nodeX, V(g)$nodeY),
                axes = FALSE,
                asp = 0,
                main = main,
                edge.arrow.size = 0.15,
                edge.arrow.width = 1,
                vertex.label.color = "grey15")
}


#'
#' Plots pathway with colored significant paths
#'
#' Plots the layout of a pathway, coloring the significant subpathways
#' in different colors depending on whether they are significantly up- or
#' down-regulated. Nodes may be also colored providing a suitable list of
#' colors for each node. Function \code{node_color_per_de}
#' assigns colors to the nodes depending on their differential expression.
#'
#' @param comp Comparison data frame as returned by the \code{do_wilcox}
#' function.
#' @param metaginfo Pathways object.
#' @param pathway Name of the pathway to be plotted.
#' @param conf Level of significance of the comparison for the adjusted
#' p-value. Default is 0.05.
#' @param node_colors List, named by the pathway name, including the
#' color of each node for each pathway.
#' @param colors Either a character vector with 3 colors (indicating,
#' in this order, down-regulation, non-significance and up-regulation colors)
#' or a key name indicating the color scheme to be used. Options are:
#' @slot classic ColorBrewer blue, white and colorBrewer red.
#' @slot hipathia Hipathia predefined color scheme: Green, white and orange.
#' By default \code{classic} color scheme is applied.
#'
#' @return Image in which a pathway is ploted. Edges are colored so that the
#' UP- and DOWN-activated subpathways are identified.
#'
#' @examples
#' data(comp)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' pathway_comparison_plot(comp, metaginfo = pathways, pathway = "hsa03320")
#'
#' data(results)
#' data(brca)
#' sample_group <- colData(brca)[,1]
#' colors_de <- node_color_per_de(results, pathways,
#' sample_group, "Tumor", "Normal")
#' pathway_comparison_plot(comp, metaginfo = pathways, pathway = "hsa04012",
#' node_colors = colors_de)
#'
#' @export
#'
pathway_comparison_plot <- function(comp, metaginfo, pathway, conf=0.05,
                                    node_colors = NULL, colors = "classic"){

    if(length(colors) == 1){
        if(colors == "hipathia"){ 
            colors <- c("#50b7ae", "darkgrey", "#f16a34")
        }else if(colors == "classic"){ 
            colors <- c("#0571b0", "darkgrey","#ca0020")}
    }
    down_col <- colors[1]
    no_col <- colors[2]
    up_col <- colors[3]

    pathigraph <- metaginfo$pathigraphs[[pathway]]
    if(all(grepl(" - ", rownames(comp)))){
        effector = FALSE
    }else{
        effector = TRUE}

    paths <- sapply(strsplit(rownames(comp), "-"), "[[", 2)
    comp <- comp[paths == pathigraph$path.id, ]

    # Find edge colors
    g <- add_edge_colors(pathigraph,
                         comp,
                         effector,
                         up_col = up_col,
                         down_col = down_col,
                         no_col = no_col)

    edge_lty <- (E(g)$relation*-1 +1)/2 + 1
    id <- metaginfo$pathigraphs[[pathway]]$path.id
    name <- metaginfo$pathigraphs[[pathway]]$path.name
    title <- paste(id, "-", name)
    plot_pathigraph(g,
                    node_color = node_colors$colors[[pathway]],
                    edge_lty = edge_lty,
                    main = title )
}



add_edge_colors <- function(pathigraph, pcomp, effector, up_col = "#ca0020",
                            down_col = "#0571b0", no_col = "darkgrey",
                            conf = 0.05){

    if(effector == TRUE){
        subgraphs <- pathigraph$effector.subgraphs
    }else{
        subgraphs <- pathigraph$subgraphs}

    g <- pathigraph$graph
    elg <- apply(get.edgelist(g), 1, paste, collapse = "_")
    states <- matrix(NA,
                     nrow = nrow(pcomp),
                     ncol= length(elg),
                     dimnames = list(rownames(pcomp), elg))
    for( path_name in rownames(pcomp)){
        if(pcomp[path_name,"FDRp.value"] <= conf){
            up_down <- pcomp[path_name,"UP/DOWN"]
        }else{
            up_down <- "N"
        }
        subgraph <- subgraphs[[path_name]]
        els <- apply(get.edgelist(subgraph), 1, paste, collapse = "_")
        states[path_name, els] <- up_down
    }

    E(g)$color <- rep(no_col, length(E(g)))

    colors <- c(up_col, no_col, down_col)
    names(colors) <- c("UP", "N", "DOWN")

    for(i in which(!grepl("_func", elg))){
        edge_states <- names(table(states[,i]))
        if(is.null(edge_states)){
            warning("Edge", elg[i], "is not present in any subpath")
        }else{
            E(g)$color[i] <- colors[edge_states[1]]
            if(length(edge_states) > 1)
                for(j in 2:length(edge_states)){
                    g <- g + edge(unlist(strsplit(elg[i], split = "\\_")))
                    E(g)$color[length(E(g))] <- colors[edge_states[j]]
                    E(g)$relation[length(E(g))] <- E(g)$relation[i]
                }
        }
    }

    return(g)
}


#'
#' Colors of the nodes by its differential expression
#'
#' Performs a differential expression on the nodes and computes the colors
#' of the nodes depending on it_ Significant up- and down-regulated nodes
#' are depicted with the selected color, with a gradient towards the
#' non-significant color depending on the value of the p-value.
#' Smaller p-values give rise to purer colors than higher p-values.
#'
#' @param results Object of results as provided by the \code{hipathia}
#' function_
#' @param metaginfo Object of pathways_
#' @param group Character indicating the column in which the group variable is 
#' stored, in case the object provided to \code{hipathia} was a 
#' SummarizedExperiment, or a vector with the class to which each sample 
#' belongs. Samples must be ordered as in \code{results}.
#' @param g1 String, label of the first group to be compared
#' @param g2 String, label of the second group to be compared
#' @param group_by How to group the subpathways to be visualized. By default
#' they are grouped by the pathway to which they belong. Available groupings
#' include "uniprot", to group subpathways by their annotated Uniprot functions,
#' "GO", to group subpathways by their annotated GO terms, and "genes", to group
#' subpathways by the genes they include. Default is set to "pathway".
#' @param colors Either a character vector with 3 colors (indicating,
#' in this order, down-regulation, non-significance and up-regulation colors)
#'  or a key name indicating the color scheme to be used. Options are:
#' @slot classic ColorBrewer blue, white and colorBrewer red.
#' @slot hipathia Hipathia predefined color scheme: 
#' Green, white and orange.
#' By default \code{classic} color scheme is applied.
#' @param conf Level of significance of the comparison for the adjusted p-value
#'
#' @return List of color vectors, named by the pathways to which they belong.
#' The color vectors represent the differential expression
#' of the nodes in each pathway.
#'
#' @examples
#' data(results)
#' data(brca)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' sample_group <- colData(brca)[,1]
#' colors_de <- node_color_per_de(results, pathways,
#' sample_group, "Tumor", "Normal")
#'
#' @export
#' @importFrom methods is
#'
node_color_per_de <- function(results, metaginfo, group, g1, g2, 
                              group_by = "pathway", colors = "classic", 
                              conf = 0.05){

    if(length(colors) == 1){
        if(colors == "hipathia"){
            colors <- c("#50b7ae", "white", "#f16a34")
        }else if(colors == "classic"){
            colors <- c("#1f9cda","white","#da1f1f")
        }
    }
    down_col <- colors[1]
    no_col <- colors[2]
    up_col <- colors[3]

    if(group_by != "pathway")
        metaginfo <- get_pseudo_metaginfo(metaginfo, group_by = group_by)

    if(is(group, "character") & length(group) == 1)
        if(group %in% colnames(colData(results[["nodes"]]))){
            group <- colData(results[["nodes"]])[[group]]
        }else{
            stop("Group variable must be a column in colData())")
        }
    difexp <- compute_difexp(assay(results[["nodes"]]), g1, g2, group)
    updown <- rep("both", length(difexp$statistic))
    updown[difexp$statistic < 0] <- "down"
    updown[difexp$statistic > 0] <- "up"
    node_colors <- get_colors_from_pval(updown,
                                        difexp$p.value,
                                        up_col = up_col,
                                        down_col = down_col,
                                        no_col = no_col,
                                        conf = conf)
    names(node_colors) <- rownames(results[["nodes"]])
    cols <- lapply(metaginfo$pathigraphs, function(pg){
        gen_nodes <- V(pg$graph)$name[V(pg$graph)$name %in% rownames(difexp)]
        path_colors <- node_colors[gen_nodes]
        # Add function colors
        toadd <- V(pg$graph)$name[!V(pg$graph)$name %in% rownames(difexp)]
        coltoadd <- rep("white", length(toadd))
        names(coltoadd) <- toadd
        path_colors <- c(path_colors, coltoadd)
        return(path_colors)
    })
    names(cols) <- names(metaginfo$pathigraphs)
    colors_de <- NULL
    colors_de$colors <- cols
    colors_de$group_by <- group_by
    return(colors_de)
}


#' @import grDevices
get_colors_from_pval <- function(updown, pvals, up_col = "#da1f1f",
                                 down_col = "#1f9cda", no_col = "white",
                                 both_col = "#959595",conf = 0.05){
    colors <- sapply(seq_along(updown), function(i){
        if(!is.na(pvals[i]) && pvals[i] <= conf){
            trans <- (1 - 18*pvals[i])
            if(is.na(updown[i])){
                return(no_col)
            }else if(updown[i] == "up"){
                cc <- grDevices::colorRamp(c(no_col, up_col))(trans)/255
                return(grDevices::rgb(cc[1], cc[2], cc[3]))
            }else if(updown[i] == "down" ){
                cc <- grDevices::colorRamp(c(no_col, down_col))(trans)/255
                return(grDevices::rgb(cc[1], cc[2], cc[3]))
            }else if(updown[i] == "both" ){
                cc <- grDevices::colorRamp(c(no_col, both_col))(trans)/255
                return(grDevices::rgb(cc[1], cc[2], cc[3]))
            }
        }else{
            return(no_col)
        }
    })
    return(colors)
}
