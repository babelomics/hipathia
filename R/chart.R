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
#' @param path.vals Matrix with the values of the subpathways.
#' Rows are subpathways and columns are samples.
#' @param sample.type Vector with the group to which each sample belongs.
#' The samples must be ordered as in \code{path.vals}. By default, all
#' samples will be assigned to the same class.
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
#' @param sample.clust Boolean, whether to cluster samples (columns).
#' By default TRUE.
#' @param variable.clust Boolean, whether to cluster variables (rows).
#' By default FALSE. If TRUE, rows with 0 variance are removed.
#' @param labRow,labCol Character vectors with row and column labels
#' to be used. By default rownames(path.vals) or colnames(path.vals)
#' are used, respectively.
#' @param sample.colors Named character vector of colors. The names of
#' the colors must be the classes in \code{sample.type}. Each sample will
#' be assigned the color corresponding to its class, taken from the
#' \code{sample.type} vector. By default a color will be assigned
#' automatically to each class.
#' @param scale Boolean, whether to scale each row to the interval [0,1].
#' Default is TRUE.
#' @param save.png Path to the file where the image as PNG will be saved.
#' By default, the image is not saved.
#' @param legend Boolean, whether to display a legend.
#' @param legend.xy Position for the legend, in case \code{legend} is TRUE.
#' @param pch Graphical parameter from \code{par()} function.
#' @param main Main title of the image
#'
#' @return Heatmap of the values of the subpathways
#'
#' @examples
#' data(brca_design)
#' data(path_vals)
#' sample.group <- brca_design[colnames(path_vals),"group"]
#' heatmap.plot(path_vals, sample.type = sample.group)
#' heatmap.plot(path_vals, sample.type = sample.group, colors = "hipathia",
#' variable.clust = TRUE)
#'
#' @export
#' @import grDevices graphics
#' @importFrom stats var
#' @importFrom stats heatmap
#'
heatmap.plot <- function(path.vals, sample.type = NULL, colors = "classic",
                         sample.clust = TRUE, variable.clust = FALSE,
                         labRow = NULL, labCol = NULL, sample.colors = NULL,
                         scale = TRUE, save.png = NULL, legend = TRUE,
                         legend.xy = "topright", pch = 15, main = NULL){
    if(length(colors) == 1){
        if(colors == "hipathia"){
            colors <- c("#007462", "white", "#e66430")
        }else if(colors == "classic"){
            colors <- c("blue","gray","red")
        }else if(colors == "redgreen"){
            colors <- c("green","black","red")}
    }
    if(is.null(sample.type))
        sample.type <- rep("A", ncol(path.vals))
    if(sample.clust==FALSE){
        colv <- NA
    } else {
        colv <- TRUE
    }
    if(variable.clust==FALSE){
        rowv <- NA
    } else {
        vars <- apply(path.vals, 1, stats::var)
        path.vals <- path.vals[!is.na(vars) & vars != 0,]
        rowv <- TRUE
    }
    if(is.null(labRow)){
        if(nrow(path.vals) < 50){
            labRow <- rownames(path.vals)
        }else{
            labRow <- FALSE
        }
    }
    if(is.null(labCol)){
        if(ncol(path.vals) < 50){
            labCol <- colnames(path.vals)
        }else{
            labCol <- FALSE
        }
    }
    if(is.null(sample.colors)){
        if(length(unique(sample.type)) <= 8){
            sample.colors <- c("#50b7ae", "#b6ebe7", "#e66430",
                               "#305f59", "#ffc868", "#152e2b",
                               "#a0170e", "#f9b493")[1:length(unique(
                                   sample.type))]
        }else{
            sample.colors <- c("#50b7ae", "#b6ebe7", "#e66430", "#305f59",
                               "#ffc868", "#152e2b", "#a0170e", "#f9b493",
                               grDevices::terrain.colors(length(unique(
                                   sample.type)) - 8))
        }
        names(sample.colors) <- unique(sample.type)
    }
    if(scale == TRUE){
        path.vals <- t(apply(path.vals, 1, function(x){
            (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) -
                                            min(x, na.rm = TRUE))
        }))
    }
    if(!is.null(save.png))
        grDevices::png(filename = save.png)
    if(!is.null(main))
        graphics::par(oma=c(1,0,3,1))
    stats::heatmap(path.vals,
                   margins = c(10,10),
                   labRow = labRow,
                   labCol = labCol,
                   scale = "none",
                   Rowv = rowv,
                   Colv = colv,
                   ColSideColors = sample.colors[sample.type],
                   col = grDevices::colorRampPalette(colors)(256))
    if(legend == TRUE)
        legend(legend.xy,
               legend = unique(sample.type),
               col = sample.colors[unique(sample.type)],
               pch = pch,
               xpd = TRUE,
               cex = 1,
               border = 0)
    if(!is.null(main)){
        title(main=main, outer = TRUE)
        graphics::par(oma = c(0,0,0,0))
    }
    if(!is.null(save.png))
        grDevices::dev.off()
}


#'
#' Plots two components of a PCA
#'
#' Plots two components of a PCA computed with \code{do.pca}
#'
#' @param fit princomp object as returned by \code{do.pca}
#' @param sample.type Vector with the group to which each sample belongs.
#' The samples must be ordered as in \code{path.vals}.
#' By default, all samples will be assigned to the same class.
#' @param sample.colors Named character vector of colors. The names of
#' the colors must be the classes in \code{sample.type}. Each sample will be
#' assigned the color corresponding to its class, taken from the
#' \code{sample.type} vector. By default a color will be assigned
#' automatically to each class.
#' @param cp1 Integer, number of the component in the X-axis.
#' Default is 1, the first component.
#' @param cp2 Integer, number of the component in the Y-axis.
#' Default is 2, the second component.
#' @param legend Boolean, whether to plot a legend in the plot.
#' Default is TRUE.
#' @param legend.xy Situation of the legend in the plot. Available
#' options are: "bottomright", "bottom", "bottomleft", "left",
#' "topleft", "top", "topright", "right" and "center".
#' @param cex Graphical parameter from \code{par()} function.
#' @param pch Graphical parameter from \code{par()} function.
#' @param mgp Graphical parameter from \code{par()} function.
#' @param main Title of the graphics
#' @param save.png Path to the file where the image as PNG will be saved.
#' By default, the image is not saved.
#'
#' @return Plots two components of a PCA
#'
#' @examples
#' data(path_vals)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' pca_model <- do.pca(path_vals[1:ncol(path_vals),])
#' pca.plot(pca_model, sample_group)
#'
#' @export
#' @import grDevices graphics
#'
pca.plot <- function(fit, sample.type = NULL, sample.colors = NULL, cp1 = 1,
                     cp2 = 2, legend = TRUE, legend.xy = "bottomleft", cex = 2,
                     pch = 20, mgp = c(3,1,0), main = "PCA plot",
                     save.png = NULL){
    if(is.null(sample.type)) sample.type <- rep("A", fit$n.obs)
    if(is.null(sample.colors)){
        if(length(unique(sample.type)) <= 8){
            sample.colors <- c("#50b7ae", "#b6ebe7", "#e66430",
                               "#305f59", "#ffc868", "#152e2b",
                               "#a0170e", "#f9b493")[1:length(unique(
                                   sample.type))]
        }else{
            sample.colors <- c("#50b7ae", "#b6ebe7", "#e66430", "#305f59",
                               "#ffc868", "#152e2b", "#a0170e", "#f9b493",
                               grDevices::topo.colors(length(unique(
                                   sample.type)) - 8))
        }
        names(sample.colors) <- unique(sample.type)
    }
    cpv1 <- fit$scores[,cp1]
    cpv2 <- fit$scores[,cp2]
    if(!is.null(save.png))
        grDevices::png(filename = save.png)
    graphics::plot(cpv1,
                   cpv2,
                   xlab = paste("PC", cp1),
                   ylab = paste("PC", cp2),
                   col = sample.colors[sample.type],
                   pch = pch,
                   cex = cex,
                   main = main,
                   mgp = mgp)
    if(legend == TRUE){
        legend(legend.xy,
               legend = unique(sample.type),
               col = sample.colors[unique(sample.type)],
               pch = pch,
               xpd = TRUE,
               cex = 1,
               border = 0)
    }
    if(!is.null(save.png))
        grDevices::dev.off()
}


#'
#' Plots multiple components of a PCA
#'
#' Plots multiple components of a PCA analysis computed with \code{do.pca}
#'
#' @param fit princomp object as returned by \code{do.pca}
#' @param sample.type Vector with the group to which each sample belongs.
#' The samples must be ordered as in \code{path.vals}.
#' By default, all samples will be assigned to the same class.
#' @param sample.colors Named character vector of colors. The names of the
#' colors must be the classes in \code{sample.type}. Each sample will be
#' assigned the color corresponding to its class, taken from the
#' \code{sample.type} vector. By default a color will be assigned
#' automatically to each class.
#' @param comps Vector with the components to be plot
#' @param plot.variance Logical, whether to plot the cumulative variance.
#' @param legend Boolean, whether to plot a legend in the plot.
#' Default is TRUE.
#' @param cex Graphical parameter from \code{par()} function.
#' @param pch Graphical parameter from \code{par()} function.
#' @param main Main title of the image
#' @param save.png Path to the file where the image as PNG will be saved.
#' By default, the image is not saved.
#'
#' @return Plots multiple components of a PCA
#'
#' @examples
#' data(path_vals)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' pca_model <- do.pca(path_vals[1:ncol(path_vals),])
#' multiple.pca.plot(pca_model, sample_group, cex = 3, plot.variance = TRUE)
#'
#' @export
#' @import graphics grDevices
#'
multiple.pca.plot <- function(fit, sample.type = NULL, sample.colors = NULL,
                              comps = 1:3, plot.variance = FALSE, legend = TRUE,
                              cex = 2, pch = 20, main = "Multiple PCA plot",
                              save.png = NULL){
    combs <- utils::combn(comps, 2)
    ncombs <- ncol(combs)
    nn <- ncombs
    if(!is.null(legend))
        nn <- nn + 1
    if(plot.variance==TRUE)
        nn <- nn + 1

    nr <- floor(sqrt(nn))
    nc <- ceiling((nn)/nr)
    oldmfrow <- par("mfrow")
    graphics::par(mfrow=c(nr, nc))
    oldmar <- par("mar")
    graphics::par(mar=c(4, 4, 1, 1))
    oldoma <- par("oma")
    graphics::par(oma=c(0, 0, 2, 0))
    if(!is.null(save.png))
        grDevices::png(filename = save.png)
    for(i in 1:(ncombs)){
        pca.plot(fit,
                 sample.type = sample.type,
                 sample.colors = sample.colors,
                 cp1 = combs[1,i],
                 cp2 = combs[2,i],
                 cex = cex,
                 pch = pch,
                 main = NULL,
                 legend = FALSE,
                 mgp = c(2.5,1,0))
    }
    if(legend == TRUE){
        if(is.null(sample.colors)){
            if(length(unique(sample.type)) <= 8){
                sample.colors <- c("#50b7ae", "#b6ebe7", "#e66430",
                                   "#305f59", "#ffc868", "#152e2b",
                                   "#a0170e", "#f9b493")[1:length(unique(
                                       sample.type))]
            }else{
                sample.colors <- c("#50b7ae", "#b6ebe7", "#e66430", "#305f59",
                                   "#ffc868", "#152e2b", "#a0170e", "#f9b493",
                                   grDevices::topo.colors(length(unique(
                                       sample.type)) - 8))
            }
            names(sample.colors) <- unique(sample.type)
        }
        graphics::plot(1, type="n", axes = FALSE, xlab = "", ylab = "")
        legend("center",
               legend = unique(sample.type),
               col = sample.colors[unique(sample.type)],
               pch = pch,
               lwd = 2,
               xpd = TRUE,
               cex = 1,
               border = NA,
               pt.cex = 1.2)
    }
    if(plot.variance == TRUE){
        plot.pca.variance(fit, acum = TRUE, thresh = 0.1)
    }
    title(main, outer = TRUE)
    graphics::par(mfrow = oldmfrow)
    graphics::par(oma = oldoma)
    graphics::par(mar = oldmar)
    if(!is.null(save.png))
        grDevices::dev.off()
}


plot.pca.variance <- function(fit, thresh = 0, acum = FALSE, minnum = 5){
    if(acum==FALSE){
        comptoplot <- fit$explain_var > thresh
        if(sum(comptoplot) < minnum)
            comptoplot <- 1:5
        graphics::barplot(fit$explain_var[ comptoplot ],
                          ylab = "explain variance",
                          xlab = "",
                          las = 2,
                          cex.names = 0.5,
                          ylim = c(0,1))
    } else {
        comptoplot <- fit$acum_explain_var < (1 - thresh)
        if(sum(comptoplot) < minnum) comptoplot <- 1:5
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
plot.pathigraph <- function(g, node.color = NULL, edge.lty = 1, main = "" ){
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
    if(is.null(node.color)){
        V(g)$color <- "white"}else{V(g)$color <- node.color}
    E(g)$lty <- edge.lty
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
#' colors for each node. Function \code{node.color.per.differential.expression}
#' assigns colors to the nodes depending on their differential expression.
#'
#' @param comp Comparison data frame as returned by the \code{do.wilcox}
#' function.
#' @param metaginfo Pathways object.
#' @param pathway Name of the pathway to be plotted.
#' @param conf Level of significance of the comparison for the adjusted
#' p-value. Default is 0.05.
#' @param node.colors List, named by the pathway name, including the
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
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' pathway.comparison.plot(comp, metaginfo = pathways, pathway = "hsa03320")
#'
#' data(results)
#' data(brca_design)
#' data(path_vals)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' colors_de <- node.color.per.differential.expression(results, pathways,
#' sample_group, "Tumor", "Normal")
#' pathway.comparison.plot(comp, metaginfo = pathways, pathway = "hsa04012",
#' node.colors = colors_de)
#'
#' @export
#'
pathway.comparison.plot <- function(comp, metaginfo, pathway, conf=0.05,
                                    node.colors = NULL, colors = "classic"){

    if(length(colors) == 1){
        if(colors == "hipathia"){ colors <- c("#50b7ae", "darkgrey", "#f16a34")
        }else if(colors == "classic"){ colors <- c("#0571b0", "darkgrey",
                                                   "#ca0020")}
    }
    down.col <- colors[1]
    no.col <- colors[2]
    up.col <- colors[3]

    pathigraph <- metaginfo$pathigraphs[[pathway]]
    if(all(grepl(" - ", rownames(comp)))){
        effector = FALSE
    }else{
        effector = TRUE}

    paths <- sapply(rownames(comp), function(x){
        unlist(strsplit(x, "-"))[2]
    })
    comp <- comp[paths == pathigraph$path.id, ]

    # Find edge colors
    g <- add.edge.colors(pathigraph,
                         comp,
                         effector,
                         up.col = up.col,
                         down.col = down.col,
                         no.col = no.col)

    edge.lty <- (E(g)$relation*-1 +1)/2 + 1
    id <- metaginfo$pathigraphs[[pathway]]$path.id
    name <- metaginfo$pathigraphs[[pathway]]$path.name
    title <- paste(id, "-", name)
    plot.pathigraph(g,
                    node.color = node.colors[[pathway]],
                    edge.lty = edge.lty,
                    main = title )
}



add.edge.colors <- function(pathigraph, pcomp, effector, up.col = "#ca0020",
                            down.col = "#0571b0", no.col = "darkgrey",
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
    for( path.name in rownames(pcomp)){
        if(pcomp[path.name,"FDRp.value"] <= conf){
            up.down <- pcomp[path.name,"UP/DOWN"]
        }else{
            up.down <- "N"
        }
        subgraph <- subgraphs[[path.name]]
        els <- apply(get.edgelist(subgraph), 1, paste, collapse = "_")
        states[path.name, els] <- up.down
    }

    E(g)$color <- rep(no.col, length(E(g)))

    colors <- c(up.col, no.col, down.col)
    names(colors) <- c("UP", "N", "DOWN")

    for(i in which(!grepl("_func", elg))){
        edge.states <- names(table(states[,i]))
        if(is.null(edge.states)){
            warning(paste("Edge", elg[i], "is not present in any subpath."))
        }else{
            E(g)$color[i] <- colors[edge.states[1]]
            if(length(edge.states) > 1)
                for(j in 2:length(edge.states)){
                    g <- g + edge(unlist(strsplit(elg[i], split = "\\_")))
                    E(g)$color[length(E(g))] <- colors[edge.states[j]]
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
#' of the nodes depending on it. Significant up- and down-regulated nodes
#' are depicted with the selected color, with a gradient towards the
#' non-significant color depending on the value of the p-value.
#' Smaller p-values give rise to purer colors than higher p-values.
#'
#' @param results Object of results as provided by the \code{hipathia}
#' function.
#' @param metaginfo Object of pathways.
#' @param groups Vector with the class to which each sample belongs.
#' Samples must be ordered as in \code{results}
#' @param group1.label String, label of the first group to be compared
#' @param group2.label String, label of the second group to be compared
#' @param group.by How to group the subpathways to be visualized. By default
#' they are grouped by the pathway to which they belong. Available groupings
#' include "uniprot", to group subpathways by their annotated Uniprot functions,
#' "GO", to group subpathways by their annotated GO terms, and "genes", to group
#' subpathways by the genes they include. Default is set to "pathway".
#' @param colors Either a character vector with 3 colors (indicating,
#' in this order, down-regulation, non-significance and up-regulation colors)
#'  or a key name indicating the color scheme to be used. Options are:
#' @slot classic ColorBrewer blue, white and colorBrewer red.
#' @slot hipathia Hipathia predefined color scheme: Green, white and orange.
#' By default \code{classic} color scheme is applied.
#' @param conf Level of significance of the comparison for the adjusted p-value
#'
#' @return List of color vectors, named by the pathways to which they belong.
#' The color vectors represent the differential expression
#' of the nodes in each pathway.
#'
#' @examples
#' data(results)
#' data(brca_design)
#' data(path_vals)
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' colors_de <- node.color.per.differential.expression(results, pathways,
#' sample_group, "Tumor", "Normal")
#'
#' @export
#'
node.color.per.differential.expression <- function(results, metaginfo, groups,
                                                   group1.label, group2.label,
                                                   group.by = "pathway",
                                                   colors = "classic",
                                                   conf = 0.05){

    if(length(colors) == 1){
        if(colors == "hipathia"){
            colors <- c("#50b7ae", "white", "#f16a34")
        }else if(colors == "classic"){
            colors <- c("#1f9cda","white","#da1f1f")
        }
    }
    down.col <- colors[1]
    no.col <- colors[2]
    up.col <- colors[3]

    if(group.by != "pathway")
        metaginfo <- get.pseudo.metaginfo(metaginfo, group.by = group.by)

    difexp <- compute.difexp(results$all$nodes.vals, group1.label,
                             group2.label, groups)
    updown <- sapply(difexp$statistic, function(x){
        if(x < 0){
            "down"
        }else if(x > 0){
            "up"
        }else{
            "both"
        }
    })
    node.colors <- get.colors.from.pval(updown,
                                        difexp$p.value,
                                        up.col = up.col,
                                        down.col = down.col,
                                        no.col = no.col,
                                        conf = conf)
    names(node.colors) <- rownames(results$all$nodes.vals)
    cols <- lapply(metaginfo$pathigraphs, function(pg){
        gen.nodes <- V(pg$graph)$name[V(pg$graph)$name %in% rownames(difexp)]
        path.colors <- node.colors[gen.nodes]
        # Add function colors
        toadd <- V(pg$graph)$name[!V(pg$graph)$name %in% rownames(difexp)]
        coltoadd <- rep("white", length(toadd))
        names(coltoadd) <- toadd
        path.colors <- c(path.colors, coltoadd)
        return(path.colors)
    })
    names(cols) <- names(metaginfo$pathigraphs)
    colors.de <- NULL
    colors.de$colors <- cols
    colors.de$group.by <- group.by
    return(colors.de)
}


#' @import grDevices
get.colors.from.pval <- function(updown, pvals, up.col = "#da1f1f",
                                 down.col = "#1f9cda", no.col = "white",
                                 both.col = "#959595",conf = 0.05){
    colors <- sapply(1:length(updown), function(i){
        if(!is.na(pvals[i]) && pvals[i] <= conf){
            trans <- (1 - 18*pvals[i])
            if(is.na(updown[i])){
                return(no.col)
            }else if(updown[i] == "up"){
                cc <- grDevices::colorRamp(c(no.col, up.col))(trans)/255
                return(grDevices::rgb(cc[1], cc[2], cc[3]))
            }else if(updown[i] == "down" ){
                cc <- grDevices::colorRamp(c(no.col, down.col))(trans)/255
                return(grDevices::rgb(cc[1], cc[2], cc[3]))
            }else if(updown[i] == "both" ){
                cc <- grDevices::colorRamp(c(no.col, both.col))(trans)/255
                return(grDevices::rgb(cc[1], cc[2], cc[3]))
            }
        }else{
            return(no.col)
        }
    })
    return(colors)
}
