
#' Normalize expression data to be used in \code{hipathia}
#'
#' Transforms the rank of the matrix of gene expression to [0,1] in order
#' to be processed by \code{hipathia}. The transformation may be performed
#' in two different ways. If \code{percentil = FALSE}, the transformation
#' is a re-scaling of the rank of the matrix. If \code{percentil = TRUE},
#' the transformation is performed assigning to each cell its percentil in
#' the corresponding distribution. This option is recommended for
#' distributions with very long tails.
#'
#' This transformation may be applied either to the whole matrix
#' (by setting \code{by.gene = FALSE}), which we strongly recommend, or to
#' each of the rows (by setting \code{by.gene = TRUE}), allowing each gene
#' to have its own scale.
#'
#' A previous quantiles normalization may be applied by setting
#' \code{by.quantiles = TRUE}. This is recommended for noisy data.
#'
#' For distributions with extreme outlayer values, a percentil \code{p}
#' may be given to the parameter \code{truncation.percentil}. When provided,
#' values beyond percentil p are truncated to the value of percentil p, and
#' values beyond 1-p are truncated to percentil 1-p. This step is performed
#'  before any other tranformation. By default no truncation is performed.
#'
#' @param exp.data Matrix of gene expression.
#' @param by.quantiles Boolean, whether to normalize the data by quantiles.
#' Default is FALSE.
#' @param by.gene Boolean, whether to transform the rank of each row of the
#' matrix to [0,1]. Default is FALSE.
#' @param percentil Boolean, whether to take as value the percentil of each
#' sample in the corresponding distribution.
#' @param truncation.percentil Real number p in [0,1]. When provided, values
#' beyond percentil p are truncated to the value of percentil p, and values
#' beyond 1-p are truncated to percentil 1-p. By default no truncation
#' is performed.
#'
#' @return Matrix of gene expression whose values are in [0,1].
#'
#' @examples data("brca_data")
#' trans.data <- translate.matrix(brca_data, "hsa")
#' exp.data <- normalize.data(trans.data)
#' exp.data <- normalize.data(trans.data, by.quantiles = TRUE,
#' truncation.percentil=0.95)
#'
#' @export
#' @import preprocessCore
#' @importFrom stats quantile
#' @importFrom stats ecdf
#'
normalize.data <- function(exp.data, by.quantiles = FALSE, by.gene = FALSE,
                           percentil = FALSE, truncation.percentil = NULL){

    # Normalize data matrix
    norm.data <- data.matrix(exp.data)
    if(!is.null(truncation.percentil)){
        if( truncation.percentil >= 0 & truncation.percentil <= 1){
            # Guarantees that truncation.percentil is in [0.5,1]
            if(truncation.percentil < (1-truncation.percentil))
                truncation.percentil <- 1-truncation.percentil
            # Truncates by the percentil
            norm.data <- t(apply(norm.data, 1, function(x){
                quan.inf <- stats::quantile(x, 1 - truncation.percentil,
                                            na.rm = TRUE)
                x[which(x < quan.inf)] <- quan.inf
                quan.sup <- stats::quantile(x, truncation.percentil,
                                            na.rm = TRUE)
                x[which(x > quan.sup)] <- quan.sup
                x
            }))
        }
        else{
            stop("Parameter truncation.percentil must be in [0,1]")
        }
    }
    if(by.quantiles == TRUE){
        norm.data <- preprocessCore::normalize.quantiles(norm.data)
    }
    if(by.gene == TRUE){
        if(percentil == TRUE){
            norm.data <- t(apply(norm.data,1, function(x){stats::ecdf(x)(x)}))
        }else{
            norm.data <- t(apply(norm.data, 1, function(x){
                (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) -
                                                min(x, na.rm = TRUE))
            }))
        }
    } else {
        if(percentil == TRUE){
            emp <- stats::ecdf(norm.data)
            norm.data <- t(apply(norm.data,1,emp))
        }else{
            norm.data <- (norm.data - min(norm.data, na.rm = TRUE))/
                (max(norm.data, na.rm = TRUE) - min(norm.data, na.rm = TRUE))
        }
    }

    colnames(norm.data) <- colnames(exp.data)
    rownames(norm.data) <- rownames(exp.data)

    return(norm.data)
}


#'@importFrom stats p.adjust
do.anova.test <- function(data, group, adjust = TRUE){
    tests <- apply(data, 1, anova.test.fun, group = group)
    out <- do.call("rbind", lapply(tests, function(x){
        c(summary(x$fit)[[1]][[1, "Pr(>F)"]])
    }))
    out[is.na(out)] <- 1
    if(adjust == TRUE){
        out <- cbind(out, stats::p.adjust(out[,1], method = "fdr"))
    } else {
        out <- cbind(out, out[,1])
    }
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    colnames(out) <- c("p.value", "adj.p.value")
    return(out)
}

#' @importFrom stats aov
#' @importFrom stats TukeyHSD
anova.test.fun <- function(values, group, verbose = FALSE){
    group <- factor(group)
    fit <- stats::aov(values ~ group)
    tuckey <- stats::TukeyHSD(fit)
    if(verbose==TRUE){
        summary(fit)
    }
    return(list(fit = fit,tuckey = tuckey))
}


#'@importFrom stats p.adjust
do.cor <- function(sel.vals, design, adjust = TRUE){

    data <- sel.vals[,design$sample]

    testData <- do.call("rbind",apply(data, 1, cor.test.fun, design$value))
    if(adjust==TRUE){
        fdrData <- stats::p.adjust(testData[,1], method = "fdr")
    } else {
        fdrData <- testData[,1]
    }
    data2 <- data.frame(testData[,1:3], fdrData, stringsAsFactors = FALSE)

    colnames(data2) <- c("p.value", "UP/DOWN", "correlation", "FDRp.value")
    data2[data2$statistic>0,"UP/DOWN"] <- "UP"
    data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
    return(data2)
}


#'@importFrom stats cor.test
cor.test.fun <- function(x, values){
    r <- try(stats::cor.test(x = as.numeric(x), y = as.numeric(values)))
    result <- cor.data.frame(r)
    return(result)
}


cor.data.frame <- function(wilcox){
    if (class(wilcox) == "try-error" | is.na( wilcox$p.value )){
        pvalue <- 1
        class <- "0"
        esti <- 0
    } else{
        pvalue <- wilcox$p.value
        esti <- wilcox$estimate[1]
        if (esti < 0){
            class <- "DOWN" ## regarding DISEASE
        } else if (esti > 0){
            class <- "UP" ## regarding DISEASE
        } else if (esti == 0){
            if (wilcox$conf.int[1] == 0){
                class <- "UP"
            } else if (wilcox$conf.int[2] == 0){
                class <- "DOWN"
            } else{
                class <- 0
            }
        }
    }
    result <- data.frame(pvalue, class, esti, stringsAsFactors = FALSE)
    return(result)
}


#' Apply Wilcoxon test
#'
#' Performs a Wilcoxon test for the values in \code{sel.vals} comparing
#' conditions \code{g1} and \code{g2}
#'
#' @param sel.vals Matrix of values. Columns represent samples.
#' @param group.value Vector with the class to which each sample belongs.
#' Samples must be ordered as in \code{sel.vals}
#' @param g1 String, label of the first group to be compared
#' @param g2 String, label of the second group to be compared
#' @param paired Boolean, whether the samples to be compared are paired.
#' If TRUE, function \code{wilcoxsign_test} from package \code{coin} is
#' used. If FALSE, function \code{wilcox.test} from package \code{stats}
#' is used.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method
#'
#' @return Dataframe with the result of the comparison
#'
#' @examples
#' data(path_vals)
#' data(brca_design)
#' sample.group <- brca_design[colnames(path_vals),"group"]
#' comp <- do.wilcoxon(path_vals, sample.group, g1 = "Tumor", g2 = "Normal")
#'
#' @export
#'
do.wilcoxon <- function(sel.vals, group.value, g1, g2, paired=FALSE,
                        adjust=TRUE){

    g1_indexes <- which(group.value == g1)
    g2_indexes <- which(group.value == g2)

    stat.vals <- suppressWarnings(calculate.wilcox.test(sel.vals,
                                                        g2_indexes,
                                                        g1_indexes,
                                                        paired = paired,
                                                        adjust = adjust))
    return(stat.vals)
}


#'@importFrom stats p.adjust
calculate.wilcox.test <- function( data, control, disease, paired, adjust=TRUE){
    if(paired == TRUE){
        dat <- apply(data, 1, wilcoxsign.test.fun, control, disease)
        testData <- do.call("rbind",dat)
        if(adjust==TRUE){
            fdrData <- stats::p.adjust(testData[,1], method = "fdr")
        } else {
            fdrData <- testData[,1]
        }
        data2 <- data.frame(testData, fdrData, stringsAsFactors=FALSE)
    }else{
        testData <- do.call("rbind",
                            apply(data, 1, wilcox.test.fun,
                                  control, disease, paired))
        if(adjust==TRUE){
            fdrData <- stats::p.adjust(testData[,1], method = "fdr")
        } else {
            fdrData <- testData[,1]
        }
        # Standardize statistic
        lc <- length(control)
        ld <- length(disease)
        m <- lc*ld/2
        sigma <- sqrt(lc*ld*(lc+ld+1)/12)
        z <- (testData[,3]-m)/sigma
        data2 <- data.frame(testData[,1:2], z, fdrData, stringsAsFactors=FALSE)
        need0 <- which(data2$pvalue==1 &
                           data2$class == "0" &
                           data2$fdrData == 1)
        data2[need0, "z"] <- 0
    }
    colnames(data2) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
    data2[data2$statistic>0,"UP/DOWN"] <- "UP"
    data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
    return(data2)
}


wilcoxsign.test.fun <- function(x, control, disease){
    r <- try(coin::wilcoxsign_test(as.numeric(x[disease])~
                                       as.numeric(x[control]),
                                   showWarnings = FALSE))
    if (class(r) == "try-error"){
        pvalue <- 1
        class <- "0"
        stat <- 0
    } else{
        pvalue <- coin::pvalue(r)
        stat <- coin::statistic(r, "standardized")[1]
        if(stat > 0){
            class <- "UP"
        }else if (stat < 0){
            class <- "DOWN"
        }else{
            class <- "0"
        }
    }
    result <- data.frame(pvalue = pvalue,
                         class = class,
                         stat = stat,
                         stringsAsFactors = FALSE)
    return(result)
}


#'@importFrom stats wilcox.test
wilcox.test.fun <- function(x, control, disease, paired){
    r <- try(stats::wilcox.test(x = as.numeric(x[disease]),
                                y = as.numeric(x[control]),
                                conf.int = TRUE,
                                alternative = "two.sided",
                                paired = paired),
             silent = TRUE)
    result <- wilcox.data.frame(r)
    return(result)
}


wilcox.data.frame <- function(wilcox){
    if (class(wilcox) == "try-error"){
        pvalue <- 1
        class <- "0"
        stat <- 0
    } else{
        pvalue <- wilcox$p.value
        esti <- wilcox$estimate
        stat <- wilcox$statistic[[1]]
        if (esti < 0){
            class <- "DOWN" ## regarding DISEASE
        } else if (esti > 0){
            class <- "UP" ## regarding DISEASE
        } else if (esti == 0){
            if (wilcox$conf.int[1] == 0){
                class <- "UP"
            } else if (wilcox$conf.int[2] == 0){
                class <- "DOWN"
            } else{
                class <- 0
            }
        }
    }
    result <- data.frame(pvalue, class, stat,stringsAsFactors=FALSE)
    return(result)
}


#'
#' Performs a Principal Components Analysis
#'
#' @param data Matrix of values to be analyzed. Samples must be represented
#' in the columns.
#' @param cor A logical value indicating whether the calculation should use
#' the correlation matrix or the covariance matrix. (The correlation matrix
#' can only be used if there are no constant variables.)
#'
#' @return \code{do.pca} returns a list with class \code{princomp}.
#'
#' @examples
#' data(path_vals)
#' pca.model <- do.pca(path_vals[1:ncol(path_vals),])
#'
#' @export
#' @importFrom stats princomp
#'
do.pca <- function(data, cor = FALSE){
    fit <- stats::princomp(t(data), cor = cor)
    fit$var <- fit$sdev^2
    fit$explain_var <- fit$var/sum(fit$var)
    fit$acum_explain_var <- cumsum(fit$explain_var)
    return(fit)
}


compute.difexp <- function(vals, group1.label, group2.label, groups){

    g1_indexes <- which(groups == group1.label)
    g2_indexes <- which(groups == group2.label)
    grupo1 <- (groups[c(g1_indexes, g2_indexes)] == group1.label) + 0
    grupo2 <- (groups[c(g1_indexes, g2_indexes)] == group2.label) + 0
    design <- cbind(grupo1, grupo2)

    fit <- limma::lmFit(vals[,c(g1_indexes, g2_indexes)], design)
    cont.matrix <- limma::makeContrasts(grupo1 - grupo2, levels = design)
    fit2 <- limma::contrasts.fit(fit, cont.matrix)
    fit2 <- limma::eBayes(fit2)
    result <- data.frame(statistic = as.numeric(fit2$t),
                         p.value = as.numeric(fit2$p.value),
                         laterality = as.factor(fit2$t>0))
    rownames(result) <- rownames(fit2)

    return(result)
}


compute.node.difexp <- function(results, groups, group1.label, group2.label,
                                verbose = FALSE){

    difexp <- list()
    for(pathway in names(results$by.path)){
        if(verbose == TRUE)
            print(pathway)
        difexp[[pathway]]<-compute.difexp(results$by.path[[pathway]]$nodes.vals,
                                          group1.label,
                                          group2.label,
                                          groups)
    }
    return(difexp)
}


#' Compute pathway summary
#'
#' Computes a summary of the results, summarizing the number and proportion
#' of up- and down-regulated subpathways in each pathway.
#'
#' @param comp Comparison data frame as returned by the \code{do.wilcoxon}
#' function.
#' @param metaginfo Pathways object
#' @param conf Level of significance of the comparison for the adjusted
#' p-value. Default is 0.05.
#'
#' @return Table with the summarized information for each of the pathways.
#' Rows are the analized pathways. Columns are:
#' * \code{num.total.paths} Number of total subpathways in which each pathway
#' is decomposed.
#' * \code{num.significant.paths} Number of significant subpathways in the
#' provided comparison.
#' * \code{percent.significant.paths} Percentage of significant subpathways
#' from the total number of subpathways in a pathway.
#' * \code{num.up.paths} Number of significant up-regulated subpathways in the
#' provided comparison.
#' * \code{percent.up.paths} Percentage of significant up-regulated subpathways
#' from the total number of subpathways in a pathway.
#' * \code{num.down.paths} Number of significant down-regulated subpathways in
#' the provided comparison.
#' * \code{percent.down.paths} Percentage of significant down-regulated
#' subpathways from the total number of subpathways in a pathway.
#'
#' @examples
#' data(comp)
#' pathways <- load.pathways(species = "hsa", pathways.list = c("hsa03320",
#' "hsa04012"))
#' get.pathways.summary(comp, pathways)
#'
#' @export
#'
get.pathways.summary <- function(comp, metaginfo, conf = 0.05){
    comp$pathways <- sapply(rownames(comp), function(n){
        unlist(strsplit(n, split = "-"))[2]
    })
    id.pathways <- unique(comp$pathways)
    name.pathways <- sapply(id.pathways, function(id){
        metaginfo$pathigraphs[[id]]$path.name
    })
    summ <- lapply(id.pathways, function(pathway){
        minicomp <- comp[comp$pathways == pathway,]
        num.total.paths <- nrow(minicomp)
        num.sig.paths <- sum(minicomp$FDRp.value < conf)
        percent.sig.paths <- round(num.sig.paths/nrow(minicomp)*100, digits = 2)
        is.up.path <- minicomp$FDRp.value < conf & minicomp$`UP/DOWN` == "UP"
        num.up.paths <- sum(is.up.path)
        percent.up.paths <- round(num.up.paths/nrow(minicomp) * 100, digits = 2)
        is.down.path <- minicomp$FDRp.value < conf &
            minicomp$`UP/DOWN` == "DOWN"
        num.down.paths <- sum(is.down.path)
        percent.down.paths <- round(num.down.paths/nrow(minicomp) * 100,
                                    digits = 2)
        data.frame(num.total.paths = num.total.paths,
                   num.significant.paths = num.sig.paths,
                   percent.significant.paths = percent.sig.paths,
                   num.up.paths = num.up.paths,
                   percent.up.paths = percent.up.paths,
                   num.down.paths = num.down.paths,
                   percent.down.paths = percent.down.paths)
    })
    summ <- do.call("rbind", summ)
    rownames(summ) <- name.pathways
    summ <- cbind(id.pathways, summ)
    summ <- summ[order(summ$percent.significant.paths, decreasing = TRUE),]
    return(summ)
}


get.pathways.pvalues <- function(comp, conf = 0.05){
    comp$pathways <- sapply(rownames(comp), function(n){
        unlist(strsplit(n, split="\\_"))[2]
    })
    name.pathways <- unique(comp$pathways)
    pvals <- lapply(name.pathways, function(pathway){
        comp[comp$pathways == pathway, "FDRp.value"]
    })
    names(pvals) <- name.pathways
    return(pvals)
}

