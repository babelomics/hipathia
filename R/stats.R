##
## stats.R
## Statistics functions
##
## Written by Marta R. Hidalgo, marta.hidalgo@outlook.es
##
## Code style by Hadley Wickham (http://r-pkgs.had.co.nz/style.html)
## https://www.bioconductor.org/developers/how-to/coding-style/
##

#' Normalize expression data from a SummarizedExperiment or matrix to be used 
#' in 
#' \code{hipathia}
#'
#' Transforms the rank of the SummarizedExperiment or matrix of gene expression 
#' to [0,1] in order
#' to be processed by \code{hipathia}. The transformation may be performed
#' in two different ways. If \code{percentil = FALSE}, the transformation
#' is a re-scaling of the rank of the matrix. If \code{percentil = TRUE},
#' the transformation is performed assigning to each cell its percentil in
#' the corresponding distribution. This option is recommended for
#' distributions with very long tails.
#'
#' This transformation may be applied either to the whole matrix
#' (by setting \code{by_gene = FALSE}), which we strongly recommend, or to
#' each of the rows (by setting \code{by_gene = TRUE}), allowing each gene
#' to have its own scale.
#'
#' A previous quantiles normalization may be applied by setting
#' \code{by_quantiles = TRUE}. This is recommended for noisy data.
#'
#' For distributions with extreme outlayer values, a percentil \code{p}
#' may be given to the parameter \code{truncation_percentil}. When provided,
#' values beyond percentil p are truncated to the value of percentil p, and
#' values beyond 1-p are truncated to percentil 1-p. This step is performed
#' before any other tranformation. By default no truncation is performed.
#'
#' @param data Either a SummarizedExperiment or a matrix of gene expression.
#' @param sel_assay Character or integer, indicating the assay to be normalized 
#' in the SummarizedExperiment. Default is 1.
#' @param by_quantiles Boolean, whether to normalize the data by quantiles.
#' Default is FALSE.
#' @param by_gene Boolean, whether to transform the rank of each row of the
#' matrix to [0,1]. Default is FALSE.
#' @param percentil Boolean, whether to take as value the percentil of each
#' sample in the corresponding distribution.
#' @param truncation_percentil Real number p in [0,1]. When provided, values
#' beyond percentil p are truncated to the value of percentil p, and values
#' beyond 1-p are truncated to percentil 1-p. By default no truncation
#' is performed.
#'
#' @return Matrix of gene expression whose values are in [0,1].
#'
#' @examples data("brca_data")
#' trans_data <- translate_data(brca_data, "hsa")
#' exp_data <- normalize_data(trans_data)
#' exp_data <- normalize_data(trans_data, by_quantiles = TRUE,
#' truncation_percentil=0.95)
#'
#' @export
#' @import preprocessCore
#' @import SummarizedExperiment
#' @importFrom stats quantile
#' @importFrom stats ecdf
#' @importFrom methods is
#'
normalize_data <- function(data, sel_assay = 1, by_quantiles = FALSE, 
                         by_gene = FALSE, percentil = FALSE, 
                         truncation_percentil = NULL){
    
    if(is(data, "SummarizedExperiment")){
        se_flag <- TRUE
        mat <- assay(data, sel_assay)
    }else if(is(data, "matrix")){
        se_flag <- FALSE
        mat <- data
    }else{
        stop("Only SummarizedExperiment or matrix classes accepted as data")
    }
    norm_mat <- normalize_matrix(mat, by_quantiles = by_quantiles, 
                                 by_gene = by_gene, percentil = percentil, 
                                 truncation_percentil=truncation_percentil)
    if(se_flag == TRUE){
        norm_data <- SummarizedExperiment(list(norm = norm_mat), 
                                          colData = colData(data))
    }else{
        norm_data <- norm_mat
    }
    return(norm_data)
}


#' @import preprocessCore
#' @importFrom stats quantile
#' @importFrom stats ecdf
#' @importFrom DelayedArray rowMins
#' @importFrom DelayedArray rowMaxs
normalize_matrix <- function(mat, sel_assay = 1, by_quantiles = FALSE, 
                           by_gene = FALSE,
                           percentil = FALSE, truncation_percentil = NULL){
    
    # Normalize data matrix
    norm_data <- data.matrix(mat)
    if(!is.null(truncation_percentil)){
        if( truncation_percentil >= 0 & truncation_percentil <= 1){
            # Guarantees that truncation_percentil is in [0.5,1]
            if(truncation_percentil < (1-truncation_percentil))
                truncation_percentil <- 1-truncation_percentil
            # Truncates by the percentil
            norm_data <- t(apply(norm_data, 1, function(x){
                quan_inf <- stats::quantile(x, 1 - truncation_percentil,
                                            na.rm = TRUE)
                x[x < quan_inf] <- quan_inf
                quan_sup <- stats::quantile(x, truncation_percentil,
                                            na.rm = TRUE)
                x[x > quan_sup] <- quan_sup
                x
            }))
        }
        else{
            stop("Parameter truncation_percentil must be in [0,1]")
        }
    }
    if(by_quantiles == TRUE){
        norm_data <- preprocessCore::normalize.quantiles(norm_data)
    }
    if(by_gene == TRUE){
        if(percentil == TRUE){
            norm_data <- t(apply(norm_data, 1, function(x){stats::ecdf(x)(x)}))
        }else{
            min <- rowMins(norm_data, na.rm = TRUE)
            max <- rowMaxs(norm_data, na.rm = TRUE)
            norm_data <- (norm_data - min)/(max - min)
        }
    } else {
        if(percentil == TRUE){
            emp <- stats::ecdf(norm_data)
            norm_data <- t(apply(norm_data, 1, emp))
        }else{
            norm_data <- (norm_data - min(norm_data, na.rm = TRUE))/
                (max(norm_data, na.rm = TRUE) - min(norm_data, na.rm = TRUE))
        }
    }
    colnames(norm_data) <- colnames(mat)
    rownames(norm_data) <- rownames(mat)

    return(norm_data)
}


#'@importFrom stats p.adjust
do_anova_test <- function(data, group, adjust = TRUE){
    tests <- apply(data, 1, anova_test_fun, group = group)
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
anova_test_fun <- function(values, group, verbose = FALSE){
    group <- factor(group)
    fit <- stats::aov(values ~ group)
    tuckey <- stats::TukeyHSD(fit)
    if(verbose==TRUE){
        summary(fit)
    }
    return(list(fit = fit,tuckey = tuckey))
}


#'@importFrom stats p.adjust
do_cor <- function(sel_vals, design, adjust = TRUE){

    data <- sel_vals[,design$sample]

    testData <- do.call("rbind",apply(data, 1, cor_test_fun, design$value))
    if(adjust==TRUE){
        fdrData <- stats::p.adjust(testData[,1], method = "fdr")
    } else {
        fdrData <- testData[,1]
    }
    data2 <- data.frame(testData[,seq_len(3)], fdrData, 
                        stringsAsFactors = FALSE)

    colnames(data2) <- c("p.value", "UP/DOWN", "correlation", "FDRp.value")
    data2[data2$statistic>0,"UP/DOWN"] <- "UP"
    data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
    return(data2)
}


#'@importFrom stats cor.test
cor_test_fun <- function(x, values){
    r <- try(stats::cor.test(x = as.numeric(x), y = as.numeric(values)))
    result <- cor_data_frame(r)
    return(result)
}


cor_data_frame <- function(wilcox){
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
#' Performs a Wilcoxon test for the values in \code{sel_vals} comparing
#' conditions \code{g1} and \code{g2}
#'
#' @param data Either a SummarizedExperiment object or a matrix, containing the
#' values. Columns represent samples.
#' @param group Either a character indicating the name of the column in colData 
#' including the classes to compare, or a character vector with the class to 
#' which each sample belongs. 
#' Samples must be ordered as in \code{data}
#' @param g1 String, label of the first group to be compared
#' @param g2 String, label of the second group to be compared
#' @param paired Boolean, whether the samples to be compared are paired.
#' If TRUE, function \code{wilcoxsign_test} from package \code{coin} is
#' used. If FALSE, function \code{wilcox.test} from package \code{stats}
#' is used.
#' @param adjust Boolean, whether to adjust the p.value with
#' Benjamini-Hochberg FDR method
#' @param sel_assay Character or integer, indicating the assay to be normalized 
#' in the SummarizedExperiment. Default is 1.
#' @param order Boolean, whether to order the results table by the 
#' \code{FDRp.value} column. Default is FALSE. 
#'
#' @return Dataframe with the result of the comparison
#'
#' @examples
#' data(path_vals)
#' data(brca_design)
#' sample_group <- brca_design[colnames(path_vals),"group"]
#' comp <- do_wilcoxon(path_vals, sample_group, g1 = "Tumor", g2 = "Normal")
#'
#' @export
#' @import SummarizedExperiment
#' @importFrom methods is
#'
do_wilcoxon <- function(data, group, g1, g2, paired = FALSE,
                        adjust = TRUE, sel_assay = 1, order = FALSE){

    if(is(data, "SummarizedExperiment")){
        se_flag <- TRUE
        if(is(group, "character") & length(group) == 1)
            if(group %in% colnames(colData(data))){
                group <- colData(data)[[group]]
            }else{
                stop("Group variable must be a column in colData(data)")
            }
        vals <- assay(data, sel_assay)
    }else if(is(data, "matrix")){
        se_flag <- FALSE
        vals <- data
    }else{
        stop("Only SummarizedExperiment or matrix classes accepted as data")
    }

    g1_indexes <- which(group == g1)
    g2_indexes <- which(group == g2)

    stat_vals <- suppressWarnings(
        calculate_wilcox_test(vals, g2_indexes, g1_indexes, paired = paired,
                              adjust = adjust))
    if(se_flag == TRUE && "subpath.name" %in% colnames(rowData(data)))
        stat_vals <- cbind(name = rowData(data)[["subpath.name"]], stat_vals)
    if(order == TRUE)
        stat_vals <- stat_vals[order(stat_vals$FDRp.value, decreasing = FALSE),]
    return(stat_vals)
}




#'@importFrom stats p.adjust
calculate_wilcox_test <- function(data, control, disease, paired, adjust=TRUE){
    if(paired == TRUE){
        dat <- apply(data, 1, wilcoxsign_test_fun, control, disease)
        testData <- do.call("rbind", dat)
        if(adjust == TRUE){
            fdrData <- stats::p.adjust(testData[,1], method = "fdr")
        } else {
            fdrData <- testData[,1]
        }
        data2 <- data.frame(testData[,c(2,3,1)], fdrData, 
                            stringsAsFactors = FALSE)
    }else{
        testData <- do.call("rbind",
                            apply(data, 1, wilcox_test_fun,
                                  control, disease, paired))
        if(adjust == TRUE){
            fdrData <- stats::p.adjust(testData[,1], method = "fdr")
        } else {
            fdrData <- testData[,1]
        }
        # Standardize statistic
        lc <- length(control)
        ld <- length(disease)
        m <- lc*ld/2
        sigma <- sqrt(lc * ld * (lc + ld + 1)/12)
        z <- (testData[,3] - m)/sigma
        data2 <- data.frame(testData[,2], z, fdrData, testData[,1],
                            stringsAsFactors = FALSE)
        rownames(data2) <- rownames(testData)
        need0 <- which(data2$pvalue == 1 &
                           data2$class == "0" &
                           data2$fdrData == 1)
        data2[need0, "z"] <- 0
    }
    colnames(data2) <- c("UP/DOWN", "statistic", "p.value", "FDRp.value")
    data2[data2$statistic>0,"UP/DOWN"] <- "UP"
    data2[data2$statistic<0,"UP/DOWN"] <- "DOWN"
    return(data2)
}


wilcoxsign_test_fun <- function(x, control, disease){
    r <- try(coin::wilcoxsign_test(
        as.numeric(x[disease]) ~ as.numeric(x[control]), showWarnings = FALSE))
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
wilcox_test_fun <- function(x, control, disease, paired){
    r <- try(stats::wilcox.test(x = as.numeric(x[disease]),
                                y = as.numeric(x[control]),
                                conf.int = TRUE,
                                alternative = "two.sided",
                                paired = paired),
             silent = TRUE)
    result <- wilcox_data_frame(r)
    return(result)
}


wilcox_data_frame <- function(wilcox){
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
#' @param data SummarizedExperiment or matrix of values to be analyzed. Samples 
#' must be represented in the columns.
#' @param sel_assay Character or integer, indicating the assay to be normalized 
#' in the SummarizedExperiment. Default is 1.
#' @param cor A logical value indicating whether the calculation should use
#' the correlation matrix or the covariance matrix. (The correlation matrix
#' can only be used if there are no constant variables.)
#'
#' @return \code{do_pca} returns a list with class \code{princomp}.
#'
#' @examples
#' data(path_vals)
#' pca_model <- do_pca(path_vals[seq_len(ncol(path_vals)),])
#'
#' @export
#' @importFrom stats princomp
#' @importFrom methods is
#'
do_pca <- function(data, sel_assay = 1, cor = FALSE){
    if(is(data, "SummarizedExperiment")){
        data <- assay(data, sel_assay)
    }else if(is(data, "matrix")){
        data <- data
    }else{
        stop("Only SummarizedExperiment or matrix classes accepted as data")
    }
    fit <- stats::princomp(t(data), cor = cor)
    fit$var <- fit$sdev^2
    fit$explain_var <- fit$var/sum(fit$var)
    fit$acum_explain_var <- cumsum(fit$explain_var)
    return(fit)
}


compute_difexp <- function(vals, group1_label, group2_label, groups){

    g1_indexes <- which(groups == group1_label)
    g2_indexes <- which(groups == group2_label)
    grupo1 <- (groups[c(g1_indexes, g2_indexes)] == group1_label) + 0
    grupo2 <- (groups[c(g1_indexes, g2_indexes)] == group2_label) + 0
    design <- cbind(grupo1, grupo2)

    fit <- limma::lmFit(vals[,c(g1_indexes, g2_indexes)], design)
    cont_matrix <- limma::makeContrasts(grupo1 - grupo2, levels = design)
    fit2 <- limma::contrasts.fit(fit, cont_matrix)
    fit2 <- limma::eBayes(fit2)
    result <- data.frame(statistic = as.numeric(fit2$t),
                         p.value = as.numeric(fit2$p.value),
                         laterality = as.factor(fit2$t>0))
    rownames(result) <- rownames(fit2)

    return(result)
}



#' Compute pathway summary
#'
#' Computes a summary of the results, summarizing the number and proportion
#' of up- and down-regulated subpathways in each pathway.
#'
#' @param comp Comparison data frame as returned by the \code{do_wilcoxon}
#' function.
#' @param metaginfo Pathways object
#' @param conf Level of significance of the comparison for the adjusted
#' p-value. Default is 0.05.
#'
#' @return Table with the summarized information for each of the pathways.
#' Rows are the analized pathways. Columns are:
#' * \code{num_total_paths} Number of total subpathways in which each pathway
#' is decomposed.
#' * \code{num_significant_paths} Number of significant subpathways in the
#' provided comparison.
#' * \code{percent_significant_paths} Percentage of significant subpathways
#' from the total number of subpathways in a pathway.
#' * \code{num_up_paths} Number of significant up-regulated subpathways in the
#' provided comparison.
#' * \code{percent_up_paths} Percentage of significant up-regulated subpathways
#' from the total number of subpathways in a pathway.
#' * \code{num_down_paths} Number of significant down-regulated subpathways in
#' the provided comparison.
#' * \code{percent_down_paths} Percentage of significant down-regulated
#' subpathways from the total number of subpathways in a pathway.
#'
#' @examples
#' data(comp)
#' pathways <- load_pathways(species = "hsa", pathways_list = c("hsa03320",
#' "hsa04012"))
#' get_pathways_summary(comp, pathways)
#'
#' @export
#'
get_pathways_summary <- function(comp, metaginfo, conf = 0.05){
    comp$pathways <- sapply(strsplit(rownames(comp), split = "-"), "[[", 2)
    id_pathways <- unique(comp$pathways)
    name_pathways <- sapply(metaginfo$pathigraphs, "[[", "path.name")
    summ <- lapply(id_pathways, function(pathway){
        minicomp <- comp[comp$pathways == pathway,]
        num_total_paths <- nrow(minicomp)
        num_sig_paths <- sum(minicomp$FDRp.value < conf)
        percent_sig_paths <- round(num_sig_paths/nrow(minicomp)*100, digits = 2)
        is_up_path <- minicomp$FDRp.value < conf & minicomp$`UP/DOWN` == "UP"
        num_up_paths <- sum(is_up_path)
        percent_up_paths <- round(num_up_paths/nrow(minicomp) * 100, digits = 2)
        is_down_path <- minicomp$FDRp.value < conf &
            minicomp$`UP/DOWN` == "DOWN"
        num_down_paths <- sum(is_down_path)
        percent_down_paths <- round(num_down_paths/nrow(minicomp) * 100,
                                    digits = 2)
        data.frame(num_total_paths = num_total_paths,
                   num_significant_paths = num_sig_paths,
                   percent_significant_paths = percent_sig_paths,
                   num_up_paths = num_up_paths,
                   percent_up_paths = percent_up_paths,
                   num_down_paths = num_down_paths,
                   percent_down_paths = percent_down_paths)
    })
    summ <- do.call("rbind", summ)
    rownames(summ) <- name_pathways
    summ <- cbind(id_pathways, summ)
    summ <- summ[order(summ$percent_significant_paths, decreasing = TRUE),]
    return(summ)
}


#' Computes pathway significance
#'
#' Performs a test for each pathway checking if the number of significant
#' paths is significant, compared to not having any of the paths as significant.
#'
#' @param comp Comparison data frame as returned by the \code{do_wilcoxon}
#' function.
#'
#' @return 
#' Table with the names of the pathways and their p-value for the Fisher test
#' comparing the proportion of significant subpaths vs. 0.
#'
#' @examples
#' data(comp)
#' top_pathways(comp)
#'
#' @export
#'
top_pathways <- function(comp){
    
    path_names <- as.character(comp$path_names)
    comp$pathways <- sapply(strsplit(path_names, split = ":"), "[[", 1)
    pathways <- unique(comp$pathways)
    
    tests <- do.call(rbind, lapply(pathways, function(path) {
        t1 <- table(comp[comp$pathways == path,"FDRp.value"] < 0.05)
        t2 <- c(sum(t1), 0)
        
        ft <- fisher.test(rbind(t1, t2))
        data.frame(pathway = path, pval = ft$p.value, stringsAsFactors = F)
    }))
    
    tests <- tests[order(tests$pval),]
    
    return(tests)
}
