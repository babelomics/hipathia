
#' BRCA gene expression dataset
#' 
#' Gene expression of 40 samples from the BRCA-US project from 
#' The Cancer Genome Atlas (TCGA). 
#' 
#' Gene expression matrix with 40 samples taken from the BRCA-US project from
#' The Cancer Genome Atlas (TCGA). The data has been log-transformed and 
#' normalized with TMM. 
#' 
#' @format Matrix with 40 columns and 18638 rows. Row names are Entrez IDs 
#' and column names are the  TCGA identifyers of the samples.
#' 
#' @return Matrix with 40 columns and 18638 rows. Row names are Entrez IDs 
#' and column names are the  TCGA identifyers of the samples.

#' @usage data(brca_data)
#' 
#' @source \url{https://cancergenome.nih.gov/} 
#' 
"brca_data"


#' BRCA experimental design
#' 
#' Experimental design of the gene expression matrix \code{brca_data} with
#' 40 samples taken from the BRCA-US project from The Cancer Genome Atlas
#' (TCGA). 20 samples are "Tumor" samples and 20 samples are "Normal" samples.
#' 
#' @format Dataframe with 1 column and 40 rows, including the experimental
#' design of the 40 samples from the BRCA-US project from TCGA. Field 
#' \code{group} is the type of sample, either "Tumor" or "Normal". 
#' 
#' @return Dataframe with 1 column and 40 rows, including the experimental
#' design of the 40 samples from the BRCA-US project from TCGA. Field 
#' \code{group} is the type of sample, either "Tumor" or "Normal". 
#' 
#' @usage data(brca_design)
#' 
#' @source \url{https://cancergenome.nih.gov/} 
#' 
"brca_design"
