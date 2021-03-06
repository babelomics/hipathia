% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{load_pathways}
\alias{load_pathways}
\title{Loads the pathways object.}
\usage{
load_pathways(species, pathways_list = NULL)
}
\arguments{
\item{species}{Species of the samples.}

\item{pathways_list}{Vector of the IDs of the pathways to load. By default
all available pathways are load.}
}
\value{
An pathways object including
* \code{species} Species to which the pathways are related.
* \code{pathigraphs} List of Pathigraph objects. Each Pathigraph contains
the necessary information of a pathway for it to be analyzed
with \code{Hipathia}.
* \code{all_genes} List of all the genes included in the selection of
pathways stored in \code{pathigraphs}.
* \code{eff_norm} Vector of normalization values for effector subpathways.
* \code{path_norm} Vector of normalization values for decomposed
subpathways.
}
\description{
Loads the pathways object, which includes information about the pathways
to be analyzed.
}
\details{
The object of pathways includes information about the pathways and the
subpathways which will be analyzed. This object must be provided to some
of the functions (like \code{hipathia} or \code{quantify_terms}) in the
package. These functions will analyze all the pathways included in this
object. By default, all available pathways are load. In order to restrict
the analysis to a predefined set of pathways, specify the set of pathways
to load with the parameter \code{pathways_list}.
}
\examples{
\dontrun{pathways <- load_pathways("hsa")   # Loads all pathways for human}
pathways <- load_pathways("mmu", c("mmu03320", "mmu04024", "mmu05200"))
   # Loads pathways 03320, 04024 and 05200 for mouse

}
