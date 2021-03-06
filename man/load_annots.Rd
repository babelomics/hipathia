% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{load_annots}
\alias{load_annots}
\title{Loads functional annotations to genes}
\usage{
load_annots(db, species)
}
\arguments{
\item{db}{Database to be used. Either "GO" or "uniprot".}

\item{species}{Species of the samples.

#@examples
#load_annots("GO", "hsa")}
}
\value{
Functional annotations from HGNC to the selected database.
}
\description{
Loads functional annotations from HGNC to the selected database.
}
