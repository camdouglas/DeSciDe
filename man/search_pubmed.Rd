% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pubmed_functions.R
\name{search_pubmed}
\alias{search_pubmed}
\title{Search PubMed with Multiple Genes and Terms}
\usage{
search_pubmed(genes_list, terms_list, rank_method = "weighted")
}
\arguments{
\item{genes_list}{A list of gene IDs.}

\item{terms_list}{A list of search terms.}

\item{rank_method}{The method to rank results, either "weighted" or "total". Defaults to "weighted".}
}
\value{
A data frame with search results.
}
\description{
Perform a PubMed search for multiple genes and terms.
}
