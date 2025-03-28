#' @importFrom rentrez entrez_search
#' @importFrom dplyr mutate select arrange desc all_of row_number group_by summarise left_join across
#' @importFrom tidyr pivot_wider
#' @importFrom utils write.csv write.table
#' @importFrom magrittr %>%
NULL

# Declare global variables
utils::globalVariables(c(
  "Gene", "PubMed_Rank", "Connectivity_Rank", "Category", "Gene_Symbol", "Degree",
  "Clustering_Coefficient_Percent", "Total", "Term", "Count", "gene", "."
))

#' Search PubMed
#'
#' Perform a PubMed search for a given gene and term.
#'
#' @param gene A character string representing the gene symbol.
#' @param term A character string representing the search term.
#' @return An integer representing the number of PubMed articles found.
#' @export
single_pubmed_search <- function(gene, term) {
  query <- paste0(gene, "[TIAB] AND ", term, "[TIAB]")
  single_search_results <- entrez_search(db = "pubmed", term = query)
  return(single_search_results$count)
}

#' Rank Search Results
#'
#' Rank search results based on a chosen method.
#'
#' @param data A data frame containing search results.
#' @param terms_list A list of search terms.
#' @param rank_method The method to rank pubmed results, either "weighted" or "total". Weighted ranks results based on order of terms inputed. Total ranks results on total sum of publications across all search term combinations. Defaults to "weighted".
#' @return A data frame with ranked search results.
#' @export
rank_search_results <- function(data, terms_list, rank_method = "weighted") {
  if (rank_method == "weighted") {
    data <- data %>%
      mutate(Total = rowSums(select(., -Gene))) %>%
      select(Gene, all_of(terms_list), Total) %>%
      arrange(desc(across(all_of(terms_list)))) %>%
      mutate(PubMed_Rank = row_number())
  } else if (rank_method == "total") {
    data <- data %>%
      mutate(Total = rowSums(select(., -Gene))) %>%
      arrange(desc(Total)) %>%
      mutate(PubMed_Rank = row_number())
  } else {
    stop("Invalid rank_method. Choose either 'weighted' or 'total'.")
  }
  return(data)
}

#' Search PubMed with Multiple Genes and Terms
#'
#' Perform a PubMed search for multiple genes and terms.
#'
#' @param genes_list A list of gene IDs.
#' @param terms_list A list of search terms.
#' @param rank_method The method to rank results, either "weighted" or "total". Defaults to "weighted".
#' @return A data frame with search results.
#' @export
search_pubmed <- function(genes_list, terms_list, rank_method = "weighted") {
  single_search_results <- data.frame(Gene = character(), Term = character(), Count = integer())

  for (gene in genes_list) {
    for (term in terms_list) {
      count <- single_pubmed_search(gene, term)
      single_search_results <- rbind(single_search_results, data.frame(Gene = gene, Term = term, Count = count))
    }
  }

  aggregated_results <- single_search_results %>%
    group_by(Gene, Term) %>%
    summarise(Count = sum(Count), .groups = 'drop')

  pubmed_search_results <- aggregated_results %>%
    pivot_wider(names_from = Term, values_from = Count, values_fill = list(Count = 0))

  pubmed_search_results <- rank_search_results(pubmed_search_results, terms_list, rank_method)

  return(pubmed_search_results)
}
