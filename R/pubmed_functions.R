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
#' @return An integer representing the number of PubMed articles found from the search query in PubMed.
#' @examples
#' # Perform a PubMed search for gene 'TP53' with term 'cancer'
#' result <- single_pubmed_search("TP53", "cancer")
#' print(result)
#' @export
single_pubmed_search <- function(gene, term) {
  query <- paste0('"', gene, "[TIAB]\"", " AND ", '"', term, "[TIAB]\"")
  single_search_results <- entrez_search(db = "pubmed", term = query)
  return(single_search_results$count)
}

#' Rank Search Results
#'
#' Rank search results based on a chosen method.
#'
#' @param data A data frame containing search results.
#' @param terms_list A list of search terms.
#' @param rank_method The method to rank pubmed results, either "weighted" or "total". Weighted ranks results based on order of terms inputted. Total ranks results on total sum of publications across all search term combinations. Defaults to "weighted".
#' @return A data frame with ranked search results, which includes the genes and their corresponding ranks based on the search method.
#' @examples
#' # Example data frame
#' data <- data.frame(Gene = c("Gene1", "Gene2"),
#'                    Term1 = c(10, 20),
#'                    Term2 = c(5, 15))
#' terms_list <- c("Term1", "Term2")
#' ranked_results <- rank_search_results(data, terms_list, rank_method = "weighted")
#' print(ranked_results)
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
#' @param verbose Logical flag indicating whether to display messages. Default is TRUE.
#' @return A data frame with search results, including genes, terms, and their corresponding publication counts and ranks.
#' @examples
#' genes <- c("TP53", "BRCA1")
#' terms <- c("cancer", "tumor")
#' search_results <- search_pubmed(genes, terms, rank_method = "weighted", verbose = FALSE)
#' print(search_results)
#' @export
search_pubmed <- function(genes_list, terms_list, rank_method = "weighted", verbose = TRUE) {
  single_search_results <- data.frame(Gene = character(), Term = character(), Count = integer())

  for (gene in genes_list) {
    for (term in terms_list) {
      if (verbose) {
        message(paste("Searching PubMed for gene:", gene, "and term:", term))
      }
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
