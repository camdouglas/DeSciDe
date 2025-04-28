#' @importFrom utils head str
NULL

#' Run DeSciDe pipeline
#'
#' Run the entire analysis pipeline including PubMed search, STRING database search, and plotting.
#'
#' @param genes_list A list of gene IDs.
#' @param terms_list A list of search terms.
#' @param rank_method The method to rank pubmed results, either "weighted" or "total". Weighted ranks results based on order of terms inputted. Total ranks results on total sum of publications across all search term combinations. Defaults to "weighted".
#' @param species The NCBI taxon ID of the species. Defaults to 9606 (Homo sapiens).
#' @param network_type The type of string network to use, either "full" or "physical". Defaults to "full".
#' @param score_threshold The minimum score threshold for string interactions. Defaults to 400.
#' @param threshold_percentage Percentage threshold for ranking (default is 20%).
#' @param export Logical indicating whether to export the results. Defaults to FALSE.
#' @param file_directory Directory for saving the output files. Defaults to NULL.
#' @param export_format Format for export, either "csv", "tsv", or "excel".
#' @return A list containing the PubMed search results, STRING results, and summary results.
#' @examples
#' \donttest{
#' genes <- c("TP53", "BRCA1")
#' terms <- c("cancer", "tumor")
#' results <- descide(genes, terms, export = FALSE)
#' str(results)
#' }
#' @export
descide <- function(genes_list, terms_list,
                    rank_method = "weighted",
                    species = 9606,
                    network_type = "full",
                    score_threshold = 400,
                    threshold_percentage = 20,
                    export = FALSE,
                    file_directory = NULL,
                    export_format = "csv") {

  log_message <- function(message) {
    message(paste0(Sys.time(), ": ", message))
  }

  log_message("Starting analysis pipeline")

  current_date <- Sys.Date()
  formatted_date <- format(current_date, "%m.%d.%Y")

  if (export && is.null(file_directory)) {
    stop("Export is set to TRUE, but file_directory is not provided.")
  }

  # Step 1: Perform PubMed search
  log_message("Performing PubMed search")
  pubmed_search_results <- search_pubmed(genes_list, terms_list, rank_method)

  log_message("PubMed search completed. Results:")
  message(head(pubmed_search_results))

  # Step 2: Plot heatmap of PubMed search results
  log_message("Plotting heatmap of PubMed search results")
  plot_heatmap(pubmed_search_results, file_directory, export)

  # Step 3: Perform STRING database search
  log_message("Performing STRING database search")
  string_db_results <- search_string_db(genes_list, species, network_type, score_threshold)

  string_results <- string_db_results$string_results
  string_db <- string_db_results$string_db
  string_ids <- string_db_results$string_ids

  log_message("STRING database search completed. Results:")
  message(head(string_results))

  # Step 4: Plot STRING network interactions
  log_message("Plotting STRING network interactions")
  plot_string_network(string_db, string_ids, file_directory, export)

  # Step 5: Combine PubMed and STRING summaries
  log_message("Combining summaries")
  combined_summary <- combine_summary(pubmed_search_results, string_results, file_directory, export_format, export, threshold_percentage)

  # Step 6: Plot degree vs. clustering coefficient
  log_message("Plotting clustering coefficient")
  plot_clustering(string_results, file_directory, export)

  # Step 7: Categorize and plot genes
  log_message("Categorizing and plotting genes")
  plot_connectivity_precedence(combined_summary, file_directory, export)

  log_message("FINISHED: Analysis pipeline completed")

  combined_summary_filename <- paste(formatted_date, "Combined_Summary", sep = "_")
  combined_summary_filename <- switch(
    export_format,
    "csv" = paste0(combined_summary_filename, ".csv"),
    "tsv" = paste0(combined_summary_filename, ".tsv"),
    "excel" = paste0(combined_summary_filename, ".xlsx")
  )

  full_combined_summary_path <- if (export) file.path(file_directory, combined_summary_filename) else NULL

  return(list(
    pubmed_results = pubmed_search_results,
    string_results = string_results,
    summary_results = if (export) full_combined_summary_path else combined_summary
  ))
}
