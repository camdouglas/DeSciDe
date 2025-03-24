#' Run the Entire Analysis
#'
#' Run the entire analysis pipeline including PubMed search, STRING database search, and plotting.
#'
#' @param genes_list A list of genes.
#' @param terms_list A list of terms.
#' @param file_directory Directory for saving the output files.
#' @param export_format Format for export, either "csv", "tsv", or "excel".
#' @param threshold_percentage Percentage threshold for ranking (default is 20%).
#' @param species The NCBI taxon ID of the species. Defaults to 9606 (Homo sapiens).
#' @param network_type The type of network to use, either "full" or "physical". Defaults to "full.
#' @param score_threshold The minimum score threshold for interactions. Defaults to 400.
#' @param rank_method The method to rank results, either "weighted" or "total". Defaults to "weighted".
#' @export
descide <- function(genes_list, terms_list, file_directory,
                    export_format = "csv",
                    threshold_percentage = 20,
                    species = 9606,
                    network_type = "full",
                    score_threshold = 400,
                    rank_method = "weighted") {

  log_message <- function(message) {
    cat(paste0(Sys.time(), ": ", message, "\n"))
  }

  log_message("Starting analysis pipeline")

  # Automatically set current_date to the system date and format it
  current_date <- Sys.Date()
  formatted_date <- format(current_date, "%m.%d.%Y")

  # Step 1: Perform PubMed search
  log_message("Performing PubMed search")
  pubmed_search_results <- search_pubmed(genes_list, terms_list, rank_method)

  log_message("PubMed search completed. Results:")
  print(pubmed_search_results)

  # Step 2: Plot heatmap of PubMed search results
  log_message("Plotting heatmap of PubMed search results")
  plot_heatmap(pubmed_search_results, file_directory)

  # Step 3: Perform STRING database search
  log_message("Performing STRING database search")
  string_db_results <- search_string_db(genes_list, species, network_type, score_threshold)

  string_results <- string_db_results$string_results
  string_db <- string_db_results$string_db
  string_ids <- string_db_results$string_ids

  log_message("STRING database search completed. Results:")
  print(head(string_results))

  # Step 4: Plot STRING network interactions
  log_message("Plotting STRING network interactions")
  plot_string_network(string_db, string_ids, file_directory)

  # Step 5: Combine PubMed and STRING summaries
  log_message("Combining summaries")
  combine_summary(pubmed_search_results, string_results, file_directory, export_format)

  # Step 6: Plot degree vs. clustering coefficient
  log_message("Plotting clustering coefficient")
  print(paste("Data passed to plot_clustering function"))
  print(str(string_results))
  plot_clustering(string_results, file_directory)

  # Step 7: Categorize and plot genes
  log_message("Categorizing and plotting genes")
  categorize_and_plot_genes(string_results, pubmed_search_results, file_directory, threshold_percentage)

  log_message("FINISHED: Analysis pipeline completed")

  combined_summary_filename <- paste(formatted_date, "Combined_Summary", sep = "_")
  combined_summary_filename <- switch(
    export_format,
    "csv" = paste0(combined_summary_filename, ".csv"),
    "tsv" = paste0(combined_summary_filename, ".tsv"),
    "excel" = paste0(combined_summary_filename, ".xlsx")
  )

  full_combined_summary_path <- file.path(file_directory, combined_summary_filename)
  return(full_combined_summary_path)
}
