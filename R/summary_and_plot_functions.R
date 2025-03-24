# Import required functions and operators
#' @importFrom dplyr left_join mutate arrange desc row_number select filter case_when all_of across
#' @importFrom ggplot2 aes geom_point geom_text scale_color_manual theme_minimal labs theme element_text element_blank element_line element_rect ggplot
#' @importFrom grDevices dev.off png
#' @importFrom stats setNames
#' @importFrom utils write.csv write.table
#' @importFrom magrittr %>%
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
NULL

# Declare global variables
utils::globalVariables(c(
  "Gene", "PubMed_Rank", "Connectivity_Rank", "Category", "Gene_Symbol", "Degree",
  "Clustering_Coefficient_Percent", "Total", "Term", "Count", "gene", "."
))

#' Combine PubMed and STRING Metrics
#'
#' Combine PubMed search summary and STRING gene metrics.
#'
#' @param pubmed_search_results Data frame with PubMed search results.
#' @param string_results Data frame with STRING metrics.
#' @param file_directory Directory for saving the output summary.
#' @param export_format Format for export, either "csv", "tsv", or "excel".
#' @param export Logical indicating whether to export the summary. Defaults to FALSE.
#' @return A data frame with combined summary.
#' @export
combine_summary <- function(pubmed_search_results, string_results, file_directory, export_format = "csv", export = FALSE) {
  current_date <- Sys.Date()
  formatted_date <- format(current_date, "%m.%d.%Y")

  colnames(pubmed_search_results)[1] <- "Gene"
  colnames(string_results)[1] <- "Gene"

  combined_summary <- pubmed_search_results %>%
    left_join(string_results, by = "Gene")

  if (export) {
    if (export_format == "csv") {
      summary_filename <- paste(formatted_date, "Combined_Summary.csv", sep = "_")
      full_summary_path <- file.path(file_directory, summary_filename)
      write.csv(combined_summary, file = full_summary_path, row.names = FALSE)
    } else if (export_format == "tsv") {
      summary_filename <- paste(formatted_date, "Combined_Summary.tsv", sep = "_")
      full_summary_path <- file.path(file_directory, summary_filename)
      write.table(combined_summary, file = full_summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
    } else if (export_format == "excel") {
      summary_filename <- paste(formatted_date, "Combined_Summary.xlsx", sep = "_")
      full_summary_path <- file.path(file_directory, summary_filename)
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Summary")
      openxlsx::writeData(wb, "Summary", combined_summary)
      openxlsx::saveWorkbook(wb, full_summary_path, overwrite = TRUE)
    } else {
      stop("Invalid export format. Choose 'csv', 'tsv', or 'excel'.")
    }
    print(paste("Summary table exported to:", full_summary_path))
  }

  return(invisible(combined_summary))
}

#' Categorize and Plot Genes
#'
#' Categorize genes and create a scatter plot of Connectivity Rank vs. PubMed Rank.
#'
#' @param string_results Data frame with STRING metrics.
#' @param pubmed_search_results Data frame with PubMed search results.
#' @param file_directory Directory for saving the output plot.
#' @param export Logical indicating whether to export the plot. Defaults to FALSE.
#' @param threshold_percentage Percentage threshold for ranking (default is 20%).
#' @export
categorize_and_plot_genes <- function(string_results, pubmed_search_results, file_directory, export = FALSE, threshold_percentage = 20) {
  current_date <- Sys.Date()
  formatted_date <- format(current_date, "%m.%d.%Y")

  combined_string_results <- string_results %>%
    left_join(select(pubmed_search_results, Gene, PubMed_Rank = PubMed_Rank), by = c("Gene_Symbol" = "Gene"))

  top_threshold <- ceiling(nrow(combined_string_results) * (threshold_percentage / 100))
  bottom_threshold <- nrow(combined_string_results) - top_threshold + 1

  combined_string_results <- combined_string_results %>%
    mutate(Category = case_when(
      Connectivity_Rank <= top_threshold & PubMed_Rank <= top_threshold ~ "High Connectivity - High Precedence",
      Connectivity_Rank <= top_threshold & PubMed_Rank >= bottom_threshold ~ "High Connectivity - Low Precedence",
      TRUE ~ "Other"
    ))

  print(paste("Number of High Connectivity - High Precedence genes:", sum(combined_string_results$Category == "High Connectivity - High Precedence")))
  print(paste("Number of High Connectivity - Low Precedence genes:", sum(combined_string_results$Category == "High Connectivity - Low Precedence")))

  plot <- ggplot(combined_string_results, aes(x = Connectivity_Rank, y = PubMed_Rank, color = Category)) +
    geom_point(alpha = 1) +
    scale_color_manual(values = c("High Connectivity - High Precedence" = "navy", "High Connectivity - Low Precedence" = "red", "Other" = "black")) +
    geom_text(data = combined_string_results %>% filter(Category != "Other"), aes(label = Gene_Symbol, color = Category), vjust = -0.5, hjust = -0.5, size = 3, show.legend = FALSE) +
    theme_minimal() +
    labs(title = "Connectivity vs. Precedence",
         x = "Connectivity Rank",
         y = "PubMed Rank",
         color = "Legend") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )

  if (export) {
    rank_scatter_output_filename <- paste(formatted_date, "Connectivity_vs_Precedence.png", sep = "_")
    full_rank_scatter_output_path <- file.path(file_directory, rank_scatter_output_filename)
    png(filename = full_rank_scatter_output_path, width = 1000, height = 800, res = 150)
    print(plot)
    dev.off()
    print("Plot exported and device closed.")
  } else {
    print(plot)
  }

  return(invisible(plot))
}
