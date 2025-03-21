#' @importFrom ComplexHeatmap Heatmap HeatmapList draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom grDevices dev.off png
#' @importFrom tibble column_to_rownames
NULL

# Declare global variables
utils::globalVariables(c(
  "Gene", "PubMed_Rank", "Connectivity_Rank", "Category", "Gene_Symbol", "Degree",
  "Clustering_Coefficient_Percent", "Total", "Term", "Count", "gene", "."
))

#' Plot Heatmap
#'
#' Create and save a heatmap of the PubMed search results.
#'
#' @param pubmed_search_results A data frame containing raw search results with genes and terms.
#' @param file_directory Directory for saving the output plot.
#' @param current_date Current date for file naming.
#' @export
plot_heatmap <- function(pubmed_search_results, file_directory, current_date) {
  # Format the data for the heatmap
  heatmap_data <- pubmed_search_results %>%
    select(-Total, -PubMed_Rank) %>%
    column_to_rownames("Gene")

  column_min <- apply(heatmap_data, 2, min, na.rm = TRUE)
  column_max <- apply(heatmap_data, 2, max, na.rm = TRUE)
  color_scales <- lapply(seq_len(ncol(heatmap_data)), function(i) {
    colorRamp2(c(column_min[i], column_max[i]), c("white", "navy"))
  })

  output_filename <- paste(current_date, "PubMed_Heatmap.png", sep = "_")
  full_output_path <- file.path(file_directory, output_filename)

  png(filename = full_output_path, width = 800, height = 1200)
  heatmap_list <- HeatmapList()

  for (i in seq_len(ncol(heatmap_data))) {
    heatmap_list <- heatmap_list +
      Heatmap(heatmap_data[, i, drop = FALSE],
              col = color_scales[[i]],
              name = colnames(heatmap_data)[i],
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_names_rot = 0,
              row_names_side = "left",
              column_title_gp = gpar(fontface = "bold"),
              column_names_gp = gpar(fontface = "bold", just = "center"))
  }

  draw(heatmap_list, column_title = "PubMed Search Results", column_title_gp = gpar(fontface = "bold", fontsize = 20))
  dev.off()
  print(paste("Heatmap plot exported to:", full_output_path))
}
