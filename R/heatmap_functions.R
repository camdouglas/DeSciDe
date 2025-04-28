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
#' Create and optionally save a heatmap of the PubMed search results.
#'
#' @param pubmed_search_results A data frame containing raw search results with genes and terms.
#' @param file_directory Directory for saving the output plot. Defaults to NULL.
#' @param export Logical indicating whether to export the plot. Defaults to FALSE.
#' @return Invisibly returns a \code{HeatmapList} object.
#' @examples
#' # Example data frame
#' data <- data.frame(Gene = c("Gene1", "Gene2"),
#'                    Term1 = c(10, 20),
#'                    Term2 = c(5, 15),
#'                    Total = c(15, 35),
#'                    PubMed_Rank = c(1, 2))
#' plot_heatmap(data, file_directory = tempdir(), export = FALSE)
#' @export
plot_heatmap <- function(pubmed_search_results, file_directory = NULL, export = FALSE) {
  if (nrow(pubmed_search_results) == 0) {
    warning("No data available to plot heatmap.")
    return(invisible(NULL))
  }

  current_date <- Sys.Date()
  formatted_date <- format(current_date, "%m.%d.%Y")

  # Prepare heatmap data
  heatmap_data <- pubmed_search_results %>%
    select(-Total, -PubMed_Rank) %>%
    column_to_rownames("Gene")

  if (nrow(heatmap_data) < 2) {
    warning("Not enough data to create a meaningful heatmap.")
    return(invisible(NULL))
  }

  # Convert data frame to matrix
  heatmap_data <- as.matrix(heatmap_data)

  if (all(is.na(heatmap_data))) {
    warning("No valid data points in heatmap data.")
    return(invisible(NULL))
  }

  column_max <- apply(heatmap_data, 2, max, na.rm = TRUE)

  if (all(column_max == 0)) {
    warning("heatmap_data has zero values in each column; unable to create meaningful color scale.")
    return(invisible(NULL))
  }

  color_scales <- lapply(seq_len(ncol(heatmap_data)), function(i) {
    colorRamp2(c(0, column_max[i]), c("white", "navy"))
  })

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

  if (export && !is.null(file_directory)) {
    output_filename <- paste(formatted_date, "PubMed_Heatmap.png", sep = "_")
    full_output_path <- file.path(file_directory, output_filename)
    png(filename = full_output_path, width = 800, height = 1200)
    draw(heatmap_list, column_title = "PubMed Search Results", column_title_gp = gpar(fontface = "bold", fontsize = 20))
    dev.off()
    message("Heatmap plot exported to:", full_output_path)
  } else {
    draw(heatmap_list, column_title = "PubMed Search Results", column_title_gp = gpar(fontface = "bold", fontsize = 20))
  }

  return(invisible(heatmap_list))
}
