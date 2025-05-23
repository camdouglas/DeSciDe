#' @importFrom STRINGdb STRINGdb
#' @importFrom data.table data.table
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom dplyr group_by mutate arrange desc slice
#' @importFrom stats setNames
#' @importFrom grDevices dev.off pdf
#' @importFrom magrittr %>%
#' @importFrom utils head str
NULL

# Declare global variables
utils::globalVariables(c(
  "Gene", "PubMed_Rank", "Connectivity_Rank", "Category", "Gene_Symbol", "Degree",
  "Clustering_Coefficient_Percent", "Total", "Term", "Count", "gene", "."
))

#' Search STRING Database
#'
#' Search the STRING database for protein interactions.
#'
#' @param genes_list A list of gene IDs.
#' @param species The NCBI taxon ID of the species. Defaults to 9606 (Homo sapiens).
#' @param network_type The type of network to use, either "full" or "physical". Defaults to "full".
#' @param score_threshold The minimum score threshold for string interactions. Defaults to 400.
#' @return A list containing the following elements:
#'   \describe{
#'     \item{string_results}{A data frame with STRING interaction metrics.}
#'     \item{string_db}{The STRINGdb object used.}
#'     \item{string_ids}{The STRING IDs for the input genes.}
#'   }
#' @examples
#' \dontrun{
#' library(STRINGdb)
#' genes <- c("TP53", "BRCA1")
#' results <- search_string_db(genes)
#' print(results)
#' }
#' @export
search_string_db <- function(genes_list, species = 9606, network_type = "full", score_threshold = 400) {
  if (length(genes_list) == 0) {
    warning("No genes provided for STRING database search.")
    return(list(string_results = data.frame(), string_db = NULL, string_ids = NULL))
  }

  string_db <- STRINGdb$new(species = species, score_threshold = score_threshold, input_directory = "", network_type = network_type, version = "12.0")

  mapped_genes <- string_db$map(data.frame(gene = genes_list), "gene")
  if (nrow(mapped_genes) == 0) {
    warning("No valid genes found in STRING database for provided genes_list.")
    return(list(string_results = data.frame(), string_db = string_db, string_ids = NULL))
  }

  unique_mapped_genes <- mapped_genes %>% group_by(gene) %>% slice(1)
  string_ids <- unique_mapped_genes$STRING_id

  interactions <- string_db$get_interactions(string_ids)
  interaction_pairs <- data.table(proteinA = pmin(interactions$from, interactions$to), proteinB = pmax(interactions$from, interactions$to))
  interaction_pairs <- unique(interaction_pairs)
  interaction_pairs <- interaction_pairs[interaction_pairs$proteinA != interaction_pairs$proteinB, ]

  if (nrow(interaction_pairs) == 0) {
    warning("No interactions found for the provided genes in STRING database.")
    return(list(string_results = data.frame(), string_db = string_db, string_ids = string_ids))
  }

  nodes <- unique(c(interaction_pairs$proteinA, interaction_pairs$proteinB, string_ids))
  adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))

  for (i in 1:nrow(interaction_pairs)) {
    adjacency_matrix[interaction_pairs$proteinA[i], interaction_pairs$proteinB[i]] <- 1
    adjacency_matrix[interaction_pairs$proteinB[i], interaction_pairs$proteinA[i]] <- 1
  }

  degree <- rowSums(adjacency_matrix)
  clustering_coefficients <- numeric(length(nodes))
  for (i in seq_along(nodes)) {
    neighborhood <- which(adjacency_matrix[i, ] == 1)
    k <- length(neighborhood)
    if (k >= 2) {
      subgraph <- adjacency_matrix[neighborhood, neighborhood]
      triangles <- sum(subgraph) / 2
      clustering_coefficients[i] <- 2 * triangles / (k * (k -1))
    }
  }

  graph_obj <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
  components <- components(graph_obj)

  gene_symbol_lookup <- setNames(unique_mapped_genes$gene, unique_mapped_genes$STRING_id)
  nodes_symbols <- sapply(nodes, function(node) {
    if(node %in% names(gene_symbol_lookup)) {
      gene_symbol_lookup[node]
    } else {
      # Return the STRING ID instead of NA
      node
    }
  })

  string_results <- data.frame(
    Gene_Symbol = nodes_symbols,
    Degree = degree,
    Clustering_Coefficient_Percent = clustering_coefficients * 100,
    Clustering_Coefficient_Fraction = sapply(1:length(nodes), function(i) {
      k <- degree[i]
      if (k < 2) return("0 / 0")
      else {
        max_possible_edges <- k * (k - 1) / 2
        actual_edges <- clustering_coefficients[i] * max_possible_edges
        return(paste(round(actual_edges), "/", max_possible_edges, sep = ""))
      }
    }),
    Connected_Component_id = as.numeric(components$membership),
    Nodes_in_Connected_Component = as.numeric(components$csize[components$membership]),
    total_number_of_connected_components = components$no,
    row.names = NULL  # Add this line
  )

  string_results <- string_results %>%
    arrange(desc(Degree), desc(Clustering_Coefficient_Percent)) %>%
    mutate(Connectivity_Rank = row_number())

  return(list(string_results = string_results, string_db = string_db, string_ids = string_ids))
}

#' Plot STRING Network
#'
#' Plot STRING network interactions using STRINGdb.
#'
#' @param string_db A STRINGdb object.
#' @param string_ids A list of STRING IDs.
#' @param file_directory Directory for saving the output plot. Defaults to NULL.
#' @param export Logical indicating whether to export the plot. Defaults to FALSE.
#' @return Invisibly returns NULL.
#' @examples
#' library(STRINGdb)
#' string_db <- STRINGdb$new(species = 9606)
#' string_ids <- c("9606.ENSP00000269305", "9606.ENSP00000357940")
#' plot_string_network(string_db, string_ids, file_directory = tempdir(), export = FALSE)
#' @export
plot_string_network <- function(string_db, string_ids, file_directory = NULL, export = FALSE) {
  if (is.null(string_db) || is.null(string_ids) || length(string_ids) == 0) {
    warning("No valid STRING data available to plot network.")
    return(invisible(NULL))
  }

  current_date <- Sys.Date()
  formatted_date <- format(current_date, "%m.%d.%Y")

  if (export && !is.null(file_directory)) {
    network_output_filename <- paste(formatted_date, "STRINGdb_Network.pdf", sep = "_")
    full_network_output_path <- file.path(file_directory, network_output_filename)
    pdf(file = full_network_output_path, width = 12, height = 12)
    string_db$plot_network(string_ids)
    dev.off()
    message(paste("Network plot exported to:", full_network_output_path))
  } else {
    string_db$plot_network(string_ids)
  }

  return(invisible(NULL))
}

#' Plot STRING Interactions
#'
#' Plot STRING interactions degree vs. clustering.
#'
#' @param string_results Data frame with STRING metrics.
#' @param file_directory Directory for saving the output plot. Defaults to NULL.
#' @param export Logical indicating whether to export the plot. Defaults to FALSE.
#' @return Invisibly returns the ggplot object.
#' @examples
#' # Example data frame
#' string_results <- data.frame(Degree = c(10, 5), Clustering_Coefficient_Percent = c(20, 10))
#' plot_clustering(string_results, file_directory = tempdir(), export = FALSE)
#' @export
plot_clustering <- function(string_results, file_directory = NULL, export = FALSE) {
  current_date <- Sys.Date()
  formatted_date <- format(current_date, "%m.%d.%Y")

  if (is.null(string_results) || !all(c("Degree", "Clustering_Coefficient_Percent") %in% names(string_results))) {
    warning("Essential columns missing in string_results")
    return(invisible(NULL))
  }

  if (nrow(string_results) == 0) {
    warning("No data available for clustering plot.")
    return(invisible(NULL))
  }

  string_results$Degree <- as.numeric(string_results$Degree)
  string_results$Clustering_Coefficient_Percent <- as.numeric(string_results$Clustering_Coefficient_Percent)

  plot <- ggplot(string_results, aes(x = Degree, y = Clustering_Coefficient_Percent)) +
    geom_point(color = "black", alpha = 1) +
    theme_minimal() +
    labs(title = "Connectivity of Query",
         x = "Degree (Number of Neighbors)",
         y = "Clustering Coefficient (%)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

  if (export && !is.null(file_directory)) {
    scatter_output_filename <- paste(formatted_date, "Degree_vs_ClusteringCoefficient.png", sep = "_")
    full_scatter_output_path <- file.path(file_directory, scatter_output_filename)
    png(filename = full_scatter_output_path, width = 1000, height = 800, res = 150)
    print(plot)
    dev.off()
    message(paste("Scatter plot file path:", full_scatter_output_path))
  } else {
    print(plot)
  }

  return(invisible(plot))
}
