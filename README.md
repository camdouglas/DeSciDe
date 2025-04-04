# DeSciDe: A Tool for unbiased gene selection by simple visual analysis.

DeSciDe (Deciphering Scientific Discoveries) is an R package designed for genomic and proteomic data analysis, PubMed search, protein interaction network visualization, and comprehensive data summarization.

Features
--------

- Perform PubMed searches for a list of genes and terms.
- Search the STRING database for protein interactions.
- Rank PubMed and STRING search results to identify high confidence and novel hits. 
- Visualize results using heatmaps and network plots.
- Summarize data from multiple sources.

Installation
------------
Before installing the `DeSciDe` package, make sure you have `BiocManager` installed:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("STRINGdb")
BiocManager::install("ComplexHeatmap")

You can install the development version of DeSciDe from GitHub using the `devtools` package:

    # Install devtools if you haven't already
    install.packages("devtools")

    # Install DeSciDe from GitHub
    devtools::install_github("camdouglas/DeSciDe")

Usage
-----

#### Basic Example

Below is a basic example of how to use DeSciDe to analyze a list of genes and terms:

    # Load the package
    library(DeSciDe)

    # Define your list of genes and terms
    genes_list <- c("JUN", "MYC", "HDAC1", "TRIM33")
    terms_list <- c("cancer", "romidepsin")

    # Run the analysis pipeline
    results <- descide(
      genes_list = genes_list, 
      terms_list = terms_list
    )

    # View PubMed search results
    head(results$summary_results)

### Parameters

- `genes_list`: A list of gene IDs.
- `terms_list`: A list of search terms.
- `rank_method`: Method to rank results ("weighted" or "total"). Defaults to "weighted".
- `species`: NCBI taxon ID of the species. Defaults to 9606 (Homo sapiens).
- `network_type`: Type of STRING network to use ("full" or "physical"). Defaults to "full".
- `score_threshold`: Minimum score threshold for STRING interactions. Defaults to 400.
- `threshold_percentage`: Percentage threshold for ranking. Defaults to 20%.
- `export`: Logical indicating whether to export the results. Defaults to FALSE.
- `file_directory`: Directory for saving the output files. Defaults to NULL.
- `export_format`: Format for export ("csv", "tsv", "excel"). Defaults to "csv".

Authors
-------

- **Cameron Douglas** - Initial work - [camerondouglas@ufl.edu](mailto:camerondouglas@ufl.edu)
