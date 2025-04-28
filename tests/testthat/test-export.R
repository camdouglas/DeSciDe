library(testthat)
library(DeSciDe)
library(withr)

test_that("descide function exports correctly", {
  skip_on_cran()  # Skipping this test on CRAN

  genes_list <- c("JUN", "MYC", "HDAC1", "TRIM33")
  terms_list <- c("cancer", "romidepsin")
  current_date <- format(Sys.Date(), "%m.%d.%Y")

  tmp_dir <- local_tempdir()

  # Capture messages and warnings to check function execution
  f_capture_output <- capture.output({
    results <- suppressWarnings({
      descide(
        genes_list = genes_list,
        terms_list = terms_list,
        export = TRUE,
        file_directory = tmp_dir,
        export_format = "csv"
      )
    })
  }, type = "message")

  # Check if the exported summary file exists
  expected_summary_file <- paste0(current_date, "_Combined_Summary.csv")
  expect_true(file.exists(file.path(tmp_dir, expected_summary_file)))

  # Check if the heatmap image exists
  expected_heatmap_file <- paste0(current_date, "_PubMed_Heatmap.png")
  expect_true(file.exists(file.path(tmp_dir, expected_heatmap_file)))

  # Check if the STRINGdb network image exists
  expected_string_network_file <- paste0(current_date, "_STRINGdb_Network.pdf")
  expect_true(file.exists(file.path(tmp_dir, expected_string_network_file)))

  # Check if the Degree vs Clustering Coefficient plot exists
  expected_clustering_plot_file <- paste0(current_date, "_Degree_vs_ClusteringCoefficient.png")
  expect_true(file.exists(file.path(tmp_dir, expected_clustering_plot_file)))

  # Check if the Connectivity vs Precedence plot exists
  expected_precedence_plot_file <- paste0(current_date, "_Connectivity_vs_Precedence.png")
  expect_true(file.exists(file.path(tmp_dir, expected_precedence_plot_file)))

  # Repeat the above tests for TSV and Excel export formats

  tmp_dir <- local_tempdir()

  # Capture messages and warnings to check function execution
  f_capture_output <- capture.output({
    results <- suppressWarnings({
      descide(
        genes_list = genes_list,
        terms_list = terms_list,
        export = TRUE,
        file_directory = tmp_dir,
        export_format = "tsv"
      )
    })
  }, type = "message")

  expected_summary_file <- paste0(current_date, "_Combined_Summary.tsv")
  expect_true(file.exists(file.path(tmp_dir, expected_summary_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_heatmap_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_string_network_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_clustering_plot_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_precedence_plot_file)))

  tmp_dir <- local_tempdir()

  # Capture messages and warnings to check function execution
  f_capture_output <- capture.output({
    results <- suppressWarnings({
      descide(
        genes_list = genes_list,
        terms_list = terms_list,
        export = TRUE,
        file_directory = tmp_dir,
        export_format = "excel"
      )
    })
  }, type = "message")

  expected_summary_file <- paste0(current_date, "_Combined_Summary.xlsx")
  expect_true(file.exists(file.path(tmp_dir, expected_summary_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_heatmap_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_string_network_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_clustering_plot_file)))
  expect_true(file.exists(file.path(tmp_dir, expected_precedence_plot_file)))
})
