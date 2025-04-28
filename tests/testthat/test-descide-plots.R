# Load necessary packages
library(testthat)
library(DeSciDe)
library(withr)

test_that("descide function works correctly without export", {
  skip_on_cran()  # Skipping this test on CRAN

  # Define the genes and terms lists
  genes_list <- c("JUN", "MYC", "HDAC1", "TRIM33")
  terms_list <- c("cancer", "romidepsin")

  # Capture messages and warnings to check function execution
  f_capture_output <- capture.output({
    f_capture_message <- capture.output({
      results <- suppressWarnings({
        descide(
          genes_list = genes_list,
          terms_list = terms_list,
          export = FALSE
        )
      })
    }, type = "message")
  })

  # Check that results is a list
  expect_type(results, "list")

  # Check that pubmed_results, string_results, and summary_results are data frames
  expect_s3_class(results$pubmed_results, "data.frame")
  expect_s3_class(results$string_results, "data.frame")
  expect_s3_class(results$summary_results, "data.frame")

  # Ensure pubmed_results has rows
  expect_gt(nrow(results$pubmed_results), 0)

  # Ensure string_results has rows
  expect_gt(nrow(results$string_results), 0)

  # Ensure summary_results has rows
  expect_gt(nrow(results$summary_results), 0)

  # Assert that certain output messages were printed during the plot functions
  expect_true(any(grepl("Starting analysis pipeline", f_capture_message)))
  expect_true(any(grepl("Performing PubMed search", f_capture_message)))
  expect_true(any(grepl("Plotting heatmap of PubMed search results", f_capture_message)))
  expect_true(any(grepl("Performing STRING database search", f_capture_message)))
  expect_true(any(grepl("Plotting STRING network interactions", f_capture_message)))
  expect_true(any(grepl("Combining summaries", f_capture_message)))
  expect_true(any(grepl("Plotting clustering coefficient", f_capture_message)))
  expect_true(any(grepl("Categorizing and plotting genes", f_capture_message)))
  expect_true(any(grepl("FINISHED: Analysis pipeline completed", f_capture_message)))
})

test_that("plotting functions execute without error", {
  skip_on_cran()  # Skipping this test on CRAN

  # Define genes and terms lists
  genes_list <- c("JUN", "MYC", "HDAC1", "TRIM33")
  terms_list <- c("cancer", "romidepsin")

  # Capture messages and warnings to check function execution
  f_capture_output <- capture.output({
    results <- suppressWarnings({
      descide(
        genes_list = genes_list,
        terms_list = terms_list,
        export = FALSE
      )
    })
  }, type = "message")

  # Ensure plot functions execute without error
  expect_no_error(plot_heatmap(results$pubmed_results))
  expect_no_error({
    if (!is.null(results$string_ids) && length(results$string_ids) > 0){
      plot_string_network(results$string_db, results$string_ids)
    }
  })
  expect_no_error(plot_clustering(results$string_results))
  expect_no_error(plot_connectivity_precedence(results$summary_results))
})
