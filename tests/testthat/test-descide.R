# Load testthat and DeSciDe packages
library(testthat)
library(DeSciDe)

# Test the descide function
test_that("descide function works correctly", {
  # Define the genes and terms lists
  genes_list <- c("JUN", "MYC", "HDAC1", "TRIM33")
  terms_list <- c("cancer", "romidepsin")

  # Run the descide function
  results <- descide(
    genes_list = genes_list,
    terms_list = terms_list,
    export = FALSE
  )

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
})
