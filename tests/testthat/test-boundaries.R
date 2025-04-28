library(testthat)
library(DeSciDe)

test_that("descide function handles minimal input", {
  skip_on_cran()  # Skipping this test on CRAN

  # Capture messages and warnings
  captured_messages_warnings <- capture.output({
    suppressWarnings({
      results <- descide(
        genes_list = c("MYC"),
        terms_list = c("cancer"),
        export = FALSE
      )
    })
  }, type = "message")

  expect_type(results, "list")
  expect_s3_class(results$pubmed_results, "data.frame")
  expect_s3_class(results$string_results, "data.frame")
  expect_gt(nrow(results$pubmed_results), 0)
  expect_equal(nrow(results$string_results), 0)

  expect_warning(res <- plot_heatmap(results$pubmed_results), "Not enough data to create a meaningful heatmap.")
  expect_null(res)

  expect_warning(res <- plot_string_network(results$string_db, results$string_ids), "No valid STRING data available to plot network.")
  expect_null(res)

  expect_warning(res <- plot_clustering(results$string_results), "Essential columns missing in string_results")
  expect_null(res)

  expect_warning(res <- plot_connectivity_precedence(results$summary_results), "Not enough data to plot genes.")
  expect_null(res)
})

test_that("descide function handles empty input", {
  skip_on_cran()  # Skipping this test on CRAN

  # Capture messages and warnings
  captured_messages_warnings <- capture.output({
    suppressWarnings({
      results <- descide(
        genes_list = character(0),
        terms_list = character(0),
        export = FALSE
      )
    })
  }, type = "message")

  expect_type(results, "list")
  expect_s3_class(results$pubmed_results, "data.frame")
  expect_s3_class(results$string_results, "data.frame")
  expect_equal(nrow(results$pubmed_results), 0)
  expect_equal(nrow(results$string_results), 0)
  expect_equal(nrow(results$summary_results), 0)

  expect_warning(res <- plot_heatmap(results$pubmed_results), "No data available to plot heatmap.")
  expect_null(res)

  expect_warning(res <- plot_string_network(results$string_db, results$string_ids), "No valid STRING data available to plot network.")
  expect_null(res)

  expect_warning(res <- plot_clustering(results$string_results), "Essential columns missing in string_results")
  expect_null(res)

  expect_warning(res <- plot_connectivity_precedence(results$summary_results), "Not enough data to plot genes.")
  expect_null(res)
})
