#' @title GCAS Package Unit Tests
#' @description Unit tests for the GCAS package functions.
#' @import testthat
#' @export

context("GCAS Package Tests")

# Test input validation functions
test_that("DEGs_analysis input validation works", {
  # Test NULL input
  expect_warning(result <- DEGs_analysis(NULL))
  expect_null(result)
  
  # Test non-dataframe input
  expect_error(DEGs_analysis("not a dataframe"))
  
  # Test dataframe without ID column
  test_df <- data.frame(Sample1 = rnorm(10), Sample2 = rnorm(10))
  expect_error(DEGs_analysis(test_df))
})

test_that("GSEA_analysis input validation works", {
  # Test NULL input
  expect_warning(result <- GSEA_analysis(NULL, gmt_file = BP_GMT_7.5.1))
  expect_null(result)
  
  # Test invalid data_type
  test_df <- data.frame(gene = c("Gene1", "Gene2"), r = c(0.5, 0.3))
  expect_error(GSEA_analysis(test_df, gmt_file = BP_GMT_7.5.1, data_type = "invalid"))
  
  # Test missing required columns
  test_df <- data.frame(gene = c("Gene1", "Gene2"))
  expect_error(GSEA_analysis(test_df, gmt_file = BP_GMT_7.5.1))
})

test_that("RRA_analysis input validation works", {
  # Test NULL input
  expect_warning(result <- RRA_analysis(NULL))
  expect_null(result)
  
  # Test non-list input
  expect_error(RRA_analysis("not a list"))
  
  # Test list with insufficient datasets
  expect_error(RRA_analysis(list("GSE1" = data.frame(gene = "Gene1", logFC = 1, P.Value = 0.01))))
})

# Test visualization utility functions
test_that("gcas_theme returns valid ggplot theme", {
  theme_result <- gcas_theme()
  expect_s3_class(theme_result, "theme")
})

test_that("gcas_palette returns valid colors", {
  # Test main palette
  colors <- gcas_palette(3)
  expect_length(colors, 3)
  expect_type(colors, "character")
  
  # Test diverging palette
  diverging <- gcas_palette(palette = "diverging")
  expect_length(diverging, 3)
  
  # Test sequential palette
  sequential <- gcas_palette(5, palette = "sequential")
  expect_length(sequential, 5)
})

# Test helper functions
test_that(".determine_sample_type works correctly", {
  test_df <- data.frame(
    subtype = c("Normal", "Tumor", "Adjacent", "LC"),
    value = rnorm(4)
  )
  
  # Test without tumor_subtype
  result <- .determine_sample_type(test_df)
  expect_equal(result$type, c("Normal", "Tumor", "Normal", "Tumor"))
  
  # Test with tumor_subtype
  result <- .determine_sample_type(test_df, tumor_subtype = "LC")
  expect_equal(result$type, c("Normal", NA, "Normal", "Tumor"))
})

# Test correlation utility function
test_that(".calculate_correlation input validation works", {
  # Test NULL input
  expect_warning(result <- .calculate_correlation(NULL, "target", c("var1")))
  expect_null(result)
  
  # Test missing columns
  test_df <- data.frame(target = rnorm(10), var1 = rnorm(10))
  expect_error(.calculate_correlation(test_df, "missing", c("var1")))
})

# Test API URL functions
test_that("API URL functions work", {
  # Test get_gcas_api_url
  url <- get_gcas_api_url()
  expect_type(url, "character")
  expect_true(grepl("https://", url))
  
  # Test set_gcas_api_url
  expect_silent(set_gcas_api_url("https://example.com/api"))
})

# Test cache directory functions
test_that("get_expr_data validates inputs", {
  expect_error(get_expr_data(123, "GAPDH"))
  expect_error(get_expr_data("GSE123", 123))
})

test_that("get_OSF_data validates inputs", {
  expect_error(get_OSF_data(123, "geo_data"))
  expect_error(get_OSF_data("GSE123", "invalid"))
})

context("Visualization Tests")

test_that("viz_TvsN handles NULL input", {
  expect_warning(result <- viz_TvsN(NULL))
  expect_null(result)
})

test_that("viz_TvsN validates df_type", {
  test_df <- data.frame(
    ID = paste0("Sample", 1:10),
    dataset = "GSE123",
    Gene1 = rnorm(10),
    subtype = rep(c("Normal", "Tumor"), 5)
  )
  
  # Test valid df_type
  expect_s3_class(viz_TvsN(test_df, df_type = "single"), "ggplot")
})

context("Data Processing Tests")

test_that("cor_gcas_TIL validates inputs", {
  expect_warning(result <- cor_gcas_TIL(NULL))
  expect_null(result)
  
  test_df <- data.frame(ID = "Sample1", dataset = "GSE123", Gene1 = 1.0, subtype = "Tumor")
  expect_silent(result <- cor_gcas_TIL(test_df))
})

test_that("cor_gcas_drug validates inputs", {
  expect_warning(result <- cor_gcas_drug(NULL))
  expect_null(result)
  
  test_df <- data.frame(ID = "Sample1", dataset = "GSE123", Gene1 = 1.0, subtype = "Tumor")
  expect_silent(result <- cor_gcas_drug(test_df))
})