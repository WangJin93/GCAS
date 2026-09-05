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
  
  # Test with tumor_subtype: rows not matching tumor_subtype/Normal/Adjacent
  # are filtered out, remaining rows are classified Tumor/Normal
  result <- .determine_sample_type(test_df, tumor_subtype = "LC")
  expect_equal(result$type, c("Normal", "Normal", "Tumor"))
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

# Regression test: a gene that cannot be measured in one group must not cause
# the whole group to disappear from the results (previously na.omit() deleted
# any group row containing an NA, discarding all other gene correlations of
# that group). Every group should be retained and uncomputable pairs returned
# as NA.
test_that(".calculate_correlation keeps all groups and returns NA per pair", {
  set.seed(123)
  test_df <- data.frame(
    ID = paste0("s", 1:21),
    dataset = rep(c("DS_A", "DS_B", "DS_C"), c(10, 8, 3)),
    target = rnorm(21),
    g1 = rnorm(21),
    g2 = c(rnorm(10), rep(NA_real_, 8), rnorm(3)),   # entirely missing in DS_B
    g3 = c(rnorm(10), rep(1, 8), rnorm(3))            # constant in DS_B
  )
  
  expect_warning(result <- .calculate_correlation(test_df, "target",
                                                  c("g1", "g2", "g3")),
                 "could not be computed")
  
  # All groups are kept as columns in every matrix
  expect_setequal(colnames(result$r), c("DS_A", "DS_B", "DS_C"))
  # All matrices share identical dimnames
  for (nm in c("p", "p_adj", "n", "t", "ci_lower", "ci_upper")) {
    expect_identical(dimnames(result$r), dimnames(result[[nm]]))
  }
  # Complete group is fully computed
  expect_true(all(is.finite(result$r[, "DS_A"])))
  # In DS_B, g2 (missing) and g3 (constant) are NA while g1 is computed
  expect_true(is.finite(result$r["g1", "DS_B"]))
  expect_true(is.na(result$r["g2", "DS_B"]))
  expect_true(is.na(result$r["g3", "DS_B"]))
  # DS_C has fewer than min_samples rows: all NA but still retained
  expect_true(all(is.na(result$r[, "DS_C"])))
  # No NaN/Inf leaks into any matrix
  expect_false(any(is.nan(result$r)))
  expect_false(any(is.nan(result$p)))
  expect_false(any(is.nan(result$ci_lower)))
})

# Regression test: confidence-interval helper must not produce NaN for the
# boundary cases that used to occur (n <= 3, perfect correlation, NA input)
test_that(".calculate_correlation_ci handles boundary cases", {
  expect_true(all(is.na(.calculate_correlation_ci(0.5, 3))))
  expect_true(all(is.na(.calculate_correlation_ci(NA_real_, 10))))
  expect_equal(.calculate_correlation_ci(1, 10), c(lower = 1, upper = 1))
  ci <- .calculate_correlation_ci(0.5, 100)
  expect_true(ci["lower"] < 0.5 && 0.5 < ci["upper"])
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