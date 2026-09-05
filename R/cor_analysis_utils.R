#' @title Internal utility function for correlation analysis
#' @description This is an internal function used by other correlation analysis functions.
#' It calculates correlations between a target variable and multiple variables across datasets.
#' @import dplyr psych stats
#' @param df A dataframe containing the data to be analyzed.
#' @param target_col The column name of the target variable for correlation.
#' @param sig_cols A vector of column names to correlate with the target variable.
#' @param group_col The column name to group the data by (default: "dataset").
#' @param cor_method The correlation method to use ("pearson", "spearman", or "kendall").
#' @param min_samples The minimum number of complete (non-missing) sample pairs required to compute a correlation for a gene in a group (default: 4). Pairs below this threshold are returned as NA.
#' @param adjust_method The method for adjusting p-values for multiple testing. Options include "none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni". Default is "BH".
#' @param conf_level The confidence level for confidence interval calculation (default: 0.95).
#' @return A list containing correlation matrix (r), p-value matrix (p), adjusted p-value matrix (p_adj), 
#'         per-pair sample size matrix (n), t-statistic matrix (t), confidence interval lower bound matrix (ci_lower), 
#'         confidence interval upper bound matrix (ci_upper), and split data (sss).
#'         All matrices keep every group x gene pair (genes in rows, groups in columns, with
#'         identical dimnames); pairs that could not be computed are returned as NA rather than
#'         removing the whole group from the results.
#' @keywords internal
.calculate_correlation <- function(df, target_col, sig_cols, 
                                   group_col = "dataset", cor_method = "pearson",
                                   min_samples = 4, adjust_method = "BH",
                                   conf_level = 0.95) {
  
  # Input validation
  if (is.null(df)) {
    warning("Input dataframe is NULL")
    return(NULL)
  }
  
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  
  if (!target_col %in% colnames(df)) {
    stop(paste("target_col", target_col, "not found in dataframe"))
  }
  
  missing_cols <- setdiff(sig_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing from dataframe:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  if (!group_col %in% colnames(df)) {
    stop(paste("group_col", group_col, "not found in dataframe"))
  }
  
  if (!cor_method %in% c("pearson", "spearman", "kendall")) {
    stop("cor_method must be 'pearson', 'spearman', or 'kendall'")
  }
  
  if (!adjust_method %in% c("none", "BH", "BY", "holm", "hochberg", "hommel", "bonferroni")) {
    stop("adjust_method must be one of: 'none', 'BH', 'BY', 'holm', 'hochberg', 'hommel', 'bonferroni'")
  }
  
  # Split the dataframe by the grouping column
  sss <- split(df, df[[group_col]])
  groups <- names(sss)
  n_groups <- length(groups)
  n_genes <- length(sig_cols)

  # Initialize full result matrices (all datasets x all genes). Cells that
  # cannot be computed are left as NA; no row/column is removed afterwards, so
  # every matrix keeps the same dimensions and dimnames. This preserves all
  # computable results (e.g. when a gene is simply not measured in a dataset,
  # the other genes of that dataset are still reported).
  new_matrix <- function() {
    matrix(NA_real_, nrow = n_groups, ncol = n_genes,
           dimnames = list(groups, sig_cols))
  }
  rvalue <- new_matrix()
  pvalue <- new_matrix()
  pvalue_adj <- new_matrix()
  nvalue <- new_matrix()
  tvalue <- new_matrix()
  ci_lower <- new_matrix()
  ci_upper <- new_matrix()

  # Human-readable notes about pairs that could not be computed
  issues <- character()

  for (i in seq_along(groups)) {
    group_data <- sss[[i]]
    n_group <- nrow(group_data)
    if (n_group < min_samples) {
      issues <- c(issues, sprintf("%s: only %d sample(s) (< %d), no correlation computed",
                                  groups[i], n_group, min_samples))
      next
    }

    x <- as.numeric(group_data[[target_col]])

    # Number of complete (target, gene) pairs per gene column. A gene may be
    # entirely missing for a dataset (e.g. not measured on its platform); such
    # columns must not prevent the remaining genes from being analysed.
    n_pair <- vapply(sig_cols, function(gc) {
      sum(!is.na(x) & !is.na(as.numeric(group_data[[gc]])))
    }, integer(1))

    computable <- n_pair >= min_samples

    if (sum(computable) > 0) {
      y_ok <- as.data.frame(lapply(group_data[, sig_cols[computable], drop = FALSE],
                                   as.numeric))
      corr_result <- tryCatch(
        suppressWarnings(psych::corr.test(x = x, y = y_ok, method = cor_method)),
        error = function(e) e
      )
      if (inherits(corr_result, "error")) {
        issues <- c(issues, sprintf("%s: correlation calculation failed (%s)",
                                    groups[i], conditionMessage(corr_result)))
        next
      }
      r_ok <- as.numeric(corr_result$r)
      p_ok <- as.numeric(corr_result$p)
      # psych::corr.test returns $n as a scalar (complete data) or as a
      # pairwise-count matrix (missing data); normalise to per-pair counts
      if (is.matrix(corr_result$n)) {
        n_ok <- as.numeric(corr_result$n)
      } else {
        n_ok <- rep(as.numeric(corr_result$n), length.out = sum(computable))
      }
      # Sample size is reported per pair (number of complete observations);
      # pairs that did not yield a coefficient (zero variance, etc.) get NA.
      n_ok[!is.finite(r_ok)] <- NA_real_
      rvalue[i, computable] <- r_ok
      pvalue[i, computable] <- p_ok
      nvalue[i, computable] <- n_ok
    }

    # Apply multiple testing correction within each group. p.adjust() keeps NA
    # entries in place and only corrects the computable p-values.
    if (adjust_method != "none") {
      pvalue_adj[i, ] <- stats::p.adjust(pvalue[i, ], method = adjust_method)
    } else {
      pvalue_adj[i, ] <- pvalue[i, ]
    }

    # t-statistic (pearson) and confidence intervals are only defined for pairs
    # where a correlation coefficient was obtained
    ok <- is.finite(rvalue[i, ]) & is.finite(nvalue[i, ])
    if (any(ok)) {
      if (cor_method == "pearson") {
        denom <- pmax(1 - rvalue[i, ok]^2, 0)   # avoid negative rounding noise
        tvalue[i, ok] <- rvalue[i, ok] * sqrt((nvalue[i, ok] - 2) / denom)
      }
      for (j in which(ok)) {
        ci <- .calculate_correlation_ci(rvalue[i, j], nvalue[i, j], conf_level)
        ci_lower[i, j] <- ci[["lower"]]
        ci_upper[i, j] <- ci[["upper"]]
      }
    }

    # Explain which pairs of this dataset remained NA
    na_genes <- sig_cols[is.na(rvalue[i, ])]
    if (length(na_genes) > 0) {
      why <- vapply(na_genes, function(gc) {
        j <- which(sig_cols == gc)
        if (n_pair[j] == 0) {
          "not measured (no expression values)"
        } else if (n_pair[j] < min_samples) {
          sprintf("fewer than %d valid pairs (%d)", min_samples, n_pair[j])
        } else {
          "no variability (constant expression) or undefined"
        }
      }, character(1))
      issues <- c(issues, sprintf("%s (n = %d): %s", groups[i], n_group,
                                  paste(paste0(na_genes, " - ", why),
                                        collapse = "; ")))
    }
  }

  if (length(issues) > 0) {
    warning("Some correlations could not be computed and are returned as NA:\n",
            paste("-", issues, collapse = "\n"), call. = FALSE)
  }

  # Transpose to genes x datasets. No NA omission: all matrices keep identical
  # dimensions and dimnames.
  plist <- list(r = t(rvalue), p = t(pvalue), p_adj = t(pvalue_adj),
                n = t(nvalue), t = t(tvalue), ci_lower = t(ci_lower),
                ci_upper = t(ci_upper), sss = sss)

  return(plist)
}

#' @title Internal utility function to extract tumor subtypes
#' @description This is an internal helper function to extract and validate tumor subtypes.
#' @param subtype_data The subtype mapping data.
#' @param tumor_subtype The tumor subtypes to extract.
#' @return A vector of extracted subtypes.
#' @keywords internal
.extract_tumor_subtype <- function(subtype_data, tumor_subtype) {
  if (is.null(tumor_subtype)) {
    return(NULL)
  }
  
  if (!is.character(tumor_subtype)) {
    stop("tumor_subtype must be a character vector")
  }
  
  extracted <- extract_subset(subtype_data, tumor_subtype)
  
  if (length(extracted) == 0) {
    warning("No matching subtypes found for the given tumor_subtype")
    return(NULL)
  }
  
  return(extracted)
}

#' @title Internal utility function to determine sample type
#' @description This is an internal helper function to determine sample types (Tumor/Normal).
#' @param df The dataframe containing sample data.
#' @param tumor_subtype The tumor subtypes to consider.
#' @param subtype_col The column name for subtype information.
#' @return The dataframe with an added "type" column.
#' @keywords internal
.determine_sample_type <- function(df, tumor_subtype = NULL, subtype_col = "subtype") {
  
  if (!subtype_col %in% colnames(df)) {
    stop(paste("subtype_col", subtype_col, "not found in dataframe"))
  }
  
  # Samples without subtype annotation cannot be assigned reliably. Without
  # tumor_subtype they would silently be classified as "Tumor" (and with
  # tumor_subtype they would be silently dropped), so flag them explicitly.
  n_missing_subtype <- sum(is.na(df[[subtype_col]]))
  
  if (is.null(tumor_subtype)) {
    df <- df %>%
      dplyr::mutate(type = ifelse(!!sym(subtype_col) %in% c("Normal", "Adjacent"), 
                                  "Normal", "Tumor"))
    if (n_missing_subtype > 0) {
      warning(n_missing_subtype, " sample(s) have missing subtype annotation and are ",
              "treated as 'Tumor' by default. If they should be excluded or assigned ",
              "differently, update/refresh the sample_subtype data.", call. = FALSE)
    }
  } else {
    df <- df %>%
      dplyr::filter(!!sym(subtype_col) %in% c(tumor_subtype, "Normal", "Adjacent")) %>%
      dplyr::mutate(type = ifelse(!!sym(subtype_col) %in% tumor_subtype, "Tumor", "Normal"))
    if (n_missing_subtype > 0) {
      warning(n_missing_subtype, " sample(s) have missing subtype annotation and were ",
              "excluded by the tumor_subtype = '", tumor_subtype, "' filter.",
              call. = FALSE)
    }
  }
  
  return(df)
}

#' @title Calculate Confidence Interval for Correlation Coefficient
#' @description Calculate confidence interval for a correlation coefficient using Fisher's z-transformation.
#' @param r The correlation coefficient.
#' @param n The sample size.
#' @param conf_level The confidence level (default: 0.95).
#' @return A vector containing the lower and upper bounds of the confidence interval.
#' @keywords internal
.calculate_correlation_ci <- function(r, n, conf_level = 0.95) {
  # Guard against undefined inputs: NA/non-finite r, sample size too small for
  # Fisher's z transformation (needs n > 3), or an invalid confidence level
  if (length(r) != 1 || is.na(r) || !is.finite(r) ||
      is.na(n) || !is.finite(n) || n <= 3 ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  
  # Perfect correlation: Fisher's z is undefined but the interval is degenerate
  if (abs(r) >= 1) {
    return(c(lower = r, upper = r))
  }
  
  # Fisher's z-transformation
  z <- 0.5 * log((1 + r) / (1 - r))
  se <- 1 / sqrt(n - 3)
  
  # Calculate confidence interval in z-space
  z_alpha <- qnorm((1 + conf_level) / 2)
  z_lower <- z - z_alpha * se
  z_upper <- z + z_alpha * se
  
  # Transform back to r-space
  r_lower <- (exp(2 * z_lower) - 1) / (exp(2 * z_lower) + 1)
  r_upper <- (exp(2 * z_upper) - 1) / (exp(2 * z_upper) + 1)
  
  return(c(lower = r_lower, upper = r_upper))
}

#' @title Calculate Effect Size for Differential Expression
#' @description Calculate various effect size metrics for differential expression analysis.
#' @param logFC The log fold change.
#' @param se The standard error of the logFC.
#' @param mean_tumor The mean expression in tumor samples.
#' @param mean_normal The mean expression in normal samples.
#' @param sd_tumor The standard deviation in tumor samples.
#' @param sd_normal The standard deviation in normal samples.
#' @return A list containing effect size metrics (Cohen's d, Hedges' g, and Glass's delta).
#' @keywords internal
.calculate_effect_size <- function(logFC, se = NULL, mean_tumor = NULL, mean_normal = NULL, 
                                   sd_tumor = NULL, sd_normal = NULL) {
  results <- list()
  
  # Cohen's d approximation from logFC and SE
  if (!is.null(logFC) && !is.null(se)) {
    results$cohens_d_approx <- abs(logFC) / se
  }
  
  # Cohen's d using means and pooled SD
  if (!is.null(mean_tumor) && !is.null(mean_normal) && !is.null(sd_tumor) && !is.null(sd_normal)) {
    pooled_sd <- sqrt((sd_tumor^2 + sd_normal^2) / 2)
    results$cohens_d <- abs(mean_tumor - mean_normal) / pooled_sd
    
    # Hedges' g (corrected for small sample bias)
    results$hedges_g <- results$cohens_d * (1 - 3 / (4 * (length(logFC) * 2) - 9))
    
    # Glass's delta (using control group SD)
    results$glass_delta <- abs(mean_tumor - mean_normal) / sd_normal
  }
  
  # Interpret effect size
  if (!is.null(results$cohens_d)) {
    results$interpretation <- .interpret_effect_size(results$cohens_d)
  } else if (!is.null(results$cohens_d_approx)) {
    results$interpretation <- .interpret_effect_size(results$cohens_d_approx)
  }
  
  return(results)
}

#' @title Interpret Effect Size
#' @description Provide interpretation of effect size based on Cohen's guidelines.
#' @param effect_size The effect size (Cohen's d).
#' @return A character string describing the effect size magnitude.
#' @keywords internal
.interpret_effect_size <- function(effect_size) {
  effect_size <- abs(effect_size)
  if (effect_size < 0.2) {
    return("negligible (d < 0.2)")
  } else if (effect_size < 0.5) {
    return("small (0.2 <= d < 0.5)")
  } else if (effect_size < 0.8) {
    return("medium (0.5 <= d < 0.8)")
  } else {
    return("large (d >= 0.8)")
  }
}