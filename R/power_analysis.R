#' @title Statistical Power Analysis for Differential Expression Studies
#' @description Perform power analysis to determine statistical power or required sample size for detecting differentially expressed genes.
#' @import stats
#' @param effect_size The expected effect size (Cohen's d or standardized mean difference). For microarray data, typical values range from 0.2 (small) to 1.0 (large).
#' @param n1 Sample size in group 1 (e.g., Tumor).
#' @param n2 Sample size in group 2 (e.g., Normal). If NULL, assumed equal to n1.
#' @param alpha Significance level (Type I error rate). Default is 0.05.
#' @param power Desired statistical power (1 - Type II error rate). Default is 0.8.
#' @param test_type Type of statistical test: "t.test" (two-sample t-test) or "anova" (one-way ANOVA). Default is "t.test".
#' @param alternative Alternative hypothesis: "two.sided" (default), "less", or "greater".
#' @return A list containing power analysis results including calculated power, sample size requirements, and effect size interpretation.
#' @examples
#' \dontrun{
#' # Calculate power given sample sizes
#' power_result <- power_analysis(effect_size = 0.8, n1 = 50, n2 = 50)
#' 
#' # Calculate required sample size for 80% power
#' sample_size_result <- power_analysis(effect_size = 0.6, power = 0.8)
#' }
#' @export
power_analysis <- function(effect_size = NULL, n1 = NULL, n2 = NULL, alpha = 0.05, power = 0.8, 
                           test_type = "t.test", alternative = "two.sided") {
  
  # Input validation
  if (!test_type %in% c("t.test", "anova")) {
    stop("test_type must be 't.test' or 'anova'")
  }
  
  if (!alternative %in% c("two.sided", "less", "greater")) {
    stop("alternative must be 'two.sided', 'less', or 'greater'")
  }
  
  if (is.null(effect_size) && is.null(n1)) {
    stop("Either effect_size or n1 must be provided")
  }
  
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  
  if (power <= 0 || power >= 1) {
    stop("power must be between 0 and 1")
  }
  
  # Set n2 equal to n1 if not provided
  if (is.null(n2)) {
    n2 <- n1
  }
  
  # Determine what to calculate
  calculate_power <- !is.null(n1) && !is.null(effect_size)
  calculate_sample_size <- is.null(n1) && !is.null(effect_size)
  
  results <- list()
  
  if (test_type == "t.test") {
    if (calculate_power) {
      # Calculate power using non-central t-distribution
      df <- n1 + n2 - 2
      ncp <- effect_size * sqrt((n1 * n2) / (n1 + n2))
      
      # Calculate critical value
      if (alternative == "two.sided") {
        t_crit <- qt(1 - alpha/2, df)
        power_val <- 1 - pt(t_crit, df, ncp) + pt(-t_crit, df, ncp)
      } else if (alternative == "greater") {
        t_crit <- qt(1 - alpha, df)
        power_val <- 1 - pt(t_crit, df, ncp)
      } else {
        t_crit <- qt(alpha, df)
        power_val <- pt(t_crit, df, ncp)
      }
      
      results$power <- power_val
      results$n1 <- n1
      results$n2 <- n2
      results$effect_size <- effect_size
    } else if (calculate_sample_size) {
      # Calculate required sample size using iterative approach
      target_n <- .calculate_sample_size_t(effect_size, alpha, power, alternative)
      results$required_n_per_group <- target_n
      results$total_required_n <- target_n * 2
      results$effect_size <- effect_size
      results$target_power <- power
    }
  } else if (test_type == "anova") {
    if (calculate_power) {
      # For ANOVA with 2 groups
      df1 <- 1
      df2 <- n1 + n2 - 2
      ncp <- effect_size^2 * (n1 * n2) / (n1 + n2)
      f_crit <- qf(1 - alpha, df1, df2)
      power_val <- 1 - pf(f_crit, df1, df2, ncp)
      
      results$power <- power_val
      results$n_per_group <- n1
      results$effect_size <- effect_size
    } else if (calculate_sample_size) {
      target_n <- .calculate_sample_size_anova(effect_size, alpha, power)
      results$required_n_per_group <- target_n
      results$effect_size <- effect_size
      results$target_power <- power
    }
  }
  
  # Add effect size interpretation
  if (!is.null(effect_size)) {
    results$effect_size_interpretation <- interpret_effect_size(effect_size)
  }
  
  # Add alpha adjustment information for multiple testing
  results$alpha_info <- list(
    nominal_alpha = alpha,
    note = "Consider adjusting alpha for multiple testing. For example, with 10,000 genes tested, Bonferroni-adjusted alpha = alpha / 10000"
  )
  
  return(results)
}

#' @title Calculate Sample Size for t-test
#' @description Iteratively calculate required sample size for t-test.
#' @param d Effect size (Cohen's d).
#' @param alpha Significance level.
#' @param power Target power.
#' @param alternative Alternative hypothesis.
#' @return Required sample size per group.
#' @keywords internal
.calculate_sample_size_t <- function(d, alpha, power, alternative = "two.sided") {
  n <- 10
  tolerance <- 0.001
  max_iter <- 1000
  
  for (i in 1:max_iter) {
    df <- 2 * n - 2
    ncp <- d * sqrt(n)
    
    if (alternative == "two.sided") {
      t_crit <- qt(1 - alpha/2, df)
      current_power <- 1 - pt(t_crit, df, ncp) + pt(-t_crit, df, ncp)
    } else {
      t_crit <- qt(1 - alpha, df)
      current_power <- 1 - pt(t_crit, df, ncp)
    }
    
    if (current_power >= power) {
      return(n)
    }
    n <- n + 1
  }
  
  return(n)
}

#' @title Calculate Sample Size for ANOVA
#' @description Iteratively calculate required sample size for ANOVA.
#' @param f Effect size (Cohen's f).
#' @param alpha Significance level.
#' @param power Target power.
#' @return Required sample size per group.
#' @keywords internal
.calculate_sample_size_anova <- function(f, alpha, power) {
  n <- 10
  max_iter <- 1000
  
  for (i in 1:max_iter) {
    df1 <- 1
    df2 <- 2 * n - 2
    ncp <- f^2 * n
    f_crit <- qf(1 - alpha, df1, df2)
    current_power <- 1 - pf(f_crit, df1, df2, ncp)
    
    if (current_power >= power) {
      return(n)
    }
    n <- n + 1
  }
  
  return(n)
}

#' @title Interpret Effect Size
#' @description Provide interpretation of effect size based on Cohen's guidelines.
#' @param effect_size The effect size (Cohen's d).
#' @return A character string describing the effect size magnitude.
#' @export
interpret_effect_size <- function(effect_size) {
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

#' @title Power Analysis for Correlation Studies
#' @description Perform power analysis for correlation studies.
#' @import stats
#' @param effect_size The expected correlation coefficient (r).
#' @param n Sample size.
#' @param alpha Significance level. Default is 0.05.
#' @param power Desired statistical power. Default is 0.8.
#' @return A list containing power analysis results.
#' @examples
#' \dontrun{
#' # Calculate power for correlation study
#' cor_power <- power_analysis_cor(effect_size = 0.3, n = 100)
#' }
#' @export
power_analysis_cor <- function(effect_size = NULL, n = NULL, alpha = 0.05, power = 0.8) {
  
  if (is.null(effect_size) && is.null(n)) {
    stop("Either effect_size or n must be provided")
  }
  
  results <- list()
  
  if (!is.null(n) && !is.null(effect_size)) {
    # Calculate power for correlation
    # Using Fisher's z-transformation
    z <- 0.5 * log((1 + effect_size) / (1 - effect_size))
    se <- 1 / sqrt(n - 3)
    z_crit <- qnorm(1 - alpha/2)
    power_val <- pnorm(z * sqrt(n - 3) - z_crit)
    
    results$power <- power_val
    results$n <- n
    results$effect_size <- effect_size
  } else if (is.null(n) && !is.null(effect_size)) {
    # Calculate required sample size
    z <- 0.5 * log((1 + effect_size) / (1 - effect_size))
    z_crit <- qnorm(1 - alpha/2)
    z_power <- qnorm(power)
    required_n <- ((z_crit + z_power) / z)^2 + 3
    results$required_n <- ceiling(required_n)
    results$effect_size <- effect_size
    results$target_power <- power
  }
  
  # Add interpretation
  if (!is.null(effect_size)) {
    results$effect_size_interpretation <- interpret_correlation(effect_size)
  }
  
  return(results)
}

#' @title Interpret Correlation Coefficient
#' @description Provide interpretation of correlation coefficient magnitude.
#' @param r The correlation coefficient.
#' @return A character string describing the correlation strength.
#' @export
interpret_correlation <- function(r) {
  r <- abs(r)
  if (r < 0.1) {
    return("negligible (|r| < 0.1)")
  } else if (r < 0.3) {
    return("weak (0.1 <= |r| < 0.3)")
  } else if (r < 0.5) {
    return("moderate (0.3 <= |r| < 0.5)")
  } else if (r < 0.7) {
    return("strong (0.5 <= |r| < 0.7)")
  } else if (r < 0.9) {
    return("very strong (0.7 <= |r| < 0.9)")
  } else {
    return("almost perfect (|r| >= 0.9)")
  }
}

#' @title Calculate Minimum Detectable Effect Size
#' @description Calculate the minimum effect size detectable with given sample size and power.
#' @import stats
#' @param n1 Sample size in group 1.
#' @param n2 Sample size in group 2. If NULL, assumed equal to n1.
#' @param power Desired statistical power. Default is 0.8.
#' @param alpha Significance level. Default is 0.05.
#' @return The minimum detectable effect size (Cohen's d).
#' @examples
#' \dontrun{
#' min_effect <- minimum_detectable_effect(n1 = 50, n2 = 50)
#' }
#' @export
minimum_detectable_effect <- function(n1, n2 = NULL, power = 0.8, alpha = 0.05) {
  if (is.null(n2)) {
    n2 <- n1
  }
  
  df <- n1 + n2 - 2
  t_crit <- qt(1 - alpha/2, df)
  z_power <- qnorm(power)
  
  # Minimum detectable d
  d <- (t_crit + z_power) * sqrt((n1 + n2) / (n1 * n2))
  
  return(d)
}

#' @title Calculate Statistical Power for DEG Analysis
#' @description Calculate the statistical power to detect differentially expressed genes.
#' @import stats
#' @param n_tumor Number of tumor samples.
#' @param n_normal Number of normal samples.
#' @param effect_size Expected log fold change (as Cohen's d).
#' @param alpha Significance level. Default is 0.05.
#' @param multiple_testing Logical, whether to adjust for multiple testing. Default is TRUE.
#' @param num_tests Number of tests (genes) for multiple testing adjustment. Default is 20000.
#' @return A list containing power analysis results for DEG detection.
#' @export
power_analysis_deg <- function(n_tumor, n_normal, effect_size, alpha = 0.05, 
                               multiple_testing = TRUE, num_tests = 20000) {
  
  # Adjust alpha for multiple testing if requested
  adjusted_alpha <- if (multiple_testing) alpha / num_tests else alpha
  
  # Calculate power
  df <- n_tumor + n_normal - 2
  ncp <- effect_size * sqrt((n_tumor * n_normal) / (n_tumor + n_normal))
  t_crit <- qt(1 - adjusted_alpha/2, df)
  power_val <- 1 - pt(t_crit, df, ncp) + pt(-t_crit, df, ncp)
  
  # Calculate minimum detectable effect size
  mdes <- minimum_detectable_effect(n_tumor, n_normal, power = 0.8, alpha = adjusted_alpha)
  
  results <- list(
    power = power_val,
    n_tumor = n_tumor,
    n_normal = n_normal,
    effect_size = effect_size,
    alpha = alpha,
    adjusted_alpha = adjusted_alpha,
    minimum_detectable_effect = mdes,
    effect_size_interpretation = interpret_effect_size(effect_size),
    notes = list(
      multiple_testing_adjustment = ifelse(multiple_testing, 
                                          paste("Adjusted alpha for", num_tests, "tests"),
                                          "No multiple testing adjustment"),
      interpretation = paste("With the given sample sizes, you have", 
                             round(power_val * 100, 1), "% power to detect genes with effect size",
                             round(effect_size, 2), "or larger")
    )
  )
  
  return(results)
}