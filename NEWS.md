# GCAS 0.2.1

## Bug fixes

* `.calculate_correlation()` — the shared engine behind `cor_gcas_genelist()`,
  `cor_gcas_TIL()` and `cor_gcas_drug()` — no longer removes an entire dataset
  from the results when one gene cannot be correlated in it (e.g. the gene is
  not measured on that dataset's platform or has constant expression).
  All result matrices (`r`, `p`, `p_adj`, `n`, `t`, `ci_lower`, `ci_upper`)
  now always retain every requested dataset × gene pair and have identical
  dimensions and dimnames. Pairs that cannot be computed are returned as `NA`
  (with a warning explaining the reason) instead of discarding the whole
  dataset together with all of its still-valid correlations.
* The sample-size matrix `n` now reports, per gene, the number of complete
  (non-missing) sample pairs actually used for that correlation, rather than
  the total number of samples of the dataset.
* Multiple-testing correction is now applied only to the computable p-values of
  each dataset (`NA` cells are left in place).
* `.calculate_correlation_ci()` no longer returns `NaN` for boundary cases
  (`|r| = 1`, `n <= 3`, missing `r`): perfect correlations return a degenerate
  interval and undefined cases return `NA`.
* `.determine_sample_type()` now warns when samples carry a missing `subtype`
  annotation and would silently be treated as "Tumor" (or dropped by a
  `tumor_subtype` filter), e.g. when the bundled `sample_subtype` data is out
  of sync with the expression database.
* A stale unit-test expectation in `.determine_sample_type` was corrected.

## Compatibility

* For every dataset–gene pair that was reported by v0.2.0, the correlation
  results are numerically identical (verified across a 15-dataset
  gene–gene correlation workflow); v0.2.1 only adds results that v0.2.0
  previously discarded.

## Tests

* Added regression tests ensuring that (i) datasets with partially missing or
  non-variable genes stay in the result matrices with per-pair `NA`s and
  (ii) the confidence-interval helper handles boundary cases without `NaN`.
