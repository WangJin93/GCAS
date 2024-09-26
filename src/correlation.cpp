#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Helper function to find the index of the gene of interest
int find_gene_index(CharacterVector rownames, std::string gene_of_interest) {
  for (int i = 0; i < rownames.size(); i++) {
    if (rownames[i] == gene_of_interest) {
      return i;
    }
  }
  return -1; // Gene not found
}

// [[Rcpp::export]]
DataFrame cal_correlation(NumericMatrix expression_matrix, CharacterVector rownames, std::string gene_of_interest, std::string method) {
  int n_genes = expression_matrix.nrow();
  int n_samples = expression_matrix.ncol();

  // Find the index of the gene_of_interest
  int gene_idx = find_gene_index(rownames, gene_of_interest);

  if (gene_idx == -1) {
    stop("Gene of interest not found in the expression matrix");
  }

  NumericVector target_gene_expression = expression_matrix(gene_idx, _);

  NumericVector correlations(n_genes);
  NumericVector p_values(n_genes);

  for (int i = 0; i < n_genes; i++) {
    if (i == gene_idx) {
      correlations[i] = NA_REAL;
      p_values[i] = NA_REAL;
      continue;
    }

    NumericVector gene_expression = expression_matrix(i, _);

    if (method == "pearson") {
      // Pearson correlation calculation
      double mean_target = mean(target_gene_expression);
      double mean_gene = mean(gene_expression);
      double numerator = sum((target_gene_expression - mean_target) * (gene_expression - mean_gene));
      double denominator = sqrt(sum(pow(target_gene_expression - mean_target, 2)) * sum(pow(gene_expression - mean_gene, 2)));
      double correlation = numerator / denominator;

      // Calculate p-value using R's pnorm function
      double t_stat = correlation * sqrt((n_samples - 2) / (1 - pow(correlation, 2)));
      // Calculate p-value
      double p_value = 2 * R::pt(-fabs(t_stat), n_samples - 2, true, false);

      correlations[i] = correlation;
      p_values[i] = p_value;
    } else {
      stop("Unsupported method");
    }
  }

  DataFrame result = DataFrame::create(
    _["gene"] = rownames,
    _["r"] = correlations,
    _["p"] = p_values
  );

  return result;
}
