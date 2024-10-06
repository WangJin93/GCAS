![Graphic abstract](abstract.tif)

1. <a name="_toc179051521"></a>**Introduction**

**Title:** Visualization of Functional Enrichment Result

**Version:** 1.0.0

**Author:** Jin wang (Jin.wang93@outlook.com)

**Maintainer**: Jin wang (Jin.wang93@outlook.com)

**Description:** The GEO Cancer Analysis Suite (GCAS) is a versatile R package designed for analyzing and visualizing gene expression data in cancer research. GCAS allows for the comparison of gene expression between normal and tumor samples, correlation analysis, immune infiltration analysis, differential expression analysis, co-expression analysis, and enrichment analysis. It includes a Shiny app for interactive visualization and can also be used directly within the R environment for advanced scripting. GCAS is ideal for researchers, clinicians, and bioinformaticians seeking to explore cancer genomics data efficiently and effectively.

**Depends:** R (>= 3.5.0)

**Imports:** RobustRankAggreg, VennDiagram, digest, dplyr, ggpubr, ggrepel, httr, jsonlite, meta, psych, shiny, stringr, sva, tibble, RColorBrewer, clusterProfiler, dplyr, ggrepel, grid, ggplot2

**Encoding:** UTF-8

**URL:** https://github.com/WangJin93/GCAS

**Bug Reports:** https://github.com/WangJin93/GCAS/issues

**License:** MIT License

2. <a name="_toc179051522"></a>**Installation**
```R
remotes::install_github("WangJin93/GCAS")
```
3. <a name="_toc179051523"></a>**Function Reference**

3.1. <a name="_toc179051524"></a>**datasets_summary**

**Description**

This function summarizes the sample counts and pairing status for datasets based on a specified tumor subtype.

**Usage**
```R

datasets_summary(tumor_subtype = "NSCLC")
```
**Arguments**

|tumor_subtype|A character string specifying the tumor subtype to filter the datasets. Default is "NSCLC".|
| :- | :- |

**Value**

A data frame summarizing the number of Normal, Tumor, and Premalignant samples, as well as the pairing status for each dataset.

**Examples**
```R

dt_sum <- datasets_summary("NSCLC")
```
3.2. <a name="_toc179051525"></a>**get_expr_data**

**Description**

Retrieve expression data for specified genes from given datasets.

**Usage**
```R

get_expr_data(datasets, genes)
```
**Arguments**

|datasets|A character vector of dataset identifiers.|
| :- | :- |
|genes|A character vector of gene identifiers.|

**Value**

A dataframe containing expression data for the specified genes from the given datasets.

**Examples**
```R

results <- get_expr_data(datasets = "GSE74706", genes = c("GAPDH","TNS1"))
results <- get_expr_data(datasets = c("GSE62113","GSE74706"), genes = "GAPDH")
results <- get_expr_data(datasets = c("GSE62113","GSE74706"), genes = c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA"))
```
3.3. <a name="_toc179051526"></a>**viz_TvsN**

**Description**

Visualizing the different expression of mRNA expression data between Tumor and Normal tissues in GEO database.

**Usage**
```R

viz_TvsN(
     df,
     df_type = c("single", "multi_gene", "multi_set"),
     tumor_subtype = NULL,
     Show.P.Value = TRUE,
     Show.P.label = TRUE,
     Method = "t.test",
     Values = c("#00AFBB", "#FC4E07"),
     Show.n = TRUE,
     Show.n.location = "default"
)
```
**Arguments**

|df|Gene expression data obtained from get_expr_data().|
| :- | :- |
|df_type|The type of gene expression data, one Value of "single","multi_gene", and "multi_set".|
|Show.P.Value|Whether to display the results of differential analysis, default TRUE.|
|Show.P.label|Whether to display significance markers for differential analysis, default TRUE.|
|Method|Methods of differential analysis, "t.test" or "limma", default "t.test".|
|Values|Color palette for normal and tumor groups. Default c("#00AFBB", "#FC4E07").|
|Show.n|Display sample size.|
|Show.n.location|Y-axis position displayed for sample size.|

**Examples**
```R

df_single <- get_expr_data(datasets = "GSE27262",genes = c("TP53"))
df_multi_gene <- get_expr_data(datasets = "GSE27262",genes = c("TP53","TNS1"))
df_multi_set <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
viz_TvsN(df_single,df_type = "single")
viz_TvsN(df_multi_gene,df_type = "multi_gene",tumor_subtype ="LC")
viz_TvsN(df_multi_set,df_type = "multi_set")
```
3.4. <a name="_toc179051527"></a>**data_summary**

**Description**

Compute summary statistics (mean, standard deviation, etc.) and perform hypothesis testing (t-test or Wilcoxon test) for gene expression data across different datasets.

**Usage**
```R

data_summary(df, tumor_subtype = NULL, method = "t.test")
```
**Arguments**

|df|A dataframe containing gene expression data, with columns: 'dataset', 'subtype', and the gene expression Values.|
| :- | :- |
|tumor_subtype|A character string specifying the tumor subtype to be analyzed. Default is NULL, which means all tumor subtypes will be included.|
|method|A character string specifying the method for hypothesis testing. Options are "t.test" for t-test and "wilcox" for Wilcoxon test. Default is "t.test".|

**Value**

A dataframe with summary statistics and p-Values for each dataset.

**Examples**
```R

df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
```
3.5. <a name="_toc179051528"></a>**plot_meta_forest**

**Description**

Plotting volcano plot for DEGs between tumor and normal samples in CPTAC datasets.

This function performs a meta-analysis on multiple datasets and generates a forest plot. It also tests for publication bias.

**Usage**
```R
plot_meta_forest(results, method = "wilcox", k.min = 10)
```
**Arguments**

|results|Data frame. The results data frame containing columns for dataset, n_Tumor, mean_Tumor, sd_Tumor, n_Normal, mean_Normal, and sd_Normal.|
| :- | :- |
|method|Character. The statistical method to use (default is "wilcox").|
|k.min|Integer. Minimum number of studies for bias test (default is 7).|
|cohort|Data cohort, for example, "LUAD_APOLLO", "LUAD_CPTAC".|
|data_input|Expression data obtained from get_expr_data() function.|

**Value**

A forest plot object.

**Examples**
```R

df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
plot_meta_forest(results)
```
3.6. <a name="_toc179051529"></a>**plot_logFC_heatmap**

**Description**

This function generates a heatmap of log fold changes (logFC) for a gene across different datasets. It includes significance annotations based on p-Values.

**Usage**
```R

plot_logFC_heatmap(results, direction = "horizontal")
```
**Arguments**

|results|Data frame. The results data frame containing columns for dataset, gene, logFC, and p.Value.|
| :- | :- |
|direction|Ploting direction, horizontal or vertical.|

**Value**

A ggplot object representing the heatmap.

**Examples**
```R

df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
heatmap <- plot_logFC_heatmap(results)
print(heatmap)
```
3.7. <a name="_toc179051530"></a>**plot_logFC_scatter**

**Description**

This function generates a scatter of log fold changes (logFC) for a gene across different datasets. It includes significance annotations based on p-Values.

**Usage**
```R

plot_logFC_scatter(
     results,
     p.cut = 0.05,
     logFC.cut = 1,
     colors = c("blue", "grey20", "red")
)
```
**Arguments**

|results|Data frame. The results data frame containing columns for dataset, gene, logFC, and p.Value.|
| :- | :- |
|p.cut|Numeric. The cutoff for adjusted p-Value to determine significance. Default is 0.05.|
|logFC.cut|Numeric. The cutoff for log fold change to determine significance. Default is 1.|
|colors|A vector of color panel, default c("blue", "grey20", "red").|

**Value**

A ggplot object representing the scatter.

**Examples**
```R

df <- get_expr_data(datasets = c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113"), genes = "GAPDH")
results <- data_summary(df, tumor_subtype = "LUAD")
scatter <- plot_logFC_scatter(results, logFC.cut = 0.5, colors = c("blue","grey20", "red"))
print(scatter)

```
3.8. <a name="_toc179051531"></a>**cor_cancer_genelist**

**Description**

Perform correlation analysis of the mRNA/protein expression data in CPTAC database.

**Usage**
```R

cor_cancer_genelist(
     dataset = "GSE62113",
     id1 = "STAT3",
     id2 = c("TNS1", "TP53"),
     tumor_subtype = NULL,
     sample_type = c("Tumor", "Normal"),
     cor_method = "pearson"
)
```
**Arguments**

|dataset|Dataset name. Use 'dataset$Abbre' to get all datasets.|
| :- | :- |
|id1|Gene symbol, you can input one gene symbol.|
|id2|Gene symbols, you can input one or multiple symbols.|
|tumor_subtype|Tumor subtype used for correlation analysis, default is NULL.|
|sample_type|Sample type used for correlation analysis, default all types: c("Tumor", "Normal").|
|cor_method|Method for correlation analysis, default "pearson".|

**Value**

A list containing the correlation results and the merged data.

**Examples**
```R

results <- cor_cancer_genelist(dataset = "GSE62113",
     id1 = "STAT3",tumor_subtype = "LC",
     id2 = c("TNS1", "TP53"),
     sample_type = c("Tumor", "Normal"),
     cor_method = "pearson")

```
3.9. <a name="_toc179051532"></a>**cor_gcas_drug**

**Description**

Calculate the correlation between target gene expression and anti-tumor drug sensitivity in multiple datasets.

**Usage**
```R

cor_gcas_drug(df, cor_method = "pearson", Target.pathway = c("Cell cycle"))
```
**Arguments**

|df|The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.|
| :- | :- |
|cor_method|Method for correlation analysis, default "pearson".|
|Target.pathway|The signaling pathways of anti-tumor drug targets, default "Cell cycle". Use "drug_info" to get the detailed information of these drugs.|

**Examples**
```R

dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210",
     "GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072",
     "GSE74706","GSE18842","GSE62113")
df <- get_expr_data(genes = "TNS1", datasets = dataset)
result <- cor_gcas_drug(df, Target.pathway = c("Cell cycle"))

```
3.10. <a name="_toc179051533"></a>**cor_gcas_genelist**

**Description**

Perform correlation analysis of the expression data in multiple datasets.

**Usage**
```R

cor_gcas_genelist(
     df,
     geneset_data,
     tumor_subtype = NULL,
     sample_type = c("Tumor", "Normal"),
     cor_method = "pearson"
)
```
**Arguments**

|df|The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.|
| :- | :- |
|geneset_data|The expression data of a genelist in multiple datasets, obtained by the get_expr_data() function.|
|tumor_subtype|Tumor subtype used for correlation analysis, default is NULL.|
|sample_type|Sample type used for correlation analysis, default all types: c("Tumor", "Normal").|
|cor_method|Method for correlation analysis, default "pearson".|

**Examples**
```R

genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113")
df <- get_expr_data(genes = "TNS1",datasets = dataset)
geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
result <- cor_gcas_genelist(df, geneset_data, sample_type = c("Tumor"))

```
3.11. <a name="_toc179051534"></a>**cor_gcas_TIL**

**Description**

Calculate the correlation between target gene expression and immune cells infiltration in multiple datasets.

**Usage**
```R

cor_gcas_TIL(df, cor_method = "spearman", TIL_type = "TIMER")
```
**Arguments**

|df|The expression data of the target gene in multiple datasets, obtained by the get_expr_data() function.|
| :- | :- |
|cor_method|Method for correlation analysis, default "spearman".|
|TIL_type|Algorithm for calculating immune cell infiltration, default "TIMER".|

**Examples**
```R

dataset <- c("GSE27262", "GSE7670", "GSE19188", "GSE19804", "GSE30219",
     "GSE31210", "GSE32665", "GSE32863", "GSE43458", "GSE46539",
     "GSE75037", "GSE10072", "GSE74706", "GSE18842", "GSE62113")
df <- get_expr_data(genes = "TNS1", datasets = dataset)
result <- cor_gcas_TIL(df, cor_method = "spearman", TIL_type = "TIMER")
```
3.12. <a name="_toc179051535"></a>**viz_cor_heatmap**

**Description**

Presenting correlation analysis results using heat maps based on ggplot2.

**Usage**
```R

viz_cor_heatmap(r, p)
```
**Arguments**

|r|The correlation coefficient matrix r of the correlation analysis results obtained from the functions cor_pancancer_genelist(), cor_pancancer_TIL(), and cor_pancancer_drug().|
| :- | :- |
|p|The P-Value matrix p of the correlation analysis results obtained from the functions cor_pancancer_genelist(), cor_pancancer_TIL(), and cor_pancancer_drug().|

**Examples**
```R

genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
dataset <-  c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113")
df <- get_expr_data(genes = "TNS1", datasets = dataset)
geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
result <- cor_gcas_genelist(df, geneset_data, sample_type = c("Tumor"))
viz_cor_heatmap(result$r, result$p)

```
3.13. <a name="_toc179051536"></a>**viz_cor_volcano**

**Description**

Plotting volcano plot for the correlation analysis of a specific gene from the result of cor_gcas_genelist().

**Usage**
```R

viz_cor_volcano(
     cor_result,
     item,
     p.cut = 0.05,
     r.cut = 0.3,
     colors = c("blue", "grey20", "red")

)
```
**Arguments**

|cor_result|DataFrame. The results from cor_gcas_genelist analysis.|
| :- | :- |
|item|Character. Gene symbol or drug name or immune cell.|
|p.cut|Numeric. The cutoff for adjusted p-Value to determine significance. Default is 0.05.|
|r.cut|Numeric. The cutoff for log fold change to determine significance. Default is 0.3.|
|colors|A vector of color panel, default c("blue", "grey20", "red").|

**Value**

A ggplot2 object representing the volcano plot.

**Examples**
```R

genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
dataset <- c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113")
df <- get_expr_data(genes = "TNS1",datasets = dataset)
geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
cor_result <- cor_gcas_genelist(df, geneset_data, sample_type = c("Tumor"))
viz_cor_volcano(cor_result, "LILRB4", p.cut = 0.5, r.cut = 0.1,colors = c("blue" ,"grey20", "red"))

```
3.14. <a name="_toc179051537"></a>**viz_corplot**

**Description**

Scatter plot with sample size (n), correlation coefficient (r) and p Value (p.Value).

**Usage**
```R

viz_corplot(
     data,
     a,
     b,
     method = "pearson",
     x_lab = " relative expression",
     y_lab = " relative expression"

)
```
**Arguments**

|data|A gene expression dataset with at least two genes included, rows represent samples, and columns represent gene expression in the matrix.|
| :- | :- |
|a|Gene A|
|b|Gene B|
|method|Method for correlation analysis, "pearson" or "spearman".|
|x_lab|X-axis label.|
|y_lab|Y-axis label.|

3.15. <a name="_toc179051538"></a>**get_OSF_data**

**Description**

Retrieve GEO expression datasets and sample information from the OSF repository.

**Usage**
```R

get_OSF_data(table = "GSE19188", action = "geo_data")
```
**Arguments**

|table|A character string specifying the GEO dataset identifier (e.g., "GSE19188").|
| :- | :- |
|action|A character string specifying the action, either "geo_data" to retrieve the expression data or "sample_info" to retrieve the sample information.|

**Value**

A data frame containing the requested data.

**Examples**
```R

df <- get_OSF_data(table = "GSE74706", action = "sample_info")
df2 <- get_OSF_data(table = "GSE74706", action = "geo_data")

```
3.16. <a name="_toc179051539"></a>**DEGs_analysis**

**Description**

Perform differential expression gene analysis on a given dataset.

**Usage**
```R

DEGs_analysis(df, tumor_subtype = NULL, ...)
```
**Arguments**

|df|A dataframe containing gene expression data with sample IDs as columns.|
| :- | :- |
|tumor_subtype|A character vector specifying the tumor subtypes to be analyzed. Default is NULL, which means all tumor subtypes will be included.|
|...|Additional Arguments passed to 'lmFit', 'contrasts.fit', and 'eBayes'.|

**Value**

A dataframe with DEG analysis results, including log fold changes and p-Values.

**Examples**
```R

df <- get_OSF_data(table = "GSE74706", action = "geo_data")
results <- DEGs_analysis(df, tumor_subtype = c("NSCLC"))

```
3.17. <a name="_toc179051540"></a>**plot_volcano**

**Description**

Plotting volcano plot for DEGs between tumor and normal samples in CPTAC datasets.

**Usage**
```R

plot_volcano(
     results,
     p.cut = 0.05,
     logFC.cut = 1,
     show.top = FALSE,
     show.labels = NULL,
     colors = c("blue", "grey20", "red")

)
```
**Arguments**

|results|DataFrame. The results from DEGs analysis containing columns 'adj.P.Val', 'P.Value', 'logFC', and 'gene'.|
| :- | :- |
|p.cut|Numeric. The cutoff for adjusted p-Value to determine significance. Default is 0.05.|
|logFC.cut|Numeric. The cutoff for log fold change to determine significance. Default is 1.|
|show.top|Logical. If TRUE, labels the top 5 up- and downregulated genes. Default is FALSE.|
|show.labels|Character vector. Specific gene labels to show. Default is NULL.|
|colors|A vector of color panel, default c("blue", "grey20", "red").|

**Value**

A ggplot2 object representing the volcano plot.

**Examples**
```R

df <- get_OSF_data(table = "GSE74706", action = "geo_data")
results <- DEGs_analysis(df)
plot_volcano(results)

```
3.18. <a name="_toc179051541"></a>**coexpression_analysis**

**Description**

This function calculates the correlation between a given gene and all other genes in the provided expression matrix. It also provides the corresponding p-Values.

**Usage**
```R

coexpression_analysis(expression_matrix, gene, method = "pearson")
```
**Arguments**

|expression_matrix|A numeric matrix where rows represent genes and columns represent samples.|
| :- | :- |
|gene|A character string representing the gene for which the correlations will be calculated.|
|method|A character string specifying the correlation method to be used. Default is "pearson".|

**Value**

A data frame containing gene names, correlation coefficients, and p-Values.

**Examples**
```R

expression_matrix <- get_OSF_data(table = "GSE74706", action = "geo_data")
     results <- coexpression_analysis(expression_matrix, "RPN1")
     print(results)

```
3.19. <a name="_toc179051542"></a>**GSEA_analysis**

**Description**

This function performs Gene Set Enrichment Analysis (GSEA) based on either correlation results or limma differential analysis results.

**Usage**
```R

GSEA_analysis(data, gmt_file, pValue_cutoff = 0.05, data_type = "correlation")
```
**Arguments**

|data|A data frame containing gene names and corresponding Values. For correlation results, the columns should be named 'gene' and 'r'. For limma results, the columns should be named 'gene' and 'logFC'.|
| :- | :- |
|gmt_file|Path to the GMT file containing gene sets, or directly pass GO/KEGG/Reactome datasets.|
|pValue_cutoff|Numeric, the p-Value threshold for significance. Default is 0.05.|
|data_type|Character, type of the input data. Either "correlation" for correlation analysis results or "limma" for limma differential analysis results.|

**Value**

A GSEA analysis result object.

**Examples**
```R

df <- get_OSF_data(table = "GSE74706", action = "geo_data")
results <- DEGs_analysis(df,tumor_subtype =c("NSCLC"))
gsea_result <- GSEA_analysis(results,  gmt_file = BP_GMT_7.5.1, data_type = "limma")
results <- coexpression_analysis(df,"RPN1")
gsea_result <- GSEA_analysis(results,  gmt_file = BP_GMT_7.5.1)

```
3.20. <a name="_toc179051543"></a>**get_DEGs_list**

**Description**

Extract significantly upregulated and downregulated genes from multiple DEG analysis results.

**Usage**
```R

get_DEGs_list(DEGs_lists, logFC_cut = 1, p_cut = 0.05)
```
**Arguments**

|DEGs_lists|A list of dataframes containing DEG analysis results.|
| :- | :- |
|logFC_cut|A numeric Value specifying the log fold change cutoff for significant DEGs. Default is 1.|
|p_cut|A numeric Value specifying the p-Value cutoff for significant DEGs. Default is 0.05.|

**Value**

A list containing two lists: one for upregulated genes and one for downregulated genes across the provided datasets.

**Examples**
```R

df1 <- get_OSF_data(table = "GSE31210", action = "geo_data")
results1 <- DEGs_analysis(df1)
df2 <- get_OSF_data(table = "GSE19188", action = "geo_data")
results2 <- DEGs_analysis(df2)
DEGs_lists <- list("GSE31210" = results1, "GSE19188" = results2)
results <- get_DEGs_list(DEGs_lists)

```
3.21. <a name="_toc179051544"></a>**plot_venn**

**Description**

This function plots a Venn diagram for lists of differentially expressed genes (DEGs) across multiple datasets.

**Usage**
```R

plot_venn(results, fill_colors = NULL, palette = "Set1", lty = 2, ...)
```
**Arguments**

|results|List of character vectors. Each vector contains DEGs for a specific dataset.|
| :- | :- |
|fill_colors|Character vector. Colors to fill the Venn diagram circles. Default is NULL, which uses a palette.|
|palette|Character. Name of the RColorBrewer palette to use if fill_colors is not specified. Default is "Set1".|
|lty|Numeric. Line type for the circles in the Venn diagram. Default is 2 (dashed line).|
|...|Additional Arguments passed to venn.diagram function.|

**Value**

A list of intersected DEGs.

**Examples**
```R

df1 <- get_OSF_data(table = "GSE31210", action = "geo_data")
results1 <- DEGs_analysis(df1)
df2 <- get_OSF_data(table = "GSE19188", action = "geo_data")
results2 <- DEGs_analysis(df2)
DEGs_lists <- list("GSE31210" = results1, "GSE19188" = results2)
results <- get_DEGs_list(DEGs_lists)
plot_venn(results$DEG_up, palette = "Set1")
plot_venn(results$DEG_up, fill_colors = c("red", "green", "blue"), alpha = 0.5, cex = 1.5)

```
3.22. <a name="_toc179051545"></a>**RRA_analysis**

**Description**

This function performs RRA analysis on differentially expressed genes (DEGs) lists obtained from various studies. It ranks genes based on their differential expression and aggregates the ranks to identify consistently regulated genes across studies.

**Usage**
```R

RRA_analysis(
     DEGs_lists,
     top.num = 0,
     rra.p = 0.05,
     logFC_cut = 1,
     p_cut = 0.05
)
```
**Arguments**

|DEGs_lists|A list of DEGs data frames. Each data frame should contain at least a 'gene' column and a 'logFC' column.|
| :- | :- |
|top.num|Numeric, the number of top genes to select based on their ranks. Default is 0, which selects all genes passing the thresholds.|
|rra.p|Numeric, the p-Value threshold for RRA. Default is 0.05.|
|logFC_cut|Numeric, the log fold change threshold for filtering genes. Default is 1.|
|p_cut|Numeric, the p-Value threshold for filtering genes. Default is 0.05.|

**Value**

A list containing the number of up- and down-regulated genes and a data matrix of aggregated log fold changes.

**Examples**
```R

df1 <- get_OSF_data(table = "GSE31210", action = "geo_data")
results1 <- DEGs_analysis(df1)
df2 <- get_OSF_data(table = "GSE19188", action = "geo_data")
results2 <- DEGs_analysis(df2)
DEGs_lists <- list("GSE31210" = results1, "GSE19188" = results2)
RRA_results <- RRA_analysis(DEGs_lists)
ComplexHeatmap::pheatmap(RRA_results$RRA_results)

```
3.23. <a name="_toc179051546"></a>**combat_datasets**

**Description**

This function performs batch correction on multiple datasets using the ComBat function from the sva package.

**Usage**
```R
combat_datasets(tables, tumor_subtype = NULL)
```
**Arguments**

|tables|A character vector of table names to be processed.|
| :- | :- |
|tumor_subtype|A character string specifying the tumor subtype to filter the datasets. If NULL, all subtypes are included.|

**Value**

A list containing the combined and batch-corrected data matrix and the sample information.

**Examples**
```R

tables <- c("GSE31210", "GSE74706")
     result <- combat_datasets(tables, tumor_subtype = "LC")
     combined_data <- result$combined_data
     sample_info <- result$sample_info
```
3.24. <a name="_toc179051547"></a>**merge_clinic_data**

**Description**

Get sample_info data and merge it with expression data.

**Usage**
```R
merge_clinic_data(table = "GSE19188", data_input)
```
**Arguments**

|table|Character. The name of the dataset table to retrieve sample information from. Default is "GSE19188".|
| :- | :- |
|data_input|Data frame. Expression data obtained from get_expr_data() function.|

**Value**

Data frame. Merged data containing both expression data and sample information.

**Examples**
```R

data_input <- get_expr_data("GSE19188", "TP53")
results <- merge_clinic_data("GSE19188",data_input)

```
3.25. <a name="_toc179051548"></a>**extract_subset**

**Description**

This function searches for specified names within a nested list structure and extracts the names of found subsets.

**Usage**
```R

extract_subset(lst, names_to_find)
```
**Arguments**

|lst|A list which may contain nested lists.|
| :- | :- |
|names_to_find|A character vector of names to search for within the list.|

**Value**

A character vector of unique names from the found subsets.

**Examples**
```R

nested_list <- list(
     a = list(
     b = 1,
     c = list(d = 2)
     ),
     e = 3
)
names_to_search <- c("b", "d", "e")
result <- extract_subset(nested_list, names_to_search)
print(result)
```
4. <a name="_toc179051549"></a>**Datasets**
   
4.1. <a name="_toc179051550"></a><a name="_toc179051557"></a>**dataset_info**

**Description**

Summary of the general informations of the GEO datasets in this tool.

**Usage**
```R

data("dataset_info")

```
**Format**

A data frame with 288 observations on the following 6 variables.

4.2. **abbr_full**

**Description**

The full name of tumor abbreviation

**Usage**
```R

data("abbr_full")

```
**Format**

A data frame with 132 observations on the following 3 variables.

4.3. <a name="_toc179051551"></a><a name="_toc179051552"></a>**Subtype**

**Description**

Subtype list of al cancer types

**Usage**
```R

data("subtype")

```
**Format**

The format is: List of 21

4.4. **TIL_map**

**Description**

Mapping of immune cell types and algorithms

**Usage**
```R

data("TIL_map")

```
**Format**

A data frame with 137 observations on the following 2 variables.

4.5. <a name="_toc179051553"></a>**sample_subtype**

**Description**

Subtype information of all samples

**Usage**
```R

data("sample_subtype")

```
**Format**

A data frame with 28032 observations on the following 5 variables.

4.6. <a name="_toc179051554"></a>**GCAS_TIL**

**Description**

Immune cell infiltration score of all samples calculated based on IOBR package.

**Usage**
```R

data("GCAS_TIL")

```
**Format**

A data frame with 19538 observations on the following 138 variables.

4.7. <a name="_toc179051555"></a>**GCAS_drug**

**Description**

Anti-tumor drug sensitivity of all samples calculated based on oncoPredict package and GDSC database.

**Usage**
```R

data("GCAS_drug")

```
**Format**

A data frame with 19816 observations on the following 199 variables.

4.8. <a name="_toc179051556"></a>**drug_info**

**Description**

Anti-tumor drugs informations obtained from GDSC2.0 datasets.

**Usage**
```R

data("drug_info")

```
**Format**

A data frame with 198 observations on the following 6 variables.

