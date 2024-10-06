\name{GCAS_TIL}
\alias{GCAS_TIL}
\docType{data}
\title{
Immune cell infiltration score of all samples
}
\description{
Immune cell infiltration score of all samples calculated based on IOBR package.
}
\usage{data("GCAS_TIL")}
\format{
  A data frame with 19538 observations on the following 138 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{T_cells_MCPcounter}}{a numeric vector}
    \item{\code{Cytotoxic_lymphocytes_MCPcounter}}{a numeric vector}
    \item{\code{B_lineage_MCPcounter}}{a numeric vector}
    \item{\code{NK_cells_MCPcounter}}{a numeric vector}
    \item{\code{Monocytic_lineage_MCPcounter}}{a numeric vector}
    \item{\code{Myeloid_dendritic_cells_MCPcounter}}{a numeric vector}
    \item{\code{Neutrophils_MCPcounter}}{a numeric vector}
    \item{\code{Endothelial_cells_MCPcounter}}{a numeric vector}
    \item{\code{Fibroblasts_MCPcounter}}{a numeric vector}
    \item{\code{Bcells_EPIC}}{a numeric vector}
    \item{\code{CAFs_EPIC}}{a numeric vector}
    \item{\code{CD4_Tcells_EPIC}}{a numeric vector}
    \item{\code{CD8_Tcells_EPIC}}{a numeric vector}
    \item{\code{Endothelial_EPIC}}{a numeric vector}
    \item{\code{Macrophages_EPIC}}{a numeric vector}
    \item{\code{NKcells_EPIC}}{a numeric vector}
    \item{\code{otherCells_EPIC}}{a numeric vector}
    \item{\code{aDC_xCell}}{a numeric vector}
    \item{\code{Adipocytes_xCell}}{a numeric vector}
    \item{\code{Astrocytes_xCell}}{a numeric vector}
    \item{\samp{B-cells_xCell}}{a numeric vector}
    \item{\code{Basophils_xCell}}{a numeric vector}
    \item{\samp{CD4+_memory_T-cells_xCell}}{a numeric vector}
    \item{\samp{CD4+_naive_T-cells_xCell}}{a numeric vector}
    \item{\samp{CD4+_T-cells_xCell}}{a numeric vector}
    \item{\samp{CD4+_Tcm_xCell}}{a numeric vector}
    \item{\samp{CD4+_Tem_xCell}}{a numeric vector}
    \item{\samp{CD8+_naive_T-cells_xCell}}{a numeric vector}
    \item{\samp{CD8+_T-cells_xCell}}{a numeric vector}
    \item{\samp{CD8+_Tcm_xCell}}{a numeric vector}
    \item{\samp{CD8+_Tem_xCell}}{a numeric vector}
    \item{\code{cDC_xCell}}{a numeric vector}
    \item{\code{Chondrocytes_xCell}}{a numeric vector}
    \item{\samp{Class-switched_memory_B-cells_xCell}}{a numeric vector}
    \item{\code{CLP_xCell}}{a numeric vector}
    \item{\code{CMP_xCell}}{a numeric vector}
    \item{\code{DC_xCell}}{a numeric vector}
    \item{\code{Endothelial_cells_xCell}}{a numeric vector}
    \item{\code{Eosinophils_xCell}}{a numeric vector}
    \item{\code{Epithelial_cells_xCell}}{a numeric vector}
    \item{\code{Erythrocytes_xCell}}{a numeric vector}
    \item{\code{Fibroblasts_xCell}}{a numeric vector}
    \item{\code{GMP_xCell}}{a numeric vector}
    \item{\code{Hepatocytes_xCell}}{a numeric vector}
    \item{\code{HSC_xCell}}{a numeric vector}
    \item{\code{iDC_xCell}}{a numeric vector}
    \item{\code{Keratinocytes_xCell}}{a numeric vector}
    \item{\code{ly_Endothelial_cells_xCell}}{a numeric vector}
    \item{\code{Macrophages_xCell}}{a numeric vector}
    \item{\code{Macrophages_M1_xCell}}{a numeric vector}
    \item{\code{Macrophages_M2_xCell}}{a numeric vector}
    \item{\code{Mast_cells_xCell}}{a numeric vector}
    \item{\code{Megakaryocytes_xCell}}{a numeric vector}
    \item{\code{Melanocytes_xCell}}{a numeric vector}
    \item{\samp{Memory_B-cells_xCell}}{a numeric vector}
    \item{\code{MEP_xCell}}{a numeric vector}
    \item{\code{Mesangial_cells_xCell}}{a numeric vector}
    \item{\code{Monocytes_xCell}}{a numeric vector}
    \item{\code{MPP_xCell}}{a numeric vector}
    \item{\code{MSC_xCell}}{a numeric vector}
    \item{\code{mv_Endothelial_cells_xCell}}{a numeric vector}
    \item{\code{Myocytes_xCell}}{a numeric vector}
    \item{\samp{naive_B-cells_xCell}}{a numeric vector}
    \item{\code{Neurons_xCell}}{a numeric vector}
    \item{\code{Neutrophils_xCell}}{a numeric vector}
    \item{\code{NK_cells_xCell}}{a numeric vector}
    \item{\code{NKT_xCell}}{a numeric vector}
    \item{\code{Osteoblast_xCell}}{a numeric vector}
    \item{\code{pDC_xCell}}{a numeric vector}
    \item{\code{Pericytes_xCell}}{a numeric vector}
    \item{\code{Plasma_cells_xCell}}{a numeric vector}
    \item{\code{Platelets_xCell}}{a numeric vector}
    \item{\code{Preadipocytes_xCell}}{a numeric vector}
    \item{\samp{pro_B-cells_xCell}}{a numeric vector}
    \item{\code{Sebocytes_xCell}}{a numeric vector}
    \item{\code{Skeletal_muscle_xCell}}{a numeric vector}
    \item{\code{Smooth_muscle_xCell}}{a numeric vector}
    \item{\code{Tgd_cells_xCell}}{a numeric vector}
    \item{\code{Th1_cells_xCell}}{a numeric vector}
    \item{\code{Th2_cells_xCell}}{a numeric vector}
    \item{\code{Tregs_xCell}}{a numeric vector}
    \item{\code{ImmuneScore_xCell}}{a numeric vector}
    \item{\code{StromaScore_xCell}}{a numeric vector}
    \item{\code{MicroenvironmentScore_xCell}}{a numeric vector}
    \item{\code{B_cells_naive_CIBERSORT}}{a numeric vector}
    \item{\code{B_cells_memory_CIBERSORT}}{a numeric vector}
    \item{\code{Plasma_cells_CIBERSORT}}{a numeric vector}
    \item{\code{T_cells_CD8_CIBERSORT}}{a numeric vector}
    \item{\code{T_cells_CD4_naive_CIBERSORT}}{a numeric vector}
    \item{\code{T_cells_CD4_memory_resting_CIBERSORT}}{a numeric vector}
    \item{\code{T_cells_CD4_memory_activated_CIBERSORT}}{a numeric vector}
    \item{\code{T_cells_follicular_helper_CIBERSORT}}{a numeric vector}
    \item{\samp{T_cells_regulatory_(Tregs)_CIBERSORT}}{a numeric vector}
    \item{\code{T_cells_gamma_delta_CIBERSORT}}{a numeric vector}
    \item{\code{NK_cells_resting_CIBERSORT}}{a numeric vector}
    \item{\code{NK_cells_activated_CIBERSORT}}{a numeric vector}
    \item{\code{Monocytes_CIBERSORT}}{a numeric vector}
    \item{\code{Macrophages_M0_CIBERSORT}}{a numeric vector}
    \item{\code{Macrophages_M1_CIBERSORT}}{a numeric vector}
    \item{\code{Macrophages_M2_CIBERSORT}}{a numeric vector}
    \item{\code{Dendritic_cells_resting_CIBERSORT}}{a numeric vector}
    \item{\code{Dendritic_cells_activated_CIBERSORT}}{a numeric vector}
    \item{\code{Mast_cells_resting_CIBERSORT}}{a numeric vector}
    \item{\code{Mast_cells_activated_CIBERSORT}}{a numeric vector}
    \item{\code{Eosinophils_CIBERSORT}}{a numeric vector}
    \item{\code{Neutrophils_CIBERSORT}}{a numeric vector}
    \item{\samp{P-value_CIBERSORT}}{a numeric vector}
    \item{\code{Correlation_CIBERSORT}}{a numeric vector}
    \item{\code{RMSE_CIBERSORT}}{a numeric vector}
    \item{\code{MHC_IPS}}{a numeric vector}
    \item{\code{EC_IPS}}{a numeric vector}
    \item{\code{SC_IPS}}{a numeric vector}
    \item{\code{CP_IPS}}{a numeric vector}
    \item{\code{AZ_IPS}}{a numeric vector}
    \item{\code{IPS_IPS}}{a numeric vector}
    \item{\code{B_cells_quantiseq}}{a numeric vector}
    \item{\code{Macrophages_M1_quantiseq}}{a numeric vector}
    \item{\code{Macrophages_M2_quantiseq}}{a numeric vector}
    \item{\code{Monocytes_quantiseq}}{a numeric vector}
    \item{\code{Neutrophils_quantiseq}}{a numeric vector}
    \item{\code{NK_cells_quantiseq}}{a numeric vector}
    \item{\code{T_cells_CD4_quantiseq}}{a numeric vector}
    \item{\code{T_cells_CD8_quantiseq}}{a numeric vector}
    \item{\code{Tregs_quantiseq}}{a numeric vector}
    \item{\code{Dendritic_cells_quantiseq}}{a numeric vector}
    \item{\code{Other_quantiseq}}{a numeric vector}
    \item{\code{StromalScore_estimate}}{a numeric vector}
    \item{\code{ImmuneScore_estimate}}{a numeric vector}
    \item{\code{ESTIMATEScore_estimate}}{a numeric vector}
    \item{\code{TumorPurity_estimate}}{a numeric vector}
    \item{\code{B_cell_TIMER}}{a numeric vector}
    \item{\code{T_cell_CD4_TIMER}}{a numeric vector}
    \item{\code{T_cell_CD8_TIMER}}{a numeric vector}
    \item{\code{Neutrophil_TIMER}}{a numeric vector}
    \item{\code{Macrophage_TIMER}}{a numeric vector}
    \item{\code{DC_TIMER}}{a numeric vector}
    \item{\code{CD8_T_cells_MCPcounter}}{a numeric vector}
  }
}
\references{
Zeng, D., Ye, Z., Shen, R., Yu, G., Wu, J., Xiong, Y., Zhou, R., Qiu, W., Huang, N., Sun, L., Li, X., Bin, J., Liao, Y., Shi, M., & Liao, W. (2021). IOBR: Multi-Omics Immuno-Oncology Biological Research to Decode Tumor Microenvironment and Signatures. Frontiers in immunology, 12, 687975. https://doi.org/10.3389/fimmu.2021.687975
}
\examples{
data(GCAS_TIL)
str(GCAS_TIL)
}
