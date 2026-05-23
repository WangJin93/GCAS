## GCAS: An Integrated R Package and Shiny App for Comprehensive Cancer Data Analysis

**GCAS (GEO Cancer Analysis Suite)** is an R package with an integrated Shiny application designed to provide a unified platform for cancer transcriptomic analysis.  
It addresses common challenges in cancer research, such as data complexity, dispersed tools, and lack of reproducible, integrated workflows.

GCAS is primarily **GEO‑centered**, but also supports user‑supplied expression matrices and associated clinical information, making it suitable for both public and in‑house datasets.

---

### Data Sources

- Public gene expression datasets (primarily **GEO**)
- User‑provided expression matrices (microarray or RNA‑seq)
- Optional clinical/phenotypic annotations supplied by the user

---

### Core Functional Modules

GCAS currently includes four main analysis modules:

1. **Differential Gene Expression Analysis**
   - Tumor vs. normal or user‑defined group comparisons
   - Supports multi‑dataset integration (e.g., batch‑corrected merged analysis)

2. **Correlation Studies**
   - Gene–gene co‑expression analysis (e.g., GAPDH with m6A regulators such as IGF2BP3)
   - Provides correlation coefficients, FDR‑adjusted p‑values, and confidence intervals

3. **Pan‑Cancer Analysis**
   - Evaluates gene expression across multiple cancer types and datasets
   - Facilitates systematic assessment of candidate biomarkers in different tumors

4. **Immune Infiltration and Drug Sensitivity Analysis**
   - Assesses associations between gene expression and immune cell infiltration scores  
   - Integrates drug response prediction (e.g., using OncoPredict) to explore links between gene expression (such as GAPDH) and sensitivity to anticancer drugs (e.g., EGFR‑targeting agents like Erlotinib)

In a representative application, GCAS revealed that **GAPDH** is upregulated in multiple lung and breast cancer datasets and is positively correlated with the m6A reader **IGF2BP3**. Follow‑up in vitro experiments suggested that IGF2BP3 regulates GAPDH mRNA stability. GCAS also indicated a negative association between GAPDH expression and CD4 T cell infiltration, and a negative correlation between GAPDH expression and sensitivity to EGFR‑targeting drugs.

---

### GCAS Overview

> *Figure X. Overview of the GCAS workflow and main functional modules.*  
> ![GCAS overview](https://www.jingege.wang/wp-content/uploads/2025/03/e-post-phd-study-gcas-manuscript-abstract.png)

---

### Availability and Documentation

- **Source code & R package:**  
  <https://github.com/WangJin93/GCAS>

- **Bug reports / Issues:**  
  <https://github.com/WangJin93/GCAS/issues>

- **User guide / Documentation:**  
  <https://wangjin93.github.io/gcas.html>

The online documentation provides installation instructions, module descriptions, and example workflows using built‑in and public datasets.

---

### Installation

You can install GCAS from GitHub using the `remotes` package:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install GCAS from GitHub
remotes::install_github("WangJin93/GCAS")
Then load GCAS and (optionally) launch the Shiny app:

library(GCAS)

# Launch the GCAS Shiny interface
GCAS::run_GCAS()
Users can either:

interact with GCAS via the Shiny graphical interface, or
script analyses directly in R using GCAS functions for fully reproducible pipelines.
