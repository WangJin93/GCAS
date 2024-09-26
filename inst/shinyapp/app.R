#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(GCAS)
# 安装并加载CRAN包
if (!require("shiny")) install.packages("shiny", update = FALSE, ask = FALSE); library(shiny)
if (!require("waiter")) install.packages("waiter", update = FALSE, ask = FALSE); library(waiter)
if (!require("ggplot2")) install.packages("ggplot2", update = FALSE, ask = FALSE); library(ggplot2)
if (!require("ggpubr")) install.packages("ggpubr", update = FALSE, ask = FALSE); library(ggpubr)
if (!require("dplyr")) install.packages("dplyr", update = FALSE, ask = FALSE); library(dplyr)
if (!require("DT")) install.packages("DT", update = FALSE, ask = FALSE); library(DT)
if (!require("shinyalert")) install.packages("shinyalert", update = FALSE, ask = FALSE); library(shinyalert)
if (!require("stringr")) install.packages("stringr", update = FALSE, ask = FALSE); library(stringr)
if (!require("shinyWidgets")) install.packages("shinyWidgets", update = FALSE, ask = FALSE); library(shinyWidgets)
if (!require("psych")) install.packages("psych", update = FALSE, ask = FALSE); library(psych)
if (!require("ggthemes")) install.packages("ggthemes", update = FALSE, ask = FALSE); library(ggthemes)
if (!require("shinyjs")) install.packages("shinyjs", update = FALSE, ask = FALSE); library(shinyjs)
if (!require("bs4Dash")) install.packages("bs4Dash", update = FALSE, ask = FALSE); library(bs4Dash)
if (!require("ggrepel")) install.packages("ggrepel", update = FALSE, ask = FALSE); library(ggrepel)
if (!require("tibble")) install.packages("tibble", update = FALSE, ask = FALSE); library(tibble)
if (!require("tidyr")) install.packages("tidyr", update = FALSE, ask = FALSE); library(tidyr)
if (!require("limma")) install.packages("limma", update = FALSE, ask = FALSE); library(tidyr)
if(!require("VennDiagram")) install.packages("VennDiagram",update = F,ask = F)
library(shinyTree)
library(htmlwidgets)
# 安装并加载Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", update = FALSE, ask = FALSE)

if (!require("aplot")) BiocManager::install("aplot", update = FALSE, ask = FALSE); library(aplot)
if (!require("ggtree")) BiocManager::install("ggtree", update = FALSE, ask = FALSE); library(ggtree)
if (!require("enrichplot")) BiocManager::install("enrichplot", update = FALSE, ask = FALSE); library(enrichplot)

source("apps/visualize.R")
source("apps/Dashboard.R")
source("apps/modules_GEO_datasets.R")
source("apps/modules_Cancer_expression.R")
source("apps/modules_Cancer_dataset_info.R")
source("apps/modules_Cancer_correlation.R")
source("apps/modules-gcas-expression.R")
source("apps/modules_Cancer_multiple.R")
source("apps/modules-gcas-til.R")
source("apps/modules-gcas-drug.R")
source("apps/modules_DEGs.R")
source("apps/modules_GSEA.R")
source("apps/modules_co-expression.R")
source("apps/modules-gcas-corr.R")
source("apps/modules_RRA.R")
source("apps/modules_venn.R")
source("apps/modules_combat.R")
source("apps/mod_feedback.R")
# Define UI for application that draws a histogram
ui <- bs4Dash::dashboardPage(
  bs4Dash::dashboardHeader(
    title = tags$a(href='https://github.com/WangJin93/GCAS',
                   tags$img(src='logo.png',width="100%"),
                   style = "position: relative; top: 5px;",align="center")
  ),

  bs4Dash::dashboardSidebar(
    bs4Dash::sidebarMenu(
      bs4Dash::menuItem("Introduction",icon = icon("info-circle") , tabName = "Welcome"),
      bs4Dash::menuItem("Single Dataset Analysis",icon = icon('palette'),
                        menuSubItem("Datasets overview", tabName = "dataset"),
                        menuSubItem("Single Gene expression", tabName = "Single"),
                        menuSubItem("Multi-Gene expression", tabName = "Multiple"),
                        menuSubItem("Correlation Analysis", tabName = "correlation"),
                        menuSubItem("Sample information", tabName = "info")
                        # menuSubItem("DEPs/DEGs Analysis", tabName = "DEGs")
      ),
      bs4Dash::menuItem("Multi-Datasets Analysis",icon = icon("deezer"),
                        menuSubItem("Multi-datasets expression", tabName = "multidata"),
                        menuSubItem("Correlation analysis", tabName = "Genelist_corr"),
                        menuSubItem("Immune infiltration", tabName = "Immune_infiltration"),
                        menuSubItem("Drug sensitivity", tabName = "drug_sensitivity")
      ),
      bs4Dash::menuItem("DEG & GSEA Analysis",icon = icon("list"),
                        menuSubItem("DEG analysis", tabName = "DEG"),
                        menuSubItem("Co-expression", tabName = "co_expression"),
                        menuSubItem("GSEA enrichment", tabName = "GSEA")
      ),
      bs4Dash::menuItem("Integrative analysis",icon = icon("list"),
                        menuSubItem("Venn diagram", tabName = "venn"),
                        menuSubItem("RobustRankAggreg", tabName = "RRA"),
                        menuSubItem("ComBat datasets", tabName = "combat")
      ),
      bs4Dash::menuItem("Feedback",icon = icon('envelope') , tabName ="feedback")
    )
  ),

  bs4Dash::dashboardBody(
    bs4Dash::tabItems(
      tabItem("Welcome",ui.modules_dash("Welcome")),
      tabItem("dataset",ui.modules_GEO_datasets("dataset")),
      tabItem("Single",ui.modules_Cancer_expression("Single")),
      tabItem("correlation",ui.modules_Cancer_corr("correlation")),
      tabItem("Multiple",ui.modules_multi_gene("Multiple")),
      tabItem("info",ui.modules_dataset_info("info")),
      tabItem("Immune_infiltration",ui.modules_gcas_til("Immune_infiltration")),
      tabItem("Genelist_corr",ui.modules_gcas_corr("Genelist_corr")),
      tabItem("multidata",ui.modules_multidata_dist("multidata")),
      tabItem("drug_sensitivity",ui.modules_gcas_drug("drug_sensitivity")),
      tabItem("DEG",ui.modules_DEG("DEG")),
      tabItem("co_expression",ui.modules_co_expression("co_expression")),
      tabItem("GSEA",ui.modules_GSEA("GSEA")),
      tabItem("venn",ui.modules_venn("venn")),
      tabItem("RRA",ui.modules_RRA("RRA")),
      tabItem("combat",ui.modules_combat("combat")),
      tabItem("feedback",mod_feedback_ui("feedback"))

    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # 创建一个全局的 reactiveValues 对象来存储共享变量
  shared_values <- reactiveValues(
    type_select = NULL,
    subtypes = NULL,
    datasets_text = NULL,
    datasets_select = NULL
  )
  callModule(server.modules_dash, "Welcome")
  callModule(server.modules_GEO_datasets, "dataset", shared_values)
  callModule(server.modules_Cancer_expression, "Single", shared_values)
  callModule(server.modules_Cancer_corr, "correlation", shared_values)
  callModule(server.modules_dataset_info, "info", shared_values)
  callModule(server.modules_multi_gene, "Multiple", shared_values)
  #  callModule(server.modules_cptac_gene_cor, "Pancan_correlation")
  callModule(server.modules_multidata_dist, "multidata")
  callModule(server.modules_gcas_til, "Immune_infiltration")
  callModule(server.modules_gcas_corr, "Genelist_corr")
  callModule(server.modules_gcas_drug, "drug_sensitivity")
  callModule(server.modules_DEG, "DEG")
  callModule(server.modules_GSEA, "GSEA")
  callModule(server.modules_venn, "venn")
  callModule(server.modules_RRA, "RRA")
  callModule(server.modules_combat, "combat")
  callModule(server.modules_co_expression, "co_expression")
  callModule(mod_feedback_server, "feedback")
}

# Run the application
shinyApp(ui = ui, server = server)
