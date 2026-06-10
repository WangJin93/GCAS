ui.modules_gcas_corr <- function(id) {
  ns <- NS(id)
  fluidPage(
        fluidRow(
          column(
            3,
            wellPanel(
              selectInput(
                inputId = ns("Type"),
                label = "Cancer type:",
                choices = c("Pancancer",unique(dataset_info %>%  .["type"])),
                selected = "Pancancer",
                multiple = F
              ),
              conditionalPanel(
                condition = "input.Type == 'Pancancer'",
                prettyRadioButtons(
                  inputId = ns("pancancer_datasets"),
                  label = "Pancancer datasets:",
                  choices = c("Panel 1", "Panel 2", "User defined"),
                  animation = "pulse",inline = T,
                  status = "info"
                ),ns=ns
              ),
              shinyTree(ns("subtype"), theme="proton",stripes=T, themeIcons = F,multiple=T, themeDots = T),
              br(),
              selectInput(
                inputId = ns("datasets_text"),
                label = "Select dataset::",
                choices = "",
                multiple = T
              ),
              hr(),
              selectizeInput(
                inputId = ns("Pancan_search"),
                label = "Input a gene symbol",
                choices = "GAPDH",
                width = "100%",
                options = list(
                  create = TRUE,
                  maxOptions = 5,
                  plugins = list("restore_on_backspace")
                )
              ),

              selectizeInput(
                inputId = ns("geneset"),
                label = "Gene sets:",
                selected = c("Immune checkpoint"),
                choices = c("Immune checkpoint","m6A methylation","Cell senescence","EMT","Cell senescence","Ferroptosis","Cuproptosis","User defined"),
                multiple = F
              ),
              textAreaInput(
                inputId = ns("genelist"), # molecule identifier
                label = "Input gene list:",
                value = "",
                placeholder = "Paste your genelist!\neg.\nTNS1\nTP53\nTXNIP",
                height = "200"
              ),
              hr(),
              wellPanel(
                selectInput(
                  inputId = ns("sample_type"),
                  label = "Only use normal or tumor",
                  choices = c("Normal","Tumor"),
                  selected = "Tumor", multiple = TRUE
                ),

                selectInput(
                  inputId = ns("cor_method"),
                  label = "Select Correlation method",
                  choices = c("spearman", "pearson"),
                  selected = "pearson"
                ),
                shinyWidgets::actionBttn(
                  inputId = ns("search_bttn"),
                  label = "Go!",
                  style = "gradient",
                  icon = icon("search"),
                  color = "primary",
                  block = TRUE,
                  size = "sm"
                )
              )

            )
          ),
      column(9,
             shinycssloaders::withSpinner(plotOutput(ns("genelist_cor"), height = "auto")),
        hr(),
        fluidRow(
          column(2,
                 br(),
                 shinyWidgets::actionBttn(
                   inputId = ns("show.single"),
                   label = "Individual plot!",
                   style = "gradient",
                   color = "default",
                   icon = icon("sliders"),
                   block = TRUE,
                   size = "sm"
                 )
          ),
          column(2, numericInput(inputId = ns("height_scatter"), label = "Height", value = 500, max = 600)),
          column(2, numericInput(inputId = ns("width_scatter"), label = "Width", value = 800)),
          column(2,
                 br(),
                 downloadBttn(
                   outputId = ns("download"),
                   style = "gradient",
                   color = "default",
                   block = TRUE,
                   size = "sm"
                 ))
        ),
        hr(),

        # h5("NOTEs:"),
        # p("1. ", tags$a("The immune checkpoint molecules (ICMs) were collected from the study conducted by Charoentong et al (43) who reported 24 immunoinhibitory genes and 45 immunostimulatory genes. ")),
        # p("2. ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/28052254/", "Charoentong P, Finotello F, Angelova M, et al. Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Rep. 2017;18(1):248-262. doi:10.1016/j.celrep.2016.12.019")),
        DT::DTOutput(outputId = ns("tbl")),
      )
    )
  )
}

server.modules_gcas_corr <- function(input, output, session) {
  ns <- session$ns
  output$subtype <- renderTree({
    if (input$Type != "Pancancer"){
      subtype[[input$Type]]
    }
  })
  observeEvent(c(input$Type,input$pancancer_datasets),{
    if (input$Type != "Pancancer"){
      ddd <-  dataset_info %>%
        dplyr::filter(type %in% input$Type)%>% .[,"dataset"] %>%  intersect(.,sample_subtype[sample_subtype$subtype %in% c(extract_subset(subtype,input$Type)),][,"dataset"])

      updateSelectInput(session,"datasets_text",
                        label = "Select dataset:",
                        choices = ddd,
                        selected = ddd
      )
    }else{
      ddd <- case_when(
        input$pancancer_datasets == "Panel 1"~ c("GSE41258", "GSE68848", "GSE76250", "GSE66229", "GSE31210", "GSE6919", "GSE25097", "GSE71729", "GSE26712", "GSE30784", "GSE13507", "GSE53624", "GSE63514", "GSE46517", "GSE40435", "GSE35570", "GSE17025","GSE50006","GSE12195","GSE99671"),
        input$pancancer_datasets == "Panel 2"~ c("GSE87211", "GSE50161", "GSE54002", "GSE29272", "GSE19188", "GSE70768", "GSE14520", "GSE62452", "GSE66957", "GSE25099", "GSE3167", "GSE23400", "GSE39001", "GSE15605", "GSE53757", "GSE33630", "GSE115810","GSE31048","GSE56315","GSE126209"),
        input$pancancer_datasets == "User defined"~ c("GSE31210")
      )
      updateSelectInput(session,"datasets_text",
                        label = "Select dataset:",
                        choices = dataset_info$dataset,
                        selected = ddd
      )
    }
  })
  subtype_selected <- reactive({
    tree <- input$subtype
    if (!is.null(tree)){
      req(tree)
      nn<- sapply(get_selected(tree, format = "classid"), function(x) x[1])
    }else{
      nn <- subtype[[input$Type]][[1]]
    }
    print(nn)
    nn
  })

  observeEvent(input$subtype,{

    nn<- subtype_selected()
    ddd <-  dataset_info %>%
      dplyr::filter(type %in% input$Type)%>% .[,"dataset"] %>%  intersect(.,sample_subtype[sample_subtype$subtype %in% c(extract_subset(subtype,nn)),][,"dataset"])

    updateSelectInput(session,"datasets_text",
                      label = "Select dataset:",
                      choices = ddd,
                      selected = ddd
    )

  })
  observe({

    updateSelectizeInput(
      session,
      "Pancan_search",
      choices = genelist,
      selected = "GAPDH",
      server = TRUE
    )
  })

  observe({
    genelist <- case_when(
      input$geneset == "Immune checkpoint" ~ c("SIRPA\nCTLA4\nTIGIT\nLAG3\nVSIR\nLILRB2\nSIGLEC7\nHAVCR2\nLILRB4\nPDCD1\nBTLA"),
      input$geneset == "m6A methylation" ~ c("METTL3\nMETTL14\nWTAP\nKIAA1429\nZC3H13\nRBM15\nRBM15B\nFTO\nALKBH5\nYTHDC1\nYTHDC2\nYTHDF1\nYTHDF2\nYTHDF3\nHNRNPC\nHNRNPA2B1\nIGF2BP1\nIGF2BP2\nIGF2BP3"),
      input$geneset == "Cell senescence" ~ c("GAPDH\nCDKN1A\nCDKN2A\nLMNB1\nMKI67"),
      input$geneset == "EMT" ~ c("CDH1\nCDH2\nSNAIL\nVIM"),
      input$geneset == "Ferroptosis" ~ c("ACSL4\nLPCAT3\nTFR1\nSLC7A11\nGPX4\nFTH1"),
      input$geneset == "Cuproptosis" ~ c("FDX1\nLIAS\nLIPT1\nDLD\nDLAT\nPDHA1\nPDHB\nMTF1\nGLS\nCDKN2A\nATP7B\nSLC31A1"),
      input$geneset == "User defined" ~ c("")
    )
    updateTextAreaInput(session, "genelist",
                        label = "Input gene list:",
                        value = genelist)
  })

  genelists <- reactive({strsplit(input$genelist,"\n")[[1]]})
  geneset_data <- eventReactive(input$search_bttn, {
    ga_ids <-  genelists()
    # print(ga_ids)
    if (length(input$datasets_text) <2){
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "Please input at least 2 datasets!"
      ))

      return(NULL)
    }
    dd <- get_expr_data(genes  = ga_ids,datasets = input$datasets_text)
    if (is.null(dd)) {
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "No results were returned, because the current dataset does not contain the gene symbols you entered."
      ))
      return(NULL)

    } else {
      # 验证必要列是否存在
      required_cols <- c("ID", "subtype", "dataset", "tissue", "Patient.ID")
      if (!all(required_cols %in% colnames(dd))) {
        missing <- setdiff(required_cols, colnames(dd))
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          paste0("Missing required columns: ", paste(missing, collapse = " ,"))
        ))
        return(NULL)
      }

      # 获取实际返回的基因列
      non_expr_cols <- c(required_cols, "type")
      actual_gene_cols <- setdiff(colnames(dd), non_expr_cols)

      if (length(actual_gene_cols) < length(ga_ids)) {
        missing_genes <- setdiff(ga_ids, actual_gene_cols)
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          paste0("The following gene symbols are not found in the selected datasets: ", paste(missing_genes, collapse = " ,"))
        ))
      }

      if (length(actual_gene_cols) == 0) {
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "No gene symbols were found in the selected datasets."
        ))
        return(NULL)
      }
    }
    return(dd)
  })

  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("genelist_cor"), html = waiter::spin_hexdots(), color = "white")
  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      Pancan_search <- input$Pancan_search
      df <- get_expr_data(genes = Pancan_search,datasets = input$datasets_text)
      if (is.null(df)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "No results were returned, please comfirm your inputed gene symbol."
        ))
        return(NULL)
      }
      nn<- subtype_selected()
      gd <- geneset_data()
      if (is.null(gd)){
        return(NULL)
      }

      # 验证 geneset_data 返回的数据格式
      required_cols <- c("ID", "subtype", "dataset", "tissue", "Patient.ID")
      gene_cols <- setdiff(colnames(gd), required_cols)
      if (length(gene_cols) < 2) {
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "Geneset must contain at least 2 genes for correlation analysis."
        ))
        return(NULL)
      }

      p <- cor_gcas_genelist(df, gd,
                             tumor_subtype = nn,
                             sample_type = input$sample_type,
                             cor_method = input$cor_method)
    }
    return(p)
  })

#
  plot_data <- eventReactive(input$search_bttn, {
      p <- plot_func()

    if (is.null(p)) return(NULL)

    # 使用 plist 返回格式：转换矩阵为长格式数据框
    # plist 包含: r, p, p_adj, n, t, ci_lower, ci_upper, sss
    # 矩阵结构: 行 = 数据集, 列 = 基因

    # 获取矩阵维度
    n_datasets <- nrow(p$r)
    n_genes <- ncol(p$r)
    total_rows <- n_datasets * n_genes

    # 创建 Dataset 和 Symbol 列
    Dataset <- rep(rownames(p$r), each = n_genes)
    Symbol <- rep(colnames(p$r), n_datasets)

    # 提取统计量，确保长度一致
    extract_stat <- function(mat) {
      vec <- as.numeric(mat)
      if (length(vec) != total_rows) {
        warning(paste("Matrix dimension mismatch, expected", total_rows, "got", length(vec)))
        vec <- rep(NA, total_rows)
      }
      return(vec)
    }

    n <- extract_stat(p$n)
    r <- extract_stat(p$r)
    t_val <- extract_stat(p$t)
    p_val <- extract_stat(p$p)
    p_adj <- extract_stat(p$p_adj)
    ci_lower <- extract_stat(p$ci_lower)
    ci_upper <- extract_stat(p$ci_upper)

    result_df <- data.frame(
      Symbol = Symbol,
      Dataset = Dataset,
      n = n,
      r = r,
      t = t_val,
      p = p_val,
      p.adj = p_adj,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )

    # 添加显著性标记
    result_df$pstar <- ifelse(result_df$p.adj < 0.05,
                             ifelse(result_df$p.adj < 0.001, "***", ifelse(result_df$p.adj < 0.01, "**", "*")),
                             "")

    # 重命名列以便显示
    colnames(result_df) <- c( "Dataset","Gene/Symbol", "n", "r", "t", "p", "p.adj", "ci_lower", "ci_upper", "Significance")

    return(result_df)
  })
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  output$genelist_cor <- renderPlot(width = width_scatter,
                                          height = height_scatter,{
    w$show() # Waiter add-ins
                                            if (!is.null(plot_func())){
                                              # 使用 plist 返回格式
                                              viz_cor_heatmap(plot_func()$r, plot_func()$p)
                                            }
                                          })

  corr_func <- eventReactive(input$tbl_rows_selected, {
    s <- input$tbl_rows_selected
    if (length(s)) {
      selected <- plot_data()[s,]
      #data<- df[which(df$type == type),]
      checkpoint <- selected$`Gene/Symbol`
      # 使用 plist 返回格式：从 sss 获取合并数据
      df <- plot_func()$sss[[selected$Dataset]]
      # 找到包含目标基因的列
      target_col_idx <- which(colnames(df) == checkpoint)
      if (length(target_col_idx) == 0) {
        return(NULL)
      }
      df <- df[,c("ID","subtype","dataset","type",colnames(df)[6],checkpoint)]
      print(df$dataset %>% unique())
      return(df)
    }
    return(df)
  })


#
  observeEvent(eventExpr = input$tbl_rows_selected,{
    s <- input$tbl_rows_selected

    if (length(s)) {
      showModal(
        modalDialog(
          title = paste("Correlation scatter plot"),
          size = "l",
          fluidPage(
            column(
              12,
              plotOutput(ns("scatter_plot"),width = "100%",height = "auto",),
              DT::DTOutput(outputId = ns("scatter_corr")),

            )
          )
        )
      )

      output$scatter_plot <- renderPlot(width = 450,
                                        height = 300,{
                                          # 使用新的返回格式
                                          cor_result <- corr_func()
                                          viz_corplot(cor_result, colnames(cor_result)[5], colnames(cor_result)[6],
                                                     method = input$cor_method, x_lab = " expression", y_lab = " expression")

                                        })
      output$scatter_corr <- DT::renderDataTable(server = TRUE,{
        DT::datatable(
          corr_func(),
          rownames = T,
          extensions = c("Buttons"),
          options = list(
            pageLength = 5,
            dom = "Bfrtip",
            buttons = list(
              list(
                extend = "csv", text = "Download table", filename =  paste0(names(corr_func())[5],"_",names(corr_func())[6],"_",corr_func()[1,3]),
                exportOptions = list(
                  modifier = list(page = "all")
                )
              )
            )
          )
        )

      }
      )
      ## downloadTable
    }
  })

#
#   observeEvent(input$search_bttn, {
#     if (nchar(input$Pancan_search) >= 1) {
#       shinyjs::show(id = "save_csv")
#     } else {
#       shinyjs::hide(id = "save_csv")
#     }
#   })
#
#
  output$tbl <- DT::renderDataTable(server = FALSE, {
    if (!is.null( plot_data())){
      DT::datatable(
        plot_data(),
        rownames = T,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,"_",input$Type,"_cor"),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )
    }
  })

  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search,"_",input$Type,"_cor",".pdf")
    },
    content = function(file) {
        pdf(file,  width =  input$width_scatter/70 ,height = input$height_scatter/70)
      # 使用 plist 返回格式
      print(viz_cor_heatmap(plot_func()$r, plot_func()$p))
      dev.off()
    }
  )
  observeEvent(eventExpr = input$show.single,{

    showModal(
      modalDialog(
        title = paste(input$Pancan_search," correlation volcano plot in multi-datasets" ),
        size = "l",
        fluidPage(
          div(style="display:flex",
              selectInput(
                inputId = ns("gene_select"), label = "Choose a specific gene:",
                choices = genelists()
              ),
              numericInput(inputId = ns("xinter"), label = "r cutpoint:", value = 0.2,max = 0.6,width = 200),

              selectInput(
                inputId = ns("yinter"),
                label = "P cutpoint:",
                choices = c("0.05", "0.005","0.001"),
                selected = "0.05",width = 200
              )
              ),

            div(style="display:flex",
                textInput(ns('values'), "colors (',' seperated)", value = 'green,black,red',width = 200),
                numericInput(inputId = ns("height_scatter1"), label = "Height", value = 500,max = 800,width = 200),
                numericInput(inputId = ns("width_scatter1"), label = "Width", value = 600,width = 200)
            ),


            plotOutput(ns("scatter_plot1"),width = "100%",height = "auto",),
            DT::DTOutput(outputId = ns("scatter_corr1")),



        )
      )
    )



    width_scatter1 <- reactive ({ input$width_scatter1 })
    height_scatter1 <- reactive ({ input$height_scatter1 })
    output$scatter_plot1 <- renderPlot(width = width_scatter1,
                                       height = height_scatter1,{
                                         # 使用新的返回格式
                                         viz_cor_volcano(plot_func(), item = input$gene_select,
                                                         r.cut = input$xinter,
                                                         p.cut = as.numeric(input$yinter),
                                                         colors = unlist(strsplit(input$values, ',')))
                                       })

    output$scatter_corr1 <- DT::renderDataTable(server = F, {
      # 获取 plot_func 的结果
      plot_result <- plot_func()
      if (is.null(plot_result)) {
        return(NULL)
      }
      
      gene <- input$gene_select
      r.cut = input$xinter
      p.cut = as.numeric(input$yinter)

      # 从矩阵格式转换为数据框格式（与 viz_cor_volcano 一致）
      if (gene %in% rownames(plot_result$r)) {
        # 只使用 r 和 p 列，然后进行 na.omit
        cor_result_df <- data.frame(
          Symbol = colnames(plot_result$r),
          r = plot_result$r[gene, ],
          p = plot_result$p[gene, ],
          p.adj = plot_result$p_adj[gene, ]
        ) %>% na.omit()
      } else {
        return(NULL)
      }

      # Determine change status
      cor_result_df$change <- ifelse(cor_result_df$p.adj < p.cut,
                                     ifelse(cor_result_df$r < -r.cut, "negative",
                                            ifelse(cor_result_df$r > r.cut, "positive", "no")),
                                     "no")
      # Arrange results by r
      results <- dplyr::arrange(cor_result_df, r)
      DT::datatable(
        results,
        rownames = F,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,  "_",gene,"_",input$Type),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

    }
    )

  })

}
