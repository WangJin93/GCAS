ui.modules_gcas_til <- function(id) {
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
                label = "Select dataset:",
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
                inputId = ns("immune_type"),
                label = "Immune infiltration algorithm:",
                selected = c("TIMER"),
                choices = TIL_map$algorithm %>% unique(),
                multiple = F
              ),
              selectInput(
          inputId = ns("cor_method"),
          label = "Select Correlation method",
          choices = c("spearman", "pearson"),
          selected = "pearson"
        ),
        tags$hr(style = "border:none; border-top:2px solid #5E81AC;"),
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
      ),
      column(9,
             shinycssloaders::withSpinner(plotOutput(ns("hm_gene_immune_cor"), height = "auto")),
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
        # p("1. ", tags$a(href = "http://timer.cistrome.org/", "TIL data source")),
        # p("2. ", tags$a(href = "https://pancanatlas.xenahubs.net/", "Genomic profile data source")),
        DT::DTOutput(outputId = ns("tbl")),
      )
    )
  )
}

server.modules_gcas_til <- function(input, output, session) {
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
      nn <- input$Type
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
      choices = get_genelist(),
      selected = "GAPDH",
      server = TRUE
    )
  })


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("hm_gene_immune_cor"), html = waiter::spin_hexdots(), color = "white")

  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      Pancan_search <- input$Pancan_search
      if (length(input$datasets_text) <2){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "Please input at least 2 datasets!"
        ))

        return(NULL)
      }
      df <- get_expr_data(genes = Pancan_search,datasets = input$datasets_text)
      if (is.null(df)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "No results were returned, please comfirm your inputed gene symbol."
        ))
        return(NULL)
      }

      p <- cor_gcas_TIL(df,
                              cor_method = input$cor_method,
                              TIL_type = input$immune_type
      )
    }
    return(p)
  })
  plot_data <- eventReactive(input$search_bttn, {
      p <- plot_func()
    if (is.null(p)) return(NULL)
    
    # 使用 plist 返回格式：转换矩阵为长格式数据框
    # plist 包含: r, p, p_adj, n, t, ci_lower, ci_upper, sss
    # 矩阵结构: 行 = 数据集, 列 = 免疫细胞
    
    # 获取矩阵维度
    n_datasets <- nrow(p$r)
    n_cells <- ncol(p$r)
    total_rows <- n_datasets * n_cells
    
    # 创建 Dataset 和 Cell_Type 列
    Dataset <- rep(rownames(p$r), each = n_cells)
    Cell_Type <- rep(colnames(p$r), n_datasets)
    
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
      Cell_Type = Cell_Type,
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
    
    # 格式化数字精度为小数点后3位
    format_num <- function(x, digits = 3) {
      ifelse(is.na(x), NA, sprintf(paste0("%.", digits, "f"), x))
    }
    
    # 格式化 p 值和 p.adj（小于 0.001 显示为 <0.001）
    format_pval <- function(x) {
      ifelse(is.na(x), NA, ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
    }
    
    result_df$n <- format_num(result_df$n)
    result_df$r <- format_num(result_df$r)
    result_df$t <- format_num(result_df$t)
    result_df$p <- format_pval(result_df$p)
    result_df$p.adj <- format_pval(result_df$p.adj)
    result_df$ci_lower <- format_num(result_df$ci_lower)
    result_df$ci_upper <- format_num(result_df$ci_upper)
    
    # 重命名列以便显示
    colnames(result_df) <- c("Dataset","Cell_Type",  "n", "r", "t", "p", "p.adj", "ci_lower", "ci_upper", "Significance")
    
    return(result_df)
  })
  TIL_corr_func <- eventReactive(input$tbl_rows_selected, {
    s <- input$tbl_rows_selected
    if (length(s)) {
      selected <- plot_data()[s,]
      checkpoint <- selected$Cell_Type
      # 使用 plist 返回格式：根据选中的数据集获取数据
      df <- plot_func()$sss[[selected$Dataset]]
      # 找到包含目标细胞类型的列
      target_col_idx <- which(colnames(df) == checkpoint)
      if (length(target_col_idx) == 0) {
        return(NULL)
      }
      print(df)
      df <- df[,c("ID","tissue","subtype","dataset",colnames(df)[3],checkpoint)]
    }
    rownames(df)<-NULL
    return(df)
  })
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$hm_gene_immune_cor <- renderPlot(width = width_scatter,
                                          height = height_scatter,{
    w$show() # Waiter add-ins
                                            if (!is.null(plot_func())){
                                              # 使用 plist 返回格式
                                              viz_cor_heatmap(plot_func()$r, plot_func()$p)

                                            }

  })

  observeEvent(eventExpr = input$tbl_rows_selected,{
    s <- input$tbl_rows_selected

    if (length(s)) {
      showModal(
        modalDialog(
          title = paste("Gene-Immunue infitration scatter plot", input$Pancan_search),
          size = "l",
          fluidPage(
            column(
              12,align = "center",
              plotOutput(ns("scatter_plot"),width = "100%",height = "auto",),
              DT::DTOutput(outputId = ns("scatter_corr")),
            )
          )
        )
      )

      output$scatter_plot <- renderPlot(width = 450,
                                        height = 300,{
                                          # 使用新的返回格式
                                          cor_result <- TIL_corr_func()
                                          viz_corplot(cor_result, colnames(cor_result)[5], colnames(cor_result)[6],
                                                     method = input$cor_method, x_lab = " expression", y_lab = "")

                                        })
      # output$scatter_corr <- DT::renderDataTable(
      #   TIL_corr_func(), server = TRUE,selection = 'single'
      # )
      output$scatter_corr <- DT::renderDataTable(server = FALSE, {
        DT::datatable(
          TIL_corr_func(),
          rownames = T,
          extensions = c("Buttons"),
          options = list(
            pageLength = 5,
            dom = "Bfrtip",
            buttons = list(
              list(
                extend = "csv", text = "Download table", filename =  paste0(names(TIL_corr_func())[4],"_",names(TIL_corr_func())[5],"_",TIL_corr_func()[1,3]),
                exportOptions = list(
                  modifier = list(page = "all")
                )
              )
            )
          )
        )
      })

      ## downloadTable
    }
  })



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
              extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,"_",input$Type,"_cor_TIL"),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )
    }
  })


  ## downloadTable
  # output$downloadTable <- downloadHandler(
  #   filename = function() {
  #     paste0(input$Pancan_search, "_", input$profile, "_pancan_TIL.csv")
  #   },
  #   content = function(file) {
  #     write.csv( plot_data(), file, row.names = FALSE)
  #   }
  # )

  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search,"_",input$Type,"TIL_cor.pdf")
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
                inputId = ns("gene_select"), label = "Choose a specific cell:",
                choices = TIL_map[TIL_map$algorithm == input$immune_type,]$cell_type
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
