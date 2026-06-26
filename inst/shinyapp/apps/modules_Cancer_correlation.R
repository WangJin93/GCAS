ui.modules_Cancer_corr <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "",
    icon = icon("database"),
    useShinyjs(),  # Set up shinyjs

    sidebarLayout(
      sidebarPanel(
        width = 3,
        HTML("<p><strong>Cancer types:</strong></p>"),
        textOutput(ns("type_selected")),
        br(),
        uiOutput(ns("datasets_text")),
        hr(),
        h4("Variable A:"),
        selectizeInput(
          inputId = ns("ga_id"), # molecule identifier
          label = "Input gene symbol:",
          choices = c("GAPDH"),
          options = list(
            create = TRUE,
            maxOptions = 5,
            placeholder = "e.g. GAPDH",
            plugins = list("restore_on_backspace")
          )
        ),
        hr(),
        h4("Variable B:"),
        selectizeInput(
          inputId = ns("ga_cor_id"), # molecule identifier
          label = "Input gene B symbol:",
          choices = c("GAPDH","TNS1"),
          multiple = TRUE,
          options = list(
            create = TRUE,
            placeholder = "You can enter multiple genes",
            plugins = list("restore_on_backspace")
          )
        )
      ),
      mainPanel(
        fluidPage(
          fluidRow(
            column(
              4,
              wellPanel(
                h4("Parameters:"),

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
                actionBttn(
                  inputId = ns("ga_cor_submit"),
                  label = "Submit",
                  style = "gradient",
                  icon = icon("check"),
                  color = "default",
                  block = TRUE,
                  size = "sm"
                )
              ),
              tags$br(),
              wellPanel(
                h4("Download figure:"),
                numericInput(inputId = ns("cor_height"), label = "Height", value = 5),
                numericInput(inputId = ns("cor_width"), label = "Width", value = 8),
                downloadBttn(
                  outputId = ns("cor_download"),
                  style = "gradient",
                  color = "default",
                  block = TRUE,
                  size = "sm"
                )
              )
            ),
            column(
              8,

              shinycssloaders::withSpinner(DT::dataTableOutput(ns("ga_cor_data"))),
              h4("Scatter plot of selected row:"),
              column(
                11,
                shinycssloaders::withSpinner(plotOutput(ns("ga_cor_output"),width = "600",height = "auto")),
              ),
              column(
                11,
                shinycssloaders::withSpinner(DT::dataTableOutput(ns("cor_scatter_data")))
              )
            ),


          ))
      )
    )
  )
}

server.modules_Cancer_corr <- function(input, output, session, shared_values) {
  ns <- session$ns
  # df <- dataset
  # df <- dataset
  observe({

    updateSelectizeInput(
      session,
      "ga_id",
      choices = c("GAPDH","TNS1"),
      selected = "GAPDH",
      server = TRUE
    )
  })
  output$type_selected  <- renderText({
    if (!is.null(input$datasets_text) && input$datasets_text == "Integrated data (ComBat)") {
      "integrated data"
    } else {
      shared_values$subtypes
    }
  })

  output$datasets_text <- renderUI({
    # 检查是否有整合数据
    has_combat <- !is.null(shared_values$combat_data)

    # 构建数据集选项
    dataset_choices <- shared_values$datasets_text

    # 如果有整合数据，添加到选项中
    if (has_combat) {
      dataset_choices <- c("Integrated data (ComBat)" = "Integrated data (ComBat)", dataset_choices)
    }

    selectInput(
      inputId = ns("datasets_text"),
      label = "Select dataset:",
      choices = dataset_choices,
      selected = shared_values$datasets_select,
      multiple = FALSE
    )
  })

  observe({

    updateSelectizeInput(
      session,
      "ga_cor_id",
      choices = get_genelist(),
      selected = "FOXM1",
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "ga_id",
      choices = get_genelist(),
      selected = "GAPDH",
      server = TRUE
    )
  })



  #########correlation
  # Show waiter for plot
  w1 <- waiter::Waiter$new(id = ns("ga_cor_output"), html = waiter::spin_hexdots(), color = "white")
  plot_cor_func <- eventReactive(input$ga_cor_submit, {
    id1 <- input$ga_id
    id2 <- input$ga_cor_id
    datasets_text <- input$datasets_text

    # 根据选中的数据集来判断是否使用整合数据
    if (datasets_text == "Integrated data (ComBat)") {
      # 使用 ComBat 整合数据
      if (is.null(shared_values$combat_data)) {
        sendSweetAlert(
          session,
          title = "Error",
          text = "No integrated data available. Please run ComBat analysis first in 'Integrative analysis' tab.",
          type = "error"
        )
        return(NULL)
      }

      combat_data <- shared_values$combat_data
      sample_info <- shared_values$combat_sample_info

      # 恢复 sample_info 中的 ID 列
      sample_info <- tibble::rownames_to_column(sample_info, "ID")

      # 检查基因是否存在于数据中
      all_genes <- c(id1, id2)
      missing_genes <- all_genes[!all_genes %in% rownames(combat_data)]
      if (length(missing_genes) > 0) {
        sendSweetAlert(
          session,
          title = "Error",
          text = paste("Gene(s)", paste(missing_genes, collapse = ", "), "not found in the integrated data."),
          type = "error"
        )
        return(NULL)
      }

      # 准备相关性分析数据
      # 过滤样本类型
      sample_types <- input$sample_type
      if (length(sample_types) == 0) {
        sample_types <- c("Tumor", "Normal")
      }

      # 过滤样本信息
      sample_info_filtered <- sample_info[sample_info$type %in% sample_types, , drop = FALSE]

      if (nrow(sample_info_filtered) == 0) {
        sendSweetAlert(
          session,
          title = "Error",
          text = "No samples match the selected sample type(s).",
          type = "error"
        )
        return(NULL)
      }

      # 获取共同的样本
      common_samples <- intersect(sample_info_filtered$ID, colnames(combat_data))
      if (length(common_samples) < 3) {
        sendSweetAlert(
          session,
          title = "Error",
          text = "Not enough samples for correlation analysis.",
          type = "error"
        )
        return(NULL)
      }

      # 提取表达数据
      expr_data <- combat_data[all_genes, common_samples, drop = FALSE]
      expr_data <- t(expr_data)
      expr_data <- as.data.frame(expr_data) %>% tibble::rownames_to_column("ID")
      print(head(expr_data))
      print(head(sample_info))
      expr_data <- merge( sample_info[-6],expr_data,by="ID")
      expr_data <- merge( expr_data,sample_info[c(1,6)],by="ID")
      print(head(expr_data))

      # 先获取所有需要的基因对，然后批量计算调整 p 值
      all_pairs <- list()
      for (i in 1:length(id1)) {
        for (j in 1:length(id2)) {
          g1 <- id1[i]
          g2 <- id2[j]
          if (g1 != g2) {
            all_pairs[[length(all_pairs) + 1]] <- list(g1 = g1, g2 = g2)
          }
        }
      }

      # 计算每对基因的相关性
      cor_results_list <- list()
      p_values <- c()
      for (pair_idx in 1:length(all_pairs)) {
        pair <- all_pairs[[pair_idx]]
        g1 <- pair$g1
        g2 <- pair$g2

        x <- expr_data[, g1]
        y <- expr_data[, g2]

        # 移除 NA 值
        valid_idx <- !is.na(x) & !is.na(y)
        x <- x[valid_idx]
        y <- y[valid_idx]

        if (length(x) < 3) next

        # 计算相关性
        if (input$cor_method == "pearson") {
          test_result <- cor.test(x, y, method = "pearson")
        } else {
          test_result <- cor.test(x, y, method = "spearman")
        }

        # 获取置信区间
        ci <- test_result$conf.int

        cor_results_list[[length(cor_results_list) + 1]] <- data.frame(
          Symbol = g2,
          n = length(x),
          r = test_result$estimate,
          t = ifelse(input$cor_method == "pearson", test_result$statistic, NA),
          p = test_result$p.value,
          ci_lower = ci[1],
          ci_upper = ci[2],
          stringsAsFactors = FALSE
        )
        p_values <- c(p_values, test_result$p.value)
      }

      # 合并结果并计算调整后的 p 值
      cor_result <- do.call(rbind, cor_results_list)
      if (nrow(cor_result) > 0) {
        cor_result$p.adj <- p.adjust(p_values, method = "BH")
        # 调整列顺序，与 cor_cancer_genelist 一致
        cor_result <- cor_result[, c("Symbol", "n", "r", "t", "p", "p.adj", "ci_lower", "ci_upper")]
      }

      # 准备绘图数据
      return(list(
        cor_result = cor_result,
        cor_data = expr_data
      ))

    } else {
      # 使用原始单数据集
      nn<- shared_values$subtypes

      results <- cor_cancer_genelist(datasets_text,
                                      id1,
                                      id2,
                                     extract_subset(subtype,nn),
                                     input$sample_type,
                                      input$cor_method)
      return(results)
    }
  })

corr_func <- eventReactive(input$ga_cor_data_rows_selected , {
    s <- input$ga_cor_data_rows_selected
    if (length(s)) {
      # 使用新的返回格式：cor_data
      result <- plot_cor_func()
      if (is.null(result)) {
        return(NULL)
      }

      cor_data <- result$cor_data
      if (is.null(cor_data) || nrow(cor_data) == 0) {
        return(NULL)
      }

      # 使用列名而不是硬编码索引
      # 获取ID列和type列的位置
      id_col <- which(colnames(cor_data) %in% c("ID","dataset","subtype","tissue","Patient.ID","type"))
print(head(cor_data))
      # 获取基因列（排除ID和type列）
      gene_cols <- setdiff(1:ncol(cor_data), c(id_col))
      print(gene_cols)

      # 获取选中的基因对
      # s是选中行的索引，我们需要映射到基因对
      selected_genes <- gene_cols[c(1,s+1)]
      print(selected_genes)

      # 构建新的数据框：ID + 选中的基因 + type
      df <- cor_data[, c(id_col, selected_genes)]
      print(head(df))
    }
    return(df)


  })




output$ga_cor_output <- renderPlot(width = 500,
                                   height = 300,{
                                      w1$show() # Waiter add-ins
                                     if (length(input$ga_cor_data_rows_selected)) {
                                                                        cor_data <- corr_func() %>% na.omit()
                                       viz_corplot(cor_data, colnames(cor_data)[7], colnames(cor_data)[8],
                                                  method = input$cor_method, x_lab = " expression", y_lab = " expression")
                                     }else{NULL}
                                    }
)
output$cor_download <- downloadHandler(
  filename = function() {
    paste0(colnames(corr_func())[7],"_",colnames(corr_func())[8],"_", input$datasets_text, ".pdf")
  },
  content = function(file) {
    cor_data <- corr_func() %>% na.omit()

    pdf(file = file, onefile = FALSE,height = input$cor_height,width = input$cor_width)
    print(viz_corplot(cor_data, colnames(cor_data)[7], colnames(cor_data)[8],
                      method = input$cor_method, x_lab = " expression", y_lab = " expression"))
    dev.off()

  }
)

  output$ga_cor_data <- DT::renderDataTable(server = FALSE, {
    if (input$ga_cor_submit>0){
      # 获取相关性结果
      result <- plot_cor_func()

      # 检查是否有结果
      if (is.null(result)) {
        return(NULL)
      }

      cor_result <- result$cor_result

      # 检查cor_result是否为空
      if (is.null(cor_result) || nrow(cor_result) == 0) {
        return(NULL)
      }

      # 格式化 p 值和调整后的 p 值：小于 0.001 的显示为 "<0.001"
      cor_result$p <- ifelse(cor_result$p < 0.001, "<0.001", sprintf("%.3f", cor_result$p))
      cor_result$p.adj <- ifelse(cor_result$p.adj < 0.001, "<0.001", sprintf("%.3f", cor_result$p.adj))

      # 格式化其他数值列保留三位小数
      cor_result$r <- sprintf("%.3f", cor_result$r)
      cor_result$t <- sprintf("%.3f", cor_result$t)
      cor_result$ci_lower <- sprintf("%.3f", cor_result$ci_lower)
      cor_result$ci_upper <- sprintf("%.3f", cor_result$ci_upper)

      DT::datatable(
        cor_result,
        rownames = FALSE,selection = 'single',
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename = paste(input$ga_id,"_", input$datasets_text,  "correlation results"),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

    }
  })
  output$cor_scatter_data <- DT::renderDataTable(server = FALSE, {
    s <- input$ga_cor_data_rows_selected
    if (length(s)) {
      DT::datatable(
        corr_func(),
        rownames = T,
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename = paste0(colnames(corr_func())[7],"_",colnames(corr_func())[8],"_", input$datasets_text ),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

      }
  })


  }
