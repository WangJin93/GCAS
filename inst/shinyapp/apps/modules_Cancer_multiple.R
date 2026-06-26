ui.modules_multi_gene <- function(id) {
  ns <- NS(id)
  fluidPage(
fluidRow(column(3,
                wellPanel(
                  HTML("<p><strong>Cancer types:</strong></p>"),
                  textOutput(ns("type_selected")),
                  br(),
                  uiOutput(ns("datasets_text")),
                  hr(),
                  textAreaInput(
                    inputId = ns("ga_ids"), # molecule identifier
                    label = "Input gene symbols:",
                    value = "YTHDC1\nYTHDC2\nYTHDF1\nYTHDF2\nYTHDF3\nIGF2BP1\nIGF2BP2\nIGF2BP3",
                    placeholder = "Paste your genelist!\neg.\nTNS1\nTP53\nTXNIP",
                    height = "150"
                  ) ,
                  selectInput(inputId = ns("method"), label = "Select statistic method", choices = c("wilcox.test", "t.test"), selected = "t.test"),

        materialSwitch(ns("Show.P.value"), "Show P value", inline = TRUE),
        materialSwitch(ns("Show.P.label"), "Show P label", inline = TRUE),
        materialSwitch(ns("Show.n"), "Show sample size", inline = T),
        materialSwitch(ns("x_axis45"), "X-axis 45°", inline = T),
        colourpicker::colourInput(inputId = ns("tumor_col"), "Normal sample color", "#00AFBB"),
        colourpicker::colourInput(inputId = ns("normal_col"), "Tumor sample color", "#FC4E07"),
        tags$hr(style = "border:none; border-top:2px solid #5E81AC;"),
        shinyWidgets::actionBttn(
          inputId = ns("search_bttn"),
          label = "Go!",
          style = "gradient",
          icon = icon("search"),
          color = "primary",
          block = TRUE,
          size = "sm"
        ))
      ),
      column(9,
             shinycssloaders::withSpinner(plotOutput(ns("gene_pancan_dist"), height = "auto")),
             fluidRow(
               column(2, numericInput(inputId = ns("height_scatter"), label = "Height", value = 500, max = 600)),
               column(2, numericInput(inputId = ns("width_scatter"), label = "Width", value = 800)),
               column(4,
                      br(),
                      downloadBttn(
                        outputId = ns("download"),
                        style = "gradient",
                        color = "default",
                        block = TRUE,
                        size = "md"
                      ))
             ),
             hr(),
             DT::DTOutput(outputId = ns("tbl")),
        #      shinyjs::hidden(
        #        wellPanel(
        #          id = ns("save_csv"),
        #          downloadButton(ns("downloadTable"), "Save as csv")
        #        )
        # )
      )
    )
  )
}

server.modules_multi_gene <- function(input, output, session, shared_values) {
  ns <- session$ns

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

  colors <- reactive({
    c(input$tumor_col, input$normal_col)
  })


  ga_ids <- reactive({strsplit(input$ga_ids,"\n")[[1]]})


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("gene_pancan_dist"), html = waiter::spin_hexdots(), color = "black")

  plot_func <- eventReactive(input$search_bttn, {
    if (length(ga_ids()) >= 1) {
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
        all_genes <- ga_ids()
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

        # 准备与 get_expr_data 返回格式一致的数据
        expr_data <- combat_data[all_genes, , drop = FALSE]
        expr_data <- t(expr_data)
        expr_data <- as.data.frame(expr_data)
        expr_data$ID <- rownames(expr_data)
        expr_data$dataset <- "Integrated data (ComBat)"

        # 合并样本信息（删除expr_data中的dataset列，使用sample_info中的）
        if ("dataset" %in% colnames(sample_info)) {
          sample_info <- sample_info[, !colnames(sample_info) %in% c("dataset")]
        }
        expr_data <- merge(expr_data, sample_info, by = "ID", all.x = TRUE)

        # 确保 subtype 字段有值
        if (!"subtype" %in% colnames(expr_data)) {
          expr_data$subtype <- shared_values$subtypes_combat[1]
        }

        # 使用 viz_TvsN 函数
        nn <- shared_values$subtypes_combat
        if (is.null(nn)) {
          nn <- shared_values$subtypes
        }
        nn <- extract_subset(subtype, nn)

        p <- viz_TvsN(expr_data,df_type = "multi_gene",tumor_subtype = nn,
                       Method =  input$method,
                       Show.P.value = input$Show.P.value,
                       Show.P.label = input$Show.P.label,
                       Show.n = input$Show.n,
                       values = colors())
      } else {
        # 使用原始单数据集
        df <- get_expr_data(datasets = datasets_text,genes = ga_ids())
        # 从树节点名称提取实际的subtypes
        nn <- shared_values$subtypes
        if (!is.null(nn) && length(nn) > 0) {
          # 如果是树节点名称（如"Lung cancer"），提取对应的subtypes
          nn <- extract_subset(subtype, nn)
        }
        p <- viz_TvsN(df,df_type = "multi_gene",tumor_subtype = nn,
                       Method =  input$method,
                       Show.P.value = input$Show.P.value,
                       Show.P.label = input$Show.P.label,
                       Show.n = input$Show.n,
                       values = colors())
      }

      if (input$x_axis45){
        p <- p+ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      }
    }
    return(p)
  })

  output$colorvalues <- reactive({
    c(input$tumor_col, input$normal_col)
  })

  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$gene_pancan_dist <- renderPlot(width = width_scatter,
                                        height = height_scatter,{
    w$show() # Waiter add-ins
                                          plot_func()
                                        })

  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0(input$ga_ids, "_", input$datasets_text, ".csv")
    },
    content = function(file) {
      write.csv(plot_func()$data, file, row.names = FALSE)
    }
  )

  output$download <- downloadHandler(
    filename = function() {
      paste0(input$ga_ids, "_", input$datasets_text, ".pdf")
    },
    content = function(file) {
      p <- plot_func()

        pdf(file, width =  input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    }
  )

  # observeEvent(input$search_bttn, {
  #   if (nchar(input$ga_ids) >= 1) {
  #     shinyjs::show(id = "save_csv")
  #   } else {
  #     shinyjs::hide(id = "save_csv")
  #   }
  # })

  output$tbl <- DT::renderDataTable(server = FALSE, {
    DT::datatable(
      plot_func()$data,
      rownames = T,
      extensions = c("Buttons"),
      options = list(
        pageLength = 5,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "csv", text = "Download table", filename =  paste0(input$ga_ids,"_", input$datasets_text),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        )
      )
    )
  })

}
