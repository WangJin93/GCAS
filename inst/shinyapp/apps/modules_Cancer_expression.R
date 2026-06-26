ui.modules_Cancer_expression <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "",
    icon = icon("database"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        br(),
    HTML("<p><strong>Cancer types:</strong></p>"),
    textOutput(ns("type_selected")),
    br(),
    uiOutput(ns("datasets_text")),
    hr(),
    hr(),
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
        br()
      ),
      mainPanel(
        width = 9,
        fluidPage(
          fluidRow(
            column(
              4,
              wellPanel(
                h4("Analysis Controls"),
                selectInput(inputId = ns("method"), label = "Select statistic method", choices = c("wilcox.test", "t.test"), selected = "t.test"),

                materialSwitch(ns("pdist_show_p_value"), "Show P value", inline = FALSE),
                materialSwitch(ns("pdist_show_p_label"), "Show P label", inline = FALSE),
                materialSwitch(ns("show_n"), "Show sample size", inline = FALSE),
                colourpicker::colourInput(inputId = ns("tumor_col"), "Normal sample color", "#00AFBB"),
                colourpicker::colourInput(inputId = ns("normal_col"), "Tumor sample color", "#FC4E07"),

                actionBttn(
                  inputId = ns("ga_submit"),
                  label = "Submit",
                  style = "gradient",
                  icon = icon("check"),
                  color = "default",
                  block = TRUE,
                  size = "sm"
                )),
            ),
            column(
              8,
              shinycssloaders::withSpinner(plotOutput(ns("ga_expr_output"),width = "100%",height = "auto")),
              br(),
              fluidRow(
                column(2, numericInput(inputId = ns("height_scatter"), label = "Height", value = 400, max = 600)),
                column(2, numericInput(inputId = ns("width_scatter"), label = "Width", value = 400)),
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
              DT::dataTableOutput(ns("ga_expr_data"))
            )

          ))

      )
    )
  )
}

server.modules_Cancer_expression <- function(input, output, session, shared_values) {
  ns <- session$ns
  # df <- dataset
  observe({

    updateSelectizeInput(
      session,
      "ga_id",
      choices = get_genelist(),
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

  colors <- reactive({
    c(input$tumor_col, input$normal_col)
  })


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("ga_expr_output"), html = waiter::spin_hexdots(), color = "white")
  plot_func <- eventReactive(input$ga_submit, {
    id <- input$ga_id
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
      if (!id %in% rownames(combat_data)) {
        sendSweetAlert(
          session,
          title = "Error",
          text = paste("Gene", id, "not found in the integrated data."),
          type = "error"
        )
        return(NULL)
      }

      # 准备与 get_expr_data 返回格式一致的数据
      expr_data <- combat_data[id, , drop = FALSE]
      expr_data <- t(expr_data)
      expr_data <- as.data.frame(expr_data, stringsAsFactors = FALSE)
      expr_data$ID <- rownames(expr_data)
      print(head(expr_data))
      # 将基因名列改名为value
      # 删除可能存在的重复列（例如gene_id恰好也是"value"的情况）
      expr_data$dataset <- "Integrated data (ComBat)"

      # 合并样本信息（删除expr_data中的dataset列，使用sample_info中的）
      if ("dataset" %in% colnames(sample_info)) {
        sample_info <- sample_info[, !colnames(sample_info) %in% c("dataset")]
      }
      expr_data <- merge(expr_data, sample_info, by = "ID", all.x = TRUE)
      expr_data <- expr_data[c(1,3,2,4:7)]
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
      p <- viz_TvsN(expr_data, df_type = "single", tumor_subtype = nn,
                    Method = input$method,
                    Show.P.value = input$pdist_show_p_value,
                    Show.P.label = input$pdist_show_p_label,
                    Show.n = input$show_n,
                    values = colors()) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = NULL)) +
        ylab(ifelse(stringr::str_detect(input$datasets_text, "Phospho"), input$phoso_site, input$ga_id))

      return(p)

    } else {
      # 使用原始单数据集
      data_select <- get_expr_data(datasets_text, id)
      if (is.null(data_select)) {
        sendSweetAlert(
          session,
          title = "Error",
          text = "Error to query data and plot. Please make sure you have typed corrected gene symbol.",
          type = "error"
        )
        return(NULL)
      }
      # 从树节点名称提取实际的subtypes
      nn <- shared_values$subtypes
      if (!is.null(nn) && length(nn) > 0) {
        # 如果是树节点名称（如"Lung cancer"），提取对应的subtypes
        nn <- extract_subset(subtype, nn)
      }
      p<- viz_TvsN(data_select,df_type = "single",tumor_subtype = nn,
                   Method =  input$method,
                   Show.P.value = input$pdist_show_p_value,
                   Show.P.label = input$pdist_show_p_label,
                   Show.n = input$show_n,
                   values = colors())  +
        ggplot2::guides(fill = ggplot2::guide_legend(title = NULL)) +
        ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))

      return(p)
    }
  })

  # 获取绘图对象
  get_plot <- reactive({
    plot_func()
  })

  # 获取数据对象
  get_plot_data <- reactive({
    result <- plot_func()
    # 从 ggplot 对象中提取数据
    if (!is.null(result) && ggplot2::is.ggplot(result)) {
      return(result$data)
    }
    return(NULL)
  })

  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$ga_expr_output <- renderPlot(width = width_scatter,
                                      height = height_scatter,{
                                        w$show() # Waiter add-ins
                                        p <- get_plot()
                                        if (!is.null(p)) {
                                          print(p + ylab(input$ga_id))
                                        }
                                      }


                                      )
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$ga_id, "_", input$datasets_text, "_expression.pdf")
    },
    content = function(file) {
      pdf(file = file, onefile = FALSE, width = input$width_scatter/70 ,height = input$height_scatter/70)
      p <- get_plot()
      if (!is.null(p)) {
        print(p + ylab(input$ga_id))
      }
      dev.off()
    }
  )
  br()
  br()
  output$ga_expr_data <- DT::renderDataTable(server = FALSE, {
    data <- get_plot_data()
    if (!is.null(data)) {
      DT::datatable(
        data,
        rownames = TRUE,
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table",
              filename = paste(input$ga_id, "in", input$datasets_text),
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
