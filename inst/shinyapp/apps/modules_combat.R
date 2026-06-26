ui.modules_combat <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(column(3,
                    wellPanel(
                      selectInput(
                        inputId = ns("Type"),
                        label = "Cancer type:",
                        choices = unique(dataset_info %>%  .[,"type"]),
                        selected = "Lung cancer",
                        multiple = F
                      ),
                      shinyTree(ns("subtype"), theme="proton",stripes=T, themeIcons = F,multiple=T, themeDots = T),
                      br(),
                      selectInput(
                        inputId = ns("datasets_text"),
                        label = "Select ≤ 5 datasets:",
                        choices = "",
                        multiple = T
                      ),

                      shinyWidgets::actionBttn(
                        inputId = ns("search_bttn"),
                        label = "Run comBat!",
                        style = "gradient",
                        icon = icon("search"),
                        color = "primary",
                        block = TRUE
                      ))
    ),
    column(9,
           bs4Dash::tabsetPanel(type = "tabs",id=ns("tcga_single"),
                       tabPanel("ComBat result",value = "Results",
                                shinycssloaders::withSpinner(DT::DTOutput(outputId = ns("tbl"))),
                                downloadBttn(ns("downloadData"), "Download combat dataset",
                                             style = "gradient",
                                             color = "default"),
                                br(),
                                tags$hr(),
                                h5("Analyze integrated data:"),
                                shinyWidgets::actionBttn(
                                  inputId = ns("go_to_analysis"),
                                  label = "Go to Datasets overview",
                                  style = "gradient",
                                  icon = icon("arrow-right"),
                                  color = "primary",
                                  block = TRUE
                                )
                       ),
                       tabPanel("DEG results",value = "DEG",
                                fluidRow(
                                  column(3,
                                         shinyWidgets::actionBttn(
                                           inputId = ns("DEG_analyze"),
                                           label = "Analyze!",
                                           style = "gradient",
                                           icon = icon("search"),
                                           color = "primary",
                                           block = TRUE
                                         )
                                  ),
                                  column(6,
                                         downloadBttn(
                                           outputId = ns("download_DEG"),  # 使用inputId而不是outputId
                                           label = "Download results",
                                           style = "gradient",
                                           color = "default",
                                           block = TRUE
                                         ))
                                ),
                                fluidRow(
                                  column(
                                    width = 8,  # 剩余8列用于volcano plot
                                    shinycssloaders::withSpinner(DT::dataTableOutput(ns("DEG_result")))
                                  )
                                )


                       ),
                       tabPanel("Volcano plot",value = "Volcano",
                                fluidRow(
                                  column(3,
                                         numericInput(ns("logFC_cutoff"), "|logFC| threshold:", value = 1,min = 0),
                                         selectInput(
                                           ns("p_cutoff"),
                                           "Adjust P value threshold:",
                                           c("0.05","0.01","0.001","0.0001"),
                                           selected = "0.05",
                                           multiple = FALSE,
                                           selectize = TRUE,
                                           width = NULL,
                                           size = NULL
                                         ),
                                         textInput(ns('values'), "colors (',' seperated)", value = 'green,black,red'),

                                         shinyWidgets::actionBttn(
                                           inputId = ns("plot_volcano"),
                                           label = "Plot volcano!",
                                           style = "gradient",
                                           icon = icon("search"),
                                           color = "primary",
                                           block = TRUE
                                         )
                                  ),
                                  column(9,
                                         shinycssloaders::withSpinner(plotOutput(ns("volcano_plot"), height = "auto")),
                                         fluidRow(
                                           column(3, numericInput(inputId = ns("height_scatter"), label = "Height", value = 500, max = 600)),
                                           column(3, numericInput(inputId = ns("width_scatter"), label = "Width", value = 800)),
                                           column(4,
                                                  br(),
                                                  downloadBttn(
                                                    outputId = ns("download"),
                                                    style = "gradient",
                                                    color = "default",
                                                    block = TRUE
                                                  ))
                                         )

                                  )
                                )
                       )
           )

    )
    )
  )
}

server.modules_combat <- function(input, output, session, shared_values) {
  ns <- session$ns

  # 按需加载sva包
  if (!requireNamespace("sva", quietly = TRUE)) {
    showModal(modalDialog(
      title = "Package Required",
      "The 'sva' package is required for ComBat analysis. Please install it using: ",
      code("BiocManager::install('sva')"),
      easyClose = TRUE,
      footer = NULL
    ))
    return()
  }
  library(sva)

  output$subtype <- renderTree({
    if (input$Type != "Pancancer"){
      subtype[[input$Type]]
    }
  })
  subtype_selected <- reactive({
    tree <- input$subtype
    if (!is.null(tree)){
      req(tree)
      nn<- sapply(get_selected(tree, format = "classid"), function(x) x[1])
      return(nn)
    }
  })
  observeEvent(input$subtype,{

    nn<- subtype_selected()
    ddd <-  dataset_info %>%
      dplyr::filter(type %in% input$Type)%>% .[,"dataset"] %>%  intersect(.,sample_subtype[sample_subtype$subtype %in% c(extract_subset(subtype,nn)),][,"dataset"])
    updateSelectInput(session,"datasets_text",
                      label = "Select ≤ 5 datasets:",
                      choices = ddd
    )

  })
  combat_results <- eventReactive(input$search_bttn, {
    if (length(input$datasets_text) > 4) {
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "You can only select up to 5 datasets."
      ))
      return(NULL)
    }
    df <- combat_datasets(tables = input$datasets_text,
                          tumor_subtype = subtype_selected())

    # 将整合数据存储到 shared_values
    if (!is.null(df)) {
      shared_values$combat_data <- df$combined_data
      shared_values$combat_sample_info <- df$sample_info
      shared_values$combat_datasets <- input$datasets_text
      shared_values$subtypes_combat <- subtype_selected()
    }

    return(df)
  })
  output$tbl <- DT::renderDataTable(server = FALSE, {
    if (!is.null(combat_results())){
      DT::datatable(
        combat_results()$combined_data,
        rownames = T,
        extensions = c("Buttons"),
        options = list(
          pageLength = 10,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  "combat_dataset",
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

    }
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("combat_dataset.csv")
    },
    content = function(file) {
      write.csv(combat_results()$combined_data, file, row.names = FALSE)
    }
  )

  # 跳转到 Datasets overview 页面
  observeEvent(input$go_to_analysis, {
    # 设置一个标志，表示刚完成整合分析，需要自动选择 Integrated data
    shared_values$just_completed_combat <- TRUE
    shinyjs::runjs("$('.nav-sidebar li a[data-value=\"dataset\"]').click()")
  })

  DEG_results <- eventReactive(input$DEG_analyze, {
    # 构建对比矩阵
    if (!is.null(combat_results())){
      sample_info <- combat_results()$sample_info
      combined_data <- combat_results()$combined_data
      design <- model.matrix(~ 0 + sample_info$type)
      colnames(design) <- levels(factor(sample_info$type))
      rownames(design) <- rownames(sample_info)
      fomula <- "Tumor-Normal"
      contrast.matrix<-makeContrasts(contrasts= fomula,
                                     levels = design)
      fit <- lmFit(combined_data,design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      DEG<-topTable(fit2, coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数
      DEG<-tibble::rownames_to_column(DEG,var = "Gene")
      return(DEG)
    }
  })

  output$DEG_result <- DT::renderDataTable(server = FALSE, {
    if (!is.null(DEG_results())){
      DT::datatable(
        DEG_results(),
        rownames = F,
        extensions = c("Buttons"),
        options = list(
          pageLength = 10,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  "combat_DEG_results",
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

    }
  })
  output$download_DEG <- downloadHandler(
    filename = function() {
      paste0("combat_DEG_results.csv")
    },
    content = function(file) {
      write.csv(DEG_results(), file, row.names = FALSE)
    }
  )

  w <- waiter::Waiter$new(id = "wb_plot0", html = waiter::spin_hexdots(), color = "black")
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  plot_vol <- eventReactive(input$plot_volcano,{
    df <-  DEG_results()
    p <- plot_volcano(df,
                      p.cut= input$p_cutoff %>% as.numeric(),
                      logFC.cut= input$logFC_cutoff,
                      colors = unlist(strsplit(input$values, ',')))
    return(p)
  })
  output$volcano_plot <- renderPlot(width = width_scatter,
                               height = height_scatter,{
                                 w$show() # Waiter add-ins
                                plot_vol()
                               })
  output$download <- downloadHandler(
    filename = function() {
      "DEG_Volcano_plot.pdf"
    },
    content = function(file) {
      pdf(file,width = input$width_scatter/70 ,height = input$height_scatter/70)
      plot_vol()
      dev.off()
    })

}
