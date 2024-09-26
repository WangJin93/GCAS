ui.modules_multidata_dist <- function(id) {
  ns <- NS(id)
  fluidPage(
fluidRow(column(3,
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
              label = "Input a protein symbol",
              choices = "GAPDH",
              width = "100%",
              options = list(
                create = TRUE,
                maxOptions = 5,
                plugins = list("restore_on_backspace")
              )
        ),
        selectInput(inputId = ns("method"), label = "Select statistic method", choices = c("wilcox.test", "t.test"), selected = "t.test"),
        materialSwitch(ns("pdist_show_p_value"), "Show P value", inline = TRUE),
        materialSwitch(ns("pdist_show_p_label"), "Show P label", inline = TRUE),
        materialSwitch(ns("show_n"), "Show sample size", inline = FALSE),
        numericInput(inputId = ns("show_n_position"), label = "Show n position", value = -3),
        colourpicker::colourInput(inputId = ns("tumor_col"), "Tumor sample color", "#00AFBB"),
        colourpicker::colourInput(inputId = ns("normal_col"), "Normal sample color", "#FC4E07"),
        tags$hr(style = "border:none; border-top:2px solid #5E81AC;"),
        shinyWidgets::actionBttn(
          inputId = ns("search_bttn"),
          label = "Go!",
          style = "gradient",
          icon = icon("search"),
          color = "primary",
          block = TRUE,
          size = "sm"
        )),
        wellPanel(
          )
      ),
      column(9,
             tabsetPanel(type = "tabs",id=ns("tcga_single"),
                         tabPanel("Box plot",value = "Results",
                                  shinycssloaders::withSpinner(plotOutput(ns("gene_pancan_dist"), height = "auto")),
                                  br(),
                                  fluidRow(
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

                                  tags$br(),
                                  DT::DTOutput(outputId = ns("tbl")),
                         ),
                         tabPanel("Data summary",value = "summary",
                                  DT::dataTableOutput(ns("data_summary"))

                         ),
                         tabPanel("Forest plot",value = "forest",
                                  shinycssloaders::withSpinner(plotOutput(ns("meta_forest"), height = "auto")),
                                  br(),
                                  fluidRow(
                                    column(2, numericInput(inputId = ns("height_forest"), label = "Height", value = 400)),
                                    column(2, numericInput(inputId = ns("width_forest"), label = "Width", value = 800)),
                                    column(2,
                                           br(),
                                           downloadBttn(
                                             outputId = ns("download_forest"),
                                             style = "gradient",
                                             color = "default",
                                             block = TRUE,
                                             size = "sm"
                                           ))
                                  )
                         )
             )

      )
    )
  )
}

server.modules_multidata_dist <- function(input, output, session) {
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
      return(nn)
    #   # 检查 nn 是否为空或长度为零
    #   if (length(nn) == 0) {
    #     nn <- input$Type
    #   } else {
    #     # 使用 any() 检查 nn 是否在 unlist(subtype[[input$Type]]) 中
    #     if (!any(nn %in% unlist(subtype[[input$Type]]))) {
    #       nn <- input$Type
    #     }
    #   }
    # }else{
    #   nn <- input$Type
    # }
    # print(nn)
    }
  })

  observeEvent(input$subtype,{

      nn<- subtype_selected()
      print(nn)
      ddd <-  dataset_info %>%
        dplyr::filter(type %in% input$Type)%>% .[,"dataset"] %>%  intersect(.,sample_subtype[sample_subtype$subtype %in% extract_subset(subtype[[input$Type]] ,nn),][,"dataset"])
      # ddd <-  subtype[[input$Type]] %>% extract(.,nn)
      #   dplyr::filter(type %in% )%>% .[,"dataset"] %>%  intersect(.,sample_subtype[sample_subtype$subtype %in% c(extract_subset(subtype,nn)),][,"dataset"])
      print(ddd)
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

  colors <- reactive({
    c(input$tumor_col, input$normal_col)
  })


  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("gene_pancan_dist"), html = waiter::spin_hexdots(), color = "black")
  expr_data <- eventReactive(input$search_bttn, {
    Pancan_search <- input$Pancan_search
    df <- get_expr_data(datasets=input$datasets_text,
                        genes = Pancan_search)
    return(df)
  })
  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      df <- expr_data()
      nn<- subtype_selected()
      p <- viz_TvsN(df,df_type = "multi_set",
                    tumor_subtype = nn,
        Method =  input$method,
        Show.P.value = input$pdist_show_p_value,
        Show.P.label = input$pdist_show_p_label,
        Show.n = input$show_n,
        Show.n.location =input$show_n_position,
        values = colors()
      )
        p <- p+ ylab(paste0(input$Pancan_search, " expression"))

    }
    return(p)
  })
  data_sum <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      df <- expr_data()
      result <- data_summary(df)
      return(result)
    }
  })
  output$data_summary <- DT::renderDataTable(server = FALSE, {
    DT::datatable(
      data_sum(),
      rownames = T,
      extensions = c("Buttons"),
      options = list(
        pageLength = 5,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,  "_",input$Type,"_summary"),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        )
      )
    )
  })

  output$colorvalues <- reactive({
    c(input$tumor_col, input$normal_col)
  })

  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$gene_pancan_dist <- renderPlot(width = width_scatter,
                                        height = height_scatter,{
    w$show() # Waiter add-ins
    plot_func() + ggplot2::theme(
        plot.margin = margin(t = 0,
                             r = 0,
                             b = 50,
                             l = 50)
      )
  })


  output$download <- downloadHandler(
    filename = function() {
      paste0(input$Pancan_search,  "_GCAS.pdf")
    },
    content = function(file) {
      p<- plot_func() + ggplot2::theme(
        plot.margin = margin(t = 0,
                             r = 0,
                             b = 50,
                             l = 50)
      )
      pdf(file, width =  input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    }
  )

  observeEvent(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      shinyjs::show(id = "save_csv")
    } else {
      shinyjs::hide(id = "save_csv")
    }
  })

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
            extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,  "_",input$Type),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        )
      )
    )
  })
  width_forest<- reactive ({ input$width_forest})
  height_forest<- reactive ({ input$height_forest})
  output$meta_forest <- renderPlot(width = width_forest,
                                   height = height_forest,{
                                     w$show() # Waiter add-ins
                                     if (!is.null(data_sum())){
                                       plot_meta_forest(data_sum())
                                     }
                                   })

}
