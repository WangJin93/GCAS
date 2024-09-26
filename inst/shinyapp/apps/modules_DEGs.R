ui.modules_DEG <- function(id) {
  ns <- NS(id)
  fluidPage(
fluidRow(column(3,
                wellPanel(
                  selectInput(
                    inputId = ns("datasets_text"),
                    label = "Select dataset:",
                    choices = dataset_info$dataset,
                    selected = "GSE31210",
                    multiple = F
                  ),
                  selectInput(
                    inputId = ns("subtype"),
                    label = "Subtype:",
                    choices = "",
                    multiple = T
                  ),
                  shinyWidgets::actionBttn(
                    inputId = ns("search_bttn"),
                    label = "Analyze!",
                    style = "gradient",
                    icon = icon("search"),
                    color = "primary",
                    block = TRUE,
                    size = "sm"
                  )
                  ),
                hr(),
                  wellPanel(
                    h5("Plotting parameter:"),
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
                    selectizeInput(
                      ns("show.lable"),
                      "Symbols to show:",
                      choices = NULL,
                      multiple = TRUE,
                      options = list(
                        create = TRUE,
                        placeholder = "Enter a Symbol symbol, e.g. GAPDH",
                        plugins = list("restore_on_backspace")
                      )
                    ),

                    selectInput(
                      ns("vol.theme"),
                      "Theme:",
                      c("Default","classic","minimal"),
                      selected = "Default",
                      multiple = FALSE,
                      selectize = TRUE,
                      width = NULL,
                      size = NULL
                    ),
                    shinyWidgets::actionBttn(
                      inputId = ns("plot"),
                      label = "Plot!",
                      style = "gradient",
                      icon = icon("search"),
                      block = TRUE,
                      size = "sm"
                    ),
                    tags$hr(style = "border:none; border-top:2px solid #5E81AC;")
                  )
      ),
      column(9,
             tabBox(
               id = ns("survtab"), width = 12,

               # DEG results tab
               tabPanel(
                 title = "DEG results",
                 shinycssloaders::withSpinner(DT::DTOutput(outputId = ns("tbl"))),
                 shinyWidgets::downloadBttn(ns("downloadData"), "Download DEG results")

               ),

               # Volcano plot tab
               tabPanel(
                 title = "Volcano plot",
                 # 让wellPanel和volcano并排显示
                 fluidRow(
                   column(
                     width = 4,  # 设置适当的宽度，wellPanel占4列
                     wellPanel(
                       numericInput(inputId = ns("height_scatter"), label = "Height", value = 600),
                       numericInput(inputId = ns("width_scatter"), label = "Width", value = 600),
                       downloadBttn(
                         outputId = ns("download"),  # 使用inputId而不是outputId
                         style = "gradient",
                         color = "default",
                         block = TRUE,
                         size = "sm"
                       )
                     )
                   ),
                   column(
                     width = 8,  # 剩余8列用于volcano plot
                     shinycssloaders::withSpinner(
                       plotOutput(outputId = ns("volcano"), width = "100%", height =  "auto")
                     )
                   )
                 )
               )
             )
      )
)
  )
}

server.modules_DEG <- function(input, output, session) {
  ns <- session$ns
  observe({

    updateSelectizeInput(
      session,
      "show.lable",
      choices = genelist,
      selected = "GAPDH",
      server = TRUE
    )
  })
  observeEvent(input$datasets_text,{
    groups <- sample_subtype %>% dplyr::filter(dataset == input$datasets_text) %>%
      .[,"subtype"] %>% unique()
    updateSelectizeInput(
      session,
      "subtype",
      choices = c(groups[which(groups != "Normal")]),
      selected = c(groups[which(groups != "Normal")])
    )
  })



  DEG_results <- eventReactive(input$search_bttn, {
    #创建数据库连接
    showNotification("Downloading data .... ",duration = 3)
    expr <- get_OSF_data(table = input$datasets_text, action = "geo_data")
    showNotification("Download completed, analyzing DEGs....",duration = 2)
    results <- DEGs_analysis(expr,tumor_subtype = input$subtype)
    # results <- viz_DEGs_volcano(expr,p.cut=as.numeric(input$p_cutoff) ,logFC.cut=input$logFC_cutoff,show.labels=input$show.lable)

    return(results)

  })
  # Show waiter for plot

  w <- waiter::Waiter$new(id = "wb_plot0", html = waiter::spin_hexdots(), color = "black")
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  plot_vol <- eventReactive(input$plot,{
    df <-  DEG_results()
    p <- plot_volcano(df,
                      p.cut= input$p_cutoff %>% as.numeric(),
                      logFC.cut= input$logFC_cutoff,
                      show.top=FALSE,
                      show.labels=input$show.lable,
                      colors = unlist(strsplit(input$values, ',')))
    return(p)
  })
  output$volcano <- renderPlot(width = width_scatter,
                               height = height_scatter,{
    w$show() # Waiter add-ins
      p <- plot_vol()
      #remove legend (if selected)
      if (input$vol.theme=="classic"){
        p<-p+theme_classic()
      } else if(input$vol.theme=="minimal"){
        p<-p+theme_minimal()
      }
      p
  })
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$datasets_text,"_",input$method,"_Volcano_plot.pdf")
    },
    content = function(file) {
      p <- plot_vol()
      #remove legend (if selected)
      if (input$vol.theme=="classic"){
        p<-p+theme_classic()
      } else if(input$vol.theme=="minimal"){
        p<-p+theme_minimal()
      }
        pdf(file,width = input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    })
  output$tbl <- DT::renderDataTable(server = T, {
    if (input$search_bttn){
      DT::datatable(
        DEG_results(),
        rownames = FALSE,
        options = list(
          pageLength = 5,
          dom = "Bfrtip"
        )
      )

    }
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$datasets_text, "_DEGs_results.csv")
    },
    content = function(file) {
      write.csv(DEG_results(), file, row.names = FALSE)
    }
  )
}
