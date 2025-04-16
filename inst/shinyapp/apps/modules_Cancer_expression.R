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
                colourpicker::colourInput(inputId = ns("normal_col"), "Normal sample color", "#00AFBB"),
                colourpicker::colourInput(inputId = ns("tumor_col"), "Tumor sample color", "#FC4E07"),

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
      choices = genelist,
      selected = "GAPDH",
      server = TRUE
    )
  })


  output$type_selected  <- renderText({
    shared_values$subtypes
  })

  output$datasets_text <- renderUI({
       selectInput(
         inputId = ns("datasets_text"),
         label = "Select dataset:",
         choices = shared_values$datasets_text,
         selected = shared_values$datasets_select,
         multiple = FALSE
       )
  })




  colors <- reactive({
    c( input$normal_col,input$tumor_col)
  })



  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("ga_expr_output"), html = waiter::spin_hexdots(), color = "white")
  plot_func <- eventReactive(input$ga_submit, {
    datasets_text <- input$datasets_text
    id <- input$ga_id
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
    nn<- shared_values$subtypes
    p<- viz_TvsN(data_select,df_type = "single",tumor_subtype = nn,
                   Method =  input$method,
                   Show.P.value = input$pdist_show_p_value,
                   Show.P.label = input$pdist_show_p_label,
                   Show.n = input$show_n,
                   values = colors())  +
      ggplot2::guides(fill = ggplot2::guide_legend(title = NULL)) +
      ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))

    return(p)
  })
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$ga_expr_output <- renderPlot(width = width_scatter,
                                      height = height_scatter,{
                                        w$show() # Waiter add-ins
                                          plot_func()+
                                            ylab(input$ga_id)
                                        }


                                      )
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$ga_id, "_",input$datasets_text,".pdf")
    },
    content = function(file) {
        pdf(file = file, onefile = FALSE, width = input$width_scatter/70 ,height = input$height_scatter/70)
        print(plot_func()+
                ylab(input$ga_id)
        )
        dev.off()

      # ggplot2::ggsave(
      #   filename = file, plot = print(p_surv(), newpage = F), device = input$device,
      #   units = "cm", width = input$width, height = input$height, dpi = 600
      # )
    }
  )
  br()
  br()
  output$ga_expr_data <- DT::renderDataTable(server = FALSE, {
      DT::datatable(
        plot_func()$data,
        rownames = T,
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename = paste(input$ga_id, "in",input$datasets_text),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )
  })


  }
