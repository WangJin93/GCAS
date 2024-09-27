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


  observe({

    updateSelectizeInput(
      session,
      "ga_cor_id",
      choices = genelist,
      selected = "FOXM1",
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "ga_id",
      choices = genelist,
      selected = "GAPDH",
      server = TRUE
    )
  })



  #########correlation
  # Show waiter for plot
  w1 <- waiter::Waiter$new(id = ns("ga_cor_output"), html = waiter::spin_hexdots(), color = "white")
  plot_cor_func <- eventReactive(input$ga_cor_submit, {
    dataset <- input$datasets_text
    id1 <- input$ga_id
    id2 <- input$ga_cor_id
###########################
    nn<- shared_values$subtypes

    results <- cor_cancer_genelist(dataset,
                                    id1,
                                    id2,
                                   extract_subset(subtype,nn),
                                   input$sample_type,
                                    input$cor_method)
    return(results)
  })

corr_func <- eventReactive(input$ga_cor_data_rows_selected , {
    s <- input$ga_cor_data_rows_selected
    if (length(s)) {
      df<-plot_cor_func()$cor_data[c(1,6,s+6)]
    }
    return(df)


  })




output$ga_cor_output <- renderPlot(width = 500,
                                   height = 300,{
                                      w1$show() # Waiter add-ins
                                     if (length(input$ga_cor_data_rows_selected)) {
                                       df <- corr_func()%>% na.omit()
                                       viz_corplot(df,colnames(df)[2],colnames(df)[3],method=input$cor_method,x_lab = " expression",y_lab = " expression")
                                     }else{NULL}
                                    }
)
output$cor_download <- downloadHandler(
  filename = function() {
    paste0(colnames(corr_func())[2],"_",colnames(corr_func())[3],"_", input$datasets_text, ".pdf")
  },
  content = function(file) {
    df <- corr_func()%>% na.omit()

      pdf(file = file, onefile = FALSE,height = input$cor_height,width = input$cor_width)
    print(viz_corplot(df,colnames(df)[2],colnames(df)[3],method=input$cor_method,
                      x_lab = " expression",y_lab = " expression"))
      dev.off()

  }
)

  output$ga_cor_data <- DT::renderDataTable(server = FALSE, {
    if (input$ga_cor_submit>0){
      DT::datatable(
        plot_cor_func()$cor_result,
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
              extend = "csv", text = "Download table", filename = paste0(colnames(corr_func())[2],"_",colnames(corr_func())[3],"_", input$datasets_text ),
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
