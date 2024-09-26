ui.modules_dataset_info <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "",
    icon = icon("database"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        br(),
                      HTML("<p><strong>Cancer types:</strong></p>"),
                      textOutput(ns("type_selected")),
                      br(),
                      uiOutput(ns("datasets_text")),
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

                    wellPanel(
                      actionBttn(
                        inputId = ns("ga_stage_submit"),
                        label = "Fetch clinic data!",
                        style = "gradient",
                        icon = icon("check"),
                        color = "default",
                        block = TRUE,
                        size = "sm"
                      ),
                      hr(),
                      selectInput(inputId = ns("feather"), label = "Select clinic feather to plot", choices = "", selected = ""),
                      prettyCheckbox(ns("iscont"),"Convert continuous variables to categorical variables?",F),
                      conditionalPanel(ns=ns,
                                       condition = "input.iscont == true",
                                       numericInput(ns("cont.cut"),"Cut point:",value = 0)
                      ),
                      hr(),
                    ),
                    tags$br(),
                    wellPanel(
                      p("Download figure:"),
                      numericInput(inputId = ns("height_stage"), label = "Height", value = 400),
                      numericInput(inputId = ns("width_stage"), label = "Width", value = 600),
                      hr(),
                      downloadBttn(
                        outputId = ns("download_stage"),
                        style = "gradient",
                        color = "default",
                        block = TRUE,
                        size = "sm"
                      )
                    )
    ),

mainPanel(
             tabBox(
               id = ns("survtab"), width = 12,

               tabPanel( title = "Clinic data",
                         shinycssloaders::withSpinner(DT::dataTableOutput(ns("ga_stage_data")))

               ),
               tabPanel( title = "Stage plot",
                         plotOutput(ns("ga_stage_output"),width = "100%",height = "auto")

               )
             )

      )
  )

)
}

server.modules_dataset_info <- function(input, output, session, shared_values) {
  ns <- session$ns
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

  clinic_merge <- eventReactive(input$ga_stage_submit,{
    data_select <- get_expr_data(input$datasets_text, input$ga_id)
    if (is.null(data_select)) {
      sendSweetAlert(
        session,
        title = "Error",
        text = "Error to query data and plot. Please make sure you have typed corrected gene symbol.",
        type = "error"
      )
      return(NULL)
    }
    dd <- merge_clinic_data(input$datasets_text,data_select)
    return(dd)
  })

  output$ga_stage_data <- DT::renderDataTable(server = FALSE, {
    DT::datatable(
      clinic_merge(),
      rownames = FALSE,
      extensions = c("Buttons"),
      options = list(
        pageLength = 10,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "csv", text = "Download table", filename = paste(input$ga_id, "_clinic_data_",input$datasets_text),
            exportOptions = list(
              modifier = list(page = "all")
            )
          )
        ),
        scrollX = TRUE
      )
    )
  })
  observeEvent(input$ga_stage_submit,{
    updateSelectInput(session, "feather",
                      label = "Select clinic feather to plot",
                      choices =colnames(clinic_merge())[c(-1:-6)],selected = "Gender|gender"
    )
  })

  observeEvent({input$feather
    input$iscont},{
      if ( input$iscont ==T){
        vv <- clinic_merge() %>% .[,input$feather] %>% as.numeric() %>% na.omit()
        print(vv)
        updateNumericInput(session,
                           "cont.cut",
                           value = median(vv))
      }

    })
  width_stage <- reactive ({ input$width_stage })
  height_stage <- reactive ({ input$height_stage })
  w <- waiter::Waiter$new(id = ns("ga_expr_output"), html = waiter::spin_hexdots(), color = "white")

  output$ga_stage_output <- renderPlot(width = width_stage,
                                       height = height_stage,
                                       {
                                         w$show() # Waiter add-ins

                                         stage_plot(clinic_merge(),input$feather,is.cont = input$iscont,cont.cut = input$cont.cut)[["p"]]+
                                           ggplot2::theme(
                                             axis.text.x = element_text(size = 20),
                                             axis.text.y = element_text(size = 20),
                                             axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20)
                                           )+
                                           ggplot2::theme(axis.text.x = element_text(hjust = .5, vjust = .5, size = 20),
                                                          axis.text.y = element_text( size = 20))+
                                           ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))



                                       })
  output$download_stage <- downloadHandler(
    filename = function() {
      paste(input$ga_id, "_",input$feather,"_",input$datasets_text,".pdf")
    },
    content = function(file) {
      p <-   stage_plot(clinic_merge(),input$feather,is.cont = input$iscont,cont.cut = input$cont.cut)[["p"]]+
        ggplot2::theme(
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20)
        )+
        ggplot2::theme(axis.text.x = element_text(hjust = .5, vjust = .5, size = 20),
                       axis.text.y = element_text( size = 20))+
        ylab(ifelse(stringr::str_detect(input$datasets_text,"Phospho"),input$phoso_site,input$ga_id))

      pdf(file,onefile = T, width = input$width_stage/70 ,height = input$height_stage/70)
      print(p)
      dev.off()
    })


}
