ui.modules_RRA <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(3,
             tabBox(
               id = ns("venntab"), height = "100%", width = 12,  collapsible = F,

               tabPanel("Upload",
                        h4(".csv file (comma delimited)"),
                        HTML("<hr> <a href='GSE31210_DEGs_results.csv'> <i class='fa fa-download'> </i> Download example data</a>"),
                        br(),
                        fileInput(ns('datafile3'), "Choose CSV File (>1)",
                                  accept = c('comma-separated-values', 'text/plain', '.csv'),
                                  multiple = TRUE
                        ),
                        br(),
                        numericInput(ns("logFC_cutoff"), "Threshold of |log2FC|:", value = 1, min = 0),
                        selectInput(
                          ns("p_cutoff"),
                          "P value cutoff:",
                          choices = c("0.1", "0.05", "0.01", "0.001", "0.0001"),
                          selected = "0.05",
                          multiple = FALSE,
                          selectize = TRUE
                        ),
                        numericInput(ns("top.n"), "Numbers of top DEGs to show: ", value = 20,min = 0),
                        p("Note: If u need show all genes, input 0"),

               ),
               tabPanel("Settings",
                        checkboxInput(
                          ns("show.value"),
                          "Show values:", F),
                        checkboxInput(
                          ns("show.gene"),
                          "Show gene:", TRUE)
               )
             ),
             hr(),
             shinyWidgets::actionBttn(
               inputId = ns("RRA_bttn"),
               label = "Submit!",
               style = "gradient",
               icon = icon("circle"),
               color = "primary",
               block = TRUE,
               size = "md"

             )
      ),
      column(9,align="center",
             tabsetPanel(id = ns("tablist"),
                         tabPanel('RRA results',
                                  shinycssloaders::withSpinner(dataTableOutput(ns("plotdata")))
                         ),
                         tabPanel('Hearmap',
                                  shinycssloaders::withSpinner(plotOutput(ns("heatmap"), height = "auto")),
                                  fluidRow(
                                    column(2),
                                    column(2, numericInput(inputId = ns("height_scatter"), label = "Height", value = 500, max = 600)),
                                    column(2, numericInput(inputId = ns("width_scatter"), label = "Width", value = 800)),
                                    column(2,
                                           br(),
                                           downloadBttn(
                                             outputId = ns("download"),label = "Download Figure",
                                             style = "gradient",
                                             color = "default",
                                             block = TRUE,
                                             size = "sm"
                                           ))
                                  )

                         )
             ),

             tags$head(tags$style(".mybutton{background-color:aliceblue;} .mybutton2{background-color:antiquewhite;} .skin-black .sidebar .mybutton{color: green;}"))
      ) # column
    ) # fluidRow
  ) # fluidPage
}

server.modules_RRA <- function(input, output, session) {
  ns <- session$ns
  RRA_data <- reactive({
    inFiles <- input$datafile3
    if (length(inFiles)<2) {
      return(NULL)
    }else{
      file_extension =  unlist(strsplit(inFiles$datapath, '[.]'))[length(unlist(strsplit(inFiles$datapath, '[.]')))]
      if(!file_extension %in% "csv"){
        print(paste("wrong file format", file_extension))
        return(NULL)

      }
      DEGs_lists<-list()
      for (i in 1:nrow(inFiles)) {
        gse.name<- stringr::str_remove(inFiles[i,1],"_DEGs_results.csv")
        dd<- list(gse.name = read.csv(inFiles[i,4]))
        names(dd) <- gse.name
        DEGs_lists<-c(DEGs_lists,dd)
      }

    }
    return(DEGs_lists)
  })

  RRA <- eventReactive(input$RRA_bttn, {
    if (is.null(RRA_data())) {
      return(NULL)
    }
    results <- RRA_analysis(RRA_data(),
                            top.num = input$top.n,
                            logFC_cut = input$logFC_cutoff,
                            p_cut =  input$p_cutoff %>% as.numeric())
    return(results)

  })



  # Show waiter for plot
  w <- waiter::Waiter$new(id = "heatmapw", html = waiter::spin_hexdots(), color = "black")
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })
  output$heatmap <- renderPlot(width =  width_scatter ,height =  height_scatter ,{
    if (is.null(RRA_data())) {
      NULL
    } else{
      rra_data <- RRA()$RRA_results
      ComplexHeatmap::pheatmap(rra_data,display_numbers=input$show.value,show_rownames=input$show.gene,
                               color = circlize::colorRamp2(seq(-max(abs(min(rra_data)), abs(max(rra_data))),
                                                                max(abs(min(rra_data)), abs(max(rra_data))), length = 3), c("blue", "#EEEEEE", "red")),
                               name = NULL,
                               angle_col = "45")
    }

  })
  output$plotdata <- renderDataTable({
    if (is.null(RRA_data())) {
      NULL
    } else{
      DT::datatable(
        RRA()$RRA_results,
        rownames = T,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  "RRA_results",
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )
    }

  })

  output$download <- downloadHandler(
    filename = function() {
      "RRA_heatmap.pdf"
    },
    content = function(file) {
      rra_data <- RRA()$RRA_results
      p <- ComplexHeatmap::pheatmap(rra_data,display_numbers=input$show.value,show_rownames=input$show.gene,
                               color = circlize::colorRamp2(seq(-max(abs(min(rra_data)), abs(max(rra_data))),
                                                                max(abs(min(rra_data)), abs(max(rra_data))), length = 3), c("blue", "#EEEEEE", "red")),
                               name = NULL,
                               angle_col = "45")
      pdf(file, width =  input$width_scatter/70 ,height = input$height_scatter/70)
      print(p)
      dev.off()

    })
}
