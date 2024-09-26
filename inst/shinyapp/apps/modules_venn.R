ui.modules_venn <- function(id) {
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
                        )
               ),
               tabPanel("Settings",
                        sliderInput(ns('marginsize'), 'Select margin size', min = 0, max = 0.5, value = 0.1, step = 0.01),
                        sliderInput(ns('linewidth'), 'Select line width', min = 0, max = 5, value = 1, step = 0.1),
                        radioButtons(ns("linetype"), "Select line type:",
                                     choices = c("solid" = 1, "dashed" = 2, "dotted" = 3, "dot-dashed" = 4, "longdashed" = 5, "blank" = 0),
                                     inline = TRUE
                        )
               ),
               tabPanel("Font & Color",
                        textInput(ns('venncolors'), "Comma separated list of colors", value = "dodgerblue,darkorange1,seagreen3,orchid3,goldenrod1"),
                        a(href = "http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf", "R-colors"),
                        hr(),
                        sliderInput(ns('num.fontsize'), 'Select number font size', min = 0, max = 5, value = 2, step = 0.1),
                        sliderInput(ns('cat.fontsize'), 'Select category font size', min = 0, max = 5, value = 2, step = 0.1),
                        selectInput(ns("cat.font"), "Select Font:", choices = c('mono','Courier', 'sans','Helvetica',
                                                                            'serif','Times','AvantGarde','Bookman', 'Helvetica-Narrow',
                                                                            'NewCenturySchoolbook', 'Palatino', 'URWGothic', 'URWBookman',
                                                                            'NimbusMon', 'URWHelvetica','NimbusSan', 'NimbusSanCond',
                                                                            'CenturySch', 'URWPalladio', 'URWTimes','NimbusRom')),
                        radioButtons(ns("cat.face"), "Select Face:", choices = c("plain","bold","italic","bold.italic"), inline = TRUE)

               )
             ),
             hr(),
             shinyWidgets::actionBttn(
               inputId = ns("venn_bttn"),
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
                         tabPanel('Up',
                                  shinycssloaders::withSpinner(plotOutput(ns("venn_diagram"), height = "600px", width = "600px")),
                                  shinyjs::hidden(
                                    downloadButton(ns("download_up"), "Download Figure", class = "mybutton")
                                  ),
                                  HTML("<hr>"),

                                  shinycssloaders::withSpinner(dataTableOutput(outputId = ns("venn_list")))
                         ),
                         tabPanel('Down',
                                  shinycssloaders::withSpinner(plotOutput(ns("venn_diagram2"), height = "600px", width = "600px")),
                                  shinyjs::hidden(
                                    downloadButton(ns("download_down"), "Download Figure", class = "mybutton")
                                  ),
                                  HTML("<hr>"),
                                  shinycssloaders::withSpinner(dataTableOutput(outputId = ns("venn_list2"))),
                                  shinyjs::hidden(
                                    downloadButton(ns("download_venn_down.csv"), "Download csv table", class = "mybutton")
                                  )
                         )
             ),

             tags$head(tags$style(".mybutton{background-color:aliceblue;} .mybutton2{background-color:antiquewhite;} .skin-black .sidebar .mybutton{color: green;}"))
      ) # column
    ) # fluidRow
  ) # fluidPage
}

server.modules_venn <- function(input, output, session) {
  ns <- session$ns
  Venn_data <- reactive({
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
  observeEvent(input$venn_bttn, {
    if (is.null(Venn_data())) {
      shinyjs::hide(id = "download_up")
      shinyjs::hide(id = "download_down")
    } else {
      shinyjs::show(id = "download_up")
      shinyjs::show(id = "download_down")}
  })

  Venn <- eventReactive(input$venn_bttn, {
    if (is.null(Venn_data())) {
      return(NULL)
    }

    return(get_DEGs_list(Venn_data(),logFC_cut = input$logFC_cutoff,
                         p_cut = input$p_cutoff %>% as.numeric()))

  })
  venn_list <- eventReactive(input$venn_bttn, {
    if (is.null(Venn_data())) {
      return(NULL)
    }
    gs <- Reduce(intersect,Venn()$DEG_up)
    dat<- data.frame(gene = gs)
    for (i in names(Venn_data())) {
      dat<-merge(dat,Venn_data()[[i]] %>% dplyr::filter(gene %in% gs) %>% dplyr::select(c("gene","logFC")),by = "gene")
    }
    colnames(dat)<- c("up_DEGs",names(Venn_data()))
    return(dat)

  })


  plot.color <- reactive({
    venncolors <- input$venncolors
    color.tmp <- gsub(' ', '', venncolors)
    color.tmp <- unlist(strsplit(venncolors, ','))[1:length(Venn()$DEG_up)]
  })

  output$venn_diagram <- renderPlot({
    if (is.null(Venn_data())) {
      NULL
    } else{
      VD <- venn.diagram(Venn()$DEG_up, filename = NULL, fill = plot.color(), cex = input$num.fontsize,
                         margin = input$marginsize, cat.cex = input$cat.fontsize,
                         cat.fontface = input$cat.face, fontface = input$cat.face,
                         cat.fontfamily = input$cat.font, fontfamily = input$cat.font, lwd = input$linewidth,
                         lty = rep(as.numeric(input$linetype), length(Venn()$DEG_up)))
      grid.draw(VD)
    }

  })
  output$venn_list <- renderDataTable({
    if (is.null(Venn_data())) {
      NULL
    } else{
      DT::datatable(
        venn_list(),
        rownames = T,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  "up_DEGs",
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )
    }

  })

  output$download_up <- downloadHandler(
    filename = function() {
      "up_DEGs.pdf"
    },
    content = function(file) {
      pdf(file, width = 12 ,height = 12)
      VD <- venn.diagram(Venn()$DEG_up, filename = NULL, fill = plot.color(), cex = input$num.fontsize,
                         margin = input$marginsize, cat.cex = input$cat.fontsize,
                         cat.fontface = input$cat.face, fontface = input$cat.face,
                         cat.fontfamily = input$cat.font, fontfamily = input$cat.font, lwd = input$linewidth,
                         lty = rep(as.numeric(input$linetype), length(Venn()$DEG_up)))
      grid.draw(VD)
      dev.off()

    })

  venn_list2 <- eventReactive(input$venn_bttn, {
    if (is.null(Venn_data())) {
      return(NULL)
    }
    gs <- Reduce(intersect,Venn()$DEG_down)
    dat<- data.frame(gene = gs)
    for (i in names(Venn_data())) {
      dat<-merge(dat,Venn_data()[[i]] %>% dplyr::filter(gene %in% gs) %>% dplyr::select(c("gene","logFC")),by = "gene")
    }
    colnames(dat)<- c("down_DEGs",names(Venn_data()))
    return(dat)

  })

  output$venn_diagram2 <- renderPlot({
    if (is.null(Venn_data())) {
      NULL
    } else{
      VD <- venn.diagram(Venn()$DEG_down, filename = NULL, fill = plot.color(), cex = input$num.fontsize,
                         margin = input$marginsize, cat.cex = input$cat.fontsize,
                         cat.fontface = input$cat.face, fontface = input$cat.face,
                         cat.fontfamily = input$cat.font, fontfamily = input$cat.font, lwd = input$linewidth,
                         lty = rep(as.numeric(input$linetype), length(Venn()$DEG_down)))
      grid.draw(VD)
    }

  })
  output$download_down <- downloadHandler(
    filename = function() {
      "down_DEGs.pdf"
    },
    content = function(file) {
      pdf(file, width = 12 ,height = 12)
      VD <- venn.diagram(Venn()$DEG_down, filename = NULL, fill = plot.color(), cex = input$num.fontsize,
                         margin = input$marginsize, cat.cex = input$cat.fontsize,
                         cat.fontface = input$cat.face, fontface = input$cat.face,
                         cat.fontfamily = input$cat.font, fontfamily = input$cat.font, lwd = input$linewidth,
                         lty = rep(as.numeric(input$linetype), length(Venn()$DEG_down)))
      grid.draw(VD)
      dev.off()

    })

  output$venn_list2 <- renderDataTable({
    if (is.null(Venn_data())) {
      NULL
    } else{
      DT::datatable(
        venn_list2(),
        rownames = T,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  "down_DEGs",
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
