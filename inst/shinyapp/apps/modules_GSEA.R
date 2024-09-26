ui.modules_GSEA <- function(id) {
  ns <- NS(id)
  fluidPage(
fluidRow(column(3,
                wellPanel(
                  tabPanel("Upload",value = "upload",
                           h4(".csv/.txt file (tab or comma delimited)")
                           ,
                           HTML("<hr> <a href='volcano.csv'> <i class='fa fa-download'> </i> Download example data</a>  "),
                           br(),
                           shinyWidgets::prettyRadioButtons(
                             inputId = ns("data_type"),
                             label = "Data type",
                             icon = icon("check"),inline = T,
                             choices = c("co-expression"="correlation", "DEG"="limma"),
                             animation = "tada",
                             status = "info"
                           ),

                           fileInput(ns('datafile'), "Choose text/csv File (co-expression/DEG result)",
                                     accept=c('.csv'),multiple = FALSE),

                  ),
                  shinyWidgets::prettyRadioButtons(
                    inputId = ns("gmt"),
                    label = "Select gmt file:",
                    icon = icon("check"),inline = T,
                    choices = c("BP_GMT_7.5.1", "KEGG_GMT_7.5.1", "Upload"),
                    animation = "tada",
                    status = "default"
                  ),
                  conditionalPanel(ns=ns,
                                   condition = "input.gmt == 'Upload'",
                                   fileInput('gmtfile', "Choose .GMT File",
                                             accept=c(".gmt"),multiple = FALSE)
                  ),
                  selectInput(
                    ns("p_cutoff"),
                    "P value cutoff:",
                    c("0.1","0.05","0.01","0.001","0.0001"),
                    selected = "0.05",
                    multiple = FALSE,
                    selectize = TRUE,
                    width = NULL,
                    size = NULL
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
                p("Note: This will take several tens of seconds!"),
      ),
      column(9,
             tabBox(
               id = ns("survtab"), width = 12,

               # DEG results tab
               tabPanel(
                 title = "GSEA results",
                 shinycssloaders::withSpinner(DT::DTOutput(outputId = ns("tbl"))),
                 shinyWidgets::downloadBttn(ns("downloadData"), "Download GSEA results")
               ),

               # Volcano plot tab
               tabPanel(
                 title = "GSEA plot",
                 # 让wellPanel和volcano并排显示
                 fluidRow(
                   column(
                     width = 3,  # 设置适当的宽度，wellPanel占4列
                     wellPanel(
                       selectizeInput(
                         inputId = ns("item"),
                         label = "Input a ID",
                         choices = "",
                         width = "100%",
                         options = list(
                           create = F,
                           maxOptions = 5
                         )
                       ),
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
                       plotOutput(outputId = ns("GSEA_plot"), width = "100%", height =  "auto")
                     )
                   )
                 )
               )
             )
      )
)
  )
}

server.modules_GSEA <- function(input, output, session) {
  ns <- session$ns
  data_input <- eventReactive(input$search_bttn, {
    inFile <- input$datafile
    if (is.null(inFile)){
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "Please upload your .csv formate file first!"
      ))
      return(NULL)
    }
    else{
      expr <-read.csv(inFile$datapath)
      print(head(expr))
    return(expr)
    }
  })
  gmt_input <- eventReactive(input$search_bttn, {
    if (input$gmt == "Upload"){
      inFile <- input$gmtfile
      if (is.null(inFile)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "Please upload your GMT file first!"
        ))
        return(NULL)
      }
      else{
        file_extension =  unlist(strsplit(inFile$datapath, '[.]'))[length(unlist(strsplit(inFile$datapath, '[.]')))]
        if(file_extension %in% c("gmt","GMT")){
          gmt <- read.gmt(inFile$datapath)
          return(gmt)
        } else {
          print(paste("wrong file format", file_extension))
          return(NULL)
        }
      }
    }

  })

  gmt <- reactive({
    if (input$gmt == "BP_GMT_7.5.1") {
      BP_GMT_7.5.1
    } else if (input$gmt == "KEGG_GMT_7.5.1") {
      KEGG_GMT_7.5.1
    } else if (input$gmt == "Upload") {
      gmt_input()
    } else {
      NULL
    }
  })

  GSEA_results <- eventReactive(input$search_bttn, {
    if (!is.null(data_input())){
      results <- GSEA_analysis(data_input(),
                               gmt_file = gmt(),
                               data_type = input$data_type, pvalue_cutoff = as.numeric(input$p_cutoff))
      # results <- viz_DEGs_volcano(expr,p.cut= ,logFC.cut=input$logFC_cutoff,show.labels=input$show.lable)
      return(results)
    }else{
      return(NULL)
    }

  })
  # Show waiter for plot
  observe({
    updateSelectizeInput(
      session,
      "item",
      choices = GSEA_results()@result[["ID"]],
      selected = GSEA_results()@result[["ID"]][1],
      server = TRUE
    )
  })

  w <- waiter::Waiter$new(id = "wb_plot0", html = waiter::spin_hexdots(), color = "black")
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  output$GSEA_plot <- renderPlot(width = width_scatter,
                               height = height_scatter,{
                                 if (!is.null(GSEA_results()) & !is.null(input$item)){
                                   w$show() # Waiter add-ins
                                   enrichplot::gseaplot(GSEA_results(),input$item)
                                 }
  })
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$item,"_GSEA_plot.pdf")
    },
    content = function(file) {
      pdf(file,width = input$width_scatter/70 ,height = input$height_scatter/70)
      enrichplot::gseaplot(GSEA_results(),input$item)
      dev.off()
    })
  output$tbl <- DT::renderDataTable(server = T, {
    if (input$search_bttn){
      DT::datatable(
        GSEA_results()@result[-2],
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
      paste0("GSEA_results.csv")
    },
    content = function(file) {
      write.csv(GSEA_results()@result[-2], file, row.names = T)
    }
  )

}
