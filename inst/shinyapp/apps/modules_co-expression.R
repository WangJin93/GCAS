ui.modules_co_expression <- function(id) {
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
                  selectizeInput(
                    inputId = ns("Pancan_search"),
                    label = "Input a gene symbol",
                    choices = "GAPDH",
                    width = "100%",
                    options = list(
                      create = TRUE,
                      maxOptions = 5,
                      plugins = list("restore_on_backspace")
                    )
                  ),
                  selectInput(
                    inputId = ns("cor_method"),
                    label = "Select Correlation method",
                    choices = c("spearman", "pearson"),
                    selected = "pearson"
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
                    numericInput(ns("r_cutoff"), "|r| threshold:", value = 0.3,min = 0),
                    selectInput(
                      ns("p_cutoff"),
                      " P value threshold:",
                      c("0.05","0.01","0.001","0.0001"),
                      selected = "0.05",
                      multiple = FALSE,
                      selectize = TRUE,
                      width = NULL,
                      size = NULL
                    ),
                    textInput(ns('values'), "colors (, seperated)", value = 'green,black,red'),

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
                 title = "Co-expression results",
                 shinycssloaders::withSpinner(DT::DTOutput(outputId = ns("tbl"))),
                 shinyWidgets::downloadBttn(ns("downloadData"), "Download co-expression results")
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

server.modules_co_expression <- function(input, output, session) {
  ns <- session$ns
  observeEvent(input$datasets_text,{
    groups <- sample_subtype %>% dplyr::filter(dataset == input$datasets_text) %>%
      .[,"subtype"] %>% unique()
    updateSelectizeInput(
      session,
      "subtype",
      choices = groups,
      selected = groups
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


  corr_results <- eventReactive(input$search_bttn, {
    #创建数据库连接
    showNotification("Downloading data .... ",duration = 3)
    expr <- get_OSF_data(table = input$datasets_text, action = "geo_data")
    showNotification("Download completed, analyzing DEGs....",duration = 2)
    sample_info <- sample_subtype[sample_subtype$ID %in% colnames(expr), ]
    sample_info <- sample_info %>% dplyr::filter(subtype %in% input$subtype)
    expr <- expr[, c("ID", sample_info$ID)]
    results <- coexpression_analysis(expr,
                                     gene = input$Pancan_search,
                                     method = input$cor_method)

    return(results)

  })
  # Show waiter for plot

  w <- waiter::Waiter$new(id = "wb_plot0", html = waiter::spin_hexdots(), color = "black")
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  plot_vol <- eventReactive(input$plot,{
    results <-  corr_results()
    results$r<- as.numeric(results$r)
    results$p<- as.numeric(results$p)
    results$correlation<-ifelse(results$p< as.numeric(input$p_cutoff) ,ifelse(results$r < -input$r_cutoff,"Negative",ifelse(results$r > input$r_cutoff,"Positive","No")),"No")
    print(table(results$correlation))
    colors <- unlist(strsplit(input$values, ','))
    color_map <- c(
      "Negative" = colors[1],
      "No" = colors[2],
      "Positive" = colors[3]
    )

    # 为 unique_values 分配颜色
    colors <- color_map[unique(results$correlation)]

    results<-arrange(results,r)
    print(results)
      top5 <-c(results[c(1:10),1], results[c((nrow(results)-9):nrow(results)),1])
      results$label <- ifelse(results$gene %in% top5,results$gene,"")

    p <-  ggplot(data = results) +
      aes(x=r %>% as.numeric(),y=-log10(p %>% as.numeric())) +
      # geom_point(alpha = input$alphaInput, size = input$pointSize, shape = 16) +
      geom_point(size = 2, shape = 16,alpha=0.8) +


      # This needs to go here (before annotations)
      theme_light(base_size = 20)

    p <- p + aes(color=correlation) + scale_color_manual(values=colors)

    p <- p + geom_hline(yintercept = c(-log10(as.numeric(input$p_cutoff))), linetype="dashed", color="grey30",size=0.8)+
      geom_vline(xintercept = c(input$r_cutoff,-input$r_cutoff), linetype="dashed", color="grey30",size=0.8)

    p <- p+ theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())+
      labs(x=bquote('Correlation coefficient'), y=bquote(''*-Log[10]*' (P value)'))

    ########## User defined labeling
      p <-  p + geom_text_repel(
        data = results,
        aes(label = label),
        size = 5,
        color="black",
        nudge_x = 0.2,
        nudge_y=0.2,max.overlaps=100000,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.3+5*0.1, "lines"),show.legend=F
      )



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
      p
        pdf(file,width = input$width_scatter/70 ,height = input$height_scatter/70)
        print(p)
        dev.off()
    })
  output$tbl <- DT::renderDataTable(server = F, {
    if (input$search_bttn){
      DT::datatable(
        corr_results(),
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
      paste0(input$datasets_text,"_",input$Pancan_search,"_coexpression.csv")
    },
    content = function(file) {
      write.csv(corr_results(), file, row.names = FALSE)
    }
  )

}
