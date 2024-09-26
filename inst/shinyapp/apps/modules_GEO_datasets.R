ui.modules_GEO_datasets <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "",
    icon = icon("database"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        br(),
    HTML("<p><strong>Cancer types:</strong></p>"),
    shinyTree(ns("subtype"), theme="proton",stripes=T, themeIcons = F,multiple=T, themeDots = T),
    br(),
        uiOutput(ns("datasets_text")),
        hr()
      ),
      mainPanel(
        width = 9,
        DT::dataTableOutput(ns("geo_table")),

      )
    )
  )
}

server.modules_GEO_datasets <- function(input, output, session, shared_values) {
  ns <- session$ns
  # df <- dataset

  output$subtype <- renderTree({
    subtype
  })
  # observeEvent(input$subtype,{
  #   tree <- input$subtype
  #   req(tree)
  #   nn<- get_selected(tree, format = "classid")[[1]][1]
  #   print(extract_subset(subtype,nn))
  # })
  selected_types <- reactive({
    tree <- input$subtype
    req(tree)
    nn<- sapply(get_selected(tree, format = "classid"), function(x) x[1])
    shared_values$subtypes <- nn

    if (length(nn) != 0){

      dd <- sample_subtype[sample_subtype$subtype %in% (extract_subset(subtype,nn) %>% strsplit(.,"/")%>% unlist()),][,"dataset"] %>% unique()
    }else{
      dd <-  NULL
    }
    shared_values$datasets_text <- dd
    # shared_values$type_select <- dataset_info %>% dplyr::filter(dataset %in% dd) %>% .[,"type"] %>% unique()
    print(dd)
    return(dd)
  })


  output$datasets_text <- renderUI({
     if (!is.null(input$subtype)){
       selectInput(
         inputId = ns("datasets_text"),
         label = "Select dataset:",
         choices = selected_types(),
         multiple = FALSE
       )

       }else{
         selectInput(
           inputId = ns("datasets_text"),
           label = "Select dataset:",
           choices = "",
           multiple = FALSE
         )

     }

  })

  df_select <- reactive({
    if (!is.null(input$subtype)){
      tree <- input$subtype
      req(tree)
      datasets_summary(sapply(get_selected(tree, format = "classid"), function(x) x[1]))
    }else{
      dataset_info
    }
  })

  output$geo_table <- DT::renderDataTable({
    df_select()%>%
      DT::datatable(
        rownames = FALSE,selection = 'single',
        options = list(
          language = list(search = "Filter with keyword:")
        )
      )

  })



  # Keep selected database
  observeEvent(input$geo_table_rows_selected,{
    s <- input$geo_table_rows_selected
    if (length(s)) {
      ddd <- df_select()%>% .[s,]
    }

    if (!is.null(input$subtype)){
      updateSelectInput(session, "datasets_text",
                        label = "Select dataset:",
                        selected = ddd$dataset
      )
    }else{
      updateSelectInput(session, "datasets_text",
                        label = "Select dataset:",
                        choices = dataset_info$dataset,
                        selected = ddd$dataset
      )
      shared_values$datasets_text <- dataset_info$dataset

    }
  })
  observe({
    shared_values$datasets_select <- input$datasets_text

  })


  }
