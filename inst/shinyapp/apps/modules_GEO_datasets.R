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
    selectInput(
      inputId = ns("cancer_type_select"),
      label = NULL,
      choices = c("Select from tree", "Integrated data"),
      selected = "Select from tree",
      multiple = FALSE
    ),
    shinyTree(ns("subtype"), theme="proton",stripes=T, themeIcons = F,multiple=T, themeDots = T),
    br(),
        uiOutput(ns("datasets_text")),
        hr()
      ),
      mainPanel(
        width = 9,
        bs4Dash::tabsetPanel(type = "tabs",id=ns("tcga_single"),
                    tabPanel("Dataset info",value = "info",
                             DT::dataTableOutput(ns("geo_table"))
                    ),
                    tabPanel("Abbreviate",value = "abbre",
                             selectInput(
                               inputId = ns("Type"),
                               label = "Cancer type:",
                               choices = unique(dataset_info %>%  .["type"]),
                               selected = "Lung cancer",
                               multiple = F
                             ),

                             DT::dataTableOutput(ns("abbreviate"))

                    )
        )

      )
    )
  )
}

server.modules_GEO_datasets <- function(input, output, session, shared_values) {
  ns <- session$ns
  # df <- dataset
  
  # 检查是否刚完成整合分析，如果是，自动选择 Integrated data
  observe({
    if (!is.null(shared_values$just_completed_combat) && shared_values$just_completed_combat) {
      # 自动选择 Integrated data
      updateSelectInput(session, "cancer_type_select", selected = "Integrated data")
      # 清除标志
      shared_values$just_completed_combat <- FALSE
    }
  })
  
  # 控制 shinyTree 的显示
  observe({
    if (input$cancer_type_select == "Integrated data") {
      shinyjs::hide(ns("subtype"))
    } else {
      shinyjs::show(ns("subtype"))
    }
  })
  
  # 监听 shinyTree 的选择变化，如果用户选择了tree中的项目，自动切换回 "Select from tree"
  observeEvent(input$subtype, {
    tree <- input$subtype
    if (!is.null(tree)) {
      selected <- get_selected(tree, format = "classid")
      # 检查是否有选中的项目
      has_selection <- length(selected) > 0 && any(sapply(selected, function(x) length(x) > 0))
      if (has_selection) {
        # 自动切换回 "Select from tree"
        updateSelectInput(session, "cancer_type_select", selected = "Select from tree")
      }
    }
  })

  output$subtype <- renderTree({
    # 只在用户选择"Select from tree"时才加载subtype数据
    req(input$cancer_type_select == "Select from tree")
    subtype
  })
  # observeEvent(input$subtype,{
  #   tree <- input$subtype
  #   req(tree)
  #   nn<- get_selected(tree, format = "classid")[[1]][1]
  #   print(extract_subset(subtype,nn))
  # })
  selected_types <- reactive({
    # 只在用户选择"Select from tree"且有树选择时才执行
    req(input$cancer_type_select == "Select from tree")

    # 如果选择了 Integrated data，返回 NULL，因为我们会在 output$datasets_text 中单独处理
    if (input$cancer_type_select == "Integrated data") {
      # 设置 subtypes 为整合数据的 subtypes
      if (!is.null(shared_values$subtypes_combat)) {
        shared_values$subtypes <- shared_values$subtypes_combat
      } else {
        shared_values$subtypes <- "Integrated"
      }
      return(NULL)
    }

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
     # 如果选择了 Integrated data
     if (input$cancer_type_select == "Integrated data") {
       has_combat <- !is.null(shared_values$combat_data)
       
       if (has_combat) {
         dataset_choices <- c("Integrated data (ComBat)" = "Integrated data (ComBat)")
         return(selectInput(
           inputId = ns("datasets_text"),
           label = "Select dataset:",
           choices = dataset_choices,
           selected = "Integrated data (ComBat)",
           multiple = FALSE
         ))
       } else {
         return(HTML("<p><strong>Please run ComBat analysis first.</strong></p>"))
       }
     }
     
     # 检查是否有整合数据
     has_combat <- !is.null(shared_values$combat_data)
     
     # 构建数据集选项
     dataset_choices <- c()
     
     if (!is.null(input$subtype)){
       dataset_choices <- selected_types()
     }
     
     # 如果有整合数据，添加到选项中
     if (has_combat) {
       dataset_choices <- c("Integrated data (ComBat)" = "Integrated data (ComBat)", dataset_choices)
     }
     
     selectInput(
       inputId = ns("datasets_text"),
       label = "Select dataset:",
       choices = dataset_choices,
       multiple = FALSE
     )
  })

  df_select <- reactive({
    # 如果选择了 Integrated data
    if (input$cancer_type_select == "Integrated data") {
      has_combat <- !is.null(shared_values$combat_data)
      
      if (has_combat) {
        # 创建一个临时的数据集信息表格，显示整合数据的信息
        combat_datasets <- shared_values$combat_datasets
        combat_subtypes <- shared_values$subtypes_combat
        
        # 计算样本数量
        sample_info <- shared_values$combat_sample_info
        tumor_count <- sum(sample_info$type == "Tumor")
        normal_count <- sum(sample_info$type == "Normal")
        
        return(data.frame(
          type = "Integrated data",
          dataset = "Integrated data (ComBat)",
          Normal_Adjacent = normal_count,
          Tumor = tumor_count,
          Premalignant = 0,
          Paired = "F",
          stringsAsFactors = FALSE
        ))
      } else {
        return(dataset_info)
      }
    }
    
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
    
    # 如果选择了 Integrated data
    if (input$cancer_type_select == "Integrated data") {
      if (!is.null(ddd)) {
        updateSelectInput(session, "datasets_text",
                          label = "Select dataset:",
                          selected = ddd$dataset
        )
      }
      return()
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
  output$abbreviate <- DT::renderDataTable({
    abbr_full %>%
      dplyr::filter(Type ==  input$Type) %>%
      DT::datatable(
        rownames = FALSE,selection = 'single'
      )

  })


  }
