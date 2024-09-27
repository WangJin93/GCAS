ui.modules_gcas_corr <- function(id) {
  ns <- NS(id)
  fluidPage(
        fluidRow(
          column(
            3,
            wellPanel(
              selectInput(
                inputId = ns("Type"),
                label = "Cancer type:",
                choices = c("Pancancer",unique(dataset_info %>%  .["type"])),
                selected = "Pancancer",
                multiple = F
              ),
              conditionalPanel(
                condition = "input.Type == 'Pancancer'",
                prettyRadioButtons(
                  inputId = ns("pancancer_datasets"),
                  label = "Pancancer datasets:",
                  choices = c("Panel 1", "Panel 2", "User defined"),
                  animation = "pulse",inline = T,
                  status = "info"
                ),ns=ns
              ),
              shinyTree(ns("subtype"), theme="proton",stripes=T, themeIcons = F,multiple=T, themeDots = T),
              br(),
              selectInput(
                inputId = ns("datasets_text"),
                label = "Select dataset::",
                choices = "",
                multiple = T
              ),
              hr(),
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

              selectizeInput(
                inputId = ns("geneset"),
                label = "Gene sets:",
                selected = c("Immune checkpoint"),
                choices = c("Immune checkpoint","m6A methylation","Cell senescence","EMT","Cell senescence","Ferroptosis","Cuproptosis","User defined"),
                multiple = F
              ),
              textAreaInput(
                inputId = ns("genelist"), # molecule identifier
                label = "Input gene list:",
                value = "",
                placeholder = "Paste your genelist!\neg.\nTNS1\nTP53\nTXNIP",
                height = "200"
              ),
              hr(),
              wellPanel(
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
                shinyWidgets::actionBttn(
                  inputId = ns("search_bttn"),
                  label = "Go!",
                  style = "gradient",
                  icon = icon("search"),
                  color = "primary",
                  block = TRUE,
                  size = "sm"
                )
              )

            )
          ),
      column(9,
             shinycssloaders::withSpinner(plotOutput(ns("genelist_cor"), height = "auto")),
        hr(),
        fluidRow(
          column(2,
                 br(),
                 shinyWidgets::actionBttn(
                   inputId = ns("show.single"),
                   label = "Individual plot!",
                   style = "gradient",
                   color = "default",
                   icon = icon("sliders"),
                   block = TRUE,
                   size = "sm"
                 )
          ),
          column(2, numericInput(inputId = ns("height_scatter"), label = "Height", value = 500, max = 600)),
          column(2, numericInput(inputId = ns("width_scatter"), label = "Width", value = 800)),
          column(2,
                 br(),
                 downloadBttn(
                   outputId = ns("download"),
                   style = "gradient",
                   color = "default",
                   block = TRUE,
                   size = "sm"
                 ))
        ),
        hr(),

        # h5("NOTEs:"),
        # p("1. ", tags$a("The immune checkpoint molecules (ICMs) were collected from the study conducted by Charoentong et al (43) who reported 24 immunoinhibitory genes and 45 immunostimulatory genes. ")),
        # p("2. ", tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/28052254/", "Charoentong P, Finotello F, Angelova M, et al. Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Rep. 2017;18(1):248-262. doi:10.1016/j.celrep.2016.12.019")),
        DT::DTOutput(outputId = ns("tbl")),
      )
    )
  )
}

server.modules_gcas_corr <- function(input, output, session) {
  ns <- session$ns
  output$subtype <- renderTree({
    if (input$Type != "Pancancer"){
      subtype[[input$Type]]
    }
  })
  observeEvent(c(input$Type,input$pancancer_datasets),{
    if (input$Type != "Pancancer"){
      ddd <-  dataset_info %>%
        dplyr::filter(type %in% input$Type)%>% .[,"dataset"] %>%  intersect(.,sample_subtype[sample_subtype$subtype %in% c(extract_subset(subtype,input$Type)),][,"dataset"])

      updateSelectInput(session,"datasets_text",
                        label = "Select dataset:",
                        choices = ddd,
                        selected = ddd
      )
    }else{
      ddd <- case_when(
        input$pancancer_datasets == "Panel 1"~ c("GSE41258", "GSE68848", "GSE76250", "GSE66229", "GSE31210", "GSE6919", "GSE25097", "GSE71729", "GSE26712", "GSE30784", "GSE13507", "GSE53624", "GSE63514", "GSE46517", "GSE40435", "GSE35570", "GSE17025","GSE50006","GSE12195","GSE99671"),
        input$pancancer_datasets == "Panel 2"~ c("GSE87211", "GSE50161", "GSE54002", "GSE29272", "GSE19188", "GSE70768", "GSE14520", "GSE62452", "GSE66957", "GSE25099", "GSE3167", "GSE23400", "GSE39001", "GSE15605", "GSE53757", "GSE33630", "GSE115810","GSE31048","GSE56315","GSE126209"),
        input$pancancer_datasets == "User defined"~ c("GSE31210")
      )
      updateSelectInput(session,"datasets_text",
                        label = "Select dataset:",
                        choices = dataset_info$dataset,
                        selected = ddd
      )
    }
  })
  subtype_selected <- reactive({
    tree <- input$subtype
    if (!is.null(tree)){
      req(tree)
      nn<- sapply(get_selected(tree, format = "classid"), function(x) x[1])
    }else{
      nn <- subtype[[input$Type]][[1]]
    }
    print(nn)
    nn
  })

  observeEvent(input$subtype,{

    nn<- subtype_selected()
    ddd <-  dataset_info %>%
      dplyr::filter(type %in% input$Type)%>% .[,"dataset"] %>%  intersect(.,sample_subtype[sample_subtype$subtype %in% c(extract_subset(subtype,nn)),][,"dataset"])

    updateSelectInput(session,"datasets_text",
                      label = "Select dataset:",
                      choices = ddd,
                      selected = ddd
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

  observe({
    genelist <- case_when(
      input$geneset == "Immune checkpoint" ~ c("SIRPA\nCTLA4\nTIGIT\nLAG3\nVSIR\nLILRB2\nSIGLEC7\nHAVCR2\nLILRB4\nPDCD1\nBTLA"),
      input$geneset == "m6A methylation" ~ c("METTL3\nMETTL14\nWTAP\nKIAA1429\nZC3H13\nRBM15\nRBM15B\nFTO\nALKBH5\nYTHDC1\nYTHDC2\nYTHDF1\nYTHDF2\nYTHDF3\nHNRNPC\nHNRNPA2B1\nIGF2BP1\nIGF2BP2\nIGF2BP3"),
      input$geneset == "Cell senescence" ~ c("GAPDH\nCDKN1A\nCDKN2A\nLMNB1\nMKI67"),
      input$geneset == "EMT" ~ c("CDH1\nCDH2\nSNAIL\nVIM"),
      input$geneset == "Ferroptosis" ~ c("ACSL4\nLPCAT3\nTFR1\nSLC7A11\nGPX4\nFTH1"),
      input$geneset == "Cuproptosis" ~ c("FDX1\nLIAS\nLIPT1\nDLD\nDLAT\nPDHA1\nPDHB\nMTF1\nGLS\nCDKN2A\nATP7B\nSLC31A1"),
      input$geneset == "User defined" ~ c("")
    )
    updateTextAreaInput(session, "genelist",
                        label = "Input gene list:",
                        value = genelist)
  })

  genelists <- reactive({strsplit(input$genelist,"\n")[[1]]})
  geneset_data <- eventReactive(input$search_bttn, {
    ga_ids <-  genelists()
    # print(ga_ids)
    if (length(input$datasets_text) <2){
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "Please input at least 2 datasets!"
      ))

      return(NULL)
    }
    dd <- get_expr_data(genes  = ga_ids,datasets = input$datasets_text)
    if (is.null(dd)) {
      showModal(modalDialog(
        title = "Message",   easyClose = TRUE,
        "No results were returned, because the current dataset does not contain the gene symbols you entered."
      ))
      return(NULL)

    }else{
      if (length(colnames(dd))-3 < length(ga_ids)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          paste0("The selected datasets only contain the following gene symbols: ", paste(colnames(dd)[6:(ncol(dd))],collapse = " ,") )
        ))}
    }
    return(dd)
  })

  # Show waiter for plot
  w <- waiter::Waiter$new(id = ns("genelist_cor"), html = waiter::spin_hexdots(), color = "white")
  plot_func <- eventReactive(input$search_bttn, {
    if (nchar(input$Pancan_search) >= 1) {
      Pancan_search <- input$Pancan_search
      df <- get_expr_data(genes = Pancan_search,datasets = input$datasets_text)
      if (is.null(df)){
        showModal(modalDialog(
          title = "Message",   easyClose = TRUE,
          "No results were returned, please comfirm your inputed gene symbol."
        ))
        return(NULL)
      }
      nn<- subtype_selected()
      if (is.null(geneset_data())){
        return(NULL)
      }

      p <- cor_gcas_genelist(df,geneset_data(),
                             tumor_subtype = nn,
                             sample_type = input$sample_type,
                             cor_method = input$cor_method)
    }
    return(p)
  })

#
  plot_data <- eventReactive(input$search_bttn, {
      p <- plot_func()

    if (is.null(p)) return(NULL)

    rvalue_T<-as.data.frame(p$r) %>%  rownames_to_column("Gene")
    pvalue_T<-as.data.frame(p$p) %>%  rownames_to_column("Gene")
    pdata<-gather(as.data.frame(rvalue_T),Dataset,corr,-Gene)
    pvalue_T<-gather(as.data.frame(pvalue_T),Dataset,PValue,-Gene)
    pdata <- merge(pdata,pvalue_T,by=c("Dataset","Gene"))
    pdata$pstar <- ifelse(pdata$PValue < 0.05,
                          ifelse(pdata$PValue < 0.001, "***", ifelse(pdata$PValue < 0.01, "**", "*")),
                          ""
    )
    return(pdata)
  })
  width_scatter <- reactive ({ input$width_scatter })
  height_scatter <- reactive ({ input$height_scatter })

  output$genelist_cor <- renderPlot(width = width_scatter,
                                          height = height_scatter,{
    w$show() # Waiter add-ins
                                            if (!is.null(plot_func())){
                                              viz_cor_heatmap(plot_func()$r,plot_func()$p)
                                            }
                                          })

  corr_func <- eventReactive(input$tbl_rows_selected, {
    s <- input$tbl_rows_selected
    if (length(s)) {
      selected <- plot_data()[s,]
      #data<- df[which(df$type == type),]
      checkpoint <- selected$Gene
      cancer <- selected$Dataset
      df <- plot_func()$sss[[cancer]]
      df <- df[,c("ID","subtype","dataset","type",colnames(df)[6],checkpoint)]
      return(df)
    }
    return(df)
  })


#
  observeEvent(eventExpr = input$tbl_rows_selected,{
    s <- input$tbl_rows_selected

    if (length(s)) {
      showModal(
        modalDialog(
          title = paste("Correlation scatter plot"),
          size = "l",
          fluidPage(
            column(
              12,
              plotOutput(ns("scatter_plot"),width = "100%",height = "auto",),
              DT::DTOutput(outputId = ns("scatter_corr")),

            )
          )
        )
      )

      output$scatter_plot <- renderPlot(width = 450,
                                        height = 300,{
                                          df<-corr_func() %>% na.omit()
                                          viz_corplot(df,colnames(df)[5],colnames(df)[6],method=input$cor_method,x_lab= " exppression",y_lab=" exppression")

                                        })
      output$scatter_corr <- DT::renderDataTable(server = TRUE,{
        DT::datatable(
          corr_func(),
          rownames = T,
          extensions = c("Buttons"),
          options = list(
            pageLength = 5,
            dom = "Bfrtip",
            buttons = list(
              list(
                extend = "csv", text = "Download table", filename =  paste0(names(corr_func())[5],"_",names(corr_func())[6],"_",corr_func()[1,3]),
                exportOptions = list(
                  modifier = list(page = "all")
                )
              )
            )
          )
        )

      }
      )
      ## downloadTable
    }
  })

#
#   observeEvent(input$search_bttn, {
#     if (nchar(input$Pancan_search) >= 1) {
#       shinyjs::show(id = "save_csv")
#     } else {
#       shinyjs::hide(id = "save_csv")
#     }
#   })
#
#
  output$tbl <- DT::renderDataTable(server = FALSE, {
    if (!is.null( plot_data())){
      DT::datatable(
        plot_data(),
        rownames = T,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,"_",input$Type,"_cor"),
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
      paste0(input$Pancan_search,"_",input$Type,"_cor",".pdf")
    },
    content = function(file) {
        pdf(file,  width =  input$width_scatter/70 ,height = input$height_scatter/70)
      print(viz_cor_heatmap(plot_func()$r,plot_func()$p))
      dev.off()
    }
  )
  observeEvent(eventExpr = input$show.single,{

    showModal(
      modalDialog(
        title = paste(input$Pancan_search," correlation volcano plot in multi-datasets" ),
        size = "l",
        fluidPage(
          div(style="display:flex",
              selectInput(
                inputId = ns("gene_select"), label = "Choose a specific gene:",
                choices = genelists()
              ),
              numericInput(inputId = ns("xinter"), label = "r cutpoint:", value = 0.2,max = 0.6,width = 200),

              selectInput(
                inputId = ns("yinter"),
                label = "P cutpoint:",
                choices = c("0.05", "0.005","0.001"),
                selected = "0.05",width = 200
              )
              ),

            div(style="display:flex",
                textInput(ns('values'), "colors (',' seperated)", value = 'green,black,red',width = 200),
                numericInput(inputId = ns("height_scatter1"), label = "Height", value = 500,max = 800,width = 200),
                numericInput(inputId = ns("width_scatter1"), label = "Width", value = 600,width = 200)
            ),


            plotOutput(ns("scatter_plot1"),width = "100%",height = "auto",),
            DT::DTOutput(outputId = ns("scatter_corr1")),



        )
      )
    )



    width_scatter1 <- reactive ({ input$width_scatter1 })
    height_scatter1 <- reactive ({ input$height_scatter1 })
    output$scatter_plot1 <- renderPlot(width = width_scatter1,
                                       height = height_scatter1,{
                                         p <- plot_func()
                                         viz_cor_volcano(p,item = input$gene_select,
                                                         r.cut = input$xinter,
                                                         p.cut = as.numeric(input$yinter),
                                                         colors = unlist(strsplit(input$values, ',')))
                                       })

    output$scatter_corr1 <- DT::renderDataTable(server = F, {
      cor_result <- plot_func()
      gene <- input$gene_select
      r.cut = input$xinter
      p.cut = as.numeric(input$yinter)
      results <- data.frame(r = cor_result$r[gene,],
                            p = cor_result$p[gene,])

      # Determine change status
      results$change <- ifelse(results$p < p.cut,
                               ifelse(results$r < -r.cut, "negative",
                                      ifelse(results$r > r.cut, "positive", "no")),
                               "no")
      # Arrange results by r
      results <- dplyr::arrange(results, r)
      DT::datatable(
        results,
        rownames = T,
        selection =  "single",
        extensions = c("Buttons"),
        options = list(
          pageLength = 5,
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "csv", text = "Download table", filename =  paste0(input$Pancan_search,  "_",gene,"_",input$Type),
              exportOptions = list(
                modifier = list(page = "all")
              )
            )
          )
        )
      )

    }
    )

  })

}
