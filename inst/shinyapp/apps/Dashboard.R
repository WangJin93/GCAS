ui.modules_dash <- function(id) {
  ns <- NS(id)
  tagList(
    bs4Dash::bs4Jumbotron(
      title = "Welcome to the GEO Cancer Analysis Suite (GCAS)!",
      lead = HTML("The GEO Cancer Analysis Suite (GCAS) is a versatile R package designed for analyzing and visualizing gene expression data in cancer research. GCAS allows for the comparison of gene expression between normal and tumor samples, correlation analysis, immune infiltration analysis, differential expression analysis, co-expression analysis, and enrichment analysis. It includes a Shiny app for interactive visualization and can also be used directly within the R environment for advanced scripting. GCAS is ideal for researchers, clinicians, and bioinformaticians seeking to explore cancer genomics data efficiently and effectively.</br><b>Citation: </b>**"),
      status = "info",
      btnName = "Github source code",
      href = "https://github.com/WangJin93/GCAS"
    ),
    fluidRow(
      bs4Dash::column(6,
                      bs4Dash::bs4UserCard(
               title = bs4UserDescription(
                 title = "Application Developer",
                 subtitle = "Jin wang",
                 image = "P1102513.jpg",
                 type = 1
               ),
               status = "info",
               width = 12,

               bs4Dash::bs4ListGroup(
                 width = 12,
                 type = "action",
                 bs4ListGroupItem(
                   h4("Email: Jinwang93@suda.edu.cn")
                 ),
                 bs4ListGroupItem(
                   h4("Affiliation: Soochow University")

                 ),
                 bs4ListGroupItem(
                   h4("Personal website: https://www.jingege.wang"),
                   href = "https://www.jingege.wang"
                 )
               )
             )
      ),
      column(6,
             tags$img(src = 'abstract.png', width = '100%')
             #
             # shinycssloaders::withSpinner(DTOutput(outputId = ns("tf_info"))),
             #   downloadButton(ns("download_tf.csv"), "Download csv table",class = "mybutton")

      )

    )
  )
}
server.modules_dash <- function(input, output, session) {
  ns <- session$ns

  # output$tf_info <- renderDT({
  #
  #   DT::datatable(
  #     tf_list,rownames = F,
  #     options = list(pageLength = 5)
  #   )
  #
  # })
  #
  # output$download_tf.csv <- downloadHandler(
  #   filename = function() {
  #     "TF_list.csv"
  #   },
  #   content = function(file) {
  #     write.csv(tf_list, file, row.names = FALSE)
  #   }
  # )

}

