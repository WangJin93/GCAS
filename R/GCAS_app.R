#' Run GCAS Shiny App
#' @import shiny
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' GCAS_app()
#' }
GCAS_app <- function() {
  shiny::shinyAppFile(system.file("shinyapp", "app.R", package = "GCAS"))
}
