#' @title Get GEO expression data
#' @description
#' Get the GEO data by using the API. All results saved in MySQL database.
#' @import jsonlite httr
#' @param table For action = expression, use dataset$Abbre to get all tables; For action = clinic, remove _protein/_mRNA/_Phospho from dataset$Abbre.
#' @param action "expression", "sample_info" or "geo_data".
#' @param genes Gene symbols, you can input one or multiple symbols.
#' @param api_url Optional custom API URL. If NULL, uses the default GCAS API URL.
#' @examples
#' \dontrun{
#' results <- get_data(table = "GSE74706", action = "expression", genes = c("GAPDH","TNS1"))
#' }
#' @export
get_data <- function(table = "GSE74706",
                     action = "expression",
                     genes = c("GAPDH", "TNS1"),
                     api_url = NULL) {
  
  # Input validation
  if (!is.character(table) || length(table) != 1) {
    stop("table must be a single character string")
  }
  
  if (!action %in% c("expression", "sample_info", "geo_data")) {
    stop("action must be 'expression', 'sample_info', or 'geo_data'")
  }
  
  # Use custom API URL or default
  base_url <- if (is.null(api_url)) {
    "https://www.jingege.wang/bioinformatics/GCAS/api.php"
  } else {
    api_url
  }
  
  # Build query parameters
  query_params <- list(action = action, table = table)
  
  if (action == "expression" && length(genes) > 0) {
    query_params$genes <- paste0(genes, collapse = ",")
  }
  
  # Handle table names with underscores for expression action
  if (action == "expression" && stringr::str_detect(table, "_")) {
    query_params$table <- strsplit(table, "_")[[1]][1]
  }
  
  # Build URL with query parameters
  url <- httr::modify_url(base_url, query = query_params)
  
  # Make API request
  tryCatch({
    response <- httr::GET(url)
    httr::stop_for_status(response)
    res <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    
    # Post-process results for expression action with underscore tables
    if (action == "expression" && stringr::str_detect(table, "_")) {
      sample_ids <- sample_subtype[sample_subtype$dataset == table, ]$ID
      res <- res[, c("row_names", sample_ids), drop = FALSE]
    }
    
    return(res)
  }, error = function(e) {
    warning(paste("Failed to retrieve data from API:", e$message))
    return(NULL)
  })
}

#' @title Get GCAS API URL
#' @description Get the current GCAS API URL.
#' @return The current API URL as a character string.
#' @export
get_gcas_api_url <- function() {
  "https://www.jingege.wang/bioinformatics/GCAS/api.php"
}

#' @title Set GCAS API URL
#' @description Set a custom GCAS API URL for the current session.
#' @param url The new API URL.
#' @return NULL
#' @export
set_gcas_api_url <- function(url) {
  if (!is.character(url) || length(url) != 1) {
    stop("url must be a single character string")
  }
  options(GCAS_api_url = url)
  invisible()
}