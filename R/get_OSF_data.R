#' @title Get GEO expression dataset and sample information in OSF repository
#' @description Retrieve GEO expression datasets and sample information from the OSF repository.
#' @param table A character string specifying the GEO dataset identifier (e.g., "GSE19188").
#' @param action A character string specifying the action, either "geo_data" to retrieve the expression data or "sample_info" to retrieve the sample information.
#' @return A data frame containing the requested data.
#' @import httr stringr
#' @examples
#' \dontrun{
#' df <- get_OSF_data(table = "GSE74706", action = "sample_info")
#' df2 <- get_OSF_data(table = "GSE74706", action = "geo_data")
#' }
#' @export
get_OSF_data <- function(table = "GSE19188", action = "geo_data", retries = 2, timeout = 60) {
  # If the table name contains an underscore, use only the part before the underscore
  if (str_detect(table, "_")) {
    table <- strsplit(table, "_")[[1]][1]
  }

  # Get the file URL
  file_url <- get_data(table = table, action = action)$link

  # If the URL is NULL, return an error message and exit
  if (is.null(file_url)) {
    cat(table, "was not found in GCAS database!\n")
    return(NULL)
  }

  # Create necessary directories
  base_dir <- "data_files"
  action_dir <- file.path(base_dir, action)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir)
  }
  if (!dir.exists(action_dir)) {
    dir.create(action_dir)
  }

  # Set the save path
  save_path <- file.path(action_dir, paste0(table, ".rds"))

  # If the file already exists, read and return it
  if (file.exists(save_path)) {
    cat("File is already in", save_path, "\n")
    dd <- readRDS(save_path)
    colnames(dd)[1] <- "ID"
    return(dd)
  } else {
    # Fetch data from the OSF database
    cat("Fetching data from OSF database.\n")

    attempt <- 1
    while (attempt <= retries) {
      response <- tryCatch(
        {
          GET(file_url, timeout(timeout), progress())
        },
        error = function(e) {
          cat("Error occurred on attempt", attempt, ":", e$message, "\n")
          return(NULL)
        }
      )

      # If the response is NULL, it means an error occurred
      if (is.null(response)) {
        if (attempt == retries) {
          stop("Max retries reached. Could not download the file. Please try again later.")
        }
      } else {
        # Check response status
        if (status_code(response) == 200) {
          cat("Download successful on attempt", attempt, ".\n")
          # Write the file
          writeBin(content(response, "raw"), save_path)
          cat("File downloaded successfully to", save_path, "\n")
          dd <- readRDS(save_path)
          colnames(dd)[1] <- "ID"
          return(dd)
        } else {
          cat("Failed to download file on attempt", attempt, ". Status code:", status_code(response), "\n")
          if (attempt == retries) {
            stop("Max retries reached. Could not download the file. Please try again later.")
          }
        }
      }
      attempt <- attempt + 1
    }
  }
}
