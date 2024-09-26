#' @title Get GEO expression dataset and sample information in OSF repository
#' @description Retrieve GEO expression datasets and sample information from the OSF repository.
#' @param table A character string specifying the GEO dataset identifier (e.g., "GSE19188").
#' @param action A character string specifying the action, either "geo_data" to retrieve the expression data or "sample_info" to retrieve the sample information.
#' @return A data frame containing the requested data.
#' @import httr
#' @examples
#' \dontrun{
#' df <- get_OSF_data(table = "GSE74706", action = "sample_info")
#' df2 <- get_OSF_data(table = "GSE74706", action = "geo_data")
#' }
#' @export
get_OSF_data <- function(table = "GSE19188", action = "geo_data") {
  # 如果表名包含下划线，则仅使用下划线前的部分
  if (stringr::str_detect(table, "_")) {
    table <- strsplit(table, "_")[[1]][1]
  }

  # 获取文件URL
  file_url <- get_data(table = table, action = action)$link

  # 如果URL为NULL，返回错误信息并退出
  if (is.null(file_url)) {
    cat(table, "was not found in GCAS database!\n")
    return(NULL)
  }

  base_dir <- "data_files"
  action_dir <- file.path(base_dir, action)
  if (!dir.exists(base_dir)) {
    dir.create(base_dir)
  }

  # 如果action目录不存在，创建action目录
  if (!dir.exists(action_dir)) {
    dir.create(action_dir)
  }
  # 设置保存路径
  save_path <- paste0("data_files/", action, "/", table, ".rds")

  # 如果文件已经存在，直接读取并返回
  if (file.exists(save_path)) {
    cat("File is already in", save_path, "\n")
    dd <- readRDS(save_path)
    colnames(dd)[1] <- "ID"
    return(dd)
  } else {
    # 从OSF数据库获取数据
    cat("Fetching data from OSF database.\n")

    # 发送GET请求下载文件
    response <- httr::GET(file_url, httr::progress())

    # 检查响应状态
    if (httr::status_code(response) == 200) {
      cat("Successful.\n")
      # 写入文件
      writeBin(httr::content(response, "raw"), save_path)
      cat("File downloaded successfully to", save_path, "\n")
      dd <- readRDS(save_path)
      colnames(dd)[1] <- "ID"
      return(dd)
    } else {
      cat("Failed to download file. Status code:", httr::status_code(response), "\n")
      return(NULL)
    }
  }
}

