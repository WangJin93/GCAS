#' @title Retrieve Gene Expression Data
#' @description Retrieve expression data for specified genes from given datasets.
#' @param datasets A character vector of dataset identifiers.
#' @param genes A character vector of gene identifiers.
#' @return A dataframe containing expression data for the specified genes from the given datasets.
#' @examples
#' \dontrun{
#' results <- get_expr_data(datasets = "GSE74706", genes = c("GAPDH","TNS1"))
#' results <- get_expr_data(datasets = c("GSE62113","GSE74706"), genes = "GAPDH")
#' results <- get_expr_data(datasets = c("GSE62113","GSE74706"), genes = c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA"))
#' }
#' @export
get_expr_data <- function(datasets, genes) {
  # 创建缓存目录（如果不存在）
  base_dir <- "data_files"
  action_dir <- file.path(base_dir, "data_temp")
  if (!dir.exists(base_dir)) {
    dir.create(base_dir)
  }
  if (!dir.exists(action_dir)) {
    dir.create(action_dir)
  }

  # 检查基因向量是否为空
  if (length(genes) == 0) {
    return(NULL)
  } else {
    message("Querying data of identifier ", paste0(genes, collapse = ", "), " from datasets ", paste0(datasets, collapse = ", "))

    # 创建一个空的数据框
    expr_data <- data.frame()

    for (x in datasets) {
      # 生成基因的哈希值作为文件名的一部分
      genes_hash <- digest::digest(genes, algo = "md5")
      cache_file <- file.path(action_dir, paste0(x, "_", genes_hash, ".RData"))

      # 检查缓存文件是否存在
      if (file.exists(cache_file)) {
        message("Loading cached data from ", cache_file)
        load(cache_file)
        expr_data <- plyr::rbind.fill(expr_data, cached_data)
        next()
      }

      data <- get_data(x, "expression", genes)

      # 检查数据检索是否成功
      if (is.null(nrow(data))) {
        message(paste0(genes, collapse = ", "), " retrieve no results in ", x)
        next()
      } else {
        row.names(data) <- NULL
        cached_data <- data %>% tibble::column_to_rownames(var = "row_names") %>% t() %>% as.data.frame()
        cached_data <- tibble::rownames_to_column(cached_data, var = "ID")
        cached_data$dataset <- x
        # 保存到缓存
        save(cached_data, file = cache_file)
        expr_data <- plyr::rbind.fill(expr_data, cached_data)
      }
    }

    # 检查是否检索到任何数据
    if (nrow(expr_data) == 0) {
      message("Retrieve no data.")
      return(NULL)
    }

    # 将基因表达列转换为数值
    if (length(genes) == 1) {
      expr_data[, 2] <- as.numeric(expr_data[, 2])
    } else {
      expr_data[,2:(ncol(expr_data)-1)] <- apply(expr_data[,2:(ncol(expr_data)-1)], 2, as.numeric)
    }

    # 子集并合并样本亚型信息
    expr_data <- expr_data[,c("ID",colnames(expr_data)[2:(ncol(expr_data)-1)], "dataset")]
    expr_data <- merge(expr_data, sample_subtype, by = c("ID", "dataset"), all.x = TRUE)

    return(expr_data)
  }
}
