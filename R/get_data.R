#' @title Get GEO expression data
#' @description
#' Get the GEO data by using the api. All results saved in MySQL database.
#' @import jsonlite
#' @param table For action = expression, use dataset$Abbre to get all tables; For action = clinic, remove _protein/_mRNA/_Phospho from dataset$Abbre.
#' @param action "expression", "sample_info" or "geo_data".
#' @param gene Gene symbols, you can input one or multiple symbols.
#' @examples
#' \dontrun{
#' results <- get_data(table = "GSE74706", action = "expression", genes = c("GAPDH","TNS1"))
#' }
#' @export
#'
get_data <- function(table = "GSE74706",
                     action = "expression",
                     genes = c("GAPDH","TNS1")){
  if (action == "expression"){
    if (length(genes)>1) genes <- paste0(genes,collapse = ",")
    if (stringr::str_detect(table,"_")){
      url <- paste0("https://www.jingege.wang/bioinformatics/GCAS/api.php?action=",action,"&table=",strsplit(table,"_")[[1]][1],"&genes=",genes)
      res <- jsonlite::fromJSON(url)
      res <- res[c("row_names",sample_subtype[which(sample_subtype$dataset == table),]$ID)]
    }else{
      url <- paste0("https://www.jingege.wang/bioinformatics/GCAS/api.php?action=",action,"&table=",table,"&genes=",genes)
      res <- jsonlite::fromJSON(url)

    }
  }else{
    url <- paste0("https://www.jingege.wang/bioinformatics/GCAS/api.php?action=",action,"&table=",table)
    res <- jsonlite::fromJSON(url)

  }


  return(res)
}
