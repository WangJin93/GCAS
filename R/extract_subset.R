#' @title Extract Names from List Subsets
#' @description This function searches for specified names within a nested list structure and extracts the names of found subsets.
#' @param lst A list which may contain nested lists.
#' @param names_to_find A character vector of names to search for within the list.
#' @return A character vector of unique names from the found subsets.
#' @export
#' @examples
#' \dontrun{
#'   nested_list <- list(
#'     a = list(
#'       b = 1,
#'       c = list(d = 2)
#'     ),
#'     e = 3
#'   )
#'   names_to_search <- c("b", "d", "e")
#'   result <- extract_subset(nested_list, names_to_search)
#'   print(result)
#' }
extract_subset <- function(lst, names_to_find) {

  # Recursive function to search for subsets
  find_subsets <- function(x, names) {
    result <- list()
    if (is.list(x)) {
      for (name in names) {
        if (name %in% names(x)) {
          result[[name]] <- x[[name]]
        }
      }
      for (elem in x) {
        res <- find_subsets(elem, names)
        if (length(res) > 0) {
          result <- c(result, res)
        }
      }
    }
    return(result)
  }

  # Call the recursive function to perform the search
  subsets <- find_subsets(lst, names_to_find)

  # Extract names of the subsets
  subset_names <- c()
  extract_names <- function(x) {
    if (is.list(x)) {
      for (name in names(x)) {
        subset_names <<- c(subset_names, name)
        extract_names(x[[name]])
      }
    }
  }

  extract_names(subsets)

  return(unique(subset_names))
}
