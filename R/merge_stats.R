#'
#' @description Merges data from the `get_stats()` function if it was used on both `First-Pass` and `MBR` stats reports.
#'
#' @param stats_tsv The stats report file. Needs to be read into R using "readr::read_tsv()" prior. This can be either the first-pass or mbr report.
#' @param firstPass_statsData A data frame containing the results from `get_stats()` where `search` was set to "First-Pass".
#' @param mbr_statsData A data frame containing the results from `get_stats()` where `search` was set to "MBR"
#' 
#' @return A data frame containing the number of IDs as reported in the stats report from both First-Pass and MBR searches. 
#' @export
#'
#' @importFrom magrittr %>%
#' 
merge_stats <- function(firstPass_statsData,
                        mbr_statsData) {
  
  first_pass <- firstPass_statsData
  mbr <- mbr_statsData
  
  combined_stats <- dplyr::bind_rows(first_pass,
                                     mbr)
  
  return(combined_stats)
  
}