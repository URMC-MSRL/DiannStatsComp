#'
#' @description Gets coverage data from a DIA-NN search stats report.
#'
#' @param stats_tsv The stats report file. Needs to be read into R using `readr::read_tsv()` prior. This can be either the First-Pass or MBR report.
#' @param experiment_metadata A data frame containing the metadata for the experiment. This should be made using the `set_metadata()` function also in this package.
#' @param researcher_name The name of the researcher (as found in the raw file name). The should be a string and formatted as follows: "Researcher_" where "Researcher" is whatever is in the raw file. If the raw file has no name, any string can be used. For example, "None" would be appropriate.
#' @param work_order_number The number of the project (as found in the raw file name). The should be a string and formatted as follows: "_WorkOrderNumber" where "WorkOrderNumber" is whatever is in the raw file; for the MSRL, this is typically something like "25-001". If the raw file has no work order number, any string can be used. For example, "None" would be appropriate.
#' @param search The type of search pass within DIA-NN the stats report was generated from. This should be a string and needs to be either "First-Pass" or "MBR". 
#' 
#' @return A data frame containing the number of IDs as reported in the stats report. To get data for both first-pass and MBR searches (if they were performed), you will need to use this function and change the `search` parameter to the appropriate argument ("First-Pass" or "MBR".
#' @export
#'
#' @importFrom magrittr %>%
#' 
get_stats <- function(stats_tsv,
                      experiment_metadata,
                      researcher_name,
                      work_order_number,
                      search) {
  
  researcher <- {{researcher_name}}
  work_order <- {{work_order_number}}
  
  search_type = {{search}}
  
  stats_data_in <- stats_tsv
  n_runs <- stats_data_in %>% 
    dplyr::distinct(File.Name) %>% 
    nrow()
  
  stats_data_out <- stats_data_in %>% 
    dplyr::select(
      c(File.Name,
        FWHM.Scans,
        FWHM.RT,
        Average.Peptide.Length,
        Average.Peptide.Charge,
        Average.Missed.Tryptic.Cleavages)
      ) %>% 
    purrr::set_names(
      'file_name',
      'scans_fwhm',
      'peak_width_fwhm_min',
      'mean_peptide_length',
      'mean_peptide_charge',
      'mean_missed_cleavage_percent'
      ) %>% 
    dplyr::mutate(
      peak_width_fwhm_sec = peak_width_fwhm_min * 60,
      scans_fw = scans_fwhm * 2.5,
      peak_width_fw_sec = peak_width_fwhm_sec * 2.5,
      mean_missed_cleavage_percent = mean_missed_cleavage_percent * 100,
      file_name = stringr::str_remove(file_name,
                                      '\\\\\\\\smdnas.*\\\\MSRL.*\\\\'),
      file_name = stringr::str_remove(file_name,
                                      '\\.raw'),
      file_name = stringr::str_remove(file_name,
                                      researcher),
      file_name = stringr::str_remove(file_name,
                                      work_order)) %>% 
    dplyr::select(-peak_width_fwhm_min) %>% 
    dplyr::relocate(peak_width_fwhm_sec:peak_width_fw_sec,
                    .after = scans_fwhm) %>% 
    dplyr::mutate(diann_pass = search_type,
                  .before = 'scans_fwhm')
  
  run_info <- do.call('rbind',
                      replicate(n_runs,
                                experiment_metadata,
                                simplify = FALSE))
  
  stats_coverage_data <- dplyr::bind_cols(run_info,
                                          stats_data_out)
  
  return(stats_coverage_data)
  
}