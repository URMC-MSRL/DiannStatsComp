#' Get data from the log of a DIA-NN search
#'
#' @description Gets coverage data from a DIA-NN search log where oxidation at M was set as the only variable modification.
#'
#' @param log_data The log file. Needs to be read into R using `readr::read_tsv()` prior.
#' @param experiment_metadata A data frame containing the metadata for the experiment. This should be made using the `set_metadata()` function also in this package.
#' @param researcher_name The name of the researcher (as found in the raw file name). The should be a string and formatted as follows: "Researcher_" where "Researcher" is whatever is in the raw file. If the raw file has no name, any string can be used. For example, "None" would be appropriate.
#' @param work_order_number The number of the project (as found in the raw file name). The should be a string and formatted as follows: "_WorkOrderNumber" where "WorkOrderNumber" is whatever is in the raw file; for the MSRL, this is typically something like "25-001". If the raw file has no work order number, any string can be used. For example, "None" would be appropriate.
#' 
#' @return A data frame containing the number of IDs as reported in the log (# precursor IDs, # protein IDs, number of precursors in the spectral library) for both First-Pass and MBR searches if they were performed.
#' @export
#'
#' @importFrom magrittr %>%

get_global_log <- function(log_data,
                           experiment_metadata,
                           researcher_name,
                           work_order_number) {
  
  researcher <- {{researcher_name}}
  work_order <- {{work_order_number}}
  
  log_file <- log_data %>% 
    purrr::set_names('line')
  
  number_files_analyzed <- log_file %>% 
    dplyr::filter(grepl('files will be processed',
                        line)) %>% 
    dplyr::mutate(
      line = as.double(str_remove(line, 
                                  ' files will be processed'))
      ) %>% 
    dplyr::pull()
  
  diann_pass <- tibble::tibble(
    diann_pass = c(rep('first pass',
                       number_files_analyzed),
                   rep('mbr',
                       number_files_analyzed))
    ) %>% 
    dplyr::mutate(
      diann_pass = stringr::str_replace(diann_pass,
                                        'first pass',
                                        'First-Pass'),
      diann_pass = stringr::str_replace(diann_pass,
                                        'mbr',
                                        'MBR')
      )
  
  run_number <- log_file %>% 
    dplyr::filter(grepl('File #',
                        line)) %>% 
    dplyr::slice(1:(number_files_analyzed*2)) %>% 
    dplyr::mutate(
      line = stringr::str_remove(line,
                                 '\\[.*\\] File #'),
      line = stringr::str_remove(line,
                                 '/.*')
      ) %>% 
    purrr::set_names('run_number')
  
  n_searches <- nrow(run_number)
  
  run_name <- log_file %>% 
    dplyr::filter(grepl('Loading run',
                        line)) %>% 
    dplyr::slice(1:(number_files_analyzed*2)) %>% 
    dplyr:: mutate(
      line = stringr::str_remove(
        line,
        '\\[.*\\] Loading run \\\\\\\\smdnas.*\\\\MSRL.*\\\\'
      ),
      line = stringr::str_remove(line,
                                 '\\.raw'),
      line = stringr::str_remove(line,
                                 researcher),
      line = stringr::str_remove(line,
                                 work_order)
      ) %>% 
    purrr::set_names('file_name')
  
  lib_size <- log_file %>% 
    dplyr::filter(grepl('library precursors are potentially detectable',
                        line)) %>% 
    dplyr::mutate(
      line = stringr::str_remove(line,
                                 '\\[.*\\] '),
      line = stringr::str_remove(line,
                                 ' library precursors.*'),
      line = as.double(line)
      ) %>% 
    purrr::set_names('precursors_in_lib') %>% 
    dplyr::distinct(precursors_in_lib) %>% 
    dplyr::bind_cols(dplyr::distinct(diann_pass)) %>% 
    dplyr::full_join(diann_pass,
                     by = 'diann_pass')
  
  precursor_ids <- log_file %>% 
    dplyr::filter(grepl('Number of IDs',
                        line)) %>% 
    dplyr::mutate(
      line = stringr::str_remove(line,
                                 '.* Number of IDs at.*:'),
      line = as.double(line)
      ) %>% 
    purrr::set_names('precursor_ids')
  
  protein_ids <- log_file %>% 
    dplyr::filter(grepl('Number of genes',
                        line)) %>% 
    dplyr::mutate(
      line = stringr::str_remove(line,
                                 '.* Number of genes.*, '),
      line = stringr::str_remove(line,
                                 ' \\(.*'),
      line = as.double(line)
      ) %>% 
    purrr::set_names('protein_ids')
  
  log_file_info <- dplyr::bind_cols(run_number,
                                    run_name) %>% 
    dplyr::bind_cols(lib_size) %>% 
    dplyr::bind_cols(precursor_ids) %>% 
    dplyr::bind_cols(protein_ids)
  
  run_info <- do.call('rbind',
                      replicate(n_searches,
                                experiment_metadata,
                                simplify = FALSE))
  
  coverage_data <- dplyr::bind_cols(run_info,
                                    log_file_info) %>% 
    dplyr::select(-run_number) %>% 
    dplyr::relocate(file_name,
                    .before = 'matrix') %>% 
    dplyr::relocate(precursors_in_lib,
                    .after = 'diann_pass')
  
  return(coverage_data)
} 