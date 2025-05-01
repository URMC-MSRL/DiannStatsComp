#'
#' @description Combines `set_metadata()`, `get_phospho_log()`, `get_stats()`, and `merge_stats()` into a single function. This is used if there is a stats report file created for both the First-Pass and MBR search results for a DIA-NN search with phosphorylation set as a variable mod of STY.
#'
#' @param log_data The log file. Needs to be read into R using `readr::read_tsv` prior.
#' @param first_pass_stats_tsv The stats report file from the First-Pass search. Needs to be read into R using "readr::read_tsv()" prior. This should just be the First-Pass report.
#' @param mbr_stats_tsv The stats report file from the MBR search. Needs to be read into R using "readr::read_tsv()" prior. This should just be the MBR report.
#' @param researcher_name The name of the researcher (as found in the raw file name). The should be a string and formatted as follows: "Researcher_" where "Researcher" is whatever is in the raw file. If the raw file has no name, any string can be used. For example, "None" would be appropriate.
#' @param work_order_number The number of the project (as found in the raw file name). The should be a string and formatted as follows: "_WorkOrderNumber" where "WorkOrderNumber" is whatever is in the raw file; for the MSRL, this is typically something like "25-001". If the raw file has no work order number, any string can be used. For example, "None" would be appropriate.
#' @param matrix The type of material analyzed (i.e. what cell line or tissue). Should be a string.
#' @param injection_amount_ng The amount of material analyzed in nanograms. Should be an integer.
#' @param experiment_type The type of experiment performed (i.e. global vs phospho). Should be a string.
#' @param injection_method The method used to inject peptides (i.e. Direct Injection vs Trap & Elute) Should be a string.
#' @param separation_column The separation column used. Should be a string.
#' @param column_temp_C The temperature of the column oven. Should be an integer. 
#' @param mass_spec The mass spec used. Should be a string.
#' @param hplc The HPLC used. Should be a string.
#' @param hplc_method The HPLC method used. Should be a string.
#' @param ms_method The mass spec method used. Should be a string.
#' @param diann_version The version of DIA-NN used for the search. Should be a string.
#' @param mz_range The precursor m/z range set in DIA-NN. Should be a string.
#' @param z_range The precursor z range set in DIA-NN. Should be a string.
#' @param peptide_length The peptide length range set in DIA-NN.Should be a string.
#' @param missed_cleavages The maximum number of missed cleavages set in DIA-NN. Should be an integer.
#' @param max_var_mods The maximum number of variable modifications set in DIA-NN. Should be an integer.
#' @param mods_enabled The modifications enabled in DIA-NN. Should be a string.
#' 
#' @return A data frame containing the experimental details along with number of IDs as reported in both the log and stats report from both First-Pass and MBR searches. 
#' @export
#'
#' @importFrom magrittr %>%
#' 

phospho_summarized <- function(log_data,
                               first_pass_stats_tsv,
                               mbr_stats_tsv,
                               researcher_name = 'None_',
                               work_order_number = '_None',
                               matrix = 'HeLa Cells',
                               injection_amount_ng = 250,
                               experiment_type = 'Phospho',
                               injection_method = 'Trap & Elute',
                               separation_column = 'IO 15cm',
                               column_temp_C = 40,
                               hplc_method = '72 SPD',
                               mass_spec = 'Astral',
                               hplc = 'Vanquish Neo',
                               ms_method = 'Variable Window 3ms ITT',
                               diann_version = '1.8.1',
                               mz_range = '380-980',
                               z_range = '2-4',
                               peptide_length = '7-30',
                               missed_cleavages = 1,
                               max_var_mods = 1,
                               mods_enabled = 'Ox(M), Phos(STY)') {
  
  log_data <-  log_data
  first_pass_stats_tsv <- first_pass_stats_tsv
  mbr_stats_tsv <- mbr_stats_tsv
  
  researcher_name <- {{researcher_name}}
  work_order_number <- {{work_order_number}}
  matrix <- {{matrix}}
  injection_amount_ng <- {{injection_amount_ng}}
  experiment_type  <- {{experiment_type}}
  injection_method  <- {{injection_method}}
  separation_column  <- {{separation_column}}
  column_temp_C  <- {{column_temp_C}}
  hplc_method  <- {{hplc_method}}
  mass_spec <- {{mass_spec}}
  hplc  <- {{hplc}}
  ms_method  <- {{ms_method}}
  diann_version  <- {{diann_version}}
  mz_range <- {{mz_range}}
  z_range  <- {{z_range}}
  peptide_length  <- {{peptide_length}}
  missed_cleavages  <- {{missed_cleavages}}
  max_var_mods <- {{max_var_mods}}
  mods_enabled <- {{mods_enabled}}
  
  metadata <- set_metadata(
    matrix = matrix,
    injection_amount_ng = injection_amount_ng,
    experiment_type  = experiment_type,
    injection_method  = injection_method,
    separation_column  = separation_column,
    column_temp_C  = column_temp_C,
    hplc_method  = hplc_method,
    mass_spec = mass_spec,
    hplc  = hplc,
    ms_method  = ms_method,
    diann_version  = diann_version,
    mz_range = mz_range,
    z_range  = z_range,
    peptide_length  = peptide_length,
    missed_cleavages  = missed_cleavages,
    max_var_mods = max_var_mods,
    mods_enabled = mods_enabled
    )
  
  phospho_log_data <- get_phospho_log(
    log_data = log_data,
    experiment_metadata = metadata,
    researcher_name = researcher_name,
    work_order_number = work_order_number
    )
  
  first_pass_stats_data <- get_stats(
    stats_tsv = first_pass_stats_tsv,
    experiment_metadata = metadata,
    researcher_name = researcher_name,
    work_order_number = work_order_number,
    search = 'First-Pass'
    )
  
  mbr_stats_data <- get_stats(
    stats_tsv = mbr_stats_tsv,
    experiment_metadata = metadata,
    researcher_name = researcher_name,
    work_order_number = work_order_number,
    search = 'MBR'
    )
  
  merge_stats_data <- merge_stats(
    first_pass_stats_data,
    mbr_stats_data
    )
  
  merged_data <- dplyr::full_join(
    phospho_log_data,
    merge_stats_data,
    by = c(colnames(metadata),
           'diann_pass',
           'file_name')
    )
  
  return(merged_data)
  
}