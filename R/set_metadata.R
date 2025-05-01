#' Define the Experimental Details
#'
#' @description Creates a data frame that contains experimental details describing the files searched by DIA-NN. This is then used combined with the run names from the DIA-NN search results.
#'
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
#' @return A data frame containing various information pertinent to the run searched. 
#' @export
#'
#' @importFrom magrittr %>%
#'

set_metadata <- function(matrix = 'HeLa Cells',
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
  
  metadata_out <- tibble::tibble(
    'matrix' = c(matrix),
    'injection_amount_ng' = c(injection_amount_ng),
    'experiment_type' = c(experiment_type),
    'injection_method' = c(injection_method),
    'separation_column' = c(separation_column),
    'column_temp_C' = c(column_temp_C),
    'hplc_method' = c(hplc_method),
    'mass_spec' = c(mass_spec),
    'hplc' = c(hplc),
    'ms_method' = c(ms_method),
    'diann_version' = c(diann_version),
    'mz_range' = c(mz_range),
    'z_range' = c(z_range),
    'peptide_length' = c(peptide_length),
    'missed_cleavages' = c(missed_cleavages),
    'max_var_mods' = c(max_var_mods),
    'mods_enabled' = c(mods_enabled)
    )
  
  return(metadata_out)
  
}