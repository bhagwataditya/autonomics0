


#===========================================
# GENERIC
#===========================================

#' Load omics
#' @param file                         path to omics data file
#' @param platform                    'exiqon', 'maxquant', 'metabolon', 'metabolonlipids', 'soma'
#' @param sheet                        excel sheet number or name if applicable
#' @param quantity                     string: which quantity to extract into exprs
#' @param design_file                  path to design file
#' @param log2_transform               logical
#' @param log2_offset                  offset in mapping x -> log2(x+offset)
#' @param infer_design_from_sampleids  logical
#' @param design_sep                   string: design separator
#' @param zero_subgroup_nas            logical: whether to zero subgroup NAs (i.e. NA value for all samples in a subgroup)
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'
#' # EXIQON
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                         package = 'subramanian.2016')
#'    file %>% autonomics.import::load_omics(
#'                platform = 'exiqon',
#'                infer_design_from_sampleids = TRUE)
#' }
#'
#' # PROTEINGROUPS
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemdiff/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    object <- file %>% autonomics.import::load_omics(
#'                          platform = 'maxquant',
#'                          quantity = 'Ratio normalized',
#'                          infer_design_from_sampleids = TRUE)
#' }
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_omics(sheet=2, platform = 'metabolon')
#' }
#'
#' # METABOLONLIPIDS
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_omics(sheet = 'Lipid Class Concentrations', platform = 'metabolonlipids')
#' }
#'
#' # RNASEQ
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    infer_design_from_sampleids <- TRUE
#' }
#' @importFrom magrittr %>%
#' @export
load_omics <- function(
   file,
   platform,
   sheet                       = 2,
   quantity                    = NULL,
   log2_transform              = TRUE,
   log2_offset                 = 0,
   design_file                 = NULL,
   infer_design_from_sampleids = FALSE,
   design_sep                  = NULL
){
   # Satisfy CHECK
   . <- NULL

   # Load components
   sdata1 <- file %>% autonomics.import::load_sdata(platform = platform, sheet = sheet, quantity = quantity)
   fdata1 <- file %>% autonomics.import::load_fdata(platform = platform, sheet = sheet)
   exprs1 <- file %>% autonomics.import::load_exprs(platform = platform, sheet = sheet, quantity = quantity)

   # Wrap into SummarizedExperiment
   object <- SummarizedExperiment::SummarizedExperiment(assays=list(exprs = exprs1))
   if (log2_transform)   autonomics.import::exprs(object) %<>% (function(x)log2(x+log2_offset))
   autonomics.import::sdata(object)  <- sdata1
   autonomics.import::fdata(object)  <- fdata1

   # Merge in design
   design_df <- autonomics.import::write_design(file, platform                    = platform,
                                                      infer_design_from_sampleids = infer_design_from_sampleids,
                                                      quantity                    = quantity,
                                                      design_sep                  = design_sep,
                                                      sheet                       = sheet) %>%
                autonomics.support::rm_empty_vars()

   object %<>% autonomics.import::merge_sdata(design_df, by = autonomics.import::sampleid_varname(platform))
   if (!is.null(design_file)){
      file_df <- autonomics.import::read_design(design_file)
      object %<>% autonomics.import::merge_sdata(file_df, by = autonomics.import::sampleid_varname(platform))
   }

   # Order on subgroup (and replicate)
   subgroup_values  <- object %>% autonomics.import::svalues('subgroup')
   replicate_values <- object %>% autonomics.import::svalues('replicate')
   if (!is.null(subgroup_values)){
      if (!is.null(replicate_values)){ object %<>% magrittr::extract(, order(subgroup_values, replicate_values))
      } else {                         object %<>% magrittr::extract(, order(subgroup_values))}
   }

   # Return
   object
}

