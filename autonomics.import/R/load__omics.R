


#===========================================
# GENERIC
#===========================================

#' Load omics
#' @param object                       SummarizedExperiment
#' @param platform                    'exiqon', 'maxquant', 'metabolon', 'soma'
#' @param sheet                        excel sheet number or name if applicable
#' @param quantity                     string: which quantity to extract into exprs
#' @param design_file                  path to design file
#' @param log2_transform               logical
#' @param log2_offset                  offset in mapping x -> log2(x+offset)
#' @param infer_design_from_sampleids  logical
#' @param design_sep                   string: design separator
#' @return sample dataframe
#' @importFrom magrittr %>%
#' @export
load_omics <- function(
   object,
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

   object %<>% autonomics.import::merge_sdata(design_df, by = autonomics.import::sampleid_varname(platform), verbose = FALSE)
   #if (!is.null(design_file)){
      #file_df <- autonomics.import::read_design(design_file)
      #object %<>% autonomics.import::merge_sdata(file_df, by = autonomics.import::sampleid_varname(platform))
   #}

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

