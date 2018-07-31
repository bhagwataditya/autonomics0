utils::globalVariables('.')


#' Write top features for specified contrast
#'
#' @param object          SummarizedExperiment
#' @param design          design matrix
#' @param contrast        see \code{analyze_eset}
#' @param direction       any value in autonomics.find::DIRECTIONS
#' @param topdef  top definition
#' @param result_dir      directory
#' @examples
#' require(magrittr)
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    result_dir <- tempdir() %T>% message()
#'    object %>% autonomics.find::write_top_features(
#'                  contrast       = c(BM_E = 'BM_E'),
#'                  direction      = 'neg', 
#'                  topdef = 'fdr < 0.5', 
#'                  result_dir     = result_dir)
#' }
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    result_dir <- tempdir() %T>% message()
#'    object %>% autonomics.find::write_top_features(
#'       contrast = billing.differentiation.data::contrasts[1], result_dir = result_dir)
#'    object %>% autonomics.find::write_top_features(
#'       contrast   = billing.differentiation.data::contrasts[1],
#'       direction  = 'both', 
#'       result_dir = result_dir)
#' }
#' @importFrom  magrittr   %<>%  %>%
#' @export
write_top_features <- function(
   object,
   design         = create_design_matrix(object),
   contrast,
   direction      = autonomics.find::DIRECTIONS[1],
   topdef = default_topdef(object),
   result_dir     = NULL
){

   # Assert valid inputs
   autonomics.import::assert_is_valid_eset(object)
   assertive.types::assert_is_matrix(design)
   assertive.types::assert_is_numeric(design)
   assertive.base::assert_is_identical_to_true(unname(ncol(object)) == nrow(design))
   assertive.strings::assert_is_a_non_empty_string(direction)
   assertive.sets::assert_is_subset(direction, autonomics.find::DIRECTIONS)

   # Keep top portion of eSet (corresponds with ora query!)
   object %<>% autonomics.find::select_fvars_and_filter_samples(design, contrast)
   object %<>% autonomics.find::filter_n_arrange_top_features(names(contrast), topdef, direction)

   # Return if eset empty
   if (nrow(object) == 0){
     autonomics.support::cmessage("\t\t%s %s 0   no top features - abort", contrast, direction_to_sign(direction))
     return(NULL)
   }

   # Create dir and file names
   if (is.null(result_dir)){
      table_file <- ""
      gene_file  <- ""
   } else {
      subdir <- autonomics.find::get_contrast_subdir(result_dir, names(contrast))
      dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
      table_file <- sprintf('%s/top_feature_table__%s__%s.txt', subdir, gsub('[\\-\\:]', '_', names(contrast)), direction)
      gene_file  <- sprintf('%s/top_gene_symbols__%s__%s.txt', subdir, gsub('[\\-\\:]', '_', names(contrast)), direction)
   }

   # Report
   autonomics.support::cmessage("\t\t%s %s 0   %s   %s",
                                contrast, 
                                autonomics.find::direction_to_sign(direction),
                                basename(table_file),
                                basename(gene_file))
   
   # Write top tables and top genes
   object %>% autonomics.import::write_fdata_to_file(table_file)
   object %>% autonomics.import::write_fvar_to_file(fvar = autonomics.import::fname_var(.), file = gene_file)
}
   
#' Write all features to file
#' @param object     eset
#' @param result_dir   dir
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::ALL
#'    result_dir <- tempdir()
#'    autonomics.find::write_all_features(object, result_dir)
#'    dir(result_dir)
#'    unlink(result_dir, recursive = TRUE)
#' }
#' @importFrom magrittr   %>%
#' @export
write_all_features <- function(object, result_dir){
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  table_file   <- sprintf('%s/all_feature_table.txt', result_dir)
  gsymbol_file <- sprintf('%s/all_gene_symbols.txt', result_dir)
  object %>% autonomics.import::write_fdata_to_file(table_file)
  object %>% autonomics.import::write_fvar_to_file(autonomics.import::fname_var(.), gsymbol_file)
  autonomics.support::cmessage("\t\tall   %s", basename(table_file))
}

#' Write feature tables to file
#' @param object         SummarizedExperiment
#' @param design         design matrix
#' @param contrasts      A character vector of model contrasts to be written.
#' @param direction     'neg', 'pos', 'both'
#' @param topdef top definition
#' @param result_dir     dir where to save results
#' @examples
#' require(magrittr)
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object  <- autonomics.data::stemcomp.proteinratios
#'    result_dir <- tempdir() %T>% message()
#'    object %>% autonomics.find::write_features(
#'                  contrasts      =  c(BM_E = 'BM_E', BM_EM = 'BM_EM', EM_E = 'EM_E'),
#'                  topdef = 'fdr < 0.05', 
#'                  result_dir     =  result_dir)
#' }
#' @importFrom magrittr  %>%
#' @export
write_features <- function(
   object,
   design    = autonomics.find::create_design_matrix(object),
   contrasts = autonomics.find::default_contrasts(object),
   direction = c('neg', 'pos'),
   topdef,
   result_dir
){
   object %>% autonomics.find::write_all_features(result_dir)
   for (i in seq_along(contrasts)){
      cur_contrast <- contrasts[i]
      for (curdirection in direction){
         object %>% autonomics.find::write_top_features(
                       design     = design, 
                       contrast   = cur_contrast, 
                       direction  = curdirection, 
                       topdef     = topdef, 
                       result_dir = result_dir)
      }
   }
}

