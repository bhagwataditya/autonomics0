utils::globalVariables('.')


#' Write top features for specified contrast
#'
#' @param object        eset
#' @param design          design matrix
#' @param contrast        see \code{analyze_eset}
#' @param direction       any value in autonomics.find::DIRECTIONS
#' @param top_definition  top definition
#' @param result_dir      directory
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::billing2016
#'    contrast <- object %>% autonomics.find::default_contrasts() %>% extract(1)
#'    object %<>% autonomics.find::add_limma_to_fdata(contrasts = contrast)
#'    direction <- 'neg'
#'    top_definition <- autonomics.find::default_top_definition(object)
#'    result_dir = tempdir()
#'    autonomics.find::write_top_features(object, contrast = contrast,
#'       direction = direction, top_definition = top_definition, result_dir = result_dir)
#' }
#' if (require(billing.differentiation.data)){
#'    require(magrittr)
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
   top_definition = default_top_definition(object),
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
   object %<>% select_fvars_and_filter_samples(design, contrast)
   object %<>% filter_n_arrange_top_features(names(contrast), top_definition, direction)

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
#' @param object       eset
#' @param design         design matrix
#' @param contrasts      A character vector of model contrasts to be written.
#' @param top_definition top definition
#' @param result_dir     dir where to save results
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object  <- autonomics.data::billing2016 %>% 
#'                 autonomics.find::add_limma_to_fdata()
#'    contrasts <- autonomics.find::default_contrasts(object)
#'    top_definition <- object %>% autonomics.find::default_top_definition()
#'    result_dir <- tempdir()
#'    autonomics.find::write_features(object, contrasts = contrasts,
#'                                    top_definition = top_definition, result_dir = result_dir)
#' }
#' @importFrom magrittr  %>%
#' @export
write_features <- function(
   object,
   design = create_design_matrix(object),
   contrasts,
   top_definition,
   result_dir
){
   object %>% write_all_features(result_dir)
   for (i in seq_along(contrasts)){
      cur_contrast <- contrasts[i]
      for (direction in c('neg', 'pos')){
         write_top_features(object, design, cur_contrast, direction = direction, top_definition, result_dir = result_dir)
      }
   }
}

