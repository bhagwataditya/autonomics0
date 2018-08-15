
#' Select relevant fvars
#' @param object  eset
#' @param contrast_names  contrast names
#' @return object
#' @examples 
#' require(magrittr)
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    contrasts <- object %>% autonomics.find::default_contrasts()
#'    object %<>% autonomics.find::add_limma_to_fdata(contrasts = contrasts)
#'    object %>% autonomics.import::fvars() %>% length()
#'    object %>% autonomics.find::select_relevant_fvars(names(contrasts[1:2])) %>% 
#'               autonomics.import::fvars() %>% length()
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% autonomics.import::fdata() %>% dim()
#'    contrasts <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    object %>% autonomics.find::select_relevant_fvars(names(contrasts)) %>% 
#'                 autonomics.import::fdata() %>% dim()
#' }
#' # This shows why it is useful to have contrast_names (only) as 
#' # inputs rather than contrasts
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::exiqon
#'    (contrast_names <- object %>% autonomics.find::infer_contrast_names())
#'    contrast_names %<>% magrittr::extract(-c(2,3))
#'    object %>% autonomics.find::select_relevant_fvars(contrast_names) %>% 
#'                 autonomics.import::fdata() %>% dim()
#' }
#' @importFrom magrittr   %>%   %<>%
#' @export
select_relevant_fvars <- function(object, contrast_names){
  analysis_fvars <- autonomics.import::fvars(object) %>% 
                    magrittr::extract(stringi::stri_detect_regex(., '^(rank|value|p|fdr|bonf|sigb|quantile)[.].*'))
  contrast_fvars <- lapply(contrast_names, 
                           function(x){
                              autonomics.import::fvars(object) %>% 
                              magrittr::extract(stringi::stri_detect_regex(., sprintf('^(rank|value|p|fdr|bonf|sigb|quantile)[.]%s$', x)))
                           }) %>% 
                    unlist()
  annotation_fvars <- autonomics.import::fvars(object) %>% setdiff('fasta_hdrs') %>% setdiff(analysis_fvars)
  autonomics.import::fdata(object) %<>% magrittr::extract(c(annotation_fvars, contrast_fvars))
  return(object)
}

#' Filter relevant samples
#' @param object  eset
#' @param design    design matrix
#' @param contrasts contrasts
#' @examples 
#' require(magrittr)
#' 
#' if (require(autonomics.data)){
#'    object <- autonomics.data::ALL
#'    unique(object$subgroup)
#'    contrasts <- c(M_F = '(TM+BM)/2 - (TF+BF)/2', B_T = '(BM+BF)/2 - (TM+TF)/2')
#'    object %<>% autonomics.find::add_limma_to_fdata(contrasts = contrasts)
#'    object %>%  autonomics.find::are_relevant_samples(contrasts = contrasts)
#'    object %>%  autonomics.find::filter_relevant_samples(contrasts = contrasts)
#' }
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    contrasts <- object$subgroup[1]
#'    object %>% autonomics.find::are_relevant_samples(contrasts = contrasts)
#'    object %>% autonomics.find::filter_relevant_samples(contrasts = contrasts)
#' }
#' 
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% dim()
#'    contrasts <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    object %>% autonomics.find::are_relevant_samples(   contrasts = contrasts)
#'    object %>% autonomics.find::filter_relevant_samples(contrasts = contrasts) %>% dim()
#' }
#' @rdname filter_relevant_samples
#' @importFrom magrittr   %<>%   %>%
#' @export
are_relevant_samples <- function(object, design = create_design_matrix(object), contrasts){
   assertive.base::assert_all_are_true(nrow(design) == ncol(object))
   contrast_mat <- design %>% create_contrast_matrix(contrasts)
   
   design_vars <- matrixStats::rowAnys(contrast_mat != 0)
   samples     <- matrixStats::rowAnys((design[, design_vars, drop = FALSE]!=0))
   samples
}

#'@rdname filter_relevant_samples
#'@importFrom magrittr %>%
#'@export
filter_relevant_samples <- function(object, design = create_design_matrix(object), contrasts){
   idx <- are_relevant_samples(object, design = design, contrasts = contrasts)
   object %>% magrittr::extract(, idx)
}

#' Select portion of SummarizedExperiment relevant to contrasts
#' @param object   eSet
#' @param design     design matrix
#' @param contrasts   contrasts of subgroup levels (named string)
#' @examples 
#' require(magrittr)
#' 
#' # STEM CELL COMPARISON   
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios %>% autonomics.find::add_limma_to_fdata()
#'    contrasts <- object %>% autonomics.find::default_contrasts() %>% magrittr::extract(1)
#'    object %>% autonomics.find::select_fvars_and_filter_samples(contrasts = contrasts)
#' }
#' 
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% dim()
#'    contrasts <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    object %>% autonomics.find::select_fvars_and_filter_samples(contrasts = contrasts) %>% 
#'                 dim()
#' }
#' @return Filtered eSet
#' @importFrom magrittr %>%
#' @export
select_fvars_and_filter_samples <- function(
   object, 
   design = create_design_matrix(object), 
   contrasts
){
  object %<>% select_relevant_fvars(names(contrasts))
  object %<>% filter_relevant_samples(design, contrasts)
  return(object)
}
