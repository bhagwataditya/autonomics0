
#' Are relevant samples for given contrast
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
#' }
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    contrasts <- object$subgroup[1]
#'    object %>% autonomics.find::are_relevant_samples(contrasts = contrasts)
#' }
#' 
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    object %>% dim()
#'    contrasts <-  c(EM0.8_0 = 'EM0.8 - EM0.0')
#'    object %>% autonomics.find::are_relevant_samples(   contrasts = contrasts)
#' }
#' @importFrom magrittr   %<>%   %>%
#' @export
are_relevant_samples <- function(object, design = autonomics.import::create_design_matrix(object), contrasts){
   assertive.base::assert_all_are_true(nrow(design) == ncol(object))
   contrast_mat <- design %>% create_contrast_matrix(contrasts)
   
   design_vars <- matrixStats::rowAnys(contrast_mat != 0)
   samples     <- matrixStats::rowAnys((design[, design_vars, drop = FALSE]!=0))
   samples
}

