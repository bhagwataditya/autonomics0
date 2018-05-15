
#' Is a valid contrast?
#' @param contrast contrast
#' @param design   design
#' @return logical
#' @examples
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    design <- autonomics.find::create_design_matrix(object)
#'    contrast <- autonomics.find::default_contrasts(object)[1]
#'    autonomics.find::is_valid_contrast(contrast, design)
#' }
#' @importFrom magrittr   %<>%   %>%
#' @export
is_valid_contrast <- function(contrast, design){
  contrast %<>% as.character()
  contrast %<>% stringi::stri_replace_all_fixed(' ', '')
  contrast %<>% stringi::stri_replace_all_regex('(/[0-9]+)', '') %>%
                stringi::stri_replace_all_fixed('(', '') %>%
                stringi::stri_replace_all_fixed(')', '')
  terms <- strsplit(contrast, '[-+ ]+') %>% unlist() %>% unname()
  all(terms %in% colnames(design))
}

#' Select valid contrasts
#' @param  contrasts vector with contrast definitions
#' @param  design    design matrix
#' @param  verbose   whether or not to report
#' @return subset of contrasts
#' @examples
#' require(magrittr)
#' @importFrom magrittr             %>%
#' @export
select_valid_contrasts <- function(contrasts, design, verbose = TRUE){
  selector <- contrasts %>% vapply(is_valid_contrast, logical(1), design)
  if (verbose){
    autonomics.support::cmessage('\t\tKeep %d valid contrasts (out of %d)', sum(selector), length(selector))
  }
  contrasts[selector]
}

#' Validy subgroups and contrasts
#' @param design    design matrix
#' @param contrasts contrasts vector
#' @return validified contrasts
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object  <- autonomics.data::billing2016
#'    design <- autonomics.find::create_design_matrix(object)
#'    contrasts <- autonomics.find::default_contrasts(object)
#'    contrasts %>% autonomics.find::validify_contrasts(design)
#' }
#' @importFrom magrittr   %>%    %<>%
#' @export
validify_contrasts <- function(contrasts, design){

  # Assert valid inputs
   assertive.types::assert_is_matrix(design)
   assertive.types::assert_is_numeric(design)

  # Validify contrast names
  if (!assertive.properties::has_names(contrasts)){
    names(contrasts) <- make.names(contrasts)
  }
  assertive.properties::assert_has_no_duplicates(names(contrasts))

  # Select valid contrasts
  contrasts %<>% select_valid_contrasts(design)
  contrasts
}
