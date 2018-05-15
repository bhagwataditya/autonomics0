has_matching_subgroup_levels <- function(contrast, an_eset, 
                                         .contrastname = assertive.base::get_name_in_parent(contrast), 
                                         .an_esetName = assertive.base::get_name_in_parent(an_eset)){
  terms <- contrast %>% strsplit('[-+() ]') %>% unlist()
  terms <- terms[terms != '']
  if (!all(terms %in% unique(autonomics.import::sdata(an_eset)$subgroup))){
    return(assertive.base::false('%s contains terms which are not among the levels of %s$subgroup', 
                 .contrastname, .an_esetName))
  }
  TRUE
}
is_valid_limma_contrast <- function(contrast, an_eset){
  if (!(ok <- assertive.types::is_character(contrast)))           return(ok)
  if (!(ok <- assertive.properties::is_scalar(contrast)))         return(ok)
  if (!(ok <- has_matching_subgroup_levels(contrast, an_eset)))   return(ok)
  TRUE
}

#' Assert validity of limma contrast
#' @param contrast limma contrast (character)
#' @param an_eset eSet
#' @export
assert_is_valid_limma_contrast <- function(contrast, an_eset){
  .contrastname = assertive.base::get_name_in_parent(contrast)
  .an_esetName  = assertive.base::get_name_in_parent(an_eset)
  assertive.base::assert_engine(
    is_valid_limma_contrast, 
    contrast = contrast, 
    an_eset = an_eset, 
    .contrastname = .contrastname, 
    .an_esetName = .an_esetName
  )
}


