#' Write fvar to file
#' @param object eset
#' @param fvar fvar
#' @param file file
#' @examples
#' if (require(autonomics.data)){
#'    object <- ALL
#'    fvar <- 'gene_symbols'
#'    file <- tempfile()
#'    write_fvar_to_file(object, fvar, file)
#'    file.exists(file)
#'    unlink(file)
#' }
#' @importFrom magrittr  %<>%   %>%
#' @export
write_fvar_to_file <- function(object, fvar, file = ""){

   # Return if no fvar
   if (length(object)==0){
      autonomics.support::cmessage('Empty eset - abort')
      return()
   }
   if (length(fvar)==0){
      return()
   }
   if (!fvar %in% autonomics.import::fvars(object)){
      warning(sprintf("'%s' not in autonomics.import::fvars(object) - abort", fvar))
      return()
   }

  # Replace separator with space
  fvaldf  <- autonomics.import::fdata(object) %>%  magrittr::extract(fvar)
  sep <- autonomics.import::sep(object)
  if (!is.null(sep)){
    fvaldf[[fvar]] %<>% stringi::stri_replace_all_fixed(sep, '\t')
  }

  # Write to file
  fvaldf %>% autonomics.support::print2txt(file, col.names = FALSE)

}

#' Write fdata to file
#' @param object    eset
#' @param file file
#' @examples
#' if (require(autonomics.data)){
#'    object <- ALL
#'    file <- tempfile()
#'    write_fdata_to_file(object, file)
#'    file.exists(file)
#'    unlink(file)
#' }
#' @importFrom magrittr  %<>%  %>%
#' @export
write_fdata_to_file <- function(object, file = ""){
  my_f   <- autonomics.import::fdata(object)
  if ('fasta_hdrs' %in% names(my_f)){
    my_f %<>% dplyr::select_(~ -fasta_hdrs)
  }
  #names(my_f) %<>% gsub(paste0('.', names(contrast)), '', .)
  my_diff <- autonomics.import::exprs(object) %>% as.data.frame()
  feature_table <- cbind(my_f, my_diff)
  feature_table %>% autonomics.support::print2txt(file)
}
