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
   if (!fvar %in% fvars(object)){
      warning(sprintf("'%s' not in fvars(object) - abort", fvar))
      return()
   }

  # Replace separator with space
  fvaldf  <- fdata(object) %>%  magrittr::extract(fvar)
  sep <- sep(object)
  if (!is.null(sep)){
    fvaldf[[fvar]] %<>% stringi::stri_replace_all_fixed(sep, '\t')
  }

  # Write to file
  fvaldf %>% autonomics.support::print2txt(file, col.names = FALSE)

}

#' Write sumexp to file
#' @param object            SummarizedExperiment
#' @param limma_quantities  character vector: which limma quantities to include
#' @return data.table       data.table
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% flatten() %>% head()
#' }
#' @importFrom magrittr  %<>%  %>%
#' @export
flatten <- function(
   object,
   limma_quantities = c('effect', 'p', 'fdr')
){
   # check
   fasta_hdrs <- F.p <- NULL

   # fdata
   dt   <- fdata(object) %>% data.table::data.table()
   if ('fasta_hdrs' %in% names(dt)){
      dt[, fasta_hdrs := NULL]
   }

   # limma
   limma_array <- limma(object)
   if (!is.null(limma_array)){
      limma_array %<>% magrittr::extract(fnames(object), , limma_quantities, drop = FALSE)
      limma_dt    <- limma_array %>%
                     matrix(nrow     = NROW(.),
                            ncol     = prod(NCOL(.), dim(.)[3]),
                            dimnames = list(rownames(.),
                                            autonomics.support::vsprintf('%s.%s', dimnames(.)[[3]], colnames(.)))) %>%
                     data.table::data.table()
      dt %<>% cbind(limma_dt)
      dt
   }

   # exprs
   exprs_dt <- exprs(object) %>% data.table::data.table()
   dt %<>% cbind(exprs_dt)

   # Order on F.p
   if (!is.null(limma_array)) dt %<>% magrittr::extract(order(F.p))

   # return
   dt
}


#' Write features to file
#' @param object SummarizedExperiment
#' @param file   txt file where to write results to
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object  <- autonomics.data::stemcomp.proteinratios
#'    file <- tempdir() %>% paste0('/features.txt') %T>% message()
#'    object %>% write_features(file)
#' }
#' @importFrom magrittr  %>%
#' @export
write_features <- function(
   object,
   file
){
   object %>% flatten() %>% autonomics.support::print2txt(file)
   autonomics.support::cmessage("\t\tall   %s", basename(file))
}

