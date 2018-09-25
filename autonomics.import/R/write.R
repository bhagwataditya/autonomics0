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

#' Write sumexp to file
#' @param object            SummarizedExperiment
#' @param limma_quantities  character vector: which limma quantities to include
#' @return data.table       data.table
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.import::flatten() %>% head()
#' }
#' @importFrom magrittr  %<>%  %>%
#' @export
flatten <- function(
   object,
   limma_quantities = c('effect', 'p', 'fdr')
){

   # fdata
   dt   <- autonomics.import::fdata(object) %>% data.table::data.table()
   if ('fasta_hdrs' %in% names(dt)){
      dt[, fasta_hdrs := NULL]
   }

   # limma
   limma_array <- autonomics.import::limma(object)
   if (!is.null(limma_array)){
      limma_array %<>% magrittr::extract(autonomics.import::fnames(object), , limma_quantities, drop = FALSE)
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
   exprs_dt <- autonomics.import::exprs(object) %>% data.table::data.table()
   dt %<>% cbind(exprs_dt)

   # Order on F.p
   dt %<>% magrittr::extract(order(F.p))

   # return
   dt
}

#' @export
write_fdata_to_file <- function(object, file = ""){
   .Deprecated('flatten')
   object %>% flatten(object) %>% autonomics.support::print2txt(file)
}
