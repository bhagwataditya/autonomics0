

#' Infer organism
#' @param keys    keys
#' @param keytype any value in names(FEATURE_IDENTIFIERS)
#' @param verbose whether to report message (logical)
#' @return organism name
#' @examples
#' require(magrittr)
#' keys <- c('ENSG00000000003',  'ENSG00000000005', 'ENSG00000000419')
#' infer_organism(keys, 'ensg')
#' @export
infer_organism <- function(keys, keytype, verbose = TRUE){

   overlaps <- FEATURE_IDENTIFIERS %>%
               extract2(keytype) %>%
               lapply(extract2, 'keys') %>%
               vapply(function(x){
                         selector <- keys %in% x
                         floor(100*sum(selector)/length(selector))
                      }, numeric(1))
   organism <- names(overlaps) %>% extract(which.max(overlaps))
   if (verbose){
      autonomics.support::cmessage('\t\t%s: %d %% %s match %s',
                                   organism,
                                   floor(overlaps[[organism]]),
                                   keytype,
                                   get_orgdb_name(organism))
   }
   return(organism)
}

