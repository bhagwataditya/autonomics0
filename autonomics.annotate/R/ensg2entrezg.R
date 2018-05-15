#' Map ensg to entrezg
#' @param x ensg vector
#' @param organism any value in ANNOTATED_ORGANISMS
#' @param verbose logical
#' @examples
#'    x <- c("ENSG00000066248", "ENSG00000066279", "ENSG99999", NA_character_)
#'    autonomics.annotate::ensg_to_entrezg(x)
#' @export
ensg_to_entrezg <- function(
   x,
   organism = infer_organism(keys = x, keytype = 'ensg'),
   verbose = TRUE
){
   entrezg <- rep(NA_character_, length(x))
   idx <- !is.na(x)
   #if (verbose) autonomics.support::cmessage('%d/%d/%d available ensg', sum(idx), length(idx))
   entrezg[idx] <- x[idx] %>% (function(y){
                      idy <- y %in% AnnotationDbi::keys(load_orgdb(organism), keytype = 'ENSEMBL')
                      if (verbose) autonomics.support::cmessage('%d ensg > %d available > %d mapped to entrezg', length(idx), length(idy), sum(idy))
                      AnnotationDbi::mapIds(load_orgdb(organism), keys = y, keytype = 'ENSEMBL', column = 'ENTREZID')
                   })
   entrezg
}

#' Map gsymbol to entrezg
#' @param x        gene symbol vector
#' @param organism any value in ANNOTATED_ORGANISMS
#' @param verbose  logical
#' @examples
#'    x <- c("NGEF", "ASPM", 'ABRACADABRA', NA_character_)
#'    organism <- 'Homo sapiens'
#'    gsymbol_to_entrezg(x, organism)
#' @export
gsymbol_to_entrezg <- function(
   x,
   organism,
   verbose = TRUE
){
   entrezg <- rep(NA_character_, length(x))
   idx <- !is.na(x)
   # if (verbose) autonomics.support::cmessage('%d/%d available gene symbols', sum(idx), length(idx))
   entrezg[idx] <- x[idx] %>% (function(y){
                      idy <- y %in% AnnotationDbi::keys(load_orgdb(organism), keytype = 'SYMBOL')
                      if (verbose) autonomics.support::cmessage('%d gsymbols > %d available > %d mapped to entrezg', length(idx), length(idy), sum(idy))
                      AnnotationDbi::mapIds(load_orgdb(organism), keys = y, keytype = 'SYMBOL', column = 'ENTREZID')
                  })
   entrezg
}
