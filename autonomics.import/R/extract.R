#' Extract features
#' @param object SummarizedExperiment
#' @param extractor logical/numeric vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    (object %<>% autonomics.import::extract_features(c(5,4)))
#'    object %>% autonomics.import::limma()
#' }
#' @importFrom magrittr %<>%
#' @export
extract_features <- function(object, extractor){
   object %<>% magrittr::extract(extractor, )
   if (!is.null(autonomics.import::limma(object))){
      autonomics.import::limma(object) %<>% magrittr::extract(autonomics.import::fnames(object), , , drop = FALSE)
   }
   object
}


#' Extract first from collapsed values
#' @param x charactervector or factorvector
#' @param sep collapsed string separator, e.g. ';'
#' @param ... to allow for S3 method dispatch
#' @return Updated x
#' @examples
#' require(magrittr)
#' x <- c('a;b;c', '1;2;3', 'alpha;beta;gamma')
#' x %>% autonomics.import::extract_first_from_collapsed(sep = ';')
#' @export
extract_first_from_collapsed <- function (x, ...) {
   UseMethod("extract_first_from_collapsed", x)
}

#' @rdname extract_first_from_collapsed
#' @importFrom magrittr %>%
#' @export
extract_first_from_collapsed.character <- function(x, sep = autonomics.import::guess_sep(x), ...){
   if (is.null(sep)) return(x)

   x %>%
   stringi::stri_split_fixed(sep) %>%
   vapply(magrittr::extract, character(1), 1)
}

#' @rdname extract_first_from_collapsed
#' @importFrom magrittr %<>%
#' @export
extract_first_from_collapsed.factor <- function(x, sep = guess_sep(x), ...){
   levels(x) %<>% extract_first_from_collapsed.character(sep=sep)
   x
}
