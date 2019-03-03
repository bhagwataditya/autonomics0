#' Switch between nondetect representations
#' @param object    SummarizedExperiment
#' @param verbose   logical(1)
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#'
#' # 0 -> NA (proteinGroups LFQ intensities)
#' #----------------------------------------
#'    if (require(graumann.lfq)){
#'       object <- 'extdata/proteinGroups.txt'            %>%
#'                  system.file(package = 'graumann.lfq') %>%
#'                  read_proteingroups()
#'       object %>% zero_to_na(verbose = TRUE)
#'    }
#'
#' # NaN -> NA (proteinGroups ratios)
#' #---------------------------------
#'    if (require(autonomics.data)){
#'       object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'                  system.file(package='autonomics.data')       %>%
#'                  read_proteingroups()
#'       object %>% nan_to_na(verbose = TRUE)
#'    }
#'
#' # -Inf -> NA (log2 transformed proteinGroups LFQ intensity)
#' #----------------------------------------------------------
#'    if (require(graumann.lfq)){
#'       object <- 'extdata/proteinGroups.txt'            %>%
#'                  system.file(package = 'graumann.lfq') %>%
#'                  read_proteingroups()
#'       object %<>% log2transform(verbose = TRUE)
#'       object %>% minusinf_to_na(verbose = TRUE)
#'    }
#'
#' # NA -> 0
#' #--------
#'    if (require(autonomics.data)){
#'       object <- 'extdata/glutaminase/glutaminase.xlsx' %>%
#'                  system.file(package = 'autonomics.data') %>%
#'                  read_metabolon()
#'       object %>% na_to_zero(verbose = TRUE)
#'    }
#'
#'
#' @importFrom magrittr %>%
#' @export
zero_to_na <- function(object, verbose = FALSE){
   selector <- autonomics.import::exprs(object) == 0
   if (any(c(selector), na.rm = TRUE)){
      if (verbose) autonomics.support::cmessage('\t\tReplace 0 -> NA for %d/%d values (in %d/%d features and %d/%d samples)',
                                                sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                                                selector %>% matrixStats::rowAnys() %>% sum(), nrow(object),
                                                selector %>% matrixStats::colAnys() %>% sum(), ncol(object))
      autonomics.import::exprs(object)[selector] <- NA_real_
   }
   object
}


#' @rdname zero_to_na
#' @importFrom magrittr %>%
#' @export
nan_to_na <- function(object, verbose = FALSE){
   selector <- is.nan(autonomics.import::exprs(object))
   if (any(c(selector), na.rm = TRUE)){
      if (verbose) autonomics.support::cmessage('\t\tReplace NaN -> NA for %d/%d values (in %d/%d features and %d/%d samples)',
                                                sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                                                selector %>% matrixStats::rowAnys() %>% sum(), nrow(object),
                                                selector %>% matrixStats::colAnys() %>% sum(), ncol(object))
      autonomics.import::exprs(object)[selector] <- NA_real_
   }
   object
}

#' @rdname zero_to_na
#' @importFrom magrittr %>%
#' @export
na_to_zero <- function(object, verbose = FALSE){
   selector <- autonomics.import::exprs(object) %>% is.na()
   if (any(selector)){
      if (verbose) autonomics.support::cmessage('\t\tReplace NA -> 0 for %d/%d values (in %d/%d features and %d/%d samples)',
                                                sum(selector), nrow(selector)*ncol(selector),
                                                selector %>% matrixStats::rowAnys() %>% sum(), nrow(object),
                                                selector %>% matrixStats::colAnys() %>% sum(), ncol(object))
      autonomics.import::exprs(object)[selector] <- 0
   }
   object
}



#' @rdname zero_to_na
#' @importFrom magrittr %>%
#' @export
minusinf_to_na <- function(object, verbose = FALSE){
   selector <- autonomics.import::exprs(object)==-Inf
   if (any(c(selector), na.rm = TRUE)){
      if (verbose) autonomics.support::cmessage('\t\tReplace -Inf -> NA for %d/%d values (in %d/%d features and %d/%d samples)',
                                                sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                                                selector %>% matrixStats::rowAnys() %>% sum(), nrow(object),
                                                selector %>% matrixStats::colAnys() %>% sum(), ncol(object))
      autonomics.import::exprs(object)[selector] <- NA_real_
   }
   object
}



