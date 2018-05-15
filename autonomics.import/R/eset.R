
# Show
#-----
.collapse_values <- function(vars){
  if (length(vars) < 5){
    paste0(vars, collapse = '  ')
  } else {
    sprintf('%s  %s  %s', paste0(vars[1:2], collapse = '  '),
            paste0('...'),
            paste0(vars[length(vars)], collapse = '  '))
  }
}

.print_feature_info <- function(object){
  sprintf( '\neSet %5d features: %s   |   %d fvars: %s', nrow(object), .collapse_values(fnames(object)),
           length(fvars(object)), .collapse_values(fvars(object)))
}

.print_sample_info <- function(object){
  sprintf( '\n     %5d  samples: %s   |   %d svars: %s', ncol(object), .collapse_values(snames(object)),
           length(svars(object)), .collapse_values(svars(object)))
}

.print_assay_info <- function(object){
  sprintf('\n          assayData: %s', .collapse_values(Biobase::assayDataElementNames(object)))
}

.show_eset <- function(object){
  feature_info <- .print_feature_info(object)
  sample_info  <- .print_sample_info(object)
  assay_info   <- .print_assay_info(object)
  rule <- paste0(rep('-', 50), collapse = '')
  rule <- paste0('\n', rule)
  cat(rule)
  cat(feature_info)
  cat(sample_info)
  cat(assay_info)
  cat(rule)
}

#' Show eSet
#' @param object eSet
#' @rdname show
#' @importFrom methods show
#' @export
setMethod("show", signature("eSet"),       function(object) {.show_eset(object)})

