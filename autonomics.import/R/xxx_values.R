
#' Split composite (s)values
#'
#' Get data.table in which composite (s)values have been split into components
#'
#' @param x    character(n), factor(n), or SummarizedExperiment
#' @param svar character(1)
#' @param sep  character(1):  separator
#' @param keep logical(1_: whether to keep original values (in first column)
#' @param ...  used for proper S3 method dispatch
#' @return data.table(x, x1, x2, ...) or data.table(x1, x2, ...)
#' @examples
#' if (require(autonomics.data)){
#'    # GLUTAMINASE
#'       require(magrittr)
#'       x <- autonomics.data::glutaminase
#'       x %>% subgroup_values() %>% split_values()
#'       x %>%                       split_values()
#' }
#' @export
split_values <- function (x, ...) {
   UseMethod("split_values", x)
}

#' @rdname split_values
#' @importFrom magrittr %>%
#' @export
split_values.character <- function(
   x,
   sep    = guess_sep(x, verbose = FALSE),
   keep   = FALSE,
   ...
){

   dt <- data.table::data.table(x = x)
   components <- if (is.null(sep)){                                                               # Single component
      data.table::data.table(x1 = x)

   } else {                                                                         # Multiple components
      dt %>% magrittr::extract(, data.table::tstrsplit(x, sep, fixed = TRUE))
   }
   components %<>% magrittr::set_names(sprintf('x%d', 1:ncol(.)))

   if (keep) cbind(dt, components) else components

   # Old approach
   # y <- values %>% stringi::stri_split_fixed(sep)
   # n.component <- length(y[[1]])
   # 1:n.component %>% lapply(function(z) vapply(y, magrittr::extract, character(1), z)) %>%
   #                   data.table::as.data.table()
}

#' @rdname split_values
#' @importFrom magrittr %>%
#' @export
split_values.factor <- function(
   x,
   sep = guess_sep(x, verbose = FALSE),
   keep = FALSE,
   ...
){
   x %>% as.character() %>% split_values(sep = sep, keep = keep)
}

#' @rdname split_values
#' @export
#' @importFrom magrittr %>%
split_values.SummarizedExperiment <- function(
   x,
   svar = 'subgroup',
   keep = FALSE,
   ...
){
   dt <- x %>%
         svalues(svar) %>%
         split_values(keep = keep)
   dt
}

#' @rdname split_values
#' @export
scomponents <- function(...){
   .Deprecated('split_values')
   split_values(..., keep_x = FALSE)
}

#' @rdname split_values
#' @export
subgroup_components <- function(...){
   .Deprecated('split_values')
   split_values(..., keep_x = FALSE)
}

#==============================================================================================

#' Matrixify datatable or dataframe
#' @param df dataframe or datatable (first column contains rownames)
#' @return matrix
#' @importFrom magrittr %>% 
#' @export
matrixify <- function(df){
  rownames1 <- df[[1]]
  df[[1]] <- NULL
  mat <- df %>% as.matrix() #data.matrix()
  rownames(mat) <- rownames1
  mat
}


#' Reshape composite values
#'
#' Reshapes a vector with composite values into a (2d) matrix:
#'    - rownames = values of component 1:(n-1)
#'    - colnames = values of component n
#'    - values can be either composite levels themselves   (fill = '', fun.aggregate = unique)
#'                 or frequency counts of composite levels (fill = 0,  fun.aggregate = length)
#' @param x               character(n) with composite values
#' @param fill            missing value placeholder (passed to data.table::dcast)
#' @param fun.aggregate   'unique' (for layout purposes) or length (for tabulate purposes)
#' @param sep             component separator
#' @examples
#' if (require(autonomics.data)){
#'    # GLUTAMINASE
#'       require(magrittr)
#'       x <- autonomics.data::glutaminase %>% subgroup_values()
#'       x %>% reshape_values(fill = '', fun.aggregate = unique)
#'       x %>% reshape_values(fill = 0,  fun.aggregate = length)
#' }
#'
#' @importFrom magrittr %>%
#' @export
reshape_values <- function(x, fill, fun.aggregate, sep = guess_sep(x)){
   n <- split_values(x, sep = sep) %>% ncol()
   formula <- if (n==1) '1' else sprintf('x%d', 1:(n-1)) %>% paste0(collapse=' + ')
   formula %<>% paste0(' ~ x', n)
   x %>% split_values(keep = TRUE) %>%
         data.table::dcast(formula = formula,
                           value.var = names(.)[1],
                           fill = fill,
                           fun.aggregate = fun.aggregate) %>%
         matrixify()
}

#============================================================================================


#' Count composite (s)values
#'
#' Reshapes a vector with composite values into a 2d matrix:
#'    + rownames = values of component 1:(n-1)
#'    + colnames = values of component n
#'    + values = frequency counts of composite levels
#'
#' @param x      character(n) factor(n), or SummarizedExperiment with composite (s)values.
#' @param svar   character(1)
#' @param ...    used for S3 method dispatch
#' @examples
#' if (require(autonomics.data)){
#'
#'    # GLUTAMINASE
#'    require(magrittr)
#'    x <- autonomics.data::glutaminase
#'    x %>% subgroup_values() %>% count_values()
#'    x %>%                       count_values()
#' }
#' @return integer
#' @export
count_values <- function(x, ...){
   UseMethod('count_values', x)
}

#' @rdname count_values
#' @importFrom magrittr %>%
#' @export
count_values.character <- function(x, ...){
   x %>%
   reshape_values(fill = 0, fun.aggregate = length)
}

#' @rdname count_values
#' @importFrom magrittr %>%
#' @export
count_values.factor <- function(x, ...){
   x %>%
   as.character() %>%
   count_values()
}

#' @rdname count_values
#' @export
count_values.SummarizedExperiment <- function(x, svar = 'subgroup', ...){
   x %>%
   svalues(svar) %>%
   count_values()
}



#==============================================================================

#' Layout composite (s)values
#' Reshapes a vector with composite values into a 2d matrix:
#'    + rownames = values of component 1:(n-1)
#'    + colnames = values of component n
#'    + values = composite levels
#'
#' @param x      character(n), factor(n), or SummarizedExperiment with composite (s)values
#' @param sep    character(1)
#' @param svar   character(1)
#' @param ...    used for proper S3 dispatching
#' @examples
#' if (require(autonomics.data)){
#'    # GLUTAMINASE
#'    require(magrittr)
#'    x <- autonomics.data::glutaminase
#'    x %>% subgroup_values() %>% layout_values()
#'    x                       %>% layout_values()
#' }
#'
#' @return integer
#' @export
layout_values <- function (x, ...) {
   UseMethod("layout_values", x)
}


#' @rdname layout_values
#' @importFrom magrittr %>%
#' @export
layout_values.character <- function(x, sep = guess_sep(x), ...){
   x  %>%
   unique() %>%
   reshape_values(fill = '', fun.aggregate = unique, sep = sep)
}

#' @rdname layout_values
#' @importFrom magrittr %>%
#' @export
layout_values.factor <- function(x, sep = guess_sep(x), ...){
   x %>%
   levels() %>%
   reshape_values(fill = '', fun.aggregate = unique, sep = sep)
}

#' @rdname layout_values
#' @export
layout_values.SummarizedExperiment <- function(x, svar = 'subgroup', sep = guess_sep(x), ...){
   x %>%
   svalues(svar) %>%
   layout_values(sep = sep)
}


# #==========================================================================================
#
# # Count no of components
# # @param x character(n) or SummarizedExperiment
# # @param svar character(1)
# # @param sep character(1)
# # @return integer(1)
# # @examples
# # require(magrittr)
# # if (require(autonomics.data)){
# #    x <- autonomics.data::glutaminase
# #    x %>% count_components()
# #    autonomics.data::glutaminase$subgroup   %>% count_components()
# # }
# # @export
# count_components <- function (x, ...) {
#    UseMethod("count_components", x)
# }
#
#
# # @rdname count_components
# # @importFrom magrittr %>%
# # @export
# count_components.SummarizedExperiment <- function(
#    x,
#    svar = 'subgroup',
#    sep  = x %>% guess_sep()
# ){
#    n <- x %>%r
#       slevels(svar) %>%
#       (function(x) if (is.null(sep)) x else x %>% stringi::stri_split_fixed(sep)) %>%
#       vapply(length,integer(1)) %>%
#       unique()
#    assertive.properties::assert_is_scalar(n)
#    n
# }
#
# # @rdname count_components
# # @importFrom magrittr %>%
# # @export
# count_components.character <- function (x,
#                                         sep = x %>% unique() %>% guess_sep(.)
# ){
#    subgroup_levels <- x %>% as.character() %>% unique()
#    n <- subgroup_levels %>% (function(x) if (is.null(sep)) x else x %>% stringi::stri_split_fixed(sep)) %>% vapply(length, integer(1)) %>% unique()
#    assertive.properties::assert_is_scalar(n)
#    n
# }
#
# n_components <- function(...){
#    .Deprecated('count_components')
#    count_components(...)
# }
#
#





