#' Convert sumexp into wide datatable
#' @param object sumexp
#' @param fid    fvar carrying feature id
#' @param fvars  additional fvars to include in table
#' @param ...    (backward compatibility)
#' @return wide data.table
#' @importFrom data.table  data.table
#' @importFrom magrittr    %>%
#' @export
sumexp_to_wide_dt <- function(
   object,
   fid = 'feature_id',
   fvars = character(0)
){
  fdata1 <- autonomics.import::fdata(object) %>% magrittr::extract(, unique(c(fid, fvars)), drop = FALSE) %>% data.table::data.table()
  exprs1 <- autonomics.import::exprs(object) %>% data.table::data.table()
  wide1  <- cbind(fdata1, exprs1)
  wide1
}

#' @rdname sumexp_to_wide_dt
eset_to_wide_table <- function(...){
   .Deprecated('sumexp_to_wide_dt')
   sumexp_to_wide_dt(...)
}

#' Convert sumexp into long datatable
#' @param object sumexp
#' @param fid    fvar carrying feature id
#' @param fvars  additional fvars to include in table
#' @param sid    svar carrying sample id
#' @param svars  additional svars to include in table
#' @param ...    (backward compatibility)
#' @return long data.table
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::protein.ratios
#'    object %>% autonomics.import::sumexp_to_long_dt()
#'    object %>% autonomics.import::sumexp_to_long_dt(svars = 'subgroup')
#' }
#' @importFrom data.table   data.table
#' @importFrom magrittr     %>%
#' @export
sumexp_to_long_dt <- function(
   object,
   fid = 'feature_id',
   fvars = character(0),
   sid = 'sample_id',
   svars = character(0)
){
  sdata1 <- autonomics.import::sdata(object) %>% magrittr::extract(, c('sample_id', svars), drop = FALSE)
  # Note: unique is to avoid duplication of same fields in fid and fvars
  object %>% autonomics.import::sumexp_to_wide_dt(fid, fvars) %>%
             data.table::melt.data.table(id.vars = unique(c(fid, fvars)), variable.name = sid, value.name = 'value') %>%
             merge(sdata1, by = sid) %>%
             magrittr::extract(, unique(c(fid, fvars, sid, svars, 'value')), with = FALSE)
}

#' @rdname sumexp_to_long_dt
eset_to_long_table <- function(...){
   .Deprecated('sumexp_to_long_dt')
   sumexp_to_long_dt(...)
}
