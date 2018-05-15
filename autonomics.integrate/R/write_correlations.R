#' Write correlations to file
#'
#' Write correlations to file
#'
#' Note that fvars are included in the file, but feature ids not
#' (except if fvars are empty).
#'
#' @param cor.dt correlation datable
#' @param eset1 eset
#' @param eset2 eset
#' @param file  file
#' @param fvars1 eset1 fvars
#' @param fvars2 eset2 fvars
#' @examples
#' if (require(subramanian.2016)){
#'    cor.dt <- subramanian.2016::top.cor.exiqon.metabolon
#'    eset1  <- subramanian.2016::exiqon
#'    eset2  <-  subramanian.2016::metabolon
#'    fvars1 <- character(0)
#'    fvars2 <- c('BIOCHEMICAL', 'SUB_PATHWAY')
#' }
#' @importFrom magrittr    %>%          %<>%
#' @importFrom data.table  data.table   :=
#' @export
write_correlations <- function(cor.dt, eset1, eset2, file, fvars1 = character(0), fvars2 = character(0)){

  # Process inputs
   name1 <- names(cor.dt)[1]
   name2 <- names(cor.dt)[2]

   # Assert
   assertive.types::assert_is_character(file)
   autonomics.import::assert_is_valid_eset(eset1)
   autonomics.import::assert_is_valid_eset(eset2)
   assertive.sets::assert_is_subset('feature_id', autonomics.import::fvars(eset1))
   assertive.sets::assert_is_subset('feature_id', autonomics.import::fvars(eset2))
   assertive.sets::assert_is_subset(cor.dt[, get(name1)], autonomics.import::fdata(eset1)$feature_id)
   assertive.sets::assert_is_subset(cor.dt[, get(name2)], autonomics.import::fdata(eset2)$feature_id)
   assertive.sets::assert_is_subset(fvars1, autonomics.import::fvars(eset1))
   assertive.sets::assert_is_subset(fvars2, autonomics.import::fvars(eset2))
   assertive.sets::assert_are_disjoint_sets(fvars1, fvars2)

   # Merge fdata
   cor.dt[, pair := paste0(get(name1), ' : ', get(name2))]
   cor.dt[, pair := autonomics.support::factorify(pair)]
   cor.dt %<>% merge(autonomics.import::fdata(eset1)[, c('feature_id', fvars1), drop = FALSE], by.x = name1, by.y = 'feature_id')
   cor.dt %<>% merge(autonomics.import::fdata(eset2)[, c('feature_id', fvars2), drop = FALSE], by.x = name2, by.y = 'feature_id')
   if (length(fvars1)==0)  fvars1 <- name1
   if (length(fvars2)==0)  fvars2 <- name2
   cor.dt <- cor.dt[, c(fvars1, fvars2, 'cor', 'pair'), with = FALSE]
   cor.dt %>% data.table::setorder(pair)
   cor.dt[, pair:=NULL]

   # Write to table
   cor.dt %>% autonomics.support::print2txt(file)
}
