utils::globalVariables(c('feature1', 'feature2', 'value1', 'value2'))

#' Correlate eset features
#'
#' Compute pairwise correlations between each featrue of eset1 and eset2.
#'
#' Both esets must contain an fvar 'feature_id' and an svar 'sample_id'.
#' Then, correlations will be computed for each pair of feature_id from eset1 and eset2.
#' Samples between the two esets will be matched on sample_id.
#'
#' @param eset1 eset
#' @param eset2 eset
#' @examples
#' if (require(subramanian.2016)){
#'    require(magrittr)
#'    exiqon <- subramanian.2016::exiqon %>%
#'              autonomics.import::filter_samples(condition==unique(.$condition)[2]) %>%
#'              magrittr::extract(1:10, )
#'    metabolon <- subramanian.2016::metabolon %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2]) %>%
#'                 magrittr::extract(rowSums(!is.na(autonomics.import::exprs(.))) == ncol(.), ) %>%
#'                 magrittr::extract(1:10, )
#'    rnaseq <- subramanian.2016::rnaseq %>%
#'              autonomics.import::filter_samples(condition==unique(.$condition)[2]) %>%
#'              magrittr::extract(1:10, )
#'    autonomics.integrate::correlate_features(exiqon, metabolon)
#'    autonomics.integrate::correlate_features(exiqon, rnaseq)
#'    autonomics.integrate::correlate_features(rnaseq, metabolon)
#' }
#' @note Avoid the use of devtools::load_all in combination with this function,
#'       as a devtools bug prevents data.table (which is used in this function) from being
#'       loaded properly:
#'       https://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
#'       Upon regular installation with devtools::install or R CMD INSTALL, this function works fine.
#' @importFrom data.table   data.table  :=
#' @importFrom magrittr     %>%   %<>%
#' @export
correlate_features <- function(eset1, eset2){

  # Process
  name1 <- assertive.base::get_name_in_parent(eset1)
  name2 <- assertive.base::get_name_in_parent(eset2)

  # Assert
  autonomics.import::assert_is_valid_eset(eset1)
  autonomics.import::assert_is_valid_eset(eset2)
  sample_id_values1 <- eset1 %>% autonomics.import::get_sample_id_values()
  sample_id_values2 <- eset2 %>% autonomics.import::get_sample_id_values()
  assertive.sets::assert_are_intersecting_sets(sample_id_values1, sample_id_values2)

  # Ensure that eset1 has the lesser number of features
  if (nrow(eset1) > nrow(eset2)){
    tmp <- eset1;  eset1 <- eset2;  eset2 <- tmp
    tmp <- name1;  name1 <- name2;  name2 <- name1
  }

  # eset1 -> data.table
  fdata1 <- data.table::data.table(feature1 = eset1 %>% autonomics.import::fid_values())
  exprs1 <- autonomics.import::exprs(eset1) %>% data.table::data.table()
  wide1  <- cbind(fdata1, exprs1)
  long1  <- wide1 %>% data.table::melt.data.table(id.vars = 'feature1', variable.name = 'sample_id', value.name = 'value1')

  # eset2 -> data.table
  fdata2 <- data.table::data.table(feature2 = eset2 %>% autonomics.import::fid_values())
  exprs2 <- autonomics.import::exprs(eset2) %>% data.table::data.table()
  wide2  <- data.table::data.table(fdata2, exprs2)
  long2  <- wide2 %>% data.table::melt.data.table(id.vars = 'feature2', variable.name = 'sample_id', value.name = 'value2')

  # Compute correlation
  corlist <- lapply(
    1:nrow(eset1),
    function(i){
      cur.wide1  <- wide1[i, ]
      cur.long1  <- cur.wide1 %>%
                    data.table::melt.data.table(id.vars = 'feature1', variable.name = 'sample_id', value.name = 'value1')
      full.dt <- merge(long2, cur.long1, by = 'sample_id')
      full.dt[,
              list(cor = stats::cor(value1, value2, use = 'pairwise.complete.obs')),
              by = list(feature1, feature2)]
    }
  )
  cor_dt <- data.table::rbindlist(corlist)
  cor_dt %>% data.table::setnames(c('feature1', 'feature2'), c(name1, name2))
  cor_dt %>% data.table::setorderv('cor')
  cor_dt

}

#' Annotate correlations
#' @param cor_dt correlation datatable returned by correlate_features
#' @param eset1 eset1
#' @param eset2 eset2
#' @param fvars1 eset1 fvars
#' @param fvars2 eset2 fvars
#' @return correlations datatable with annotations
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    eset1 <- subramanian.2016::exiqon %>%
#'             autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    eset2 <- subramanian.2016::metabolon %>%
#'             autonomics.import::filter_samples(condition==unique(.$condition)[2]) %>%
#'             magrittr::extract(rowSums(!is.na(autonomics.import::exprs(.))) == ncol(.))
#'    cor_dt <- autonomics.integrate::correlate_features(eset1, eset2)
#'    cor_dt %<>% autonomics.integrate::annotate_correlations(
#'                   eset1, eset2, character(0), c('BIOCHEMICAL', 'SUB_PATHWAY'))
#' }
#' @importFrom magrittr     %<>%         %>%
#' @importFrom data.table   data.table   :=
#' @export
annotate_correlations <- function(cor_dt, eset1, eset2, fvars1 = character(0), fvars2 = character(0)){

  # Process inputs
  name1 <- names(cor_dt)[1]
  name2 <- names(cor_dt)[2]
  assertive.base::is_not_true(name1 == name2)  # Two names should differ, as they are used

  # Assert
  autonomics.import::assert_is_valid_eset(eset1)
  autonomics.import::assert_is_valid_eset(eset2)
  assertive.types::assert_is_character(fvars1)
  assertive.types::assert_is_character(fvars2)
  assertive.sets::assert_is_subset(fvars1, autonomics.import::fvars(eset1))
  assertive.sets::assert_is_subset(fvars2, autonomics.import::fvars(eset2))
  assertive.sets::assert_are_disjoint_sets(fvars1, fvars2)
  assertive.sets::assert_are_disjoint_sets(c(name1, name2), c(fvars1, fvars2))

  # Annotate
  if (length(fvars1)>0){
     autonomics.import::fdata(eset1)$feature_id <- autonomics.import::fdata(eset1)[[autonomics.import::fid_var(eset1)]]
     annotations1 <- autonomics.import::fdata(eset1) %>% magrittr::extract(, c('feature_id', fvars1), drop = FALSE)
     cor_dt %<>% merge(annotations1, by.x = name1, by.y='feature_id')
  }
  if (length(fvars2)>0){
     autonomics.import::fdata(eset2)$feature_id <- autonomics.import::fdata(eset2)[[autonomics.import::fid_var(eset2)]]
     annotations2 <- autonomics.import::fdata(eset2) %>% magrittr::extract(, c('feature_id', fvars2))
     cor_dt %<>% merge(annotations2, by.x = name2, by.y='feature_id')
  }

  # Order columns and rows
  cor_dt %<>% magrittr::extract(, c(name1, name2, 'cor', fvars1, fvars2), with = FALSE)
  cor_dt %>% data.table::setorderv('cor')

  # Return
  cor_dt
}

#' Re-order datatable on key vars
#'
#' The columns are reshuffled such that that the key variables come first.
#' The rows are ordered on the first key variable
#'
#' @param cor_dt correlation datatable
#' @param key_vars key variables
#' @return correlation datable, ordered on variable of interest
#' @examples
#' if (require(subramanian.2016)){
#'    require(magrittr)
#'    eset1  <- subramanian.2016::exiqon
#'    eset2  <- subramanian.2016::metabolon
#'    cor_dt <- subramanian.2016::top.cor.exiqon.metabolon
#'    cor_dt %<>% autonomics.integrate::annotate_correlations(
#'                   eset1, eset2, fvars2 = c('BIOCHEMICAL', 'SUB_PATHWAY'))
#'    cor_dt %<>% autonomics.integrate::order_correlations('BIOCHEMICAL')
#' }
#' @importFrom magrittr     %>%          %<>%
#' @importFrom data.table   data.table   :=
#' @export
order_correlations <- function(cor_dt, key_vars){

  # Assert
  assertive.types::assert_is_data.table(cor_dt)
  assertive.types::assert_is_character(key_vars)
  assertive.sets::assert_is_subset(key_vars, names(cor_dt))

  # Order columns
  other_vars <- setdiff(names(cor_dt), c(key_vars, 'cor'))
  cor_dt %<>% magrittr::extract(, c(key_vars, other_vars, 'cor'), with = FALSE)

  # Order rows
  order_var <- key_vars[1]
  cor_dt %>% data.table::setorderv('cor')
  cor_dt[, (order_var) := autonomics.support::factorify(get(order_var))]
  cor_dt %>% data.table::setorderv(c(order_var, 'cor'))
  cor_dt[, (order_var) := as.character(get(order_var))]
  cor_dt
}
