


#' Split, arrange, count subgroups
#' @param subgroup_values vector with subgroup values
#' @param sep subgroup separator
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    object <- 'extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat' %>%
#'               system.file(package='atkin.2014') %>%
#'               autonomics.import::load_soma()
#'    subgroup_values <- object %>% autonomics.import::subgroup_values()
#'    subgroup_values %>% autonomics.import::sep
#'    subgroup_values %>% count_components()
#'    subgroup_values %>% split_components()
#'    subgroup_values %>% layout_subgroups()
#'    subgroup_values %>% count_subgroups()
#' }
#' @return integer

#' @importFrom magrittr %>%
dcast_subgroups <- function(subgroup_values, fill, fun.aggregate){
   n <- count_components(subgroup_values)
   formula <- sprintf('V%d', 1:(n-1)) %>% paste0(collapse=' + ') %>% paste0(' ~ V', n)
   subgroup_values %>% split_components() %>%
      data.table::dcast(formula = formula,
                        value.var = 'subgroup',
                        fill = fill,
                        fun.aggregate = fun.aggregate) %>%
      autonomics.support::matrixify()
}

#' Layout/Count subgroups
#' @param subgroup_values vector with subgroup values
#' @param sep subgroup separator
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    object <- 'extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat' %>%
#'               system.file(package='atkin.2014') %>%
#'               autonomics.import::load_soma()
#'    subgroup_values <- object %>% autonomics.import::subgroup_values()
#'    subgroup_values %>% layout_subgroups()
#'    subgroup_values %>% count_subgroups()
#' }
#' @return integer
#' @importFrom magrittr %>%
#' @export
layout_subgroups <- function(subgroup_values) subgroup_values %>% as.character() %>% unique() %>%
   dcast_subgroups(fill = '', fun.aggregate = unique)

#' @rdname layout_subgroups
#' @importFrom magrittr %>%
#' @export
count_subgroups  <- function(subgroup_values) subgroup_values %>%
   dcast_subgroups(fill = 0, fun.aggregate = length)
