

#' Extract detectome go/interpro sets
#' @param object      SummarizedExperiment
#' @param oraset_var  ora set variable
#' @return list (names = setids, values = fids)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    object %>% autonomics.ora::extract_ora_sets('goid')     %>% head(2)
#'    object %>% autonomics.ora::extract_ora_sets('interpro') %>% head(2)
#' }
#' @importFrom magrittr  %>% 
#' @export
extract_ora_sets <- function(
   object, 
   oraset_var
){
   oraid_var <- autonomics.import::oraid_var(object)
   
   object                                                                  %>% 
   autonomics.import::fdata()                                                %>% 
   magrittr::extract(, c(oraid_var, oraset_var))                           %>% 
   tidyr::separate_rows(oraset_var, sep = ';')                             %>%
   data.table::data.table()                                                %>% 
   magrittr::extract(, (oraset_var) := trimws(get(oraset_var)))            %>%
   magrittr::extract(get(oraid_var)   %>% (function(x) !is.na(x) & x!='')) %>% 
   magrittr::extract(get(oraset_var) %>% (function(x) !is.na(x) & x!=''))  %>% 
   unique()                                                                %>%  # do not double count discriminated isoforms
  (function(x) split(x[[oraid_var]], x[[oraset_var]], drop = TRUE))
}

#' Extract ora universe
#' @param object SummarizedExperiment
#' @return character vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    object %>% autonomics.ora::extract_ora_universe()       %>% head(2)
#' }
#' @importFrom magrittr %>% 
#' @export
extract_ora_universe <- function(object){
   object                             %>% 
   autonomics.import::oraid_values()  %>% 
   unique()  # do not double count discriminated isoforms
}

#' Extract ora query 
#' @param object         SummarizedExperiment
#' @param contrast_name  character
#' @param topdef character
#' @param direction      'pos' or 'neg'
#' @return chracter vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    object <- autonomics.data::billing2016
#'    object %>% autonomics.ora::extract_ora_query('BM_E', 'fdr<0.05', 'pos')          %>% head(2)
#' }
#' @importFrom magrittr %>% 
#' @export
extract_ora_query <- function(object, contrast_name, topdef, direction){
   object                                                                                    %>% 
   autonomics.find::filter_n_arrange_top_features(contrast_name, topdef, direction)  %>% 
   autonomics.import::oraid_values()                                                         %>% 
   unique()  # do not double count discriminated isoforms
}

utils::globalVariables('features')
utils::globalVariables('fvars')
utils::globalVariables('key')
utils::globalVariables('fdata')

add_fvars_to_ora <- function(ora_res, object){
   
   fdata1 <- autonomics.import::fdata(object) %>% 
             magrittr::extract(, c('feature_id', 'Uniprot accessions')) %>% 
             data.table::data.table()
   
   ora_res1 <- ora_res %>% 
               data.table::copy() %>%
               tidyr::separate_rows('features', sep = ';')
   ora_res1 %>% utils::head(3)
   
   ora_res1 %<>% merge(fdata1, by.x = 'features', by.y = 'Uniprot accessions')
   ora_res1[, features:=as.character(features)]
   ora_res1[, features   := paste0(features,   collapse = ';'), by = 'pathway']
   
   
   
   merge(fdata1 %>% magrittr::extract(, fvars, with = FALSE), by = key, all.x = TRUE)   %>%
   magrittr::extract(, c(setdiff(names(.), fvars), fvars), with = FALSE)                %>%
   data.table::setorderv(c('rank', key)) %>%
   (function(dt){
      for (curfvar in names(fdata)){
         dt %>% magrittr::extract(, (curfvar) := paste0(get(curfvar), collapse=';'), by = 'pathway')
      }
      dt
   }) %>%
   unique()
}
