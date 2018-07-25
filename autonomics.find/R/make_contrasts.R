
#' Make ref contrasts within stratum
#' 
#' Make reference contrasts for selected variable within stratum of other variables
#' @param dt         subgroup data.table for single stratum
#' @param sep        subgroup separator
#' @param component  subgroup component for which to define contrasts
#' @param ref_levels reference levels (character vector)
#' @examples
#' require(magrittr) 
#' 
#' # STEMCELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    dt  <- object %>% autonomics.import::subgroup_components()
#'    sep <- object %>% autonomics.import::subgroup_sep()
#'    autonomics.find::make_ref_contrasts_within_stratum(dt, sep, 1)
#'    
#' }
#' 
#' if (require(subramanian.2016)){
#'    dt  <- subramanian.2016::metabolon %>% autonomics.import::subgroup_components()
#'    sep <- subramanian.2016::metabolon %>% autonomics.import::subgroup_sep()
#'    dt %<>% magrittr::extract(V2 == 'w08')
#'    autonomics.find::make_ref_contrasts_within_stratum(dt, sep) %>% print()
#' }
#' @importFrom data.table   data.table
#' @importFrom magrittr     %>% 
#' @export
make_ref_contrasts_within_stratum <- function(
   dt, 
   sep,
   component  = 1,
   ref_levels = vapply(dt, magrittr::extract, character(1), 1)
){
   ref_level <- ref_levels[[component]]
   ref_dt    <- dt %>% magrittr::extract(dt[[component]] == ref_level)
   other_dt  <- dt %>% magrittr::extract(dt[[component]] != ref_level)
   paste_cols <- if (is.null(sep)){   function(dt) dt %>% unlist() %>% unname()
                 } else {             function(dt) do.call(function(...) paste(..., sep = sep), dt)
                 }
   
   if (nrow(ref_dt)>0 & nrow(other_dt)>0){
      name  <- sprintf('%s_%s',    other_dt %>% magrittr::extract2(component),  ref_dt %>% magrittr::extract2(component))
      term  <- sprintf('%s - %s',  other_dt %>% paste_cols(),                   ref_dt %>% paste_cols())
   } else {
      name <- term <- character(0)
   }
   data.table::data.table(name = name, term = term)
}

#' Aggregate strata
#' @param  x character vector
#' @return string 
#' @importFrom magrittr %>% 
#' @export
aggregate_strata <- function(x){
   sprintf('(%s)/%d', x, length(x)) %>% paste0(collapse = ' + ')
}

#' Make reference contrasts across strata
#' 
#' Make ref contrasts for selected variable across strata of other variables
#' @param component  subgroup component for which to formulate contrasts
#' @param dt         data.table with all subgroup components
#' @param sep        subgroup separator
#' @return data.table(name = contrast name, contrast = contrast formula)
#' @examples
#' require(magrittr)
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    dt  <- object %>% autonomics.import::subgroup_components()
#'    sep <- object %>% autonomics.import::subgroup_sep()
#'    autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)
#' }
#' 
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    dt  <- object %>% autonomics.import::subgroup_components()
#'    sep <- object %>% autonomics.import::subgroup_sep()
#'    autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)
#'    autonomics.find::make_ref_contrasts_across_strata(2, dt, sep)
#' }
#' 
#' if (require(subramanian.2016)){
#'    dt  <- subramanian.2016::metabolon  %>%  autonomics.import::subgroup_components()
#'    sep <- subramanian.2016::metabolon  %>%  autonomics.import::subgroup_sep()
#'    autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)
#'    autonomics.find::make_ref_contrasts_across_strata(2, dt, sep)
#' }
#' 
#' if (require(graumann.lfq)){
#'    dt <- graumann.lfq::lfq.intensities %>% autonomics.import::subgroup_components()
#'    sep <- graumann.lfq::lfq.intensities %>% autonomics.import::subgroup_sep()
#'    autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)
#'    autonomics.find::make_ref_contrasts_across_strata(2, dt, sep)
#' }
#' @importFrom data.table  data.table
#' @importFrom magrittr    %>%
#' @export
make_ref_contrasts_across_strata <- function(component, dt, sep){
   .SD <- term <- name <- NULL
   
   extract_ref_level <- function(x) as.character(x) %>% magrittr::extract(1)
   ref_levels <- dt %>% vapply(extract_ref_level, character(1))
   dt %>% magrittr::extract(, autonomics.find::make_ref_contrasts_within_stratum(.SD, sep, component, ref_levels), 
                              by = eval(names(dt)[-component]), 
                             .SDcols = names(dt)) %>% 
          magrittr::extract(, list(contrast = autonomics.find::aggregate_strata(term)), 
                              by = 'name') #%>% 
          #magrittr::extract(order(name))
}

#' Make reference contrasts
#' @param  object SummarizedExperiment, eSet, or EList
#' @return named character vector
#' @examples
#' require(magrittr)
#' 
#' # STEMCELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    object %>% autonomics.find::make_ref_contrasts() %>% data.frame()
#' }
#' 
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.find::make_ref_contrasts() %>% data.frame()
#' }
#' 
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>%  
#'    autonomics.find::make_ref_contrasts() %>% data.frame()
#' }
#' @importFrom  data.table  data.table
#' @importFrom  magrittr    %>% 
#' @export
make_ref_contrasts <- function(object){
   dt  <- object  %>%  autonomics.import::subgroup_components()
   sep <- object  %>%  autonomics.import::subgroup_sep()
   1:ncol(dt) %>% lapply(autonomics.find::make_ref_contrasts_across_strata, dt, sep) %>% 
                         data.table::rbindlist() %>% 
                        (function(x){
                           if (nrow(x)==0){  character(0)
                           } else {          x$contrast %>% magrittr::set_names(x$name)
                           }
                        })
}
