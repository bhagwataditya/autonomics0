#=============================================

#' Contrast values in character vector
#' @param x character(n)
#' @param subtractor character(1)
#' @param parenthesize logical(1): whether to parenthesize subtrahend
#' @return character(n-1)
#' @examples
#' require(magrittr)
#' letters[1:5] %>% contrast_strings()
#' letters[1:5] %>% contrast_strings(parenthesize = TRUE)
#' @importFrom magrittr %<>%
#' @export 
contrast_strings <- function(x, subtractor = ' - ', parenthesize = FALSE){
   assertive.types::assert_is_character(x)
   x %<>% magrittr::extract(.!='')
   n <- length(x)
   minuend    <- x[-1]
   subtrahend <- x[-n]
   if (parenthesize){
      minuend    %<>% paste0('(', ., ')')
      subtrahend %<>% paste0('(', ., ')')
   }
   sprintf('%s%s%s', minuend, subtractor, subtrahend)
}

#' Contrast subgroups. Contrast contrasts.
#' 
#' @param x SummarizedExperiment, matrix, or character(.)
#' @return character(.)
#' @examples
#' require(magrittr)
#' 
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    x <- autonomics.data::stemdiff.proteinratios
#'    x %>% autonomics.find::contrast_subgroups()
#' }
#' 
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#' 
#'    # SummarizedExperiment
#'    x <- autonomics.data::glutaminase
#'    x %>% contrast_subgroups()
#'    x %>% contrast_contrasts()
#'    
#'    # Matrix
#'    x <- autonomics.data::glutaminase %>% autonomics.find::layout_composite_values()
#'    x %>% contrast_subgroups()
#'    x %>% contrast_contrasts()
#'    
#'    # Vector
#'    x <- autonomics.data::glutaminase %>% autonomics.find::layout_composite_values() %>% magrittr::extract(1,)
#'    x %>% contrast_subgroups()
#' }
#' 
#' # ATKIN
#' if (require(atkin.2014)){
#' 
#'    # SummarizedExperiment
#'    x <- atkin.2014::soma.2018
#'    x %>% contrast_subgroups()
#'    x %>% contrast_contrasts()
#'            
#'  # matrix 
#'    x <- atkin.2014::soma.2018 %>% autonomics.find::layout_composite_values()
#'    x %>% contrast_subgroups()
#'    x %>% contrast_contrasts()
#'    
#'  # character vector
#'    x <- atkin.2014::soma.2018 %>% autonomics.find::layout_composite_values() %>% 
#'                        magrittr::extract(1, )
#'    x %>% contrast_subgroups()
#' }
#'
#'@export
contrast_subgroups <- function (x, ...) {
   UseMethod("contrast_subgroups", x)
}


#' @rdname contrast_subgroups
#' @importFrom magrittr %>%  %<>% 
#' @export 
contrast_subgroups.character <- function(x){
   x %<>% magrittr::extract(.!='')
   n <- length(x)
   contrast_strings(x) %>% magrittr::set_names(contrast_strings(x, '_'))
}



#' @rdname contrast_subgroups
#' @importFrom magrittr %>% 
#' @export
contrast_subgroups.matrix <- function(x){
   
   # Drop columns with missing values
   idx <- matrixStats::colAnys(x == '')
   if (any(idx)){
      autonomics.support::cmessage('\t\tDrop column(s) with missing values: %s', colnames(x)[idx] %>% paste0(collapse = ', '))
      x %<>% magrittr::extract(, !idx)
   }

   per_row_contrasts  <- x %>% apply(1, contrast_subgroups) %>% t()
   dim(per_row_contrasts) <- c(nrow(x), ncol(x)-1)
   rownames(per_row_contrasts) <- rownames(x)
   colnames(per_row_contrasts) <- colnames(x) %>% contrast_strings('_')
   
   per_col_contrasts <- x %>%  apply(2, contrast_strings)
   dim(per_col_contrasts) <- c(nrow(x) - 1, ncol(x))
   rownames(per_col_contrasts) <- rownames(x) %>% contrast_strings('_')
   colnames(per_col_contrasts) <- colnames(x)
   
   list(per_row_contrasts, per_col_contrasts)
}

#' @rdname contrast_subgroups
#' @importFrom magrittr %>% 
#' @export
contrast_subgroups.SummarizedExperiment <- function(x){
   x %>% 
   autonomics.find::layout_composite_values() %>% 
   contrast_subgroups()
}

#' @rdname contrast_subgroups
#' @export
contrast_contrasts <- function (x, ...) {
   UseMethod("contrast_contrasts", x)
}

#' @rdname contrast_subgroups
#' @importFrom magrittr %>% 
#' @export
contrast_contrasts.matrix <- function(x){
   y <- x %>% contrast_subgroups()
   per_row_contrasts <- y[[1]]
   per_col_contrasts <- y[[2]]
   
   contrast_differences <- if (nrow(per_row_contrasts)>1)   per_row_contrasts %>% apply(2, contrast_strings, parenthesize = TRUE)  else  NULL
   dim(contrast_differences) <- c(nrow(per_row_contrasts)-1, ncol(per_row_contrasts))
   rownames(contrast_differences) <- sprintf('%s_%s', rownames(per_row_contrasts)[-1], rownames(per_row_contrasts)[-nrow(per_row_contrasts)])
   colnames(contrast_differences) <- colnames(per_row_contrasts)

   contrast_differences
}

#' @rdname contrast_subgroups
#' @importFrom magrittr %>% 
#' @export
contrast_contrasts.SummarizedExperiment <- function(x){
   x %>% autonomics.find::layout_composite_values() %>% contrast_contrasts()
}




#=====


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
#'    dt  <- object %>% autonomics.import::split_values(keep = FALSE)
#'    sep <- object %>% autonomics.import::guess_sep()
#'    autonomics.find::make_ref_contrasts_within_stratum(dt, sep, 1)
#'    
#' }
#' 
#' if (require(subramanian.2016)){
#'    dt  <- subramanian.2016::metabolon %>% autonomics.import::split_values(keep = FALSE)
#'    sep <- subramanian.2016::metabolon %>% autonomics.import::guess_sep()
#'    dt %<>% magrittr::extract(x2 == 'w08')
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
#'    dt  <- object %>% autonomics.import::split_values(keep = FALSE)
#'    sep <- object %>% autonomics.import::guess_sep()
#'    autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)
#' }
#' 
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    dt  <- object %>% autonomics.import::split_values(keep = FALSE)
#'    sep <- object %>% autonomics.import::guess_sep()
#'    autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)
#'    autonomics.find::make_ref_contrasts_across_strata(2, dt, sep)
#' }
#' 
#' if (require(subramanian.2016)){
#'    dt  <- subramanian.2016::metabolon  %>%  autonomics.import::split_values(keep_original = FALSE)
#'    sep <- subramanian.2016::metabolon  %>%  autonomics.import::subgroup_sep()
#'    \dontrun{autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)}
#'    \dontrun{autonomics.find::make_ref_contrasts_across_strata(2, dt, sep)}
#' }
#' 
#' if (require(graumann.lfq)){
#'    dt <- graumann.lfq::lfq.intensities %>% autonomics.import::split_values(keep = FALSE)
#'    sep <- graumann.lfq::lfq.intensities %>% autonomics.import::guess_sep()
#'    autonomics.find::make_ref_contrasts_across_strata(1, dt, sep)
#'    \dontrun{autonomics.find::make_ref_contrasts_across_strata(2, dt, sep)}
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
   dt  <- object  %>%  autonomics.import::split_values(keep = FALSE)
   sep <- object  %>%  autonomics.import::guess_sep()
   1:ncol(dt) %>% lapply(autonomics.find::make_ref_contrasts_across_strata, dt, sep) %>% 
                         data.table::rbindlist() %>% 
                        (function(x){
                           if (nrow(x)==0){  character(0)
                           } else {          x$contrast %>% magrittr::set_names(x$name)
                           }
                        })
}
