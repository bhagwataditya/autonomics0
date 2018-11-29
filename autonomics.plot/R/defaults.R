#' Create default ggplot colors for factor levels
#' @param factor_levels character vector
#' @param show logical
#' @return color vector (character)
#' @author John Colby
#' @seealso \href{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}{stackoverflow}
#' @importFrom magrittr  %>%
#' @export
make_gg_colors <- function(factor_levels, show) {
   n <- length(factor_levels)
   hues <- seq(15, 375, length = n + 1)
   color_levels <- grDevices::hcl(h = hues, l = 65, c = 100)[1:n] %>%
                   magrittr::set_names(factor_levels)
   if (show) color_levels %>% (function(x) graphics::pie(rep(1, length(x)), names(x), col = x))
   color_levels
}

#' Make composite colors
#' @param svalues character vector
#' @param show    logical(1): whether to show colors in pie plot
#' @return named character vector (elements = colors, names = color_var levels)
#' @examples
#' require(magrittr)
#' svalues <- c('inflammation.all', 'inflammation.C2_C0', 'inflammation.D2_D0',
#'              'oxidstress.all',   'oxidstress.C2_C0',   'oxidstress.D2_D0')
#' svalues %>% autonomics.plot::make_composite_colors(show = TRUE)
#'
#' if (require(subramanian.2016)){
#'    svalues <- subramanian.2016::metabolon %>%
#'               autonomics.import::subgroup_levels()
#'    svalues %>% autonomics.plot::make_composite_colors(show = TRUE)
#' }
#'
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase %>%
#'    autonomics.import::subgroup_levels() %>%
#'    autonomics.plot::make_composite_colors(show = TRUE)
#' }
#' @importFrom magrittr %>%
#' @importFrom data.table   data.table   :=
#' @export
make_composite_colors <- function(
   svalues,
   sep  = autonomics.import::guess_sep(svalues),
   show = FALSE
){

   # Assert
   assertive.properties::assert_is_not_null(sep)

   # Satisfy CHECK
   subgroup <- V1 <- V2 <- color <- hue <- luminance <- NULL

   # Split into components
   #    * V1: first n-1 components => will be mapped to hue
   #    * V2: last component       => will be mapped to luminance
   # This approach works also when more than two components are present
   # It is therefore used instead of autonomics.import::split_values()
   V1 <- svalues %>% stringi::stri_split_fixed(sep) %>% vapply(function(x) x %>% magrittr::extract(-length(x)) %>% paste0(collapse = sep), character(1))
   V2 <- svalues %>% stringi::stri_split_fixed(sep) %>% vapply(function(x) x %>% magrittr::extract( length(x)),                            character(1))
   V1levels <- sort(unique(V1))
   V2levels <- sort(unique(V2))
   n1 <- length(V1levels)
   n2 <- length(V2levels)
   hues <- seq(15, 375, length = n1 + 1)[1:n1] %>% magrittr::set_names(V1levels)

   color_levels <- character(0)
   for (i in seq_along(hues)){
      color_levels %<>% c(colorspace::sequential_hcl(n2, h = hues[[i]], power = 1, c = c(50, 100), l = c(90, 30)) %>%
                             magrittr::set_names(paste0(V1levels[[i]], sep, V2levels)))
   }
   if (show) color_levels %>% (function(x) graphics::pie(rep(1, length(x)), names(x), col = x))

   return(color_levels)
}

#' Make fitting colors
#' @param x        colorvar levels vector
#' @param show     logical(1)
#' @param verbose  logical(1)
#' @return color vectors (values = colors, names = colorvar levels)
#' @export
make_colors <- function(
   x,
   sep = autonomics.import::guess_sep(x, verbose = FALSE),
   show = FALSE,
   verbose = FALSE
){

   # 0D colors
   if (is.null(x)){
      return(autonomics.plot::make_gg_colors('default'))
   }

   # 1D colors
   if (is.null(sep)){
      if (verbose) autonomics.support::cmessage('\t\tMake default ggplot colors')
      return(autonomics.plot::make_gg_colors(x, show = show))

   # 2D colors
   } else {
      if (verbose) autonomics.support::cmessage('\t\tMake composite colors')
      return(autonomics.plot::make_composite_colors(x, sep = sep, show = show))
   }


}

#' Default color values
#' @param object SummarizedExperiment
#' @param color_var color variable
#' @param show logical
#' @return default color values vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>%
#'    autonomics.plot::default_color_values(show = TRUE)
#' }
#'
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>%
#'    autonomics.plot::default_color_values(show=TRUE)
#' }
#'
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase %>%
#'    autonomics.plot::default_color_values(show=TRUE)
#' }
#'
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    autonomics.data::stemdiff.proteinratios %>%
#'    autonomics.plot::default_color_values(show=TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
default_color_values <- function(
   object,
   color_var = autonomics.plot::default_color_var(object),
   show = FALSE
){

   # Assert
   autonomics.import::assert_is_valid_eset(object)
   assertive.sets::assert_is_subset(color_var, autonomics.import::svars(object))

   # sep
   sep              <- object %>% autonomics.import::guess_sep(color_var)
   color_var_levels <- object %>% autonomics.import::slevels(color_var)

   # Make colors
   autonomics.plot::make_colors(color_var_levels, sep, show)

}

#' default fvars
#' @param  object eset
#' @return character vector
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::billing2016 %>% autonomics.plot::default_fvars()
#' }
#' if (require(atkin.2014)){
#'    atkin.2014::soma %>% autonomics.plot::default_fvars()
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.counts
#' }
#' @importFrom magrittr               %>%
#' @export
default_fvars <- function(object){

   if (autonomics.import::is_metabolon_eset(object)){
      return(c('BIOCHEMICAL', 'SUB_PATHWAY'))
   } else {
      fidvar   <- object %>% autonomics.import::fid_var()
      fnamevar <- object %>% autonomics.import::fname_var()
      tmp <- c(fnamevar, fidvar)
      if (assertive.properties::is_empty(tmp)) NULL else tmp
   }
}

#' default feature_plots
#' @param object eset
#' @return default value of feature_plots
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>%
#'    autonomics.plot::default_feature_plots()
#'
#'    autonomics.data::stemcomp.soma %>%
#'    autonomics.plot::default_feature_plots()
#' }
#'
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    autonomics.data::stemdiff.proteinratios %>%
#'    autonomics.plot::default_feature_plots()
#' }
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>%
#'       default_feature_plots()
#' }
#'
#' @importFrom magrittr %>%
#' @export
default_feature_plots <- function(object){

   single_subgroup            <- length(unique(autonomics.import::sdata(object)$subgroup)) == 1
   single_replicate           <- ncol(object) == 1
   upto_3_singleton_subgroups <- all(table(autonomics.import::sdata(object)$subgroup)==1)  &  (length(unique(autonomics.import::sdata(object)$subgroup)) <= 4)
   some_singletons            <- any(table(autonomics.import::sdata(object)$subgroup)==1)
   block_present              <- 'block' %in% autonomics.import::svars(object)
   each_subgroup_has_three_plus_replicates <- autonomics.import::subgroup_values(object) %>%
                                              table() %>%
                                              magrittr::is_greater_than(3) %>%
                                              all()

   if (single_subgroup & single_replicate)            return('hbar')
   if (single_subgroup | upto_3_singleton_subgroups)  return('bar')
   if (some_singletons | block_present)               return('point')
   if (each_subgroup_has_three_plus_replicates)       return('boxplot')
   else                                               return('point')
}

#' default color_var
#' @param object eset
#' @return default value of color_var
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% autonomics.plot::default_color_var()
#'    autonomics.data::stemcomp.soma          %>% autonomics.plot::default_color_var()
#'}
#'
#'# STEM CELL DIFFERENTIATION
#'if (require(billing.differentiation.data)){
#'   billing.differentiation.data::rna.voomcounts %>%
#'      default_color_var()
#'}
#' # ATKIN 2014
#'if (require(atkin.2014)){
#'   atkin.2014::soma %>% autonomics.plot::default_color_var()
#'}
#'
#' @export
default_color_var <- function(object){
   if (     'block'    %in% autonomics.import::svars(object))  'block'
   else if ('subgroup' %in% autonomics.import::svars(object))  'subgroup'
   else                                                         NULL
}

#' Default group var
#' @param object   eset
#' @return default value of group_var
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% autonomics.plot::default_group_var()
#'    autonomics.data::billing2016 %>% autonomics.plot::default_group_var()
#'}
#'
#' # ATKIN 2014
#'if (require(atkin.2014)){
#'   atkin.2014::soma %>% autonomics.plot::default_group_var()
#'}
#'
#'# STEM CELL DIFFERENTIATION
#'if (require(billing.differentiation.data)){
#'   billing.differentiation.data::rna.voomcounts %>%
#'      default_group_var()
#'}
#' @export
default_group_var <- function(object){
   if ('block' %in% autonomics.import::svars(object))   'block'
   else                                                  1
}

#' Default facet_var value
#' @export
default_facet_var <- function(){
   NULL
}

#' Default line value
#' @param object   eset
#' @return default line value
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% autonomics.plot::default_line()
#'}
#'
#'# ATKIN.2014
#'if (require(atkin.2014)){
#'   atkin.2014::soma %>% autonomics.plot::default_line()
#'}
#'
#'# STEM CELL DIFFERENTIATION
#'if (require(billing.differentiation.data)){
#'   billing.differentiation.data::rna.voomcounts %>%
#'      default_line()
#'}
#' @export
default_line <- function(object){
   if ('block' %in% autonomics.import::svars(object))  TRUE
   else                                                FALSE
}

#' default color_var
#' @param object eset
#' @return default value of color_var
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% autonomics.plot::default_shape_var()
#'}
#'
#'# STEM CELL DIFFERENTIATION
#'if (require(billing.differentiation.data)){
#'   billing.differentiation.data::rna.voomcounts %>% default_shape_var()
#'}
#' @export
default_shape_var <- function(object){
   if('replicate' %in% autonomics.import::svars(object)){
      nreplicate <- length(unique(autonomics.import::sdata(object)$replicate))
      if (nreplicate > 1  &  nreplicate <= 3){
         return('replicate')
      }
   }
   return(NULL)
}

#' Default x
#' @param object      eset
#' @param feature_plots  feature plots
#' @return default value of x
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% autonomics.plot::default_x()
#'}
#'
#'# STEM CELL DIFFERENTIATION
#'if (require(billing.differentiation.data)){
#'   billing.differentiation.data::rna.voomcounts %>% default_x()
#'}
#' @export
default_x <- function(object, feature_plots = default_feature_plots(object)){

   autonomics.import::assert_is_valid_eset(object)
   assertive.sets::assert_is_subset(feature_plots, FEATURE_PLOTS)

   sample_id_in_sdata <- 'sample_id' %in% autonomics.import::svars(object)
   subgroup_in_sdata  <- 'subgroup'  %in% autonomics.import::svars(object)

   default_x_single_plot <- function(cur_plot){

      if (cur_plot == 'bar' & sample_id_in_sdata){
         return('sample_id')

      } else if (cur_plot %in% c('violin', 'boxplot')){
         assertive.sets::assert_is_subset('subgroup', autonomics.import::svars(object))
         return('subgroup')

      } else if (cur_plot == 'point' & subgroup_in_sdata){
         return('subgroup')

      } else {
         return('snames')
      }
   }

   Map(default_x_single_plot, feature_plots) %>% unlist()

}

#' Default txt var
#' @param object eset
#' @return default value of txt_var
#' @export
default_txt_var <- function(object){
   NULL
}

#' Default zero hline
#' @param   object   eset
#' @return  logical
#' @examples
#' require(magrittr)
#'
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>%
#'    autonomics.plot::default_zero_hline()
#' }
#'
#' # STEM CELL DIFFERENTIATION
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>%
#'       autonomics.plot::default_zero_hline()
#' }
#' @export
default_zero_hline <- function(object){
   if (autonomics.import::contains_ratios(object))  TRUE
   else                                               FALSE
}
