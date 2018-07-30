#' Create default ggplot colors for factor levels
#' @param factor_levels character vector
#' @return color vector (character)
#' @author John Colby
#' @seealso \href{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}{stackoverflow}
#' @importFrom magrittr  %>%
#' @export
make_gg_colors <- function(factor_levels) {
   n <- length(factor_levels)
   hues <- seq(15, 375, length = n + 1)
   grDevices::hcl(h = hues, l = 65, c = 100)[1:n] %>%
      magrittr::set_names(factor_levels)
}

#' Make composite colors
#' @param object SummarizedExperiment
#' @param color_var svar mapped to color
#' @return named character vector (elements = colors, names = color_var levels)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.plot::make_composite_colors()
#' }
#'
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    object <- autonomics.data::glutaminase
#'    object %>% autonomics.plot::make_composite_colors()
#' }
#' @importFrom magrittr %>%
#' @importFrom data.table   data.table   :=
#' @export
make_composite_colors <- function(
   object,
   color_var = autonomics.plot::default_color_var(object)
){

   # Satisfy CHECK
   subgroup <- V1 <- V2 <- color <- hue <- luminance <- NULL

   components <- object %>% autonomics.import::scomponents(color_var)
   sep        <- object %>% autonomics.import::ssep(color_var)
   components[, subgroup:=paste0(V1,sep,V2)]
   V1levels <- unique(components$V1)
   V2levels <- unique(components$V2)
   n1 <- length(V1levels)
   n2 <- length(V2levels)
   hues       <- data.table::data.table(V1 = V1levels, hue = seq(15, 375, length = n1 + 1)[1:n1])
   luminances <- data.table::data.table(V2 = V2levels, luminance = seq(100, 25, length = n2))
   expand.grid(V1=V1levels, V2=V2levels) %>%
      data.table::data.table() %>%
      merge(luminances, by = 'V2') %>%
      merge(hues, by = 'V1') %>%
      magrittr::extract(, color := grDevices::hcl(h = hue, l = luminance, c = 100), by = c('V1', 'V2')) %>%
      merge(components, by = c('V1', 'V2')) %>%
      magrittr::extract(, color %>% magrittr::set_names(subgroup))
}


# default_color_values <- function(object){
#    NULL
# }



# make_composite_colors <- function(
#    object,
#    color_var = autonomics.plot::default_color_var(object)
# ){
#    # Assert
#    autonomics.import::assert_is_valid_eset(object)
#    assertive.sets::assert_is_subset(color_var, autonomics.import::svars(object, color_var))
#    assertive.base::assert_is_identical_to_true(autonomics.plot::fit_for_composite_coloring(object, color_var))
#
#    # Extract subgroup components
#    svar_components <- object %>% autonomics.import::scomponents(color_var)
#    sep <- autonomics.import::subgroup_sep(object)
#    svar_components[, subgroup:= paste0(V1, sep, V2)]
#
#    # Create colors
#    V1levels <- unique(svar_components$V1)
#    V2levels <- unique(svar_components$V2)
#    all_components <- expand.grid(V1 = V1levels, V2 = V2levels) %>% data.table::data.table()
#
#    palettes <- BREWERPALETTES %>% magrittr::extract(1:length(V1levels))
#    all_components %<>% merge(data.table::data.table(V1 = V1levels, palette = palettes), by = 'V1')
#    all_components %>%  magrittr::extract(,
#                                          color:= RColorBrewer::brewer.pal(length(V2)+1, unique(palette)) %>%
#                                                  magrittr::extract(length(.):2), # first is too light
#                                          by = 'V1')
#    all_components %<>% merge(svar_components, by = c('V1', 'V2'))
#    all_components %>%  magrittr::extract(, color %>% magrittr::set_names(subgroup))
# }

#' Default color values
#' @param object SummarizedExperiment
#' @param color_var color variable
#' @return default color values vector
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% autonomics.plot::default_color_values()
#' }
#'
#' # GLUTAMINASE
#' if (require(autonomics.data)){
#'    autonomics.data::glutaminase %>% autonomics.plot::default_color_values()
#' }
#'
#' # STEM CELL DIFFERENTIATION
#' if (require(autonomics.data)){
#'    autonomics.data::stemdiff.proteinratios %>% autonomics.plot::default_color_values()
#' }
#' @importFrom magrittr %>%
#' @export
default_color_values <- function(
   object,
   color_var = autonomics.plot::default_color_var(object)
){

   # Assert
   autonomics.import::assert_is_valid_eset(object)
   assertive.sets::assert_is_subset(color_var, autonomics.import::svars(object))

   # No color var
   if (is.null(color_var)){
      return(make_gg_colors('default'))
   }

   # Two component subgroups
   if (autonomics.import::svar_has_two_components(object, color_var)){
      autonomics.support::cmessage('\t\tCreating composite colors')
      return(make_composite_colors(object, color_var))
   }

   # Default ggplot colors
   autonomics.support::cmessage('\t\tCreating default ggplot colors')
   object %>% autonomics.import::slevels(color_var) %>% autonomics.plot::make_gg_colors()
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
