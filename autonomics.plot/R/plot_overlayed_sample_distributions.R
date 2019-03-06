#' Plot overlayed sample distributions
#' @param object        SummarizedExperiment
#' @param color_var     string: svar mapped to color
#' @param color_values  string vector: values = colors, names = color_var levels
#' @return ggplot2 object
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% plot_overlayed_sample_distributions()
#' }
#' @importFrom magrittr %>%
#' @export
plot_overlayed_sample_distributions <- function(
   object,
   color_var    = default_color_var(object),
   color_values = default_color_values(object, color_var)
){

   plotdt <- object %>% autonomics.import::sumexp_to_long_dt(svars = color_var)

   ggplot2::ggplot(plotdt, ggplot2::aes_string(x = 'value', group = 'sample_id', color = color_var)) +
   ggplot2::geom_density(na.rm = TRUE) +
   ggplot2::scale_color_manual(values = color_values) +
   ggplot2::geom_hline(yintercept=0, colour="white", size=1)
}

plot_overlayed_feature_distributions <- function(
   object,
   color_var    = default_color_var(object),
   color_values = default_color_values(object, color_var)
){

   plotdt <- object[1:10, ] %>% autonomics.import::sumexp_to_long_dt()
   ggplot2::ggplot(plotdt, ggplot2::aes_string(x = 'value', group = 'feature_id', color = 'feature_id')) +
   ggplot2::geom_density(na.rm=TRUE)

}
