#' Plot overlayed sample distributions
#' @param object    SummarizedExperiment
#' @param facet_var svar to facet on
#' @examples
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::exiqon
#'    facet_var <- 'subgroup'
#' }
#' @importFrom magrittr %>%
#' @export
plot_overlayed_sample_distributions <- function(
   object,
   color_var    = autonomics.plot::default_color_var(object),
   color_values = autonomics.plot::default_color_values(object, color_var)
){

   plotdt <- object %>% autonomics.import::sumexp_to_long_dt(svars = color_var)

   ggplot2::ggplot(plotdt, ggplot2::aes_string(x = 'value', group = 'sample_id', color = color_var)) +
   ggplot2::geom_density(na.rm = TRUE) +
   ggplot2::scale_color_manual(values = color_values) +
   ggplot2::geom_hline(yintercept=0, colour="white", size=1)
}
