#' Plot confounder effects
#' @param object exprs object
#' @param fvars  fvars
#' @param color_var   svar mapped to color
#' @param confounders confounding svars
#' @export
#' @importFrom magrittr %>% %<>%
#' @examples
#' require(magrittr)
plot_confounder_effects <- function(
   object,
   confounders = character(),
   fvars,
   color_var   = 'subgroup'
){

   # Limma
   object2 <- object %>% autonomics.find::eliminate_confounders(confounders = confounders)

   # Adjusted values
   nplot <- 2 + length(confounders)
   plotlist <- vector('list', nplot) %>%
               magrittr::set_names(c('original','adjusted', confounders))
   plotlist$adjusted <- object2 %>%
                        autonomics.plot::create_feature_plot_df(fvars = fvars)  %>%
                        ggplot2::ggplot() +
                        ggplot2::facet_wrap(c('feature_facet'), ncol = 1, scales = 'free') +
                        ggplot2::geom_boxplot(ggplot2::aes_string(x = 'subgroup', y = 'value', fill = 'subgroup')) +
                        ggplot2::ylab('adjusted value') +
                        ggplot2::xlab('') +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
                        ggplot2::guides(fill = FALSE)

   # Original values
   plot_df <- autonomics.plot::create_feature_plot_df(object, fvars = fvars)
   p <- ggplot2::ggplot(plot_df) +
        ggplot2::facet_wrap(c('feature_facet'), scales = 'free', ncol = 1) +
        ggplot2::ylab('value') +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
        ggplot2::guides(fill = FALSE, color = FALSE)
   plotlist$original <- p + ggplot2::geom_boxplot(ggplot2::aes_string(x = 'subgroup', y = 'value', fill = 'subgroup')) + ggplot2::xlab('')

   # Confounder effects
   for (cur_confounder in confounders){
      plotlist[[cur_confounder]] <- if (is.numeric(object[[cur_confounder]])){
                                       p + ggplot2::geom_point(  ggplot2::aes_string(x = cur_confounder, y = 'value', color = color_var))
                                    } else {
                                       p + ggplot2::geom_boxplot(ggplot2::aes_string(x = cur_confounder, y = 'value', fill  = color_var))
                                    }
   }
   autonomics.plot::combine_plots(plotlist = plotlist, cols = nplot)

}
