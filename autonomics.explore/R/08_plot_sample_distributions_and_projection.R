#' Plot sample distributions and pca
#' @param object            eset
#' @param descr               description: used as plot title and plot name
#' @param color_var           svar mapped to color
#' @param color_values        color vector (names = color_var levels, values = colors)
#' @param shape_var           svar mapped to shape
#' @param txt_var             svar mapped to txt
#' @param displayed_features  features to be displayed in the sample distributions (vector of numeric indexes or character \code{feature_id}s)
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    billing.differentiation.data::rna.voomcounts %>%
#'       autonomics.explore::plot_sample_distrs_n_pca('rna.voomcounts')
#' }
#' if (require(halama.2016)){
#'    halama.2016::cell.metabolites %>% autonomics.explore::plot_sample_distrs_n_pca(
#'       descr = 'custom colors',
#'       color_var    = 'GROUP_DESCRIPTION',
#'       color_values = c(Control = 'orange', Vehicle = 'red', `Concentration 1` = 'green',
#'                       `Concentration 2` = 'blue'))
#' }
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% 
#'    autonomics.explore::plot_sample_distrs_n_pca()
#' }
#' @return list(descr.distr = plot1, descr.pca = plot2)
#' @importFrom magrittr   %>%
#' @export
plot_sample_distrs_n_pca <- function(
   object,
   descr = '',
   color_var    = autonomics.plot::default_color_var(object),
   color_values = autonomics.plot::default_color_values(object),
   shape_var    = autonomics.plot::default_shape_var(object),
   txt_var      = autonomics.plot::default_txt_var(object),
   displayed_features = NULL
){
   plotlist <- vector(mode = 'list', length = 2) %>% magrittr::set_names(paste0(descr, c('.distr', '.pca')))
   plotlist[[1]] <- object %>% magrittr::extract(, 1:min(25, ncol(.))) %>%
                                 autonomics.plot::plot_sample_distributions(
                                    title              = descr,
                                    color_var          = color_var,
                                    color_values       = color_values,
                                    displayed_features = displayed_features)

   plotlist[[2]] <- object %>% autonomics.explore::plot_pca_samples(
                                    title        = descr,
                                    color_var    = color_var,
                                    color_values = color_values,
                                    shape_var    = shape_var,
                                    txt_var      = txt_var)
   plotlist
}
