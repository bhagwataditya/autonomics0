#' Plot detects
#'
#' Plot number of detects, partial detects, and nondetects (or imputes) per subgroup
#'
#' @param object SummarizedExperiment
#' @param svar   svar
#' @param color_values named character vector (names = colorvar levels, values = color names)
#' @return ggplot object
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'
#'    # Before imputation
#'    object <-  'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'                system.file(package = 'autonomics.data')     %>%
#'                autonomics.import::read_proteingroups()      %>%
#'                autonomics.import::invert(subgroups = c('E_EM', 'E_BM', 'EM_BM'))
#'
#'    object %>% plot_detects_per_subgroup()
#'    object %>% autonomics.import::impute_consistent_nas() %>%
#'               plot_detects_per_subgroup()
#' }
#' @importFrom magrittr %>%
#' @export
plot_detects_per_subgroup <- function(
   object,
   svar = 'subgroup',
   color_values = default_color_values(object, svar)
){
   # Initialize variables
   variable <- subgroup <- value <- NULL

   # Assert
   assertive.sets::assert_is_subset('SummarizedExperiment', class(object))
   assertive.sets::assert_is_subset(c(svar), autonomics.import::svars(object))
   assertive.sets::assert_is_subset(autonomics.import::slevels(object, svar), names(color_values))

   # Prepare datatable
   split_objects  <- object %>% autonomics.import::split_by_svar(svar)
   imputes        <- split_objects %>% vapply(function(x) x %>% autonomics.import::is_imputed() %>% matrixStats::rowAlls() %>% sum(), numeric(1))
   nondetects     <- split_objects %>% vapply(function(x) x %>% autonomics.import::is_na()      %>% matrixStats::rowAlls() %>% sum(), numeric(1))
   partialdetects <- split_objects %>% vapply(function(x) x %>% autonomics.import::is_na()      %>% matrixStats::rowAnys() %>% sum(), numeric(1)) %>%
                     magrittr::subtract(nondetects)
   fulldetects    <- nrow(object) - partialdetects - imputes - nondetects
   plot_dt <- data.table::data.table(subgroup       = names(split_objects) %>% factor(autonomics.import::slevels(object, svar)))
   if (any(fulldetects   !=0))                plot_dt %<>% magrittr::extract(,  fulldetects    := fulldetects)
   if (any(partialdetects!=0))                plot_dt %<>% magrittr::extract(,  partialdetects := partialdetects)
   if (any(nondetects!=0) | any(imputes!=0))  plot_dt %<>% magrittr::extract(,  nondetects     := nondetects + imputes)
   plot_dt %<>% data.table::melt(id.vars = 'subgroup')

   # Set order of variable levels
   variable_levels <- c('nondetects', 'partialdetects', 'fulldetects') %>%
                      magrittr::extract(.%in% plot_dt$variable)
   plot_dt$variable %<>% factor(variable_levels)

   # Plot
   title <- if        ( all(imputes==0) & !all(nondetects==0)){ 'fulldetects  |  partialdetects  |  nondetects'
            } else if (!all(imputes==0) &  all(nondetects==0)){ 'fulldetects  |  partialdetects  |  imputed nondetects'
            } else if ( all(imputes==0) &  all(nondetects==0)){ 'fulldetects  |  partialdetects'}
   ggplot2::ggplot(plot_dt, ggplot2::aes(x=forcats::fct_rev(subgroup), y = value, fill = subgroup, group = variable)) +
   ggplot2::ggtitle(title) + ggplot2::theme_bw() +
   ggplot2::geom_col(color = 'black', position = ggplot2::position_stack()) +
   ggplot2::scale_fill_manual(values = color_values) +
   ggplot2::geom_text(ggplot2::aes(label=value, x = subgroup),
                     position = ggplot2::position_stack(vjust=0.5),
                     size = ggplot2::rel(3)) +
   ggplot2::theme(#axis.text.x        = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.border       = ggplot2::element_blank(),
                 plot.title         = ggplot2::element_text(hjust = 0.5)) +
   ggplot2::xlab(NULL) +
   ggplot2::ylab(NULL) +
   ggplot2::coord_flip()
}

