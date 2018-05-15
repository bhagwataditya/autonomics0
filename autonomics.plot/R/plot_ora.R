utils::globalVariables(c('collection', 'direction', 'contrast', 'pathway', 'p'))

#' Plot ORA results
#'
#' @param ora_results    \code{\link{data.frame}} as (silently) returned from \code{autonomics.ora::run_ora_on_eset}
#' @param top_n          number of factors chosen for plotting (by rank)
#' @param collections    (sub)ontologies/collections to plot
#' @param directions     what coefficient/effect direction of a given contrast the set tested is chosen from
#' @param contrasts      contrasts to filter on
#' @param facet_var_x    variable name to \code{\link[ggplot2]{facet_grid}} on
#' @param facet_var_y    variable name to \code{\link[ggplot2]{facet_grid}} on
#' @return ggplot2 object
#' @author Aditya Bhagwat, Johannes Graumann
#' @importFrom  magrittr  %>%
#' @export
plot_ora <- function(
   ora_results,
   top_n = 5,
   collections = c('gobp', 'gocc', 'gomf', 'kegg'),
   directions = c('neg', 'pos', 'both'),
   contrasts  = NULL,
   facet_var_x = 'contrast',
   facet_var_y = 'direction')
{
# Check prerequisites & input processing ----------------------------------
   ora_results %>%
      assertive.types::assert_is_data.frame() %>%
      names() %>%
      assertive.sets::assert_are_set_equal(
         c("rank", "p", "pathway", "n_total", "n_detected", "n_selected",
           "features", "direction", "contrast", "collection"))

   top_n %>%
      assertive.types::assert_is_a_number() %>%
      assertive.numbers::assert_all_are_whole_numbers() %>%
      assertive.numbers::assert_all_are_positive()

   collections %>%
      assertive.sets::assert_are_intersecting_sets(ora_results$collection)

   directions %>%
      assertive.sets::assert_are_intersecting_sets(ora_results$direction)

   if(!is.null(contrasts))
   {
      contrasts %>%
        assertive.sets::assert_are_intersecting_sets(ora_results$contrast)
   }

   if(!is.null(facet_var_x))
   {
      facet_var_x %>%
         assertive.sets::assert_is_subset(names(ora_results))
   } else {
      facet_var_x <- '.'
   }
   if(!is.null(facet_var_y))
   {
      facet_var_y %>%
         assertive.sets::assert_is_subset(names(ora_results))
   } else {
      facet_var_y <- '.'
   }

# Processing of data set --------------------------------------------------
   ora_results %<>%
      dplyr::filter(
        rank <= top_n,
        collection %in% collections,
        direction %in% directions) %>%
      dplyr::mutate(
         direction = dplyr::case_when(
               direction == 'both' ~ 'Both',
               direction == 'neg'  ~ 'Negativ',
               direction == 'pos'  ~ 'Positive',
               TRUE ~ direction) %>%
            as.factor() %>%
            factor(levels = rev(levels(.))))
   if(!is.null(contrasts))
   {
      ora_results %<>%
         dplyr::filter(
            contrast %in% contrasts)
   }

# Plotting ----------------------------------------------------------------
   tmp_plot <- ora_results %>%
      ggplot2::ggplot(
         ggplot2::aes(x = pathway, y = -log10(p))) +
      ggplot2::geom_bar(stat = 'identity') +
      ggplot2::coord_flip() +
      ggplot2::labs(
         x = 'Pathway ID',
         y = expression(-log[10](p))
      )
   if(
      c(facet_var_x, facet_var_y) %>%
        stringi::stri_detect_regex(
          pattern = '^\\.$',
          negate  = TRUE) %>%
        any())
   {
      faceting_term <- paste(facet_var_x, facet_var_y, sep = ' ~ ')
      tmp_plot <- tmp_plot +
         ggplot2::facet_grid(
            faceting_term,
            scales = 'free_y',
            space  = 'free_y')
   }
   tmp_plot
}
