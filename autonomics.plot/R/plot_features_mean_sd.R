#' Mean/SD Plot(s)
#'
#' Mean/variance plots are helpful for the evaluation of homoscedasticity.
#' @param object              eSet, EList or SummarizedExperiment
#' @param ...                 additonal eSets, ELists, SummarizedExperiments
#' @param exponent            exponent with which the standard deviation is plotted;
#' \code{2}: $sigma^2 = 'Varianve'$; \code{1}: $sigma^1 = 'Standard Deviation'$;
#' \code{0.5}: $sigma^0.5 = sqrt(sigma)$ (default; as used by Law et al (2014));
#' @param complete_data       \code{\link{logical}} indicating whether to plot only complete data without
#' \code{\link{NA}} in any replicate;
#' @param facet_labels        Labels used for \code{object} and \code{...} by \code{\link[ggplot2]{facet_grid}}
#' @param facet_grid_labeller Labeller (see \code{\link[ggplot2]{labellers}}) for facetting labels
#' @param alpha        handed on to \code{\link[ggplot2]{geom_point}}
#' @references Law, C.W., Chen, Y., Shi, W., and Smyth, G.K. (2014). voom: precision
#' weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29.
#' @return \code{ggplot2} object
#' @export
#' @examples
#' \dontrun{
#' # Works, but is rather slow
#' if(require(autonomics.data)){
#'    require(magrittr)
#'   # Simple plot
#'   autonomics.plot::plot_features_mean_sd(autonomics.data::ALL)
#'   # Compare with the loged version of the data
#'   logALL <- autonomics.data::ALL
#'   autonomics.import::exprs(logALL) %<>% log2()
#'   plot_features_mean_sd(logALL)
#'   # Combined plot
#'   plot_features_mean_sd(autonomics.data::ALL, logALL)
#'   # Same thing, but with descriptive facet labels
#'   plot_features_mean_sd(
#'     autonomics.data::ALL,
#'     logALL,
#'     facet_labels = c("'Raw Data'", "'Logarithmized Data'~(log[2])"))
#' }
#' }
plot_features_mean_sd <- function(
   object,
   ...,
   exponent            = c(0.5, 1, 2),
   complete_data       = TRUE,
   facet_labels        = NULL,
   facet_grid_labeller = 'label_parsed',
   alpha               = 0.1)
{
   obj_list <- list(object) %>%
      c(list(...))

   for(obj in obj_list)
   {
      obj %>%
         assertive.types::assert_is_any_of(c('eSet', 'EList', 'SummarizedExperiment'))
   }

   exponent %<>%
      magrittr::extract2(1) %>%
      as.character() %>%
      match.arg(
         choices    = c(0.5, 1, 2),
         several.ok = FALSE)

   complete_data %>%
      assertive.types::assert_is_a_bool()

   if(is.null(facet_labels))
   {
      facet_labels <- paste(
         'Data',
         'Set',
         seq_along(obj_list),
         sep = '~')
   }

   facet_labels %>%
      assertive.types::assert_is_character() %>%
      assertive.properties::assert_are_same_length(obj_list)

   facet_grid_labeller %>%
      assertive.types::assert_is_a_string() %>%
      assertive.strings::assert_all_are_matching_regex('^label_') %>%
      exists(where=asNamespace('ggplot2'), mode='function') %>%
      assertive.base::assert_all_are_true()

   alpha %>%
      assertive.types::assert_is_a_number() %>%
      assertive.numbers::assert_all_are_in_closed_range(lower = 0, upper = 1)

   obj_list %<>%
      lapply(autonomics.preprocess::mean_and_sd,MARGIN = 1, by_subgroup = TRUE) %>%
      lapply(autonomics.import::fdata)

  mn_obj_list <- obj_list %>%
      lapply(dplyr::select, dplyr::ends_with('__mean')) %>%
      lapply(tidyr::gather, key = 'subgroup', value = 'mean') %>%
      lapply(dplyr::mutate, subgroup = stringi::stri_replace_all_fixed(subgroup, pattern = '__mean', replacement = ''))

  sd_obj_list <- obj_list %>%
     lapply(dplyr::select, dplyr::ends_with('__sd')) %>%
     lapply(tidyr::gather, key = 'subgroup', value = 'sd') %>%
     lapply(dplyr::mutate, subgroup = stringi::stri_replace_all_fixed(subgroup, pattern = '__sd', replacement = ''))

  cmpl_obj_list <- obj_list %>%
     lapply(dplyr::select, dplyr::ends_with('__cmpl')) %>%
     lapply(tidyr::gather, key = 'subgroup', value = 'cmpl') %>%
     lapply(dplyr::mutate, subgroup = stringi::stri_replace_all_fixed(subgroup, pattern = '__cmpl', replacement = ''))

   for(i in seq_along(mn_obj_list))
   {
      mn_obj_list[[i]] %<>%
         dplyr::mutate(
            object_id = facet_labels[i]) %>%
         dplyr::bind_cols(
            sd_obj_list[[i]] %>%
               dplyr::select(
                  sd),
            cmpl_obj_list[[i]] %>%
               dplyr::select(
                  cmpl))
   }

   mn_obj_list %<>%
      # Combine the objects
      dplyr::bind_rows() %>%
      # Ensure deterministic facet order
      dplyr::mutate(
         object_id = factor(object_id, levels = facet_labels))

   # Filter for completed (if appropriate)
   if(complete_data)
   {
      mn_obj_list %<>%
         dplyr::filter(cmpl == 1)
   }

   # Plot
   mn_obj_list %>%
      ggplot2::ggplot(
         ggplot2::aes(x = mean, y = sd^exponent)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::geom_smooth(color = 'red') +
      ggplot2::facet_grid(
         subgroup ~ object_id,
         scales   = 'free',
         labeller = facet_grid_labeller) +
      ggplot2::labs(
         x = 'Mean',
         y = switch(
            as.character(exponent),
            '0.5' = expression(sqrt(sigma)),
            '1'   = expression(sigma),
            '2'   = expression(sigma^2)))
}

utils::globalVariables(c('cmpl', 'object_id', 'sd', 'subgroup'))
