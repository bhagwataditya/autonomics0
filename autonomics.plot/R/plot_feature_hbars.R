
utils::globalVariables('.')
#' Plot feature bars
#'
#' Create bar plot of exprs(object). Annotate it with fvars. Color per 'color_var'. Print to 'file'
#'
#' @param object          \code{eSet} with \code{fvars} in fData.
#' @param fvars             fvars used to annotate in plot
#' @param color_var         svar mapped to color
#' @param color_values      color value vector (names = subgroups, contents = colours)
#' @param facet_def         facet definition string
#' @param alpha_var         svar mapped to transparancy
#' @param file              name of file to write to
#' @param width             width (inches)
#' @param height            height (inches)
#' @param xlab              xlab annotation of plot
#' @param title             plot title
#' @param legend.position   legend position
#' @param x                 NULL. Generifies interface of plot_xxx functions.
#' @param shape_var         NULL. Generifies interface of plot_xxx functions.
#' @importFrom ggplot2      aes_string     coord_flip   facet_wrap    geom_bar   geom_hline   ggplot
#' @importFrom ggplot2      theme_bw       theme        scale_alpha_discrete xlab
#' @importFrom magrittr     %>%   %<>%
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    small_eset <- autonomics.data::ALL[1:10, 1:10]
#'    plot_feature_hbars(small_eset)
#'    plot_feature_hbars(small_eset, fvars = 'gene_names')
#'    plot_feature_hbars(small_eset, color_var = 'sex')
#'    plot_feature_hbars(small_eset, facet_def = '~age+sex')
#'    small_eset$alpha <- small_eset$sex == 'M'
#'    plot_feature_hbars(small_eset, alpha_var = 'alpha')
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts %>%
#'                extract(1:10, )
#'    object %>% autonomics.plot::plot_feature_hbars(
#'       color_var = 'cell', color_values = c(EM = 'red', BM = 'green'))
#' }
#' @export
plot_feature_hbars <- function(object,
  fvars           = default_fvars(object),
  color_var       = autonomics.plot::default_color_var(object),
  color_values    = autonomics.plot::default_color_values(object, color_var),
  facet_def       = '~ sample',
  alpha_var       = NULL,
  file            = NULL,
  width           = NULL,
  height          = NULL,
  xlab            = NULL,
  title           = '',
  legend.position = 'right',
  x               = NULL,
  shape_var       = NULL
){

  # Check object, fvars, color_var, alpha_var
  if (nrow(object) == 0){
    message('\t\tno features to plot, no graph created')
    return()
  }
  assertive.sets::assert_is_subset(color_var, autonomics.import::svars(object))
  assertive.sets::assert_is_subset(fvars, autonomics.import::fvars(object))
  if (!is.null(alpha_var)){
    assertive.sets::assert_is_subset(alpha_var, autonomics.import::svars(object))
    assertive.types::assert_is_logical(autonomics.import::sdata(object)[[alpha_var]])
  }

  # Add variables 'feature', 'sample'
  autonomics.import::fdata(object)$feature <- autonomics.import::fnames(object) %>% autonomics.support::factorify(reverse = TRUE)
  autonomics.import::sdata(object)$sample  <- autonomics.import::snames(object)  %>% autonomics.support::factorify()

  # Check facet_def (after having added 'sample' var!)
  assertive.base::assert_is_identical_to_true(grepl('~', facet_def))
  facet_vars <- strsplit(facet_def, '~') %>% unlist() %>% getElement(2) %>% strsplit('\\+') %>% unlist() %>% trimws()
  assertive.sets::assert_is_subset(facet_vars, autonomics.import::svars(object))

  # Create variable plotvar
  autonomics.import::fdata(object)$plotvar <- autonomics.import::fdata(object)[, fvars, drop = FALSE]   %>%
                                      autonomics.support::characterify(char_length = 100)   %>%
                                      apply(1, paste0, collapse = '   ')                %>%
                                      unname()                                          %>%
                                      autonomics.support::uniquify()                        %>%
                                      autonomics.support::factorify(reverse = TRUE)

  # Prepare plot DF
  L2R <- data.frame(feature = autonomics.import::fdata(object)$feature, autonomics.import::exprs(object), stringsAsFactors = FALSE, check.names = FALSE)
  L2R %<>% tidyr::gather_(key_col = 'sample', value_col = 'L2R', gather_cols = setdiff(colnames(.), 'feature'))
  L2R$sample %<>% factor(levels = autonomics.import::sdata(object)$sample)
  plotDF <- L2R %>% merge(autonomics.import::fdata(object)[, c('feature', 'plotvar')], by = 'feature') %>%
                    merge(autonomics.import::sdata(object),                            by = 'sample')

  # Create plot
  # Note: position = 'identity' is needed to avoid unnecessary warnings
  #       See http://stackoverflow.com/a/19482830/2894278
  myPlot <- ggplot2::ggplot(plotDF, ggplot2::aes_string(x = 'plotvar', y = 'L2R', fill = color_var, alpha = alpha_var)) +
            ggplot2::geom_bar(stat = 'identity', position = 'identity') +
            ggplot2::geom_hline(yintercept = 0) + ggplot2::coord_flip() +
            ggplot2::facet_wrap(stats::as.formula(facet_def), nrow = 1) +
            ggplot2::xlab('') + ggplot2::ggtitle(title) + ggplot2::ylab(xlab)

  # Take care of alpha transparancy
  if (!is.null(alpha_var)){
    if (all(plotDF[[alpha_var]] == TRUE)){
       myPlot <- myPlot + scale_alpha_discrete(range = c(1,1), guide = FALSE)        # If all subgroups are selected, none should be faded out!
    } else {
       myPlot <- myPlot + scale_alpha_discrete(range = c(0.4, 1), guide = FALSE)
    }
  }

  # Customize color
  if (!is.null(color_values)){
     myPlot <- myPlot + ggplot2::scale_fill_manual(values = color_values)
  }

  # Legend position
  myPlot <- myPlot + ggplot2::theme(legend.position = legend.position)

  # Return graph itself or file to which printed
  if (is.null(file)){
     return(myPlot)
  } else {
     if (is.null(width))   width  <- 5 + 1.5 * ncol(object) + (1 + nrow(object)/40) * 2
     if (is.null(height))  height <- 1 + 0.2 * nrow(object)
     suppressWarnings(myPlot %>% autonomics.support::print2pdf(file, width = width, height = height))
     return(invisible(file))
  }
}


