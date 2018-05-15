# run_single_pca <- function(object, dims, verbose){
#    pca_res       <- pca_transform(object, dims, verbose = verbose)
#    plot_df       <- autonomics.import::sdata(object)
#    plot_df$x     <- pca_res$samples$x
#    plot_df$y     <- pca_res$samples$y
#    plot_df$dimx  <- dims[1]
#    plot_df$dimy  <- dims[2]
#    plot_df$pctx  <- pca_res$var$x
#    plot_df$pcty  <- pca_res$var$y
# #   plot_df$facet <- facet
#    plot_df
# }
# 
# make_df_single_facet <- function(object, dims, verbose){
#    plot_df <- run_single_pca(object, dims, verbose = verbose)
# }
# 
# # @importFrom magrittr   %>%   %<>%
# make_df_svar_facets <- function(object, dims, verbose, facet_var){
#    # Convert factors to characters to prevent problems when rbinding later
#    # http://stackoverflow.com/questions/2851015/convert-data-frame-columns-from-factors-to-characters
#    factor_columns <- autonomics.import::sdata(object) %>% sapply(is.factor) %>% magrittr::extract()
#    autonomics.import::sdata(object)[factor_columns][] %<>% lapply(as.character)
#    facet_values <- autonomics.import::sdata(object)[[facet_var]] %>% autonomics.support::factorify() %>% levels()
#    extract_subset <- function(facet_value){
#       object %>%
#       autonomics.import::filter_samples_(sprintf('%s == "%s"', facet_var, facet_value))
#    }
#    eset_list <- Map(extract_subset, facet_values)
#    plot_df <-  mapply(run_single_pca,
#                       object = eset_list,
#                       MoreArgs = list(dims = dims, verbose = verbose),
#                       SIMPLIFY = FALSE
#    ) %>% dplyr::bind_rows()
#    plot_df$facet <- sprintf('%s . . . %d%% : %d%%',
#                             plot_df[[facet_var]],
#                             plot_df$pctx, plot_df$pcty)
#    plot_df
# }
# 
# # @importFrom magrittr   %>%   %<>%
# make_df_dim_facets <- function(object, dims, verbose){
# 
#    # 1.2  3.4  5.6  7.8  9.10  11.12
#    modulo <- length(dims) %% 2
#    dims %<>% magrittr::extract(seq(1, length(dims) - modulo ))        # should be even number
#    dims %<>% magrittr::extract(1:min(ncol(object), length(dims)))   # not more dims than samples!
#    dims_list <- dims %>% split(rep(seq(1, length(dims)/2), each = 2))
# 
#    # 1.2  1.3  1.4  2.3  2.4  3.4
#    # comb_mat <- combn(dims, 2)
#    # dims_list <- comb_mat %>% split(rep(1:ncol(.), each = nrow(.)))
#    # dims_list %<>% set_names(sprintf('X%d.X%d', comb_mat[1,], comb_mat[2,]))
#    plot_df <- mapply(run_single_pca,
#                      dims     = dims_list,
#                      MoreArgs = list(object = object, verbose = verbose),
#                      SIMPLIFY = FALSE
#    ) %>% dplyr::bind_rows()
#    plot_df$facet <- sprintf('X%02d : X%02d . . . %d%% : %d%%',
#                             plot_df$dimx, plot_df$dimy,
#                             plot_df$pctx, plot_df$pcty)
#    plot_df
# }
# 
# plot_pca_samples_base <- function(plot_df,
#                                   color_var,
#                                   color_values,
#                                   txt_var,
#                                   shape_var, # 'fixed', 'free', or NULL
#                                   group_var,
#                                   facet_var,
#                                   nrow = NULL,
#                                   title,
#                                   xlab,
#                                   ylab,
#                                   legend.position,
#                                   scales){
#    # Plot
#    p <- ggplot2::ggplot(plot_df)
#    if (is.null(shape_var)){
#       p <- p + ggplot2::geom_point(ggplot2::aes_string(x = 'x', y = 'y', color = color_var), shape = 15, size = 3)
#    } else {
#       p <- p + ggplot2::geom_point(ggplot2::aes_string(x = 'x', y = 'y', color = color_var, shape = shape_var), size = 3)
#    }
# 
#    # Connect samples from same group with line
#    if (!is.null(group_var)){
#       p <- p + ggplot2::geom_path(ggplot2::aes_string(x = 'x', y = 'y', color = color_var, group = group_var))
#    }
# 
#    # Axes
#    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
#             ggplot2::geom_hline(yintercept = 0, linetype = 'dashed') +
#             ggplot2::xlab(xlab) +
#             ggplot2::ylab(ylab) +
#             ggplot2::ggtitle(title)
# 
#    # Annotate samples
#    if (!is.null(txt_var)){
#       p <- p + ggrepel::geom_text_repel(ggplot2::aes_string(label = txt_var, x='x', y='y', color = color_var), show.legend = FALSE)
#    }
# 
#    # Facet
#    if (!is.null(facet_var)){
#       p <- p + ggplot2::facet_wrap(stats::as.formula(sprintf('~ %s', facet_var)), scales = scales, nrow = nrow)
#    }
# 
#    # Customize colors
#    if (!is.null(color_values)){
#       p <- p + ggplot2::scale_color_manual(values = color_values)
#    }
# 
#    # Legend
#    p <- p + ggplot2::theme(legend.position = legend.position)
#    # Return
#    return(p)
# 
# }
# 
# # Plot PCA sample biplot
# #
# # Perform a double centered PCA transformation and plot the sample biplot. Samples are
# # colored based on different levels of a specified variable.
# #
# # The double data centering operation prior to principal component analysis
# # was found by Wouters et al. (2003, see reference) to substantially improve
# # the biological relevance over a traditional PCA in microarray studies.
# #
# # @param object        eset
# # @param split_var       svar on which to split data before pca analysis. Plot will be faceted on split_var too.
# # @param color_var       svar mapped to color
# # @param color_values    color vector (names = color_var levels, values = colors)
# # @param txt_var         svar mapped to txt
# # @param shape_var       svar mapped to shape
# # @param group_var       svar used to connect points with a line
# # @param facet_var       svar on which to facet plot
# # @param nrow            number of rows in faceted plot
# # @param size_fixed      size of the squares (samples) in the biplot
# # @param title           Figure title
# # @param dims            vector specifying which spm dimensions on biplot (e.g. c(1,2))
# # @param verbose         logical. If \code{TRUE}, more output is provided.
# # @param na.impute       logical
# # @param legend.position 'none', 'right', ...
# # @return ggplot object
# # @author Aditya Bhagwat
# # @importFrom magrittr  %<>%    %>%
# # @examples
# # require(magrittr)
# # plot_pca_samples(NULL)
# # if (require(autonomics.data)){
# #    object <- ALL
# #    plot_pca_samples(object, color_var = 'BT')
# #    plot_pca_samples(object, color_var = c('subgroup', 'BT', 'sex'))
# #    plot_pca_samples(object, color_var = 'BT', dims = 1:12)
# #    plot_pca_samples(object, color_var = 'BT', facet_var = 'sex')
# # }
# # if (require(billing.differentiation.data)){
# #    billing.diff.rna     %>%  plot_pca_samples(na.impute = FALSE, shape_var = NULL)
# #    protein.ratios       %>%  plot_pca_samples(na.impute = FALSE, shape_var = NULL)
# #    protein.ratios       %>%  plot_pca_samples(na.impute = TRUE,  shape_var = NULL)
# #    phospho.ratios       %>%  plot_pca_samples(na.impute = FALSE, shape_var = NULL)
# #    phospho.ratios       %>%  plot_pca_samples(na.impute = TRUE,  shape_var = NULL)
# #    phospho.occupancies  %>%  plot_pca_samples(na.impute = TRUE)
# #
# #    billing.differentiation.data::rna.voomcounts %>% plot_pca_samples()
# # }
# # if (require(subramanian.2016)){
# #    object <- subramanian.2016::metabolon
# #    object %>% autonomics.explore::plot_pca_samples(
# #                    color_var = 'condition')
# #    object %>% autonomics.explore::plot_pca_samples(
# #               color_var ='time', split_var = 'condition')
# # }
# # if (require(halama.2016)){
# #    halama.2016::cell.metabolites %>% autonomics.explore::plot_pca_samples(
# #       color_var    = 'GROUP_DESCRIPTION',
# #       color_values = c(Control = 'orange', Vehicle = 'red', `Concentration 1` = 'green',
# #                       `Concentration 2` = 'blue'))
# # }
# # @references L. Wouters et al. (2003). Graphical exploration of gene expression data:
# #             a comparative study of three multivariate methods. Biometrics 59, 1131-1140.
# # @export
# plot_pca_samples <- function(
#    object,
#    split_var       = NULL,
#    color_var       = autonomics.plot::default_color_var(object),
#    color_values    = autonomics.plot::default_color_values(object),
#    txt_var         = autonomics.plot::default_txt_var(object),
#    shape_var       = autonomics.plot::default_shape_var(object),
#    group_var       = NULL,
#    facet_var       = NULL,
#    nrow            = NULL,
#    size_fixed      = 3,
#    title           = NULL,
#    dims            = c(1, 2),
#    verbose         = TRUE,
#    na.impute       = FALSE,
#    legend.position = 'right'
# ){
#    # Return empty ggplot object if NULL eset
#    # This option is useful in loops
#    if (is.null(object)){
#       return(ggplot2::ggplot())
#    }
# 
#    # Check input args
#    autonomics.import::assert_is_valid_eset(object)
#    if (ncol(object) < 3){
#      autonomics.support::cmessage('\tExit PCA: only %d samples', ncol(object))
#      return(invisible(NULL))
#    }
#    assertive.types::assert_is_character(color_var)
#    assertive.sets::assert_is_subset(color_var, autonomics.import::svars(object))
#    if (!is.null(shape_var)){
#       assertive.types::assert_is_a_string(shape_var)
#       assertive.sets::assert_is_subset(shape_var, autonomics.import::svars(object))
#       autonomics.import::sdata(object)[[shape_var]] %<>% autonomics.support::factorify()
#    }
#    if (!is.null(facet_var)){
#       if (is.factor(autonomics.import::sdata(object)[[facet_var]])){
#          autonomics.import::sdata(object)[[facet_var]] %<>% as.character()
#       }
#    }
#    if (!is.null(title)) assertive.types::assert_is_a_string(title)
#    assertive.types::assert_is_numeric(dims)
#    dims %<>% unique()
#    assertive.numbers::assert_all_are_less_than_or_equal_to(dims, ncol(object))
#    assertive.numbers::assert_all_are_greater_than_or_equal_to(length(dims), 2)
# 
#    # Title
#    if (is.null(title)){
#       selector <- !matrixStats::rowAnys(is.na(autonomics.import::exprs(object))) &
#                   !matrixStats::rowAnys(is.infinite(autonomics.import::exprs(object)))
#       title <- if (na.impute){ sprintf('PCA on all %d features (%d imputed)', length(selector), sum(!selector))
#                } else {        sprintf('PCA on common %d features (%d total)', sum(selector),    length(selector)) }
#    }
# 
#    # Impute if required
#    if (na.impute){
#       #selector <- !matrixStats::rowAnys(is.na(exprs(object)))
#       object %<>% autonomics.import::replace_nas_with_zeros(verbose = FALSE)
#    }
# 
#    # Facet on color_var
#    if (length(color_var)> 1){
#       plot_list <- list()
#       for (i  in seq_along(color_var)){
#          title <- if (i==1) title else ''
#          plot_df <- make_df_single_facet(object, dims, verbose)
#          plot_list[[i]] <- plot_df %>%
#                            plot_pca_samples_base(color_var       = color_var[[i]],
#                                                  color_values    = color_values,
#                                                  txt_var         = txt_var,
#                                                  shape_var       = shape_var,
#                                                  group_var       = group_var,
#                                                  facet_var       = NULL,
#                                                  scales          = NULL,
#                                                  title           = title,
#                                                  xlab            = sprintf('X%d (%d%%)', dims[1], plot_df$pctx),
#                                                  ylab            = sprintf('X%d (%d%%)', dims[2], plot_df$pcty),
#                                                  legend.position = legend.position)
#       }
#       ncols <- (length(color_var) - 1) %>% sqrt() %>% floor() %>% magrittr::add(1)
#       autonomics.plot::multiplot(plotlist = plot_list, cols = ncols)
#       return(plot_list)
# 
#    # Split analysis and facet on split_var
#    } else if (!is.null(split_var)){
#       assertive.properties::assert_is_of_length(dims, 2)
#       plot_df <- make_df_svar_facets(object, dims, verbose, facet_var = split_var)
#       plot_df %>% plot_pca_samples_base(color_var       = color_var,
#                                         color_values    = color_values,
#                                         txt_var         = txt_var,
#                                         shape_var       = shape_var,
#                                         group_var       = group_var,
#                                         facet_var       = 'facet',
#                                         scales          = 'free',
#                                         title           = title,
#                                         xlab            = sprintf('X%d', dims[1]),
#                                         ylab            = sprintf('X%d', dims[2]),
#                                         legend.position = legend.position)
#    # Facet on dims
#    } else if (length(dims) > 2){
#       plot_df <- make_df_dim_facets(object, dims, verbose)
#       plot_df %>% plot_pca_samples_base(color_var       = color_var,
#                                         color_values    = color_values,
#                                         txt_var         = txt_var,
#                                         shape_var       = shape_var,
#                                         group_var       = group_var,
#                                         facet_var       = 'facet',
#                                         scales          = 'free',
#                                         title           = title,
#                                         xlab            = '',
#                                         ylab            = '',
#                                         legend.position = legend.position)
# 
#    # Facet on facet_var (or don't facet)
#    } else {
#       plot_df <- make_df_single_facet(object, dims, verbose)
#       plot_df %>% plot_pca_samples_base(color_var       = color_var,
#                                         color_values    = color_values,
#                                         txt_var         = txt_var,
#                                         shape_var       = shape_var,
#                                         group_var       = group_var,
#                                         facet_var       = facet_var,
#                                         scales          = 'fixed',
#                                         nrow            = nrow,
#                                         title           = title,
#                                         xlab            = sprintf('X%d (%d%%)', dims[1], plot_df$pctx),
#                                         ylab            = sprintf('X%d (%d%%)', dims[2], plot_df$pcty),
#                                         legend.position = legend.position)
#    }
# }
# 
# # # Plot sample biplots for multiple PC combinations
# # #
# # # Create multipanel figure with sample biplots for upto 6
# # # different principal component pairs.
# # #
# # # @param object eSet
# # # @param color_var sample variable mapped to color
# # # @param shape_var sample variable mapped to shape
# # # @param file file to print to
# # # @return gList object
# # # @author Aditya Bhagwat
# # # @examples
# # # if (require(autonomics.data)){
# # #    object <- ALL
# # #    plot_pca_samples_multi_dim(object = ALL, color_var = 'BT')
# # # }
# # # @importFrom multipanelfigure multi_panel_figure add_panel
# # # @importFrom grid grid.draw
# # # @importFrom grDevices pdf dev.off
# # # @export
# # plot_pca_samples_multi_dim <- function(object, color_var = NULL, shape_var = NULL, file = NULL){
# #
# #    run_pca_single_dim <- function(dims){
# #       pca_res       <- pca_transform(object, dims = dims, verbose = FALSE)
# #       sampleDF     <- autonomics.import::sdata(object)
# #       sampleDF$x   <- pca_res$samples$x
# #       sampleDF$y   <- pca_res$samples$y
# #       sampleDF$dim <- sprintf('X%d:%d%%  X%d:%d%%', dims[1], pca_res$var$x, dims[2], pca_res$var$y)
# #       return(sampleDF)
# #    }
# #
# #    plot_df <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)) %>%
# #               lapply(run_pca_single_dim) %>%
# #               dplyr::bind_rows()
# #
# #    p <- ggplot(sampleDF)
# #    if (is.null(shape_var)){
# #       p <- p + geom_point(aes_string(x = 'x', y = 'y', color = color_var), shape = 15, size = size_fixed)
# #    } else {
# #       p <- p + geom_point(aes_string(x = 'x', y = 'y', color = color_var, shape = shape_var), size = size_fixed)
# #    }
# #
# #    p <- p +
# #       geom_vline(xintercept = 0, linetype = 'dashed') +
# #       geom_hline(yintercept = 0, linetype = 'dashed') +
# #       xlab(xlab) +
# #       ylab(ylab) +
# #       ggtitle(title)
# #
# #    if (!is.null(txt_var)){
# #       p <- p + geom_text(aes_string(label = txt_var, x='x', y='y', color = color_var, hjust=0, vjust=0))
# #    }
# #
# #    if (!is.null(facet_var)){
# #       p <- p + facet_wrap(as.formula(sprintf('~ %s + pct', facet_var)), scales = 'free')
# #    }
# #
# #
# #
# #
# #    addSingleSpm <- function(multiFig, dim, top_panel, left_panel){
# #       p <- plot_pca_samples(object, color_var, shape_var = shape_var, dim = dim, title = '')
# #       multiFig <- add_panel(p, multiFig, top_panel = top_panel, left_panel = left_panel)
# #       return(multiFig)
# #    }
# #
# #    nGraph <- min(floor(ncol(object)/2), 6)
# #    nCol <- c('1' = 1, '2' = 2, '3' = 2, '4' = 2, '5' = 3, '6' = 3)[[nGraph]]
# #    nRow <- c('1' = 1, '2' = 1, '3' = 2, '4' = 2, '5' = 2, '6' = 2)[[nGraph]]
# #
# #
# #    figWidth  <- 5*nCol # inch
# #    figHeight <- 4*nRow # inch
# #    multiFig <- multi_panel_figure(width = figWidth, height = figHeight, columns = nCol, rows = nRow, unit = 'inch')
# #    for (iGraph in seq(1, nGraph)){
# #       iRow <- 1 + (iGraph - 1) %%  2
# #       iCol <- 1 + (iGraph - 1) %/% 2
# #       multiFig <- addSingleSpm(multiFig, dim = c(iGraph*2-1,iGraph*2),   top_panel = iRow, left_panel = iCol)
# #    }
# #
# #    if (!is.null(file)){
# #       pdf(file, width = figWidth, height = figHeight)
# #       grid.draw(multiFig)
# #       dev.off()
# #    } else {
# #       return(multiFig)
# #    }
# #
# # }
