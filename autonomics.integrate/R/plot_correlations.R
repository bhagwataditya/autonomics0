utils::globalVariables(c('R2', 'cor', 'facet', 'pair', 'sample_id', '.SD'))

#' Align platform columns
#'
#' Ensure: col1 = eset1 features, col2 = eset2 features
#' @param cor_dt correlation data.table
#' @param eset1 eset1
#' @param eset2 eset2
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    cor_dt <- subramanian.2016::top.cor.exiqon.metabolon
#'    exiqon    <- subramanian.2016::exiqon    %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    metabolon <- subramanian.2016::metabolon %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    cor_dt
#'    cor_dt %>% autonomics.integrate::align_platform_columns(exiqon, metabolon)
#'    cor_dt %>% autonomics.integrate::align_platform_columns(metabolon, exiqon)
#' }
#' align_platform_columns
#' @importFrom magrittr %>%
#' @export
align_platform_columns <- function(cor_dt, eset1, eset2){
   name1 <- names(cor_dt) %>% magrittr::extract(cor_dt[1,] %in% autonomics.import::fnames(eset1))
   name2 <- names(cor_dt) %>% magrittr::extract(cor_dt[1,] %in% autonomics.import::fnames(eset2))
   other_names <- names(cor_dt) %>% setdiff(name1) %>% setdiff(name2)
   cor_dt %<>% magrittr::extract(, c(name1, name2, other_names ), with = FALSE)
   cor_dt
}

#' Create correlation plot data.table
#' @param cor_dt correlation data.table
#' @param eset1  first eset
#' @param eset2  second eset
#' @param fvars1    eset1 fvars mapped to panel annotation
#' @param fvars2    eset2 fvars mapped to panel annotation
#' @param color_var eset1 svar mapped to color
#' @param shape_var eset1 svar mapped to shape
#' @return data.table(facet, sample_id, x, y, color)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    cor_dt    <- subramanian.2016::top.cor.exiqon.metabolon
#'    exiqon    <- subramanian.2016::exiqon    %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    metabolon <- subramanian.2016::metabolon %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    autonomics.integrate::create_cor_plot_dt(cor_dt[1:6, ], exiqon, metabolon)
#' }
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%          %<>%
#' @export
create_cor_plot_dt <- function(
   cor_dt,
   eset1,
   eset2,
   fvars1 = autonomics.integrate::default_corplot_fvars(eset1),
   fvars2 = autonomics.integrate::default_corplot_fvars(eset2),
   color_var = 'subgroup',
   shape_var = NULL
){

   # Re-order columns if required
   #cor_dt %<>% align_platform_columns(eset1, eset2)
   name1 <- names(cor_dt)[1]
   name2 <- names(cor_dt)[2]

   # eset -> long table
   long1  <- autonomics.import::sumexp_to_long_table(eset1, fvars = fvars1)  %>%
             data.table::setnames(autonomics.import::fid_var(eset1), name1)  %>%
             data.table::setnames('value', 'x')
   long2  <- autonomics.import::sumexp_to_long_table(eset2, fvars = fvars2)  %>%
             data.table::setnames(autonomics.import::fid_var(eset2), name2) %>%
             data.table::setnames('value', 'y')
   svars1 <- c('sample_id', color_var)
   if (!is.null(shape_var)) svars1 %<>% c(shape_var)
   sdata1 <- autonomics.import::sdata(eset1) %>% magrittr::extract(, svars1, drop = FALSE)

   # Merge and order
   cor_dt[, pair:= paste0(get(name1), ' : ', get(name2))]
   cor_dt[, pair:= autonomics.support::factorify(pair)]
   plotdf <- cor_dt %>% # https://stackoverflow.com/questions/23087358/why-is-allow-cartesian-required-at-times-when-when-joining-data-tables-with-dupl
             merge(long1, by = c(name1), allow.cartesian = TRUE) %>%
             merge(long2, by = c(name2, 'sample_id')) %>%
             merge(sdata1, by = c('sample_id'), all.x = TRUE, all.y = FALSE)

   # Define facet in format fid1 : fid2 : fvars1 : fvars2
   #    fid1   and fid2:  make the facets unique
   #    fvars1 and fars2: annotate plot
   plotdf[, R2   := paste0('R2 = ', as.character(signif(cor, 2)))]
   my_SDcols <- c(name1, name2) %>%
                 (function(x) if (length(fvars1)==0)   x %>% c(name1)  else  x %>% c(fvars1)) %>%
                 (function(x) if (length(fvars2)==0)   x %>% c(name2)  else  x %>% c(fvars2)) %>%
                c('R2')
   plotdf[, facet := apply(.SD, 1, paste, collapse = ' : '), .SDcols = my_SDcols]
   plotdf[, (c('R2', name1, name2, fvars1, fvars2, 'cor')) := NULL]
   plotdf %>% data.table::setorder(pair, sample_id)
   plotdf[, facet := autonomics.support::factorify(facet)]
   plotdf[, pair:= NULL]

   # Re-order columns
   plotdf %<>% magrittr::extract(, c('facet', 'sample_id', 'x', 'y', color_var, shape_var), with = FALSE)
   plotdf
}

# wrap on "facet" to have order correct
# use this labeler to annotate nicely
my_labeller <- function(labels){
   labels <- lapply(labels, as.character)
   labels$facet %<>% strsplit(' : ') %>%
      vapply(function(x){x %>%
            magrittr::extract(3:length(x)) %>%
            paste0(collapse = ' : ')
      }, character(1))

   n <- length(unlist(strsplit(labels$facet[1], ' : ')))
   for (i in 1:n){
      labels[[paste0('label', i)]] <- labels$facet %>% strsplit(' : ') %>% vapply(magrittr::extract, character(1), i)
   }
   labels$facet <- NULL
   labels
}

# wrap on "facet" to have order correct
# use this labeler to annotate nicely
my_labeller <- function(labels){
   labels <- lapply(labels, as.character)
   labels$facet %<>% strsplit(' : ') %>%
      vapply(function(x){x %>%
            magrittr::extract(3:length(x)) %>%
            paste0(collapse = ' : ')
      }, character(1))

   n <- length(unlist(strsplit(labels$facet[1], ' : ')))
   for (i in 1:n){
      labels[[paste0('label', i)]] <- labels$facet %>% strsplit(' : ') %>% vapply(magrittr::extract, character(1), i)
   }
   labels$facet <- NULL
   labels
}

#' Plot single page with correlation plots
#' @param cor_dt correlation datatable
#' @param eset1  first eset
#' @param eset2  second eset
#' @param fvars1    eset1 fvars mapped to panel annotation
#' @param fvars2    eset2 fvars mapped to panel annotation
#' @param color_var eset1 svar mapped to color
#' @param shape_var eset1 svar mapped to shape
#' @param xlab      x axis label
#' @param ylab      y axis label
#' @param ncol number of columns
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    cor_dt    <- subramanian.2016::top.cor.exiqon.metabolon
#'    exiqon    <- subramanian.2016::exiqon    %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    metabolon <- subramanian.2016::metabolon %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    autonomics.integrate::plot_correlations(cor_dt[1:10, ], exiqon, metabolon)
#' }
#' @importFrom data.table   data.table
#' @importFrom magrittr     %>%
#' @export
plot_correlations <- function(
   cor_dt,
   eset1,
   eset2,
   fvars1    = autonomics.integrate::default_corplot_fvars(eset1),
   fvars2    = autonomics.integrate::default_corplot_fvars(eset2),
   color_var = 'subgroup',
   shape_var = NULL,
   xlab      = NULL,
   ylab      = NULL,
   ncol      = 3
){
   # Validify and set defaults
   cor_dt %<>% autonomics.integrate::align_platform_columns(eset1, eset2)
   name1 <- names(cor_dt)[1]
   name2 <- names(cor_dt)[2]
   if (is.null(xlab))   xlab <- name1
   if (is.null(ylab))   ylab <- name2

   # Create plot data.table
   plotdf <- create_cor_plot_dt(cor_dt, eset1, eset2, fvars1, fvars2, color_var = color_var, shape_var = shape_var)
   title  <- sprintf('%s ~ %s', name2, name1) # cor_dt[, unique(as.character(get(name1)))] %>% paste0(collapse = ' + ') %>% sprintf('%s ~ %s', name2, .)
   p <- ggplot2::ggplot(plotdf, ggplot2::aes_string(x = 'x', y = 'y',  color = color_var, shape = shape_var)) +
        ggplot2::facet_wrap(~ facet, scales = 'free', labeller = my_labeller, ncol = ncol) +
        ggplot2::geom_point() +
        # geom_point(size = 2.5) +
        ggplot2::ggtitle(title) +
        ggplot2::ylab(ylab) +
        ggplot2::xlab(xlab) +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 9))
   print(p)
}

#' Plot correlations
#' @param cor_dt     correlation datable, as returned by autonomics.integrate::correlate_features
#' @param eset1      eset1 first eset
#' @param eset2      eset2 second eset
#' @param file       file to which to print plots
#' @param fvars1     eset1 fvars used for faceting
#' @param fvars2     eset2 fvars used for faceting
#' @param color_var  eset1 svar mapped to color
#' @param shape_var  eset1 svar mapped to shape
#' @param xlab       x axis label
#' @param ylab       y axis label
#' @param ncol       no of columns on one page
#' @param nrow       no of rows on one page
#' @param max_pages  max no of pages
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    cor_dt    <- subramanian.2016::top.cor.exiqon.metabolon
#'    exiqon    <- subramanian.2016::exiqon    %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    metabolon <- subramanian.2016::metabolon %>%
#'                 autonomics.import::filter_samples(condition==unique(.$condition)[2])
#'    file <- tempfile() %>% paste0('.pdf') %T>% message()
#'    autonomics.integrate::plot_correlations_to_file(cor_dt, exiqon, metabolon, file = file)
#'    autonomics.integrate::plot_correlations_to_file(cor_dt, exiqon, metabolon, file = file,
#'                                            fvars2 = c('BIOCHEMICAL', 'SUB_PATHWAY'))
#' }
#'@export
plot_correlations_to_file <- function(
   cor_dt,
   eset1,
   eset2,
   file,
   fvars1         = autonomics.integrate::default_corplot_fvars(eset1),
   fvars2         = autonomics.integrate::default_corplot_fvars(eset2),
   color_var      = 'subgroup',
   shape_var      = NULL,
   xlab           = NULL,
   ylab           = NULL,
   nrow           = 3,
   ncol           = 6,
   max_pages      = Inf
){

   # Assert
   assertive.types::assert_is_data.table(cor_dt)
   autonomics.import::assert_is_valid_eset(eset1)
   autonomics.import::assert_is_valid_eset(eset2)
   assertive.sets::assert_is_subset(fvars1, autonomics.import::fvars(eset1))
   assertive.sets::assert_is_subset(fvars2, autonomics.import::fvars(eset2))
   assertive.sets::assert_are_disjoint_sets(fvars1, fvars2)
   assertive.sets::assert_is_subset(color_var, autonomics.import::svars(eset1))
   assertive.types::assert_is_numeric(nrow)
   assertive.types::assert_is_numeric(ncol)

   # Plot
   npanel <- nrow*ncol
   npage  <- ceiling(nrow(cor_dt) / npanel)
   npage %<>% min(max_pages)
   #grDevices::pdf(file, paper = 'USr', width = ncol*5) # width = ncol*4, height = nrow*4,
   panel_width  <- 11.69/3   # A4 dimensions are  11.69 x 8.27
   panel_height <- 8.27/2    # A4 fits 2 x 3 panels well

   grDevices::pdf(file, height = nrow*panel_height, width = ncol*panel_width)
   for (page in 1:npage){
      cor_cur <- cor_dt[seq((page-1)*npanel + 1,  page*npanel)]
      plot_correlations(cor_cur,
                        eset1,
                        eset2,
                        fvars1        = fvars1,
                        fvars2        = fvars2,
                        shape_var     = shape_var,
                        color_var     = color_var,
                        xlab          = xlab,
                        ylab          = ylab,
                        ncol          = ncol)
   }
   grDevices::dev.off()
}

