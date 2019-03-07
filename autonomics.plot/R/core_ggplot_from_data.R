#' Centralized convenience function to shape data and provide a core ggplot2 object to build on
#'
#' @param object              SummarizedExperiment
#' @param ...                 additonal SummarizedExperiments or matrices
#' @param listed_objects      additonal SummarizedExperiments or matrices (when wrapped in a \code{\link{list}})
#' @param MARGIN              Whether the reference axis is \code{samples} (default) or \code{features}
#' @param x                   string: variable (samples or features depending on \code{MARGIN}) mapped to x (\code{horizontal == FALSE}) or y axis (\code{horizontal == TRUE})
#' @param alpha_var           string: variable (samples or features depending on \code{MARGIN}) mapped to alpha aesthetic
#' @param color_var           string: variable (samples or features depending on \code{MARGIN}) mapped to color aesthetic
#' @param fill_var            string: variable (samples or features depending on \code{MARGIN}) mapped to fill aesthetic
#' @param group_var           string: variable (samples or features depending on \code{MARGIN}) mapped to group aesthetic
#' @param linetype_var        string: variable (samples or features depending on \code{MARGIN}) mapped to linetype aesthetic
#' @param shape_var           string: variable (samples or features depending on \code{MARGIN}) mapped to shape aesthetic
#' @param size_var            string: variable (samples or features depending on \code{MARGIN}) mapped to size aesthetic
#' @param stroke_var          string: variable (samples or features depending on \code{MARGIN}) mapped to stroke aesthetic
#' @param weigth_var          string: variable (samples or features depending on \code{MARGIN}) mapped to weight aesthetic
#' @param facet1_var          string: variable (samples or features depending on \code{MARGIN}) for faceting (\code{ggplot2::facet_grid(facet1_var ~ .)})
#' @param facet2_var          string: variable (samples or features depending on \code{MARGIN}) for faceting (\code{ggplot2::facet_grid(. ~ facet2_var)})
#' @param marked_ids          vector of explicitly marked indexes in the non-\code{MARGIN} dimension (in \code{features} if \code{MARGIN == 'samples'} and vice versa)
#' @param horizontal          TRUE or FALSE: layout plot horizontally (measurements on the x, sample names on the y axis)?
#' @param scales	            "fixed" (default), "free", "free_x", or "free_y"
#' @param facet_grid_labeller labeller (see \code{\link[ggplot2]{labellers}}) for facetting labels
#' @param xlab                string: label of y axis
#' @param ylab                string: label of y axis
#' @param title               string: title
#' @details A variable (\code{color_var} etc.) may come in two forms:
#' 1) as a string and subset of column names in the annotation of the reference \code{object},
#' resulting in the use of the corresponding data or
#' 2) as a character object of \code{length(c(object, ...))}, resulting in the augmentation of the
#' data to be plotted by the variable (variable 1 associated with \code{c(object, ...)[[1]]}, etc.).
#' See the example section.
#'
#' If \code{MARGIN == 'samples'}, Objects in \code{c(object, ...)} must have identical number of
#' samples (columns). In the case of \code{MARGIN == 'features'}, an identical number of features
#' (rows) must be present).
#'
#' Annotation is fetched from the (master)\code{object} and used on all.
#' @return ggplot2 object
#' @author Johannes Graumann
#' @examples
#' require(magrittr)
#' requireNamespace('ggplot2')
#' requireNamespace('ggstance')
#' if (require(autonomics.data)){
#'
#'    # Just the basics ...
#'       core_ggplot_from_data(autonomics.data::ALL[, 1:30])
#'
#'    # Add a 'geom' layer
#'       core_ggplot_from_data(autonomics.data::ALL[, 1:30]) +
#'       ggplot2::geom_violin()
#'
#'    # (Nearly) same thing, but horizontal and filled ...
#'       core_ggplot_from_data(autonomics.data::ALL[, 1:15],
#'                             fill_var = 'sex',
#'                             horizontal = TRUE) +
#'       ggstance::geom_violinh()
#'
#'    # Combine 2 'data sets', separating them by fill (augmenting the data on-the-fly)
#'       core_ggplot_from_data(autonomics.data::ALL[, 1:15],
#'                             autonomics.data::ALL[, 16:30],
#'                             color_var = NULL,
#'                             fill_var = c('A', 'B')) +
#'       ggplot2::geom_boxplot()
#'
#'    # Combine 2 'data sets', separating them by facetting
#'       core_ggplot_from_data(autonomics.data::ALL[, 1:30],
#'                             autonomics.data::ALL[, 31:60],
#'                             facet2_var = c('A', 'B'),
#'                             horizontal = TRUE) +
#'       ggstance::geom_boxploth()
#'
#'    # Combine with further facetting ... using sdata ...
#'       core_ggplot_from_data(autonomics.data::ALL[, 1:15],
#'                             autonomics.data::ALL[, 16:30],
#'                             facet1_var = 'sex',
#'                             facet2_var = c('A', 'B'),
#'                             horizontal = TRUE) +
#'       ggstance::geom_boxploth()
#'
#'    # Combine with further facetting ... using on-the-fly augmentation (only) ...
#'       core_ggplot_from_data(autonomics.data::ALL[, 1:15],
#'                             autonomics.data::ALL[, 16:30],
#'                             facet1_var = c('C', 'D'),
#'                             facet2_var = c('A', 'B'),
#'                             horizontal = TRUE) +
#'       ggstance::geom_boxploth()
#' }
#' @export
core_ggplot_from_data <- function(
   object,
   ...,
   listed_objects      =  list(),
   MARGIN              =  c('samples', 'features'),
   x                   =  NULL,
   alpha_var           =  NULL,
   color_var           =  default_color_var(object),
   fill_var            =  NULL,
   group_var           =  NULL,
   linetype_var        =  NULL,
   shape_var           =  NULL,
   size_var            =  NULL,
   stroke_var          =  NULL,
   weigth_var          =  NULL,
   facet1_var          =  NULL,
   facet2_var          =  NULL,
   marked_ids          =  NULL,
   horizontal          =  FALSE,
   scales              =  ifelse(horizontal, 'free_x', 'free_y'),
   facet_grid_labeller = 'label_parsed',
   xlab                =  '',
   ylab                =  '',
   title               =  ''
){

# Catch args INCUDING defaults ('match.call' doesn't) without 'NULL' --------
# https://stackoverflow.com/a/14398674/2103880
   # Catch args including defaults
   arg_list <- mget(names(formals()), sys.frame(sys.nframe()))
   # Discard NULL defaults
   arg_list %<>% magrittr::extract(!sapply(arg_list, is.null))

   # Extract args corresponding to ggplot2 aes definitions, reformat etc.
   vars <- names(arg_list) %>%
           setdiff(c('', 'object', '...', 'listed_objects', 'MARGIN', 'marked_ids', 'horizontal', 'scales', 'facet_grid_labeller', 'xlab', 'ylab', 'title'))
   special_vars <- vars %>%
                   stringi::stri_replace_all_regex('_var$','') %>%
                   paste0('_plotvar')
   all_special_vars <- special_vars %>%
                       c('value_plotvar', 'feature_plotvar', 'x_plotvar', 'marked_ids_plotvar') %>%
                       unique()

# Process arguments & check prerequisites ---------------------------------
   MARGIN %<>% match.arg(choices = c('samples', 'features'),several.ok = FALSE)

   obj_list <- list(object) %>% c(list(...), listed_objects)
   obj_list %<>% magrittr::extract(!sapply(obj_list, is.null))

   n_obj <- obj_list %>% length()

   object %>% assertive.types::assert_is_any_of(c('eSet', 'EList', 'SummarizedExperiment'))

   # Check all '...' objects and ensure compatibility with master 'object'
   if (n_obj > 1){
      for (i in utils::tail(obj_list, n = -1)){
         i %>% assertive.types::assert_is_any_of(c('eSet', 'EList', 'SummarizedExperiment', 'matrix'))
         if(MARGIN == 'samples'){ assertive.base::assert_are_identical(ncol(object), ncol(i))
         } else {                 assertive.base::assert_are_identical(nrow(object), nrow(i))
         }
      }
   }

   # Extract data in 'MARGIN' dependent way
   localData      <- switch(MARGIN, samples  = autonomics.import::sdata(object),
                                    features = autonomics.import::fdata(object))
   localVars      <- switch(MARGIN, samples  = autonomics.import::svars(object),
                                    features = autonomics.import::fvars(object))
   localNames     <- switch(MARGIN, samples  = autonomics.import::snames(object),
                                    features = autonomics.import::fnames(object))
   localIdsToMark <- switch(MARGIN, samples  = autonomics.import::fdata(object) %>% magrittr::extract2('feature_id'),
                                    features = autonomics.import::sdata(object) %>% magrittr::extract2('sample_id'))

   # Step through variables corresponding to ggplot2 aes definitions
   for(va in vars){
      va <- arg_list[[va]]
      va %>% assertive.types::assert_is_character() %>%
             assertive.sets::assert_are_disjoint_sets(all_special_vars) %>%
             length() %>%
             assertive.sets::assert_is_subset(c(1, n_obj))
      if(length(va) == 1){
         va %>% assertive.sets::assert_is_subset(localVars)
         object[[va]] <- factorize_numerics(object, va, localData)
      }
   }

   if(!is.null(marked_ids)){
      marked_ids %>% assertive.numbers::assert_all_are_whole_numbers() %>%
                     assertive.numbers::assert_all_are_in_range(lower = 0,
                                                                upper = switch(MARGIN, samples  = nrow(object), features = ncol(object)))
      marked_ids <- localIdsToMark %>% magrittr::extract2(marked_ids)
   }

   horizontal %>% assertive.types::assert_is_a_bool()

   facet_grid_labeller %>% assertive.types::assert_is_a_string() %>%
                           assertive.strings::assert_all_are_matching_regex('^label_') %>%
                           exists(where=asNamespace('ggplot2'), mode='function') %>%
                           assertive.base::assert_all_are_true()

   scales %<>% match.arg(choices = c('fixed', 'free', 'free_x', 'free_y'), several.ok = FALSE)

   if (!is.null(xlab))  xlab  %>% assertive.types::assert_is_a_string()
   if (!is.null(ylab))  ylab  %>% assertive.types::assert_is_a_string()
   if (!is.null(title)) title %>% assertive.types::assert_is_a_string()

# Combine and shape the data ----------------------------------------------
   # Extract names
   xVec <- autonomics.support::factorify(localNames)
   if ('x' %in% vars) {
      if(MARGIN == 'samples'){ xVec <- object %>% magrittr::extract(,x)
      } else {                 xVec <- object %>% magrittr::extract(x,) }
   }

   # Extract/shape data matrices
   obj_list %<>% lapply(function(ob){
                           if (any(assertive.base::is2(ob, c('eSet', 'EList', 'SummarizedExperiment')))){
                              ob %<>% autonomics.import::exprs()
                           }
                           ob %>% return() })
   if(MARGIN == 'samples') obj_list %<>% lapply(t)
   obj_list %<>% lapply(data.frame, x_plotvar = xVec)

   # Deal with other variables
   for (va in stringi::stri_subset_regex(vars, pattern = '_var$')){
     plotvar <- va %>% stringi::stri_replace_all_regex('_var$', '_plotvar')
     ct <- arg_list[[va]]
     if(length(ct) == 1){
        obj_list %<>% lapply(cbind, localData[ct] %>% magrittr::set_colnames(plotvar))
     } else {
        for (i in seq_along(obj_list)){
           obj_list[[i]] %<>% cbind(rep(ct[i], times = nrow(obj_list[[i]])) %>%
                              as.data.frame() %>%
                              magrittr::set_colnames(plotvar))
           # Ensure original order
           obj_list[[i]][[plotvar]] %<>% factor(levels = ct)
        }
     }
   }

   # Combine data.frames
   prepDF <- obj_list %>% dplyr::bind_rows()

   # 'gather' into long form
   plotDF <- suppressMessages(tidyr::gather(prepDF,
                                            "feature_plotvar",
                                            "value_plotvar",
                                             setdiff(colnames(prepDF), all_special_vars)))

   # Mark IDs (as appropriate)
   if(!is.null(marked_ids)){
      plotDF %<>% dplyr::mutate(marked_plotvar = feature_plotvar %in% marked_ids)
   }

# Generate basic plot object ----------------------------------------------
   p <- ggplot2::ggplot(plotDF,
                 ggplot2::aes_string(
                    x        = ifelse(horizontal, 'value_plotvar', 'x_plotvar'    ),
                    y        = ifelse(horizontal, 'x_plotvar',     'value_plotvar'),
                    alpha    = if('alpha_var'    %in% vars) 'alpha_plotvar'     else  NULL,
                    color    = if('color_var'    %in% vars) 'color_plotvar'     else  NULL,
                    fill     = if('fill_var'     %in% vars) 'fill_plotvar'      else  NULL,
                    group    = if('group_var'    %in% vars) 'group_plotvar'     else  NULL,
                    linetype = if('linetype_var' %in% vars) 'linetype_plotvar'  else  NULL,
                    shape    = if('shape_var'    %in% vars) 'shape_plotvar'     else  NULL,
                    size     = if('size_var'     %in% vars) 'size_plotvar'      else  NULL,
                    stroke   = if('stroke_var'   %in% vars) 'stroke_plotvar'    else  NULL,
                    weight   = if('stroke_var'   %in% vars) 'stroke_plotvar'    else  NULL))

   # Add faceting info
   if (c('facet1_var', 'facet2_var') %>%
       magrittr::is_in(vars) %>%
       any())
   {
      facetting_string <- c(ifelse('facet1_var' %in% vars, 'facet1_plotvar', '.'),
                            ifelse('facet2_var' %in% vars, 'facet2_plotvar', '.')) %>%
                          paste(collapse = ' ~ ')
      p <- p + ggplot2::facet_grid(facetting_string,
                                   scales   = scales,
                                   labeller = facet_grid_labeller)
   }

   # Add labels & title and return
   p <- p + ggplot2::labs(x        = xlab,
                          y        = ylab,
                          title    = title,
                          alpha    = ifelse('alpha_var'    %in% vars && length(arg_list[['alpha_var']])    == 1, arg_list[['alpha_var']],    ''),
                          color    = ifelse('color_var'    %in% vars && length(arg_list[['color_var']])    == 1, arg_list[['color_var']],    ''),
                          fill     = ifelse('fill_var'     %in% vars && length(arg_list[['fill_var']])     == 1, arg_list[['fill_var']],     ''),
                          group    = ifelse('group_var'    %in% vars && length(arg_list[['group_var']])    == 1, arg_list[['group_var']],    ''),
                          linetype = ifelse('linetype_var' %in% vars && length(arg_list[['linetype_var']]) == 1, arg_list[['linetype_var']], ''),
                          shape    = ifelse('shape_var'    %in% vars && length(arg_list[['shape_var']])    == 1, arg_list[['shape_var']],    ''),
                          size     = ifelse('size_var'     %in% vars && length(arg_list[['size_var']])     == 1, arg_list[['size_var']],     ''),
                          stroke   = ifelse('stroke_var'   %in% vars && length(arg_list[['stroke_var']])   == 1, arg_list[['stroke_var']],   ''),
                          weight   = ifelse('weight_var'   %in% vars && length(arg_list[['weight_var']])   == 1, arg_list[['weight_var']],   ''))

   if(!horizontal) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1))
   p %>% return()
}

utils::globalVariables(c('feature_plotvar'))

factorize_numerics <- function(object, var, localData){
   retr_var <- localData %>% magrittr::extract2(var)
   if(is.numeric(retr_var)) retr_var %<>% as.character() %>% factor()
   retr_var %>% return()
}
