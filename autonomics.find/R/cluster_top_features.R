utils::globalVariables('feature_id')
utils::globalVariables('subgroup')
utils::globalVariables('value')
utils::globalVariables('.')

#' Compute median per subgroup
#' @param object    eSet
#' @examples
#' require(magrittr)
#' 
#' # STEM CELL COMPARISON
#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% autonomics.find::compute_median_per_subgroup() %>% head()
#' }
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    object %>% autonomics.find::compute_median_per_subgroup()
#' }
#' @importFrom data.table  data.table   :=
#' @importFrom magrittr    %>%
#' @export
compute_median_per_subgroup <- function(object){
  sample_id_var <- autonomics.import::get_sample_id_var(object)
  subgroup_data <- autonomics.import::sdata(object)    %>%
                   data.table::as.data.table() %>%
                   magrittr::extract(, c(sample_id_var, 'subgroup'), with = FALSE) %>%
                   data.table::setkeyv(sample_id_var)

  sg_exprs_dt <- autonomics.import::exprs(object)  %>%
                 data.table::as.data.table()       %>%
                 magrittr::extract(, feature_id := autonomics.import::fnames(object)) %>%
                 data.table::melt.data.table(id.vars = 'feature_id', variable.name = sample_id_var) %>%
                 data.table::setkeyv(sample_id_var) %>%
                 merge(subgroup_data)               %>%
                 magrittr::extract(,
                         list(median = stats::median(value, na.rm = TRUE),
                              mad    = stats::mad(value,   na.rm = TRUE)),
                         by = list(feature_id, subgroup)) %>%
                 data.table::dcast.data.table(feature_id ~ subgroup, value.var = 'median') %>% 
                 data.table::setkeyv('feature_id') %>% 
                 magrittr::extract(autonomics.import::fnames(object))

  sg_exprs_mat <- sg_exprs_dt                 %>%
                  magrittr::extract(, -1, with = FALSE) %>%
                  data.matrix()               %>%
                  magrittr::set_rownames(sg_exprs_dt[, feature_id])
                  sg_exprs_mat
}



#' Filter top significant features per contrast
#' @param object  SummarizedExperiment
#' @param n       number of top features to select.
#'                Actual number of features in returned eset will be lower,
#'                due to overlaps of top features among contrasts.
#' @examples
#' library(magrittr)
#' if (require(autonomics.data)){
#'    autonomics.data::stemcomp.proteinratios %>% 
#'    autonomics.find::filter_top_features_per_contrast()
#' }
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%
#' @export
filter_top_features_per_contrast <- function(object, n = 1000){
   
   # Satisfy CHECK
   variable <- NULL
   
   # Return if no limma in fdata
   if (!contains_limma_in_fdata(object)){
      autonomics.support::cmessage('\t\tNo limma in fdata - abort')
      return(object)
   }
   
   contains_p <- autonomics.import::fvars(object) %>%
                 stringi::stri_detect_fixed('p.') %>%
                 any()
   top_var <- ifelse(contains_p, 'p.', 'rank.')

   n_contrasts <- autonomics.import::fvars(object) %>%
                  stringi::stri_detect_fixed('rank.') %>% sum()
   n_per_contrast <- 1000 %/% n_contrasts
   top_features <- autonomics.import::fdata(object) %>%
                   magrittr::extract(, stringi::stri_detect_fixed(names(.), top_var)) %>%
                   data.table::as.data.table() %>%
                   magrittr::extract(, feature_id := autonomics.import::fnames(object)) %>%
                   data.table::melt.data.table(id.vars = 'feature_id') %>%
                   magrittr::extract(order(variable, value), ) %>%
                   magrittr::extract(, rank := seq(1, .N), by = 'variable') %>%
                   magrittr::extract(rank <= n_per_contrast) %>%
                   magrittr::extract(, feature_id) %>%
                   unique()
   object %>% magrittr::extract(top_features, )
}

#' @importFrom magrittr   %>%
print_cluster_i <- function(
   i,
   top_eset,
   apres,
   result_dir,
   geom = autonomics.plot::default_feature_plots(top_eset) %>% setdiff('bars'),
   x             = autonomics.plot::default_x(top_eset, geom[1]),
   color_var     = autonomics.plot::default_color_var(top_eset),
   shape_var     = autonomics.plot::default_shape_var(top_eset),
   group_var     = autonomics.plot::default_group_var(top_eset),
   txt_var       = autonomics.plot::default_txt_var(top_eset),
   line          = autonomics.plot::default_line(top_eset)
){
   subdir <- sprintf('%s/cluster%03d', result_dir, i)
   dir.create(subdir, showWarnings = FALSE)
   exemplar_feature <- names(apres@exemplars[i])
   other_features <- names(apres@clusters[[i]]) %>% setdiff(exemplar_feature)
   cluster_features <- c(exemplar_feature, other_features) # exemplar first
   iset <- top_eset %>% magrittr::extract(cluster_features, )

   iset %>% autonomics.import::write_fdata_to_file(file = sprintf('%s/cluster%03d_features.txt', subdir, i))

   plotargs <- list(object = iset, x = x, color_var = color_var, shape_var = shape_var, group_var = group_var, txt_var = txt_var, line = line)
   for (curplot in geom){
      curfile <- sprintf('%s/cluster%03d_%s.pdf', subdir, i, curplot)
      curargs <- c(plotargs, geom = curplot, file = curfile)
      autonomics.plot::plot_features %>% do.call(curargs)
   }
}

#' Cluster features on subgroup profiles
#' @param object      eset
#' @param result_dir    result directory
#' @param geom which types of feature plots to generate
#' @param x             svar mapped to x     in feature plots
#' @param color_var     svar mapped to color in feature plots
#' @param shape_var     svar mapped to shape in feature plots
#' @param group_var     svar mapped to group in feature plots
#' @param txt_var       svar mapped to txt   in feature plots
#' @param line          whether to connect points with a line in feature plots
#' @examples
#' \dontrun{
#'    require(magrittr)
#'    result_dir <- tempdir() %T>% message()
#'    if (require(subramanian.2016)){
#'       object <- subramanian.2016::metabolon[1:100,]
#'       object %>% cluster_features_on_subgroups(
#'                       result_dir, color_var = 'condition', group_var = 'condition')
#'    }
#' }
#' @importFrom magrittr   %>%   %<>%   %T>%
#' @export
cluster_features_on_subgroups <- function(
   object,
   result_dir,
   geom = autonomics.plot::default_feature_plots(object) %>% setdiff('bars'),
   x             = autonomics.plot::default_x(object, geom[1]),
   color_var     = autonomics.plot::default_color_var(object),
   shape_var     = autonomics.plot::default_shape_var(object),
   group_var     = autonomics.plot::default_group_var(object),
   txt_var       = autonomics.plot::default_txt_var(object),
   line          = TRUE
){

   # Return if less than 3 samples
   if (ncol(object) < 3){
      autonomics.support::cmessage('\t\tAbort clustering: %d sample(s) only', ncol(object))
      return(invisible(NULL))
   }

   # Summarize subgroup exprs
   median_exprs <- object %>% compute_median_per_subgroup()

   # Replace NA with 0
   selector <- is.na(median_exprs)
   autonomics.support::cmessage('\t\tReplace %d NA values with zero', sum(selector))
   median_exprs[selector] <- 0

   # Cluster
   autonomics.support::cmessage('\t\tCluster')
   cormat <- stats::cor(t(median_exprs))
   apres <- apcluster::apcluster(s = cormat, details = TRUE, q=0)

   # Abort if no clusters identified
   if(length(apres) == 0){
      autonomics.support::cmessage('\t\t\tNo clusters found')
      return(invisible(NULL))
   }

   # Plot heatmap
   autonomics.support::cmessage('\t\tHeatmap')
   dir.create(result_dir, showWarnings = FALSE)
   grDevices::pdf(paste0(result_dir, '/00_heatmap.pdf'))
   apcluster::heatmap(apres, cormat)
   grDevices::dev.off()

   # Plot exemplar profiles
   autonomics.support::cmessage('\t\tExemplars')
   file_name <- paste0(result_dir, '/00_exemplars.pdf')
   plotargs <- list(object  = object %>% magrittr::extract(names(apres@exemplars), ),
                    x         = x,
                    color_var = color_var,
                    shape_var = shape_var,
                    group_var = group_var,
                    txt_var   = txt_var,
                    line      = line)
   for (curplot in geom){
      curfun <- sprintf('plot_feature_%s', curplot) %>%
                utils::getFromNamespace('autonomics.plot')
      curfile <- sprintf('%s/00_exemplar_%s.pdf', result_dir, curplot)
      curargs <- plotargs %>% c(list(file = curfile))
      curfun %>% do.call(curargs)
   }

   # Print cluster features
   autonomics.support::cmessage('\t\tFeatures')
   autonomics.import::fdata(object) %<>% magrittr::extract(, !stringi::stri_detect_regex(names(.), '(rank|value|p|fdr|bonf)[.]'))
   seq(1, length(apres@clusters)) %>% plyr::l_ply(print_cluster_i,
                                                  object,
                                                  apres      = apres,
                                                  result_dir = result_dir,
                                                  x          = x,
                                                  color_var  = color_var,
                                                  shape_var  = shape_var,
                                                  group_var  = group_var,
                                                  txt_var    = txt_var,
                                                  line       = line)
}


#' Cluster top features (on subgroup profiles)
#' @param object      eSet updated with add_limma_to_fdata
#' @param design        design matrix
#' @param contrasts     contrasts definition
#' @param result_dir    result directory
#' @param n             number of top features
#' @param x             svar mapped to x in feature plots
#' @param geom which type of feature plots to generate
#' @param color_var     svar mapped to color in feature plots
#' @param shape_var     svar mapped to shape in feature plots
#' @param group_var     svar mapped to group in feature plots
#' @param txt_var       svar mapped to txt   in feature plots
#' @param line          whether to connect points in feature plot with line (logical)
#' @examples
#' \dontrun{
#' library(magrittr)
#' if (require(autonomics.data)){
#'    result_dir <- tempdir() %T>% message()
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object$subgroup
#'    contrasts <- c(BM_E = 'BM_E', BM_EM = 'BM_EM', EM_E = 'EM_E')
#'    object %>% autonomics.find::cluster_top_features_on_subgroups(
#'                    contrasts  = contrasts,
#'                    result_dir = result_dir,
#'                    n          = 500)
#' }
#'
#' if (require(billing.differentiation.data)){
#'    library(magrittr)
#'    result_dir <- tempdir() %T>% message()
#'    object <- billing.differentiation.data::protein.ratios
#'    contrasts <- c(EM01_0   = 'EM01-EM0.0',
#'                   EM02_01  = 'EM02-EM01',
#'                   EM05_02  = 'EM05-EM02',
#'                   EM15_05  = 'EM15-EM05',
#'                   EM30_15  = 'EM30-EM15')
#'    object %<>% autonomics.import::filter_samples(subgroup != 'BM0.0')
#'    object %>% cluster_top_features_on_subgroups(
#'                    contrasts  = contrasts,
#'                    result_dir = result_dir,
#'                    n          = 500)
#' }
#'
#' if (require(atkin.2014)){
#'    library(magrittr)
#'    result_dir <- tempdir() %T>% message()
#'    object <- atkin.2014::soma
#'    object$subgroup
#'    contrasts <- c(Diabetic_Control = 'D.t0 - C.t0')
#'    object %>% cluster_top_features_on_subgroups(
#'                    contrasts = contrasts,
#'                    result_dir = result_dir)
#' }
#'
#' if (require(subramanian.2016)){
#' library(magrittr)
#'    result_dir <- tempdir() %T>% message()
#'    object <- subramanian.2016::metabolon
#'    contrasts <- subramanian.2016::contrasts.metabolon
#'    object %>% cluster_top_features_on_subgroups(
#'                    contrasts = contrasts,
#'                    result_dir = result_dir,
#'                    color_var = 'condition',
#'                    group_var = 'condition')
#' }
#' }
#'
#' @importFrom  magrittr  %<>%   %>%
#' @export
cluster_top_features_on_subgroups <- function(
   object,
   design = autonomics.import::create_design_matrix(object),
   contrasts = autonomics.find::default_contrasts(object),
   result_dir,
   n             = 1000,
   geom = autonomics.plot::default_feature_plots(object) %>% setdiff('bars'),
   x             = autonomics.plot::default_x(object, geom[1]),
   color_var     = autonomics.plot::default_color_var(object),
   shape_var     = autonomics.plot::default_shape_var(object),
   group_var     = autonomics.plot::default_group_var(object),
   txt_var       = autonomics.plot::default_txt_var(object),
   line          = TRUE
){

   # Report
   autonomics.support::cmessage('\tCluster top features on subgroups')

   # Exit
   if (length(unique(object$subgroup)) < 2){
      autonomics.support::cmessage('\t\tAbort: # subgroups < 2')
      return(invisible(NULL))
   }


   # Run limma if requiredl
   cluster_dir <- paste0(result_dir, '/feature_clusters')
   if (!contains_limma_in_fdata(object)){
      object %<>% add_limma_to_fdata(contrasts, design)
   }

   # Filter significant and top ranked features
   object %<>% filter_significant_features()    %>%
                 filter_top_ranked_features()
   
   # Cluster
   object %>% cluster_features_on_subgroups(
                   cluster_dir,
                   geom     = geom,
                   x         = x,
                   color_var = color_var,
                   shape_var = shape_var,
                   group_var = group_var,
                   txt_var   = txt_var,
                   line      = line
                )

}
