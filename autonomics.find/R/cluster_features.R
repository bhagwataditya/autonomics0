devtools::load_all('../autonomics.import')

#' Create contrasts for t-based clustering
create_cluster_contrasts <- function(object){
   dt <- split_subgroup_levels(object)
   contrastdefs <- c(
      dt[, .(contrast = sprintf('%s - %s', subgroup[-1], subgroup[-.N])), by = 'V1'][['contrast']],
      dt[, .SD[1], by = 'V1'][, sprintf('%s - %s', subgroup[-1], subgroup[-.N])])
   contrastnames <- c(
      dt[, .(contrast = sprintf('%s_%s', subgroup[-1], subgroup[-.N])), by = 'V1'][['contrast']],
      dt[, .SD[1], by = 'V1'][, sprintf('%s_%s', subgroup[-1], subgroup[-.N])])
   contrastdefs %<>% magrittr::set_names(contrastnames)   
   contrastdefs
   #design <- object %>% autonomics.import::create_design_matrix(intercept = TRUE)
   #contrastdefs <- colnames(design)[-1] %>% magrittr::set_names(., .)
   #object %<>% autonomics.find::add_limma(contrastdefs = contrastdefs, design = design)
}

#' Plot heatmap
#' @param object SummarizedExperiment
#' @param scale_features TRUE or FALSE: whether to scale (i.e. z-score) features
#' @examples 
#' if (require(autonomics.data)){
#'     object <- autonomics.data::glutaminase
#'     object %<>% extract(, order(object$subgroup))
#'     feature_sample_heatmap(object)
#' }
feature_sample_heatmap <- function(object){
   
   # Reverse clutser information to have first item on top 
   features_clustered <- all(c('cluster', 'cluster_order') %in% fvars(object))
   if (features_clustered){
      object %<>% extract(rev(fdata(object)$cluster_order), )
      fdata(object)$cluster %<>% factor(rev(levels(.)))
   }
   
   # Z-score features   
   exprs(object) %<>% t() %>% scale() %>% t()
   
   # Plot image
   fields::image.plot(
      x = seq(1,ncol(object)), 
      y = seq(1,nrow(object)),
      z = t(autonomics.import::exprs(object)), 
      col = hcl.colors(12, "YlOrRd", rev = TRUE),
      xaxt= "n", yaxt= "n", xlab = "", ylab = "")
   axis( 1, at=seq(1,length.out=ncol(object) ), labels= colnames(object), las= 2, cex.axis = 0.7)
   axis( 2, at=seq(1,length.out=nrow(object) ), labels= rownames(object), las= 2, cex.axis = 0.7)
   
   # Add subgroup lines. Add cluster lines
   abline(v = cumsum(table(object$subgroup))+0.5)
   abline(v = cumsum(table(split_subgroup_values(object)$V1))+0.5, lwd = 3)
   
   if (features_clustered){
      abline(h = cumsum(table(fdata(object)$cluster))+0.5)
   }
   
}

fsplit <- function(object, fvar){
   Map(function(x) object[fdata(object)[[fvar]] == x, ], 
       flevels(object, fvar))
}

plot_cluster_contrastograms <- function(object){
   n <- length(autonomics.import::flevels(object, 'cluster'))
   par(mfrow = c(ceiling(sqrt(n)), floor(sqrt(n))),     # 2x2 layout
       oma = c(0, 0, 0, 0), # bottom, left, top, right space
       mar = c(2, 2, 2, 2), # space for one row of text at ticks and to separate plots
       mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
       xpd = NA)
   lapply(fsplit(object, 'cluster'), plot_contrastogram, sparse = TRUE)
}

plot_cluster_boxplots <- function(object){
   lapply(flevels(object, 'cluster'), 
          function(cluster_i){
             autonomics.plot::plot_features(
                object[fdata(object)$cluster == cluster_i, ], geom = 'boxplot', 
                title = cluster_i)})
}


#' Add feature clustering
#' 
#' Cluster features on t-values
#' 
#' @param object SummarizedExperiment
#' @examples 
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::exiqon
#'    object %>% autonomics.plot::plot_overlayed_sample_distributions()
#' }
#' if require(autonomics.data){
#'     object <- autonomics.data::glutaminase
#'     object$TREATMENT <- stringi::stri_split_fixed(object$subgroup, '_') %>% 
#'     vapply(extract, character(1), 1)
#'     group_var <- 'TREATMENT'
#'     object %<>% add_clustering()
#' }
#' @importFrom magrittr %>% 
#' @export
add_clustering <- function(object){
   
   # Get t matrix
   contrastdefs <- create_cluster_contrasts(object)
   object %<>% autonomics.find::add_limma(contrastdefs)
   object %<>% autonomics.import::filter_features_('F.p < 0.05')
   tmat <- autonomics.import::limma(object)[, , 't']
   # How to deal with NAs - no NAS left in exiqon data
   #idx <- tmat %>% is.na() %>% matrixStats::rowAnys()
   #tmat[idx, ]
   
   # Compute correlations
   cormat <-  propagate::bigcor(t(tmat))
   cormat %<>% magrittr::extract(1:nrow(.), 1:ncol(.))
   
   # Cluster and add to object
   apres <- apcluster::apcluster(s = cormat, details = TRUE, q=0)
   heatmapres <- apcluster::heatmap(apres, cormat, col = rev(heat.colors(12)))
   exemplar_idx <- heatmapres@exemplars %>% extract2(length(.)) %>% extract(heatmapres@order)
   exemplars  <- fid_values(object)[exemplar_idx]
   fdata(object)$cluster <- factor(autonomics.import::fid_values(object)[apres@idx], exemplars)
   fdata(object)$cluster_order <- heatmapres@clusters  %>% extract2(length(.)) %>% extract(heatmapres@order) %>% unlist()
   
   # Plot
   feature_sample_heatmap(object)
   plot_cluster_contrastograms(object)
   print(autonomics.plot::plot_features(object[exemplars, ], geom = 'boxplot'))
   
   # Return   
   object
}

# 
# par(mfrow = c(2,2),     # 2x2 layout
#    oma = c(0, 0, 0, 0), # bottom, left, top, right space
#    mar = c(2, 2, 2, 2), # space for one row of text at ticks and to separate plots
#    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
#    xpd = NA)
# plot_contrastogram(fsplit(object, 'cluster')[[2]], sparse = TRUE)
# plot_contrastogram(fsplit(object, 'cluster')[[3]], sparse = TRUE)
# plot_contrastogram(fsplit(object, 'cluster')[[5]], sparse = TRUE)
# plot_contrastogram(fsplit(object, 'cluster')[[9]], sparse = TRUE)
# 
# 
# heatmap_order <- autonomics.import::fdata(object)$cluster_order
# plotobject <- object[heatmap_order, ]
# autonomics.import::exprs(plotobject) %<>% t() %>% scale() %>% t()
# 
# ggplot2::ggplot(autonomics.import::sumexp_to_long_dt(plotobject[1:10, ])) + 
# ggplot2::geom_density(ggplot2::aes(x=value, color = feature_id))
# 
# ggplot2::ggplot(autonomics.import::sumexp_to_long_dt(plotobject)) + 
# ggplot2::geom_density(ggplot2::aes(x=value, color = sample_id)) + 
# ggplot2::guides(color = FALSE)


