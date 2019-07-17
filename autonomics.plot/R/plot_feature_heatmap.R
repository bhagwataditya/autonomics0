#' Plot feature heatmap
#' @param object          SummarizedExperiment
#' @param fvar            string
#' @param filename        string
#' @param width           number: width  in inches
#' @param height          number: height in inches
#' @param cellheight      number: heatmap cell height
#' @param fontsize        number: fontsize
#' @param treeheight_row  number: height of tree for rows
#' @param treeheight_col  number: height of tree for col
#' @param cluster_cols    TRUE or FALSE
#' @param ...             passed to pheatmap::pheatmap
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::glutaminase[1:50, ]
#'    object %>% plot_feature_heatmap_with_pheatmap(fvar='BIOCHEMICAL')
#'    filename <- tempfile() %T>% message()
#'    object %>% plot_feature_heatmap_with_pheatmap(filename = filename)
#' }
#' @importFrom magrittr %>%
#' @export
plot_feature_heatmap_with_pheatmap <- function(
   object,
   fvar       = 'feature_id',
   filename   = NULL,
   width      = 7,
   height     = 7,
   cellheight = 8,
   fontsize   = 6,
   treeheight_row = 0, #0 used for not showing dendogram for rows
   treeheight_col = 0, #0 used for not showing dendogram for cols
   cluster_cols = FALSE,
   ...
){
   #.Deprecated("plot_feature_heatmap")
   autonomics.support::cmessage('\t\tImpute NA(N) values')
   object %<>% autonomics.preprocess::impute(method = 'impute.QRILC')
   autonomics.import::fnames(object) <- object %>% autonomics.import::fvalues(fvar) %>% autonomics.support::uniquify()
   x <- autonomics.import::exprs(object)
   bwcolor = grDevices::colorRampPalette(c("yellow","grey", "blue"))
   if (!is.null(filename))   grDevices::pdf(filename, width = width, height = height, onefile = FALSE)
   pheatmap::pheatmap(
      x,
      filename   = NA,
      clustering_dist_rows = "correlation",
      scale      = 'row',
      cellheight = cellheight,
      fontsize   = fontsize,
      col        = bwcolor(50),
      treeheight_row = treeheight_row,
      treeheight_col = treeheight_col,
      cluster_cols = cluster_cols,
      border_color = NA,
      ...)
   if (!is.null(filename))   grDevices::dev.off()
   if (!is.null(filename))   autonomics.support::cmessage(filename)
}


#' Plot feature heatmap
#' @param object    SummarizedExperiment
#' @param fvar      string
#' @param filename  string
#' @param lowColor  color for downregulation (Default color yellow)
#' @param highColor color for upregulation (Default color royal blue)
#' @param midColor  color for midpoint (Default color gray)
#' @param width     number: width  in inches (Default width 2 inches)
#' @param height    number: height in inches (Default height 4 inches)
#' @param fontsize  number: font size for feature name (Default size 0 for not showing feature names)
#' @param ...       additional parameters handed through to \code{\link[ggplot2]{theme}}
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::glutaminase[1:10, ]
#'    object %>% plot_feature_heatmap(fvar='BIOCHEMICAL')
#'    filename <- tempfile() %T>% message()
#'    object %>% plot_feature_heatmap(filename = filename)
#' }
#' @importFrom magrittr %>%
#' @export
plot_feature_heatmap <- function(
   object,
   fvar      = 'feature_id',
   filename  = NULL,
   lowColor  = "yellow",
   highColor = "blue1",
   midColor  = "gray",
   width     = 2,
   height    = 4,
   fontsize  = 0,
   ...
){
   autonomics.import::fnames(object) <- object %>%
                                        autonomics.import::fvalues(fvar) %>%
                                        autonomics.support::uniquify()
   autonomics.support::cmessage('\t\tImpute NA(N) values')
   object %<>% autonomics.preprocess::impute(method = 'impute.QRILC')
   x <- autonomics.import::exprs(object)
   hres <- stats::hclust(stats::as.dist(1 - stats::cor(t(x), method="pearson"))) #compute correlation for matrix
   h.row.order <- hres$order # clustering on row
   row.order.df <- x[h.row.order,]
   scaled.x = t(scale(t(row.order.df),center=T)) #scaling dataframe
   df.molten.dat <- reshape2::melt(scaled.x) #resturcturing dataframe

   myPlot <- ggplot2::ggplot(data = df.molten.dat, ggplot2::aes_string(x = 'Var2', y = 'Var1', fill = 'value')) +
             ggplot2::geom_tile() +
             #ggplot2::facet_wrap(~ branch, scales = 'free') +
             ggplot2::scale_fill_gradient2(low = lowColor, mid = midColor, high = highColor) +
             ggplot2::theme(axis.text.x   = ggplot2::element_text(angle = 90, hjust = 1, size=width),
                            axis.text.y   = ggplot2::element_text(size=fontsize),
                            axis.ticks.length = grid::unit(.01, "cm"),
                            axis.ticks    = ggplot2::element_blank(),
                            axis.title.x  = ggplot2::element_blank(),
                            axis.title.y  = ggplot2::element_blank(),
                            legend.title  = ggplot2::element_blank(),
                            panel.background = ggplot2::element_blank(),
                            legend.justification = "top",
                            legend.key.width =   grid::unit(0.2, "cm"),
                            legend.key.height =  grid::unit(0.3,"cm"),
                            legend.text =        ggplot2::element_text(size=fontsize),
                            legend.box.margin =  ggplot2::margin(-10,-10,-10,-12),
                            ...)

   if (!is.null(filename)) myPlot %>% print2pdf(file = filename, width = width, height = height, onefile = FALSE)
   if (!is.null(filename)) autonomics.support::cmessage(filename)
   myPlot
}


#scale_fill_gradientn(colours = viridisLite::viridis(3, option = "D", alpha=0))

#df_molten_dat %<>% magrittr::set_names(stringi::stri_replace_first_fixed(names(.), 'Var1', 'feature'))
#branches_df <- cutree(hres, 2) %>% (function(x) data.frame(feature = names(x), branch = x))
#df_molten_dat %<>% merge(branches_df, 'feature', sort = FALSE)
#df_molten_dat[1,]
#levels(df_molten_dat$Var1)

#load('/home/shh2026/Desktop/project/autonomics/subramanian.bfa.lna148a/subramanian.bfa.lna148a/subramanian.2015/data/selected_genes.RData')

#select LNACTR LNA148a .. Exculte outliers .. select genes .. calculate logcpm .. generate heatmap
#gene.list <- data.frame(gene_name=c("RPS27","PPP2CA","RPS23","HRAS","EIF4A2","RPS7","FGFR3","RHOB","EIF4G2","EIF3D","RPS13","FGFR4","TSC2","RPS2","AKT1S1","EIF4B","RPS6KB1","AKT2","PLD3","DDIT4","RRAS","STK11","EIF3E","MLST8","EIF3M","RPS6","PPP2CB","PPP2R1A","RND3","EIF4A3","IRS1","RPS15A","RPSA","EIF3K","ULK1","RPS3A","RPS18","RHOT2","HIF1A","RPS4X","EIF4E","PRKAG1","EIF4EBP1","PGF","RPS28","HMOX1","AKT1","RHOT1","RHOD","PIK3C3","RPS9","RPTOR","RPS6KB2","RPS3","RHOF","GNB1L","RPS12","RPS24","RPS19","NRAS","EIF3H","RHOC","PPP2R5A","FAU","EIF3G","DGKZ","RPS4Y1","RPS15","RPS27L","EIF3I","PRR5","RPS27A","RPS6KA4","EIF3L","MVD","NSDHL","ACAT2","IDI1","MVK","MSMO1","TM7SF2","SC5D","GGPS1","FDPS","FDFT1","EBP","ACAT1","HADHB","LSS","HMGCR","HMGCS1","HADHA","CYP51A1","FURIN","NDUFA9","ATP5D","COX8A","ATP5L","NDUFB8","CYB5R3","NDUFA1","NDUFB3","NDUFS1","NDUFA5","NDUFS6","GPX4","ATP5F1","NDUFA8","ATP5J","NDUFS7","ATP5A1","COX7C","ATP5O","SDHC","NDUFB1","NDUFAF2","MAPK12","MT-ND4L","NDUFS5","GSR","ATP5C1","PRDX3","UQCRC2","NDUFB7","NDUFA12","VDAC1","SNCA","TXNRD2","NDUFA4","ATP5G1","COX7B","ATP5H","UQCRH","COX6C","RHOT2","NDUFB9","PARK7","ATP5J2","ATPAF1","OGDH","COX4I1","NDUFS4","NDUFAF1","NDUFV1","VDAC3","NDUFA13","APP","UQCRB","VPS9D1","NDUFS8","NDUFA11","CAT","COX7A2","SDHD","NDUFA3","MAOA","NDUFA4","COX7B","ATP5G1","NDUFA9","ATP5H","ATP5D","COX6C","UQCRH","COX8A","ATP5L","NDUFB8","NDUFA1","NDUFB3","NDUFS1","NDUFA5","NDUFB9","NDUFS6","ATP5J2","ATPAF1","ATP5F1","COX4I1","NDUFS4","NDUFA8","ATP5J","NDUFV1","NDUFS7","ATP5A1","COX7C","ATP5O","NDUFB1","SDHC","MT-ND4L","NDUFA13","UQCRB","VPS9D1","NDUFS5","ATP5C1","NDUFS8","NDUFA11","UQCRC2","NDUFB7","COX7A2","NDUFA12","SDHD","NDUFA3","FDFT1","EBP","NSDHL","MSMO1","LSS","TM7SF2","SC5D","CYP51A1","BDH2","ACAT2","ACAT1","HADHB","HMGCS1","HADHA","MAP2K7","SCAP","DDIT3","HSPA1A","HSPA1B","HSPA9","CANX","ATF6","HSPA2","EIF2A","HSPA8","HSPA4","SYVN1","HSP90B1","TRAF2","SREBF1","SREBF2","CEBPA","EIF2AK3"))
#gene.list <- gene.list[!duplicated(gene.list),]

#rnacounts1 <- rnacounts %>% filter_samples(subgroup %in% c('LNACTR', 'LNA148a')) %>%
#                           filter_samples(!sample_id %in% c('LNACTR.R3','LNA148a.R3')) %>%
#                           filter_features(gene_name %in% gene.list)

#%>% autonomics.import::logcpm()
