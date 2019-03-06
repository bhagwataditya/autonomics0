invwhich <- function(indices, totlength) is.element(seq_len(totlength), indices)

top_left_in_volcano <- function(effect, fdr, mlp, ntop){
   fdr_ok   <- fdr  < 0.05
   coef_ok  <- effect < -1
   coef_top <- if (any(fdr_ok)){
                  threshold <- effect %>% magrittr::extract(fdr_ok)  %>% autonomics.support::nmin(ntop+1)
                  effect < threshold
               } else {
                  rep(FALSE, length(effect))
               }
   mlp_top  <- if (any(coef_ok)){
                  threshold <- mlp %>% magrittr::extract(coef_ok) %>% autonomics.support::nmax(ntop+1)
                  mlp  > threshold
               } else {
                  rep(FALSE, length(effect))
               }
   fdr_ok & coef_ok & (coef_top | mlp_top)
}

top_right_in_volcano <- function(effect, fdr, mlp, ntop){
   fdr_ok  <- fdr  < 0.05
   coef_ok <- effect >  1
   coef_top <- if(any(fdr_ok)){
                  threshold <- effect %>% magrittr::extract(fdr_ok) %>% autonomics.support::nmax(ntop+1)
                  effect > threshold
               } else {
                  rep(FALSE, length(effect))
               }
   mlp_top <- if (any(coef_ok)){
                 threshold <- mlp %>% magrittr::extract(coef_ok)  %>% autonomics.support::nmax(ntop+1)
                 mlp > threshold
              } else {
                 rep(FALSE, length(effect))
              }
   fdr_ok & coef_ok & (coef_top | mlp_top)
}

#' Create volcano datatable
#' @param object SummarizedExperiment
#' @param ntop   number of top features (either FC wise or p wise) to be annotated
#' @return data.table
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% make_volcano_dt() %>% print()
#' }
#' @importFrom magrittr %>%
#' @importFrom data.table   data.table   :=
#' @export
make_volcano_dt <- function(object, ntop = 3){

   # Satisfy CHECK
   effect <- p <- mlp <- is.top.down <- is.top.up <- color <- fdr <- NULL

   # Extract limma datatable
   point_dt <- autonomics.import::extract_limma_dt(object)
   point_dt %<>% magrittr::extract(!is.na(effect) & !is.na(p))
   point_dt %>% magrittr::extract(, mlp  := -log10(p))

   # Prepare volcano datatable
   # Note: Using effect <= 0 (rather than effect <0) is required.
   #       Otherwise (the very few) features with effect=0 will have no effect for 'color'
   point_dt %>% magrittr::extract(effect <= 0,               color := 'effect<0')
   point_dt %>% magrittr::extract(effect >  0,               color := 'effect>0')
   point_dt %>% magrittr::extract(effect < -1 & fdr<0.05,    color := 'effect<1  &  FDR<0.05')
   point_dt %>% magrittr::extract(effect >  1 & fdr<0.05,    color := 'effect>1  &  FDR<0.05')
   point_dt %>% magrittr::extract(, is.top.down := top_left_in_volcano( effect, fdr, mlp, ntop), by = 'contrast')
   point_dt %>% magrittr::extract(, is.top.up   := top_right_in_volcano(effect, fdr, mlp, ntop), by = 'contrast')
   point_dt %>% magrittr::extract(is.top.down==TRUE, color:= 'TOP down')
   point_dt %>% magrittr::extract(is.top.up  ==TRUE, color:= 'TOP up')
   point_dt %>% magrittr::extract(, color:=factor(color, c('TOP down', 'effect<1  &  FDR<0.05', 'effect<0',
                                                           'effect>0', 'effect>1  &  FDR<0.05', 'TOP up')))
   # Return
   point_dt
}

#' Plot volcano
#' @param object SummarizedExperiment
#' @param ntop n top features to be annotated
#' @param nrow number of rows in faceted plot
#' @param legend_position handed to \code{\link[ggplot2]{theme}} as 'legend.position'
#' @return ggplot object
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    object <- autonomics.data::stemcomp.proteinratios
#'    object %>% plot_volcano()
#' }
#' @importFrom data.table   data.table   :=
#' @export
plot_volcano <- function(object, ntop = 3, nrow = NULL, legend_position = NULL){

   # Satisfy CHECK
   is.top.up <- is.top.down <- effect <- mlp <- color <- fname <- NULL

   # Set colors
   down_colors <- c(`effect<0` = 100, `effect<1  &  FDR<0.05` = 70, `TOP down` = 20) %>% vapply(function(l) grDevices::hcl(h=0,  l=l,c=100), character(1))
   up_colors   <- c(`effect>0` = 100, `effect>1  &  FDR<0.05` = 70, `TOP up`   = 20) %>% vapply(function(l) grDevices::hcl(h=120,l=l,c=100), character(1))
   color_values <- c(down_colors, up_colors)

   # Draw plot
   point_dt <- object %>% autonomics.plot::make_volcano_dt(ntop = ntop)
   txt_dt <- data.table::copy(point_dt) %>%
             magrittr::extract(is.top.up==TRUE | is.top.down==TRUE)
   tmp_plot <- ggplot2::ggplot(point_dt) +
      ggplot2::facet_wrap(~ contrast, nrow = nrow, scales = 'fixed') +
      ggplot2::geom_point(ggplot2::aes(x=effect, y=mlp, color = color), na.rm = TRUE) +
      ggrepel::geom_text_repel(data = txt_dt,
                               ggplot2::aes(x=effect, y=mlp, label=fname, color = color),
                               #hjust = 'outward',
                               na.rm = TRUE,
                               show.legend = FALSE#,
                               #direction = 'x'
      ) +
      ggplot2::scale_color_manual(values = color_values, name = NULL) +
      ggplot2::xlab(expression(log[2](FC))) +
      ggplot2::ylab(expression(-log[10](p)))
   if(!is.null(legend_position))
   {
      tmp_plot <- tmp_plot +
         ggplot2::theme(legend.position = legend_position)
   }
   return(tmp_plot)
}
