#' Create arrow diagram
#' @param object SummarizedExperiment
#' @param topdef character(1): top definition
#' @return character(1)
#' @examples

#' if (require(autonomics.data)){
#'    object <- autonomics.data::stemdiff.proteinratios
#'    object$subgroup
#'    
#'    topdef <- 'fdr < 0.05'
#'    cat(create_arrow_diagram(object, topdef))
#' }
#' @importFrom magrittr %>%
#' @export
create_arrow_diagram <- function(object, topdef){
   posdef <- paste0(topdef, '& effect > 0')
   negdef <- paste0(topdef, '& effect < 0')
   s <- '$$'
   for (contrast_name in names(autonomics.import::contrastdefs(object))){
      s <- sprintf('%s\n\\xrightarrow[%d]{%d} %s', s,
                   autonomics.find::n_contrast_features(object, contrast_name, negdef),
                   autonomics.find::n_contrast_features(object, contrast_name, posdef),
                   contrast_name)
   }
   s <- sprintf('%s\n$$', s)
   s
}



#' Create grid arrow plot
#' @examples 
#' object <- subramanian.2016::metabolon
create_grid_arrow_plot <- function(
   object, 
   color_var    = object %>% autonomics.plot::default_color_var(), 
   color_values = object %>% autonomics.plot::default_color_values(color_var)
){
   
   # Organize subgroups
   n.components <- object %>% autonomics.import::split_values(keep=TRUE) %>% ncol()
   object %>% autonomics.import::subgroup_levels()
      
   sep <- object %>% autonomics.import::subgroup_levels() %>% autonomics.import::guess_sep(.)
          object %>% autonomics.import::subgroup_levels() %>% stringi::stri_split_fixed(sep) %>% vapply(length, integer(1))
   
   
              stringi::stri_split_fixed() %>% 
              vapply(extract)
   
   object %>% autonomics.import::split_values(keep=TRUE)
   
   
   
   require(magrittr)
   subgroup_labels <- autonomics.import::subgroup_levels(object)
   ntick <- length(subgroup_labels)*2
   subgroup_x <- grid::unit(seq(1, by = 2, to = ntick)   / ntick, 'npc')
   contrast_x <- grid::unit(seq(2, by = 2, to = ntick-1) / ntick, 'npc')
   
   grid::grid.text(label = subgroup_labels, 
                   x     = subgroup_x, 
                   gp    = grid::gpar(col = color_values))
   
   grid::grid.segments(x0    = contrast_x - grid::unit(0.4/ntick, 'npc'), 
                       x1    = contrast_x + grid::unit(0.4/ntick, 'npc'), 
                       y0    = grid::unit(0.5, 'npc'), 
                       y1    = grid::unit(0.5, 'npc'), 
                       arrow = grid::arrow(length = grid::unit(0.25/ntick, 'npc'), ends = 'both', angle = 30))
   
}
