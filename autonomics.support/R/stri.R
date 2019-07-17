
#' Vectorized stri_replace_first_fixed
#' @param x character vector
#' @param pattern_vector character vector
#' @param replacement_vector character vector
#' @examples
#' require(magrittr)
#' x <- c("Intensity L", 
#'        "Intensity M", 
#'        "Intensity H",
#'        "Intensity L E(L).EM(M).BM(H).R1", 
#'        "Intensity M E(L).EM(M).BM(H).R1", 
#'        "Intensity H E(L).EM(M).BM(H).R1",
#'        "Intensity BM(L).E(M).EM(H).R2",   
#'        "Intensity L BM(L).E(M).EM(H).R2", 
#'        "Intensity M BM(L).E(M).EM(H).R2",
#'        "Intensity H BM(L).E(M).EM(H).R2", 
#'        "Intensity EM(L).BM(M).E(H).R3",   
#'        "Intensity L EM(L).BM(M).E(H).R3",
#'        "Intensity M EM(L).BM(M).E(H).R3", 
#'        "Intensity H EM(L).BM(M).E(H).R3")
#' pattern_vector <- c("E(L).EM(M).BM(H).R1", "BM(L).E(M).EM(H).R2", "EM(L).BM(M).E(H).R3")
#' replacement_vector <- ''
#' x %>% autonomics.support::vstri_replace_first_fixed(pattern_vector, replacement_vector = '')
#' @importFrom magrittr %>%
#' @export
vstri_replace_first_fixed <- function(x, pattern_vector, replacement_vector){
   if (length(replacement_vector)==1)   replacement_vector %<>% rep(length(pattern_vector))

   assertive.properties::assert_are_same_length(pattern_vector, replacement_vector)
   for (pattern in pattern_vector){
      for (replacement in replacement_vector){
         x %<>% stringi::stri_replace_first_fixed(pattern, replacement)
      }
   }
   x
}

