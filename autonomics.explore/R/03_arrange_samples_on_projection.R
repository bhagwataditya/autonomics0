# # Arrange samples on PCA
# # @param object  eset
# # @param na.impute logical
# # @examples 
# # require(magrittr)
# # if (require(billing.differentiation.data)){
# #    billing.differentiation.data::rna.voomcounts %>% 
# #       arrange_samples_on_pca() %>% colnames() %>% head()
# # }
# # @return     ordered eset
# # @importFrom magrittr   %>%   %<>%
# # @export
# arrange_samples_on_pca <- function(object, na.impute = FALSE){
#    ordered_samples <- object %>% 
#                       pca_transform(na.impute = na.impute) %>%
#                       magrittr::extract2('samples')   %>%
#                       magrittr::extract(order(.$x), ) %>%
#                       magrittr::extract2('sample')    %>% 
#                       as.character()
#    object %<>% magrittr::extract(, ordered_samples)
#    object
# }