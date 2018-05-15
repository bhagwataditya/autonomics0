require(magrittr)
object  <- autonomics.data::billing2016
contrasts <- autonomics.find::default_contrasts(object) %>% magrittr::extract(1:2)
object %<>% autonomics.import::filter_samples(subgroup %in% c('BM_E', 'BM_EM'))
design    <- autonomics.find::create_design_matrix(object)

new_limma <- object %>% autonomics.find::run_limma(    contrasts, design)
old_limma <- object %>% autonomics.find::run_limma_old(contrasts, design)

# dimensions: identical
(names(new_limma) == names(old_limma))         %>%  all()
(dim(new_limma[[1]]) == dim(old_limma[[1]]))   %>%  all()
(dim(new_limma[[2]]) == dim(old_limma[[2]]))   %>%  all()
(dim(new_limma[[3]]) == dim(old_limma[[3]]))   %>%  all()
(dim(new_limma[[4]]) == dim(old_limma[[4]]))   %>%  all()
(dim(new_limma[[5]]) == dim(old_limma[[5]]))   %>%  all()
(dim(new_limma[[6]]) == dim(old_limma[[6]]))   %>%  all()

# colnames: identical
(colnames(new_limma[[1]]) == colnames(old_limma[[1]]))  %>%  all()
(colnames(new_limma[[2]]) == colnames(old_limma[[2]]))  %>%  all()
(colnames(new_limma[[3]]) == colnames(old_limma[[3]]))  %>%  all()
(colnames(new_limma[[4]]) == colnames(old_limma[[4]]))  %>%  all()
(colnames(new_limma[[5]]) == colnames(old_limma[[5]]))  %>%  all()
(colnames(new_limma[[6]]) == colnames(old_limma[[6]]))  %>%  all()

# rownames: identical
(rownames(new_limma[[1]]) == rownames(old_limma[[1]]))  %>%  all()
(rownames(new_limma[[2]]) == rownames(old_limma[[2]]))  %>%  all()
(rownames(new_limma[[3]]) == rownames(old_limma[[3]]))  %>%  all()
(rownames(new_limma[[4]]) == rownames(old_limma[[4]]))  %>%  all()
(rownames(new_limma[[5]]) == rownames(old_limma[[5]]))  %>%  all()
(rownames(new_limma[[6]]) == rownames(old_limma[[6]]))  %>%  all()

# Updated limma improves group-specific NA handling
which.na <- list(BM_E  = autonomics.import::exprs(object)[, 1:3] %>% is.na() %>% matrixStats::rowAlls() %>% which(), 
                 BM_EM = autonomics.import::exprs(object)[, 4:6] %>% is.na() %>% matrixStats::rowAlls() %>% which())
VennDiagram::venn.diagram(which.na, filename = NULL) %>% autonomics.support::cdraw()
which.na %<>% Reduce(union, .)
autonomics.import::exprs(object) %>% magrittr::extract(which.na, )
old_limma$coef %>% magrittr::extract(which.na, )
new_limma$coef %>% magrittr::extract(which.na, )

# Non NA results remain identical
list(old = old_limma$coef %>% magrittr::extract(-which.na, ), 
     new = new_limma$coef %>% magrittr::extract(-which.na, )) %>% 
   (function(x) x$old == x$new) %>% 
   matrixStats::rowAlls() %>% 
   all()
