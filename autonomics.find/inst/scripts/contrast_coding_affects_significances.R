#################################################
# Choice of contrast coding affects significances
#################################################
eset1 <- WCQA.0109.16MLBL.CDT::metabolon
autonomics.import::fdata(eset1) %<>% magrittr::extract(, c('BIOCHEMICAL', 'feature_id'))

# contr.treatment
model.matrix(~ subgroup, data = autonomics.import::sdata(eset1)) %>% 
   magrittr::set_colnames(colnames(.) %>% stringi::stri_replace_first_fixed(., 'subgroup', '') %>% make.names()) %>% 
   limma::lmFit(eset1, .) %>% 
   limma::eBayes() %>% 
   limma::topTable(c('IPF', 'RA_ILD', 'RA_NOILD')) %>% 
   head(3)

# contr.sum
model.matrix(~ subgroup, data = autonomics.import::sdata(eset1), contrasts = list(subgroup = 'contr.sum')) %>% 
   magrittr::set_colnames(colnames(.) %>% 
                          stringi::stri_replace_first_fixed(., 'subgroup1', 'CTRL')  %>% 
                          stringi::stri_replace_first_fixed(., 'subgroup2', 'IPF')   %>% 
                          stringi::stri_replace_first_fixed(., 'subgroup3', 'RAILD') %>% 
                          make.names()) %>% 
   limma::lmFit(eset1, .) %>% 
   limma::eBayes() %>% 
   limma::topTable(c('CTRL', 'IPF', 'RAILD')) %>% 
   head(3)

# all pairwise comparisons between subgroups
model.matrix(~ subgroup, data = autonomics.import::sdata(eset1), contrasts = list(subgroup = 'contr.helmert')) %>% cbind(eset1$subgroup)
   



##########
# GENDER #
##########

# GENDER as a global variable is not always sufficient
eset1 <- WCQA.0109.16MLBL.CDT::metabolon
design <- model.matrix(~ 0 + subgroup + GENDER, data = autonomics.import::sdata(eset1)) %>% 
          magrittr::set_colnames(colnames(.) %>% stringi::stri_replace_first_fixed(., 'subgroup', ''))
limma::lmFit(eset1, design = design) %>% 
   limma::eBayes() %>% 
   limma::topTable(coef = 'GENDERM') %>% rownames() %>% extract(1) %>% 
   magrittr::extract(eset1, ., ) %>% 
   autonomics.plot::plot_feature_boxes(x = 'subgroup', color_var = 'GENDER')

# subgroup-specific GENDER
design <- model.matrix(~ 0 + subgroup / GENDER, data = autonomics.import::sdata(eset1)) %>% 
          magrittr::set_colnames(colnames(.) %>% 
                                 stringi::stri_replace_first_fixed(., 'subgroup', '') %>%
                                 make.names())
colnames(design)
limma::lmFit(eset1, design = design) %>% 
   limma::eBayes() %>% 
   limma::topTable(coef = c('CTRL.GENDERM', 'IPF.GENDERM', 'RA_ILD.GENDERM', 'RA_NOILD.GENDERM')) %>% 
   rownames() %>% extract(1) %>%    magrittr::extract(eset1, ., ) %>% 
   autonomics.plot::plot_feature_boxes(x = 'subgroup', color_var = 'GENDER')

# It also affetcs subgroupo performance
model.matrix(~ subgroup, data = autonomics.import::sdata(eset1), contrasts = list(subgroup = "contr.sum")) %>% cbind(eset1$subgroup)

design <- model.matrix(~ subgroup / GENDER, data = autonomics.import::sdata(eset1)) %>%
          magrittr::set_colnames(colnames(.) %>% 
                                 stringi::stri_replace_first_fixed(., 'subgroup', '') %>% 
                                 make.names())
colnames(design)
limma::lmFit(eset1, design = design) %>% 
   limma::eBayes() %>% 
   limma::topTable(coef = c('IPF', 'RA_ILD', 'RA_NOILD')) %>% 
   rownames() %>% extract(1) %>%    magrittr::extract(eset1, ., ) %>% 
   autonomics.plot::plot_feature_boxes(x = 'subgroup', color_var = 'GENDER')
