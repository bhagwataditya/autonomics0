
#' Load metabolon lipidomics
#' @param file path to metabolon lipidomics file
#' @param entity 'Lipid Class', 'Species', or 'Fatty Acid'
#' @param quantity 'Concentrations' or 'Compositions'
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% autonomics.import::load_metabolon_lipidomics(entity = 'Lipid Class', quantity = 'Concentrations')
#'    file %>% autonomics.import::load_metabolon_lipidomics('Lipid Class', 'Compositions')
#'    file %>% autonomics.import::load_metabolon_lipidomics('Species',     'Concentrations')
#'    file %>% autonomics.import::load_metabolon_lipidomics('Species',     'Compositions')
#'    file %>% autonomics.import::load_metabolon_lipidomics('Fatty Acid',  'Concentrations')
#'    file %>% autonomics.import::load_metabolon_lipidomics('Fatty Acid',  'Compositions')
#' }
#' @importFrom magrittr %>%
#' @export
load_metabolon_lipidomics <- function(file, entity = 'Species', quantity = 'Concentrations'){

   # Assert
   assertive.sets::assert_is_subset(entity,   c('Lipid Class', 'Species', 'Fatty Acid'))
   assertive.sets::assert_is_subset(quantity, c('Concentrations', 'Compositions'))

   # Read
   x <- file %>% readxl::read_excel(sheet = paste0(entity, ' ', quantity))

   # sdata
   row1 <- which(x[[2]]=='Client Identifier')
   coln <- which(x[row1, ] == 'Unit')
   sdata1 <- x %>% magrittr::extract((1+row1):nrow(x), 1:coln) %>%
                   data.frame() %>%
                   magrittr::set_names(x[row1, 1:coln] %>% unname() %>% unlist())

   # exprs
   fnamesrow <- if (entity == 'Lipid Class') row1 else row1-1
   exprs1 <- x %>% magrittr::extract((1+row1):nrow(x), (coln+1):ncol(x))  %>%
                   magrittr::set_names(x[fnamesrow, (coln+1):ncol(x)])    %>%
                  (function(y) {y[y=='.'] <- NA;  y})                     %>%
                   data.matrix()                                          %>%
                   magrittr::set_rownames(sdata1$`Client Identifier`)     %>%
                   t()
   # fdata
   fdata1 <- data.frame(feature_id = x %>% magrittr::extract(fnamesrow, (coln+1):ncol(x)) %>% unname() %>% unlist(),
                        lipidclass = x %>% magrittr::extract(row1,      (coln+1):ncol(x)) %>% unname() %>% unlist()) %>%
             magrittr::set_rownames(.$feature_id)

   # sumexp
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs1))
   autonomics.import::sdata(object) <- sdata1
   autonomics.import::fdata(object) <- fdata1
   object
}

