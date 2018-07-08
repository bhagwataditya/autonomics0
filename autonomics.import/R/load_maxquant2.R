#' Extract injection values from maxquant file
#' @param file 'proteinGroups.txt'
#' @return character vector: injection values
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_injections()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt',package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_injections()
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_injections()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_injections <- function(file){
   file %>%
   autonomics.support::cfread() %>%
   names() %>%
   magrittr::extract(stringi::stri_detect_fixed(., 'Razor + unique peptides ')) %>%
   stringi::stri_replace_first_fixed('Razor + unique peptides ', '')
}

#' Extract channel values from maxquant file
#' @param file string: path to proteinGroups.txt file (or other maxquant file)
#' @return character vector: maxquant channel values
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_channels()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt',package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_channels()
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_channels()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_channels <- function(file){
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   file %>% autonomics.support::cfread() %>%
            names()                      %>%
            # channel specific intensities available, except for reporter intensities
            magrittr::extract(autonomics.support::vstri_detect_fixed(., c('Intensity ', 'Reporter intensity corrected'))) %>%
            autonomics.support::vstri_replace_first_fixed(c('Intensity ', 'Reporter intensity corrected'), '')            %>%
            autonomics.support::vstri_replace_first_fixed(injections, '')                                                 %>%
            trimws()                     %>%
            unique()                     %>%
            setdiff("")
}


#' Extract maxquant snames
#' @param file full path to proteinGroups.txt
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_snames() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_snames <- function(file){
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()
   injections %>% lapply(function(x) sprintf('%s[%s]', x, channels)) %>% unlist()
}


#' Extract maxquant fnames
#' @param file full path to proteinGroups.txt
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_fnames() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_fnames <- function(file){
   file %>%
   autonomics.support::cfread() %>%
   magrittr::extract2('Majority protein IDs') %>%
   stringi::stri_split_fixed(';') %>%
   vapply(extract, character(1), 1)
}


#' Extract maxquant intensity colnames
#' @param file full path to proteinGroups.txt
#' @return character vector with intensity column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_intensity_colnames <- function(file){
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()
   injections %>% lapply(function(x) sprintf('Intensity %s %s', channels, x)) %>% unlist()
}

#' Extract maxquant ratio colnames
#' @param file full path to proteinGroups.txt
#' @param normalized logical
#' @return character vector with intensity column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_ratio_colnames()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_ratio_colnames <- function(file, normalized = FALSE){
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()
   ratios <- character(0)
   for (i in length(channels):2){
      ratios %<>% c(sprintf('%s/%s', channels[i], channels[1:(i-1)]))
   }
   autonomics.support::vsprintf('Ratio %s%s %s',
                                ratios, ifelse(normalized, ' normalized', ''), injections,
                                first_slowest = FALSE)
}



#' Extract maxquant exprs
#' @param file full path to proteinGroups.txt
#' @param what 'ratios', 'normalized ratios', 'intensities', 'lfq intensities', or 'reporter intensities'
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_exprs('intensities') %>% magrittr::extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_exprs <- function(file, what){

   fnames1    <- file %>% autonomics.import::extract_maxquant_fnames()
   snames1    <- file %>% autonomics.import::extract_maxquant_snames()
   colnamefun <- switch(what,
                       `intensities`       = autonomics.import::extract_maxquant_intensity_colnames,
                       `normalized ratios` = function(file) autonomics.import::extract_maxquant_ratio_colnames(file, normalized = TRUE),
                       `ratios`            = function(file) autonomics.import::extract_maxquant_ratio_colnames(file, normalized = FALSE))

   file %>%
   autonomics.support::cfread(select = colnamefun(.)) %>%
   data.matrix() %>%
   magrittr::set_rownames(fnames1) %>%
   magrittr::set_colnames(snames1)
}




#' Load proteingroups snames
#' @param file full path to proteinGroups.txt
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt'  %>%
#'              system.file(package = 'autonomics.data')
#'     file %>% autonomics.import::infer_maxquant_type() %>% head()
#' }
#' if (require(graumann.lfq)){
#'     file <- 'extdata/proteinGroups.txt'  %>%
#'              system.file(package = 'graumann.lfq')
#' }
#' @importFrom magrittr %>%
#' @export
infer_maxquant_value_pattern <- function(file){
   cols <- file %>% autonomics.support::cfread() %>% names()

   # Normalized Ratios
   pattern <- 'Ratio .+ normalized '
   if (any(stringi::stri_detect_regex(cols, pattern)))  return(pattern)

   # LFQ intensities
   pattern <- 'LFQ intensity '
   if (any(stringi::stri_detect_fixed(cols, pattern)))  return(pattern)

   # Intensities

}

#' Load proteingroups fnames
#' @param file full path to proteinGroups.txt
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    'extdata/stemcell.comparison/maxquant/proteinGroups.txt'  %>%
#'     system.file(package = 'autonomics.data')                 %>%
#'     autonomics.import::load_proteingroups_fnames()           %>%
#'     head()
#' }
#' @importFrom magrittr %>%
#' @export
load_proteingroups_fnames <- function(file){
   x <- file                                          %>%
        autonomics.support::cfread()                  %>%
        #magrittr::extract2('id')                     %>%
        magrittr::extract2('Majority protein IDs')    %>%
        stringi::stri_split_fixed(';')                %>%
        vapply(extract, character(1), 1)
   assertive::assert_has_no_duplicates(x)
   x
}

#' Load proteingroups exprs
#' @param file full path to proteinGroups.txt
#' @param value_type value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'
#'    # Normalized Ratios
#'    pattern <- 'Ratio (.+) normalized (.+)'
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt'  %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::load_proteingroups_exprs(pattern) %>% str()
#'
#'    # Ratios
#'    pattern <- 'Ratio .+ .'
#' }
#' @importFrom magrittr %>%
#' @export
load_proteingroups_exprs <- function(file, pattern){
   dt <- file %>% autonomics.support::cfread()
   fnames1 <- dt %>% magrittr::extract2('Majority protein IDs') %>%
                     stringi::stri_split_fixed(';') %>%
                     vapply(magrittr::extract, character(1), 1)
   dt %>% magrittr::extract(, stringi::stri_detect_regex(names(.), pattern), with = FALSE) %>%
          data.matrix() %>%
          magrittr::set_rownames(fnames1)
}

#' Load proteingroups
#' @param file full path to proteinGroups.txt
#' @examples
#' file <- system.file('extdata/stemcell.comparison/maxquant/proteinGroups.txt',
#'                      package = 'autonomics.data')
#' @importFrom magrittr %>%
#' @export
load_proteingroups2 <- function(file){
   DT <- autonomics.support::cfread(file)
}
