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
#' @param quantity 'Intensity', 'LFQ intensity', 'Reporter intensity'
#' @return character vector with intensity column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames() %>% head(3)
#'    file %>%
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames() %>% head(3)
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames('LFQ intensity') %>% head(3)
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_intensity_colnames('Reporter intensity') %>% head(3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_intensity_colnames <- function(file, quantity = 'Intensity'){

   # Asser
   assertive.files::assert_all_are_existing_files(file)
   assertive.sets::assert_is_subset(quantity, c('Intensity', 'LFQ intensity', 'Reporter intensity'))

   # Deduce injections and channels from unambiguous peptide columns
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()

   # Construct intensity colnames
   intensity_colnames <- if (length(channels)==0){                      sprintf('%s %s',    quantity,           injections)
                         } else {                  autonomics.support::vsprintf('%s %s %s', quantity, channels, injections) }

   # Ensure identical order as in actual file
   names(autonomics.support::cfread(file)) %>% magrittr::extract(. %in% intensity_colnames)
}

#' Extract maxquant ratio colnames
#' @param file full path to proteinGroups.txt
#' @param quantity 'Ratio' or 'Ratio normalized'
#' @return character vector with ratio column names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_ratio_colnames() %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_ratio_colnames <- function(file, quantity){

   assertive.sets::assert_is_subset(quantity, c('Ratio', 'Ratio normalized'))

   # Deduce injections and channels from unambiguous peptide columns
   injections <- file %>% autonomics.import::extract_maxquant_injections()
   channels   <- file %>% autonomics.import::extract_maxquant_channels()

   # Construct possible ratio colnames
   possible_ratio_columns <- autonomics.support::vsprintf('Ratio %s/%s%s %s',
                                                          channels,
                                                          channels,
                                                          if(quantity == 'Ratio normalized') ' normalized' else '',
                                                          injections)
   # Return actual colnames in correct order
   file %>%
   autonomics.support::cfread() %>%
   names() %>%
   magrittr::extract(. %in% possible_ratio_columns)
}


#' Extract maxquant exprs
#' @param file full path to proteinGroups.txt
#' @param quantity 'Ratio', 'Ratio normalized', Intensity', 'LFQ intensity', 'Reporter intensity'
#' @return character vector with sample names
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::extract_maxquant_exprs(quantity = 'Intensity') %>%
#'             magrittr::extract(1:3, 1:3)
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::extract_maxquant_exprs('Intensity')     %>%
#'             magrittr::extract(1:3, 1:3)
#'    file %>% autonomics.import::extract_maxquant_exprs('LFQ intensity') %>%
#'             magrittr::extract(1:3, 1:3)
#' }
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::extract_maxquant_exprs('Reporter intensity') %>%
#'             magrittr::extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
extract_maxquant_exprs <- function(file, quantity = 'Intensity'){

   # Assert
   assertive.sets::assert_is_subset(quantity, c('Ratio', 'Ratio normalized', 'Intensity', 'LFQ intensity', 'Reporter intensity'))

   # Construct maxquant colnames
   col_names <- if (quantity %in% c('Ratio', 'Ratio normalized')){ file %>% autonomics.import::extract_maxquant_ratio_colnames(quantity)
                } else {                                           file %>% autonomics.import::extract_maxquant_intensity_colnames(quantity)
                }

   # Extract exprs matrix
   exprs_mat <- file %>%
                autonomics.support::cfread(select = col_names) %>%
                data.matrix() %>%
                magrittr::set_rownames(autonomics.import::extract_maxquant_fnames(file))

   # Rename samples
   if (quantity %in% c('Ratio', 'Ratio normalized')){
     colnames(exprs_mat) %<>% stringi::stri_replace_first_regex('Ratio (./.) (?:normalized )?(.+)',   '$2[$1]')

   } else {
     # It is better to do it this way (rather than use a stri_replace_regex as for the ratios)
     # because ' ' separators in sample names are difficult to differentiate from ' H' constructs
     injections <- file %>% autonomics.import::extract_maxquant_injections()
     channels   <- file %>% autonomics.import::extract_maxquant_channels()
     colnames(exprs_mat) <- if (length(channels)==0) injections else autonomics.support::vsprintf('%s[%s]', injections, channels)
   }

   # Return
   exprs_mat
}

#' Extract maxquant fdata
#' @param file path to proteinGroups.txt file
#' @return dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#' }
#' @importFrom magrittr %>%
#' @export
extract_proteingroups_fdata <- function(file){
   file %>%
   autonomics.support::cfread(select = c('Majority protein IDs', 'Gene names', 'Protein names')) %>%
   data.frame(stringsAsFactors = FALSE) %>%
   magrittr::set_rownames(autonomics.import::extract_maxquant_fnames(file))
}


#' Load proteingroups
#' @param file full path to proteinGroups.txt
#' @param quantity 'Ratio normalized', 'Ratio', 'Intensity', 'LFQ intensity', 'Reporter intensity'
#' @param design_file full path to design file (created with write_maxquant_design)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/maxquant/proteinGroups.txt',
#'                         package = 'autonomics.data')
#'    quantity <- 'Ratio normalized'
#'    file %>%
#' }
#' @importFrom magrittr %>%
#' @export
load_proteingroups2 <- function(file, quantity, design_file){

   exprs_mat <- file %>% autonomics.import::extract_maxquant_exprs(quantity)
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs_mat))
   autonomics.import::fdata(object) <- file %>% autonomics.import::extract_proteingroups_fdata()
   autonomics.import::sdata(object) <- file %>% autonomics.import::read_maxquant
   object

}
