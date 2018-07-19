#=====================
# VALIDIFY SAMPLE IDS
#=====================

#' Validify sample/feature names
#' @param x character vector with sample ids
#' @return character vector
#' @importFrom magrittr %<>%
#' @export
validify_sample_ids <- function(x){
   . <- NULL
   selector <- substr(x,1,1) %in% 0:9
   x[selector] %<>% paste0('s', .)
   x
}



#=====================
# LOAD_SDATA_MAXQUANT
#=====================

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


#' Designify maxquant sampleids
#' @param sampleids character vector
#' @examples
#' require(magrittr)
#' sampleids <- c("STD(L).EM00(M).EM01(H).R1[H/L]",
#'                "STD(L).EM00(M).EM01(H).R1[M/L]")
#' sampleids %>% autonomics.import::designify_maxquant_sampleids()
#'
#' sampleids <- c("Gel 11 1", "Gel 11 2", "Gel Ctrl 1")
#' sampleids %>% autonomics.import::designify_maxquant_sampleids()
#'
#' sampleids <- c("ESC(0).NCM(1).CM(2).MV(3).EX(4).KIT(5).R1[0]",
#'                "ESC(0).NCM(1).CM(2).MV(3).EX(4).KIT(5).R1[1]")
#' sampleids %>% autonomics.import::designify_maxquant_sampleids()
#' @importFrom magrittr %>%
#' @export
designify_maxquant_sampleids <- function(sampleids){

   if (is.factor(sampleids)) sampleids %<>% as.character()

   # Break into parts         #pattern <- '(.*)\\.R\\((.+)\\)\\[(.+)\\]'
   sep <- sampleids %>% autonomics.import::infer_design_sep()
   parts <- strsplit(sampleids, split = sep, fixed = TRUE)
   n <- length(parts[[1]])

   # replicate values
   replicate_values <- parts %>% vapply(extract, character(1), n)
   is_labeled <- parts %>% vapply(extract, character(1), n) %>% stringi::stri_detect_fixed('[') %>% any()
   if (is_labeled){
      label_values <- parts %>% vapply(extract, character(1), n) %>%
         stringi::stri_extract_first_regex('(?<=\\[)(.+)(?=\\])')  %>%
         strsplit(split = '/', fixed = TRUE)
      replicate_values %<>% stringi::stri_replace_first_regex('(.+)(\\[.+\\])', '$1')
      replicate_values %<>% paste0('_', label_values %>% vapply(paste0, character(1), collapse = ''))
   }
   parts %<>% lapply(function(x)x %>% magrittr::extract(1:(length(x)-1)))

   # subgroup values
   subgroup_values <- if (is_labeled){
      parts %>% lapply(function(x){
         cursample <- x %>% stringi::stri_replace_first_regex('(.+)\\(([HML0-9]+)\\)', '$1')
         curname   <- x %>% stringi::stri_replace_first_regex('(.+)\\(([HML0-9]+)\\)', '$2')
         cursample %>% magrittr::set_names(curname)
      }) %>%
         autonomics.support::vextract(label_values) %>%
         vapply(paste0, character(1), collapse = '_')
   } else {
      parts %>% vapply(paste0, character(1), collapse = '.')
   }

   # sampleid values
   long_sampleid_values <- sprintf('%s.%s', subgroup_values, replicate_values)
   sampleid_values <- long_sampleid_values %>% stringi::stri_replace_first_regex('_[HML0-9]+', '')
   idx <- sampleid_values %>% autonomics.support::cduplicated()
   sampleid_values[idx] <- long_sampleid_values[idx]
   sampleid_values
}


#' Load maxquant snames
#' @param file path to maxquant file
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'
#'    file %>% load_snames_maxquant('Ratio')
#'    file %>% load_snames_maxquant('Ratio', clean = TRUE)
#'    file %>% load_snames_maxquant('Ratio', clean = TRUE, designify = TRUE)
#'
#'    file %>% load_snames_maxquant('Intensity')
#'    file %>% load_snames_maxquant('Intensity', clean = TRUE)
#'    file %>% load_snames_maxquant('Intensity', clean = TRUE, designify = TRUE)
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% load_snames_maxquant('Intensity')
#'    file %>% load_snames_maxquant('Intensity', clean = TRUE)
#'    file %>% load_snames_maxquant('Intensity', clean = TRUE, designify = TRUE)
#'    file %>% load_snames_maxquant('LFQ intensity')
#'    file %>% load_snames_maxquant('LFQ intensity', clean = TRUE)
#'    file %>% load_snames_maxquant('LFQ intensity', clean = TRUE, designify = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
load_snames_maxquant <- function(
   file,
   quantity  = autonomics.import::infer_maxquant_quantity(file),
   clean     = FALSE,
   designify = FALSE
){

   # Extract
   #--------
   quantity_colnames <- if (quantity %in% c('Ratio', 'Ratio normalized')){
                           file %>% autonomics.import::extract_maxquant_ratio_colnames(quantity)
                        } else {
                           file %>% autonomics.import::extract_maxquant_intensity_colnames(quantity)
                        }
   if (!clean) return(quantity_colnames)

   # Clean
   #------
   # 'Ratio normalized' is spread out into 'Ratio H/L normalized'. This requires a separate approach.
   if (quantity %in% c('Ratio', 'Ratio normalized')){
      snames1 <- quantity_colnames %>% stringi::stri_replace_first_regex('Ratio (./.) (?:normalized )?(.+)',   '$2[$1]')

   # 'Intensity', 'LFQ intensity', 'Reporter intensity' can all be treated in a similar fashion
   } else {
      injections <- file %>% autonomics.import::extract_maxquant_injections()
      channels   <- file %>% autonomics.import::extract_maxquant_channels()

      # For unlabeled data, sample names are simply injections
      if (length(channels)==0){
         snames1 <- quantity_colnames %>% stringi::stri_replace_first_fixed(paste0(quantity, ' '), '')

      # For labeled data 'stringi::stri_replace_first_xxx' does not work, because
      # ' ' separators in sample names are difficult to differentiate from ' H' constructs
      # Instead we generate them from channels and injections, but we make sure the order matches the actual order.
      } else {
         idx <- autonomics.support::vsprintf('%s %s %s', quantity, channels, injections) %>%
                match(quantity_colnames)
         snames1 <- autonomics.support::vsprintf('%s[%s]', injections, channels) %>% magrittr::extract(idx)
      }
   }
   if (!designify) return(snames1)

   # Designify
   #----------
   snames1 %>% autonomics.import::designify_maxquant_sampleids()

}


#' Load maxquant sdata
#' @param file path to maxquant file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemcell.comparison/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% load_sdata_maxquant()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% load_sdata_maxquant()
#' }
#' @importFrom magrittr %>%
#' @export
load_sdata_maxquant <- function(
   file,
   quantity = autonomics.import::infer_maxquant_quantity(file)
){
   data.frame(sample_id = file %>% autonomics.import::load_snames_maxquant(quantity, clean = TRUE)) %>%
   magrittr::set_rownames(.$sample_id)
}


#==========================================
# LOAD_SDATA_METABOLON
#==========================================

#' Load metabolon sdata
#' @param file path to metabolon file
#' @param sheet excell sheet number or name
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% load_sdata_metabolon(2) %>% extract(1:3, 1:3)
#' }
#' @importFrom magrittr %>%
#' @export
load_sdata_metabolon <- function(file, sheet){

   # Load sample data
   . <- NULL
   df <- file %>% readxl::read_excel(sheet = sheet, col_names = FALSE)
   sstart <- which(!is.na(df[1,]))[1]
   fstart <- which(!is.na(df[,1]))[1]
   df %<>% magrittr::extract(1:fstart, (sstart+1):ncol(.)) %>%
      t() %>%
      data.frame(stringsAsFactors = FALSE, check.names = FALSE)            %>%
      magrittr::set_names(df[[sstart]][1:fstart])     %>%
      magrittr::set_names(names(.) %>% stringi::stri_replace_first_fixed( 'Group   HMDB_ID', 'Group') %>% # recent metabolon files
                             stringi::stri_replace_first_fixed('Sample   HMDB_ID', 'Group'))    # older metabolon files
   df %<>% magrittr::set_rownames(.$CLIENT_IDENTIFIER)

   # Return
   df
}


#==========================================
# LOAD_SDATA_METABOLONLIPIDS
#==========================================

#' Possible sheet values for metabolonlipids data
#' @export
METABOLONLIPIDS_SHEETS <- c('Lipid Class Concentrations',
                            'Lipid Class Compositions',
                            'Species Concentrations',
                            'Species Compositions',
                            'Fatty Acid Concentrations',
                            'Fatty Acid Compositions')

#' Load metabolonlipids sdata
#'
#' Load sdata from metabolon clp (complex lipid panel) file
#' @param file  path to metabolon clp (complex lipid panel) file
#' @param sheet name of excel sheet (any value in METABOLONLIPIDS_SHEETS)
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_sdata_metabolonlipids('Lipid Class Concentrations') %>% head()
#'    file %>% load_sdata_metabolonlipids(    'Species Concentrations') %>% head()
#'    file %>% load_sdata_metabolonlipids( 'Fatty Acid Concentrations') %>% head()
#'    file %>% load_sdata_metabolonlipids('Lipid Class Compositions')   %>% head()
#'    file %>% load_sdata_metabolonlipids(    'Species Compositions')   %>% head()
#'    file %>% load_sdata_metabolonlipids( 'Fatty Acid Compositions')   %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_sdata_metabolonlipids <- function(file, sheet){
   x <- file %>% readxl::read_excel(sheet = sheet)
   row1 <- which(x[[2]]=='Client Identifier')
   coln <- which(x[row1, ] == 'Unit')

   x %>% magrittr::extract((1+row1):nrow(x), 1:coln) %>%
      data.frame() %>%
      magrittr::set_names(x[row1, 1:coln] %>% unname() %>% unlist()) %>%
      magrittr::set_rownames(.$`Client Identifier`)
}


#==========================================
# LOAD_SDATA_SOMA
#==========================================

#' Load soma sdata
#' @param file string: path to adat file
#' @return sample dataframe
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                         package = 'autonomics.data')
#'    file %>% autonomics.import::load_sdata_soma() %>% head()
#' }
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-14-130_Set_A_RPT.HybMedNormCal_20140925.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::load_sdata_soma() %>% head()
#' }
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
load_sdata_soma <- function(file){
   x <- file %>% autonomics.import::identify_soma_structure()

   file %>% autonomics.support::cfread(header = FALSE, sep = '\t', fill = TRUE) %>%
            magrittr::extract((x$row-1):nrow(.),  1:(x$col-2), with = FALSE)    %>%
            magrittr::set_names(unlist(unname(.[1,])))                          %>%
            magrittr::extract(-1, )                                             %>%
            data.frame(row.names = .$SampleId)
}


#==========================================
# LOAD_SDATA GENERIC
#==========================================

#' Load sdata
#' @param file path to omics data file
#' @param platform 'maxquant', 'metabolon', 'metabolonlipids', 'soma'
#' @param sheet excel sheet number or name if applicable
#' @param quantity string: which quantity should be extracted (only applicable for maxquant platform)
#' @return sample dataframe
#' @examples
#'  require(magrittr)
#'  if (require(autonomics.data)){
#'     file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'     file %>% load_sdata(2, 'metabolon') %>% extract(1:3, 1:3)
#'  }
#' file <- '../../datasets/WCQA-01-18MLCLP-1/WCQA-01-18MLCLP CLP  6-TAB FILE (180710).XLSX'
#' if (file.exists(file)){
#'    file %>% load_sdata('Lipid Class Concentrations', 'metabolonlipids') %>% head()
#' }
#' @importFrom magrittr %>%
#' @export
load_sdata <- function(file, platform, sheet = NULL, quantity = NULL){
   switch(platform,
          maxquant        = file %>% load_sdata_maxquant(quantity = quantity),
          metabolonlipids = file %>% load_sdata_metabolonlipids(sheet = sheet),
          metabolon       = file %>% load_sdata_metabolon(sheet = sheet),
          soma            = file %>% load_sdata_soma())
}
