# Project: comics
# 
# Author: Aditya Bhagwat
###############################################################################

# Load packages
require(magrittr)
require(UniProt.ws)
require(org.Hs.eg.db)

# Function to rm NAs in a vector
rmNA <- function(aVector){
  aVector[!is.na(aVector)]
}

# Function to start a uniprot web service
start_uniprot_ws <- function(organism){
  sprintf('\nStarting UniProt.ws') %>% message()
  dt <- system.time({
    uniprot_ws <- UniProt.ws::UniProt.ws(taxId=as.numeric(ANNOTATED_ORGANISMS[[organism]]['taxonid']))
  })
  sprintf('Done in %s min', ceiling(dt[['elapsed']]/60)) %>% message()
  uniprot_ws
}

# Function to create a map of entrez gene ids to uniprot acessions
make_eg2uniprot <- function(organism){
  sprintf('\nAssembling eg2uniprot map') %>% message()
  sprintf(
      'org.%s.%s.db',
      ANNOTATED_ORGANISMS[[organism]]['abrev'],
      ANNOTATED_ORGANISMS[[organism]]['source']) %T>%
    autonomics.annotate::install_package_if_necessary() %>%
    library(character.only = TRUE)
  sprintf(
      'org.%s.%sUNIPROT',
      ANNOTATED_ORGANISMS[[organism]]['abrev'],
      ANNOTATED_ORGANISMS[[organism]]['source']) %>%
    get() %>%
    AnnotationDbi::as.list()                                   %>% 
    AnnotationDbi::unlist2()                                   %>% 
    rmNA()
}
# eg2uniprot <- make_eg2uniprot('human')

# Function to create a map of uniprot accessions to keywords
make_accessions2keywords <- function(uniprot_accessions, uniprot_ws){
  sprintf('\nAssembling uniprot accessions to keywords map') %>% message()
  dt <- system.time({
    accessions2keywords <- UniProt.ws::select(uniprot_ws, uniprot_accessions, 'KEYWORDS')
  })
  sprintf('Done in %s min', ceiling(dt[['elapsed']]/60)) %>% message()
  accessions2keywords
}

make_keywords2accessions <- function(accessions2keywords, organism){
  
  single_keyword_to_accessions <- function(keyword){
    uniprot_accessions <- accessions2keywords$UNIPROTKB
    uniprot_keywords   <- strsplit(accessions2keywords$KEYWORDS, '; ')
    uniprot_accessions[grep(keyword, uniprot_keywords, fixed = TRUE)]
  }
  
  sprintf('\nCreating map of uniprot keywords to accessions') %>% message()
  my_time <- system.time({
    available_keywords <- accessions2keywords  %>% 
      extract2('KEYWORDS')  %>% 
      strsplit('; ')        %>% 
      unlist()              %>% 
      unique()
    keywords2accessions <- Map(single_keyword_to_accessions, available_keywords)    # 6 minutes
  }) %>% extract2('elapsed')
  
  sprintf('Done in %s min', ceiling(my_time/60)) %>% message()
  #names(accessions2keywords) <- available_uniprot_keywords
  attr(keywords2accessions, 'date') <- lubridate::today()
  mapping_package <- sprintf(
     'org.%s.%s.db',
     ANNOTATED_ORGANISMS[[organism]]['abrev'],
     ANNOTATED_ORGANISMS[[organism]]['source'])
  attr(keywords2accessions, 'source_of_keys') <- list(
     package = mapping_package,
     version = mapping_package %>%
        utils::packageVersion())
  keywords2accessions
}

# head(eg2uniprot)
# map_uniprot2eg <- setNames(names(eg2uniprot), eg2uniprot)

# Map uniprot IDs to keywords: 40 minutes (Tuesday 4 August around 18h30)
make_uniprot_mappings <- function(organism){
  
  sprintf('-----\n%s\n-----', ANNOTATED_ORGANISMS[[organism]]['abrev']) %>% message()
  uniprot_ws <- start_uniprot_ws(organism)
  file_name <- sprintf('inst/extdata/uniprot_%s_keywords_to_accessions.rds', ANNOTATED_ORGANISMS[[organism]]['kegg'])
  
  make_eg2uniprot(organism)                %>%
    make_accessions2keywords(uniprot_ws)   %>%
    dplyr::filter(!is.na(KEYWORDS))        %>% 
    make_keywords2accessions(organism = organism)             %>%
    saveRDS(file = file_name)              %>%
    message('\n')
}

Map(make_uniprot_mappings, names(ANNOTATED_ORGANISMS))

# rm uninformative uniprot keywords
rm_uninformative <- function(organism, verbose = TRUE){
  file <- sprintf('inst/extdata/uniprot_%s_keywords_to_accessions.rds', ANNOTATED_ORGANISMS[[organism]]['kegg'])
  upkw <- readRDS(file); n0 <- length(upkw)
  uninformative <- c('3D-structure', 'Direct protein sequencing', 'Complete proteome', 'Reference proteome')
  for (keyword in uninformative){
    upkw[[keyword]] <- NULL
  }
  if (verbose){
    message(sprintf('%d - %d = %d %s uniprot keywords', n0, length(uninformative), length(upkw), organism))
  }
  saveRDS(upkw, file)
}
Map(rm_uninformative, names(ANNOTATED_ORGANISMS))