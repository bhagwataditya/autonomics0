#=================================================
# GENERIC
#=================================================

nrows <- function(x, sheet=1){
   if (tools::file_ext(x) %in% c('xls', 'xlsx')){
      x %>% readxl::read_excel(sheet = sheet, .name_repair = 'minimal', col_names = FALSE) %>%  nrow()
   } else {
      x %>% readLines(warn=FALSE) %>% length()
   }
}


ncols <- function(x, sheet=1){
   if (tools::file_ext(x) %in% c('xls', 'xlsx')){
      x %>% readxl::read_excel(sheet=sheet, .name_repair = 'minimal', col_names = FALSE) %>% ncol()
   } else {
      x %>% utils::count.fields(quote = "", sep = '\t') %>% max()
   }
}

#=================================================
# extract_rectangle
#=================================================

#' Extract rectangle from omics file, data.table, or matrix
#'
#' @param x          omics datafile or datatable
#' @param rows       numeric vector
#' @param cols       numeric vector
#' @param verbose    logical
#' @param transpose  logical
#' @param drop       logical
#' @param sheet      numeric or string
#' @param ...        allow for S3 method dispatch
#' @return matrix
#' @examples
#' require(magrittr)
#' if (require('autonomics.data')){
#'
#'    # filepath
#'       x <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'       x %>% extract_rectangle(rows = 11:401, cols = 15:86, sheet = 2) %>% extract(1:3, 1:3) # exprs
#'       x %>% extract_rectangle(rows = 11:401, cols = 5,     sheet = 2) %>% extract(1:3, )    # fids
#'       x %>% extract_rectangle(rows = 2,      cols = 15:86, sheet = 2) %>% extract(,1:3)     # sids
#'       x %>% extract_rectangle(rows = 10:401, cols = 1:14,  sheet = 2) %>% extract(1:3, 1:3) # fdata
#'       x %>% extract_rectangle(rows = 1:10,   cols = 14:86, sheet = 2,
#'                               transpose = TRUE)                       %>% extract(1:3, 1:3) # sdata
#'    # matrix
#'       x <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data') %>%
#'            extract_rectangle(sheet = 2)
#'       x %>% extract_rectangle(rows = 11:401, cols = 15:86, sheet = 2) %>% extract(1:3, 1:3) # exprs
#'       x %>% extract_rectangle(rows = 11:401, cols = 5,     sheet = 2) %>% extract(1:3, )    # fids
#'       x %>% extract_rectangle(rows = 2,      cols = 15:86, sheet = 2) %>% extract(,1:3)     # sids
#'       x %>% extract_rectangle(rows = 10:401, cols = 1:14,  sheet = 2) %>% extract(1:3, 1:3) # fdata
#'       x %>% extract_rectangle(rows = 1:10,   cols = 14:86, sheet = 2,
#'                               transpose = TRUE)                       %>% extract(1:3, 1:3) # sdata
#' }
#' @importFrom magrittr %>% %<>%
#' @export
extract_rectangle <- function(x, ...){
   UseMethod('extract_rectangle')
}

#' @rdname extract_rectangle
#' @importFrom magrittr %>%
#' @export
extract_rectangle.character <- function(
   x,
   rows      = 1:nrows(x, sheet=sheet),
   cols      = 1:ncols(x, sheet=sheet),
   verbose   = FALSE,
   transpose = FALSE,
   drop      = FALSE,
   sheet     = 1,
   ...
){

   # Assert
   assertive.files::assert_all_are_existing_files(x)
   assertive.types::assert_is_numeric(rows)
   assertive.types::assert_is_numeric(cols)
   row1 <- min(rows); rown <- max(rows)
   col1 <- min(cols); coln <- max(cols)

   # Read file
   dt <- if (tools::file_ext(x) %in% c('xls', 'xlsx')){
            x %>%
            readxl::read_excel(sheet       = sheet,
                               col_names   = FALSE,
                               range       = sprintf('R%dC%d:R%dC%d', row1, col1, rown, coln),
                              .name_repair = 'minimal') %>%
            data.table::data.table()

         } else {
            x %>%
            data.table::fread(na.strings = "",
                              header     = FALSE,
                              integer64  = 'numeric',
                              skip       = row1-1,
                              nrow       = 1+rown-row1,
                              select     = col1:coln)
         }

   # Extract rectangle
   dt %>% extract_rectangle.data.table(transpose = transpose, drop = drop)

}

# Extract row
extract_dt_row <- function(dt, i) dt %>% magrittr::extract(i,) %>% as.matrix() %>% magrittr::extract(1, ) %>% unname()
extract_dt_col <- function(dt, i) dt[[i]]


#' @rdname extract_rectangle
#' @importFrom magrittr %<>% %>%
#' @export
extract_rectangle.data.table <- function(
   x,
   rows = 1:nrow(x),
   cols = 1:ncol(x),
   transpose = FALSE,
   drop = FALSE,
   ...
){
   x %>%
   magrittr::extract(rows, cols, with = FALSE) %>%
   as.matrix() %>%
   extract_rectangle.matrix(transpose = transpose, drop = drop)
}


#' @rdname extract_rectangle
#' @importFrom magrittr %<>% %>%
#' @export
extract_rectangle.matrix <- function(
   x,
   rows      = 1:nrow(x),
   cols      = 1:ncol(x),
   transpose = FALSE,
   drop      = FALSE,
   ...
){
   rectangle <- x %>% magrittr::extract(rows, cols, drop = FALSE)
   if (transpose) rectangle %<>% t()
   if (drop) if (nrow(rectangle)==1 | ncol(rectangle)==1) rectangle %<>% as.vector('character')
   rectangle
}


#' Read omics data from file
#'
#' @param file           string: name of text (txt, csv, tsv, adat) or excel (xls, xlsx) file
#' @param sheet          integer/string: only relevant for excel files
#' @param fid_rows       numeric vector: featureid rows
#' @param fid_cols       numeric vector: featureid cols
#' @param sid_rows       numeric vector: sampleid rows
#' @param sid_cols       numeric vector: sampleid cols
#' @param expr_rows      numeric vector: expr rows
#' @param expr_cols      numeric vector: expr cols
#' @param fvar_rows      numeric vector: fvar rows
#' @param fvar_cols      numeric vector: fvar cols
#' @param svar_rows      numeric vector: svar rows
#' @param svar_cols      numeric vector: svar cols
#' @param fdata_rows     numeric vector: fdata rows
#' @param fdata_cols     numeric vector: fdata cols
#' @param sdata_rows     numeric vector: sdata rows
#' @param sdata_cols     numeric vector: sdata cols
#' @param transpose      logical
#' @param verbose        logical
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # RNASEQ
#'       file <- system.file('extdata/stemcomp/rnaseq/gene_counts.txt', package = 'autonomics.data')
#'       file %>% read_omics(fid_rows   = 2:58736,   fid_cols   = 1,
#'                           sid_rows   = 1,         sid_cols   = 4:11,
#'                           expr_rows  = 2:58736,   expr_cols  = 4:11,
#'                           fvar_rows  = 1,         fvar_cols  = 1:3,
#'                           fdata_rows = 2:58736,   fdata_cols = 1:3,
#'                           transpose  = FALSE)
#'
#'    # LCMSMS PROTEINGROUPS
#'       file <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'                system.file(package = 'autonomics.data')
#'       file %>% read_omics(fid_rows   = 2:9783,  fid_cols   = 383,
#'                           sid_rows   = 1,       sid_cols   = seq(124, 316, by = 6),
#'                           expr_rows  = 2:9783,  expr_cols  = seq(124, 316, by = 6),
#'                           fvar_rows  = 1,       fvar_cols  = c(2, 6, 7, 383),
#'                           fdata_rows = 2:9783,  fdata_cols = c(2, 6, 7, 383),
#'                           transpose  = FALSE)
#'
#'    # SOMASCAN
#'       file <- system.file('extdata/stemcomp/soma/stemcomp.adat', package = 'autonomics.data')
#'       file %>% read_omics(fid_rows   = 21,       fid_cols   = 19:1146,
#'                           sid_rows   = 30:40,    sid_cols   = 4,
#'                           expr_rows  = 30:40,    expr_cols  = 19:1146,
#'                           fvar_rows  = 21:28,    fvar_cols  = 18,
#'                           svar_rows  = 29,       svar_cols  = 1:17,
#'                           fdata_rows = 21:28,    fdata_cols = 19:1146,
#'                           sdata_rows = 30:40,    sdata_cols = 1:17,
#'                           transpose  = TRUE)
#'
#'    # METABOLON
#'      file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'      file %>% read_omics(
#'                  sheet      = 2,
#'                  fid_rows   = 11:401,    fid_cols   = 5,
#'                  sid_rows   = 3,         sid_cols   = 15:86,
#'                  expr_rows  = 11:401,    expr_cols  = 15:86,
#'                  fvar_rows  = 10,        fvar_cols  = 1:14,
#'                  svar_rows  = 1:10,      svar_cols  = 14,
#'                  fdata_rows = 11:401,    fdata_cols = 1:14,
#'                  sdata_rows = 1:10,      sdata_cols = 15:86,
#'                  transpose  = FALSE)
#'
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_omics <- function(
   file,
   sheet = 1,
   fid_rows,           fid_cols,
   sid_rows,           sid_cols,
   expr_rows,          expr_cols,
   fvar_rows  = NULL,  fvar_cols  = NULL,
   svar_rows  = NULL,  svar_cols  = NULL,
   fdata_rows = NULL,  fdata_cols = NULL,
   sdata_rows = NULL,  sdata_cols = NULL,
   transpose  = FALSE,
   verbose    = TRUE
){
   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.types::assert_is_a_bool(transpose)

   # Read
   if (verbose) autonomics.support::cmessage('\t\tRead %s %s', if (sheet==1) '' else paste0("sheet '", sheet, "' of"), file)
   fixed_no_of_cols <- utils::count.fields(file, quote = "", sep = '\t') %>% (function(x)all(x==x[1]))
   x <- if (fixed_no_of_cols) extract_rectangle.character(file, sheet=sheet) else file

   # Extract exprs
   fids1  <- x %>% extract_rectangle(rows = fid_rows,  cols = fid_cols,  transpose = transpose, drop = TRUE,  sheet = sheet)
   sids1  <- x %>% extract_rectangle(rows = sid_rows,  cols = sid_cols,  transpose = transpose, drop = TRUE,  sheet = sheet)
   exprs1 <- x %>% extract_rectangle(rows = expr_rows, cols = expr_cols, transpose = transpose, drop = FALSE, sheet = sheet) %>% (function(y){class(y) <- 'numeric'; y})

   # Extract feature annotations
   #    Leave rownames(fdata1) empty: fids1 may contain non-valid values
   #    This happens in MaxQuant files, which sometimes contain missing rows (I think after opening in excell)
   fdata1 <- data.frame(feature_id = fids1, stringsAsFactors = FALSE)
   fdata_available <- !is.null(fvar_rows) & !is.null(fvar_cols)
   if (fdata_available){
      fvars1 <- x %>% extract_rectangle(fvar_rows,  fvar_cols,  transpose = transpose, drop = TRUE)
      fdata1 %<>% cbind(extract_rectangle(x, fdata_rows, fdata_cols, transpose = transpose, drop = FALSE, sheet = sheet) %>%
                        magrittr::set_colnames(fvars1) %>%
                        data.frame(stringsAsFactors = FALSE, check.names = FALSE))
   }

   # Extract sample annotations
   #    Leave rownames(sdata1) empty: sids may contain non-valid values
   #    This happens in SOMA files, where CLIENT_IDENTIFIER is not unique for calibrator and buffer samples
   sdata1 <- data.frame(sample_id = sids1, stringsAsFactors = FALSE)
   sdata_available <- !is.null(svar_rows) & !is.null(svar_cols)
   if (sdata_available){
      svars1 <- x %>% extract_rectangle(svar_rows,  svar_cols,  transpose =  transpose, drop = TRUE, sheet = sheet)
      sdata1 %<>% cbind(extract_rectangle(x, sdata_rows, sdata_cols, transpose = !transpose, drop = FALSE, sheet = sheet) %>%
                           magrittr::set_colnames(svars1) %>%
                           data.frame(stringsAsFactors = FALSE, check.names = FALSE))
   }

   # Rm features with missing ids (happens in MaxQuant files due to empty interspersed lines)
   idx <- !is.na(fids1)
   if (any(!idx)){
      if (verbose) autonomics.support::cmessage("\t\tRm %d features with missing 'feature_id' values", sum(!idx))
      fids1  %<>% magrittr::extract(idx)
      fdata1 %<>% magrittr::extract(idx,)
      exprs1 %<>% magrittr::extract(idx, )
   }

   # Rm samples with missing ids
   idx <- !is.na(sids1)
   if (any(!idx)){
      if (verbose) autonomics.support::cmessage("\t\tRm %d samples with missing 'sample_id' values", sum(!idx))
      sids1  %<>% magrittr::extract(idx)
      sdata1 %<>% magrittr::extract(idx,)
      exprs1 %<>% magrittr::extract(, idx)
   }

   # Name features and samples
   fids1 %<>% autonomics.support::uniquify('make.unique')
   sids1 %<>% autonomics.support::uniquify('make.unique')
   rownames(exprs1) <- rownames(fdata1) <- fids1
   colnames(exprs1) <- rownames(sdata1) <- sids1

   # Wrap into Sumexp
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs1))
   fdata(object) <- fdata1
   sdata(object) <- sdata1

   # Return
   object
}


#=========================================================
# RNASEQ
#=========================================================
#download_gtf
release_to_build <- function(release, organism){
   if        (organism == 'Homo sapiens'){          if (release >= 76)  'GRCh38'   else 'GRCh37'  }
   else if   (organism == 'Mus musculus'){          if (release >= 68)  'GRCm38'   else 'NCBIM37'  }
   else if   (organism == 'Rattus norvegicus'){     if (release >= 80)  'Rnor_6.0' else 'Rnor_5.0'  }
}

#' @examples
#' make_gtf_link(organism = 'Homo sapiens', release = 95)
#' make_gtf_link(organism = 'Mus musculus', release = 95)
#' @importFrom magrittr %>%
#' @export
make_gtf_link <- function(organism, release){
   sprintf('ftp://ftp.ensembl.org/pub/release-%s/gtf/%s/%s.%s.%s.gtf.gz',
           release,
           organism %>% tolower() %>% stringi::stri_replace_first_fixed(' ', '_'),
           organism %>%               stringi::stri_replace_first_fixed(' ', '_'),
           release_to_build(release, organism),
           release)
}


#' Download feature annotations
#'
#' Download feature annotations in GTF format
#' @param organism    'Homo sapiens', 'Mus musculus' or 'Rattus norvegicus'
#' @param release      GTF release. By default release 95 selected
#' @examples
#' \dontrun{
#'    download_gtf(organism = 'Homo sapiens')
#'    download_gtf(organism = 'Mus musculus')
#'    download_gtf(organism = 'Rattus norvegicus')
#' }
#' @export
download_gtf <- function(
   organism,
   release = 95,
   gtffile = sprintf("~/.autonomics/gtf/%s", basename(make_gtf_link(organism, release)))
){

   # Assert validity
   assertive.sets::assert_is_subset(organism, c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'))

   # Satisfy CHECK
   . <- NULL

   remote <- make_gtf_link(organism, release)

     if(file.exists(gtffile)){
      message(sprintf("\t\tGTF file already available at %s", gtffile))
   } else {
      message(sprintf("\t\tDownload and unzip GTF file to %s" , gtffile %>% substr(1, nchar(.)-3)))
      dir.create(dirname(gtffile), showWarnings = FALSE, recursive = TRUE)
      utils::download.file(url = remote, destfile = gtffile)
      R.utils::gunzip(gtffile,  remove = TRUE, overwrite = TRUE)
      gtffile %<>% stringi::stri_replace_last_fixed('.gz', '')
   }
   return(gtffile)
}


#get feature annotations
select_organism_database <- function(organism){
   if          (organism == 'Homo sapiens'){       'hsapiens_gene_ensembl'  }
   else if     (organism == 'Mus musculus'){       'mmusculus_gene_ensembl'  }
   else if     (organism == 'Rattus norvegicus'){  'rnorvegicus_gene_ensembl'  }
}


#' Read GTF file
#' @param gtffile       string: path to gtf file (which can be downloaded with download_gtf)
#' @param filter        Filter on feature name or feature id. By default all features are extracted.
#' @param return_value  Display feature annotation on console. By default FALSE
#' @seealso download_gtf
#' @examples
#' read_gtf(filter = 'ENSG00000198947', return_value = FALSE)
#' @importFrom magrittr %>%
#' @export
read_gtf <- function(
   gtffile,
   filter = NULL,
   return_value = FALSE
){
   # Satisfy CHECK
   . <- NULL

   # Assert validity

   # Import gtf file
   message("\t\tLoading GTF file")
   import_gtf <- rtracklayer::import(gtffile)


   if(isTRUE(return_value)){
      if (is.null(filter)){
         feature_df <- import_gtf %>%
            as.data.frame(import_gtf)
      }
      else {
         feature_df <- import_gtf %>%
            as.data.frame(import_gtf) %>%
            dplyr::filter(gene_id %in% c(filter) | gene_name %in% c(filter))
      }

      feature_df
   }
   else {
      if (is.null(filter)){
         feature_df <- import_gtf %<>%
            as.data.frame(import_gtf)
      }
      else {
         feature_df <- import_gtf %<>%
            as.data.frame(import_gtf) %>%
            dplyr::filter(gene_id %in% c(filter) | gene_name %in% c(filter))
      }

      write.table(feature_df,sprintf("~/.autonomics/annotations/%s_release%s_feature_annotation.txt", stringi::stri_replace_first_fixed(organism,' ', '_'), release), quote=FALSE, sep="\t", row.names=FALSE)
      Sys.sleep(2)
      message(sprintf("\t\t%s_release%s_feature_annotation.txt written under ~/.autonomics/annotations", stringi::stri_replace_first_fixed(organism,' ', '_'), release))

   }
}


#' Generates feature counts text file
#' @param filedir     string: path to SAM or BAM file directory (one SAM or BAM file per sample)
#' @param gtffile     string: path to (organism specific) gtffile
#' @param paired_end  TRUE or FALSE: paired end reads?
#' @param ...         passed to Rsubread::featureCounts
#' @importFrom magrittr %>%
#' @examples
#' download .zip file from "https://bitbucket.org/graumannlab/billing.stemcells/downloads/rnaseq_example_data.zip" and unzip it under ~/.autonomics/
#' get_feature_counts( filedir = "~/.autonomics/rnaseq_example_data/comparison",
#'                     organism = 'Homo sapiens',
#'                     release = 95,
#'                     paired_end = TRUE)
#' @export
get_feature_counts <- function(
   filedir,
   gtffile,
   paired_end      = FALSE,
   ...
){
   #set directory to samples with bam files
   setwd(sprintf("%s",filedir))

   #list directories with bam files
   files <- list.files(filedir, pattern = ".sam$|.bam$", full.names = TRUE, recursive = TRUE)
   filenames   <- files %>% stringi::stri_split_fixed('/') %>% vapply(function(y) y[length(y)], 'character') %>% stringi::stri_replace_first_regex('.bam|.sam', '')
   subdirnames <- files %>% stringi::stri_split_fixed('/') %>% vapply(function(y) y[length(y)-1], 'character')
   sample_names <- if (assertive.properties::has_no_duplicates(filenames)){           filenames
                  } else if (assertive.properties::has_no_duplicates(subdirnames)){  subdirnames
                  } else {                                                           paste0(subdirnames, '_', filenames) }

   message("\t\tCounting reads per feature for given samples")
   #count features for the list of samples
   gene_counts <-  Rsubread::featureCounts(files = files,
                                             annot.ext = gtffile,
                                             isGTFAnnotationFile = TRUE,
                                             GTF.attrType.extra = c("gene_name", "gene_biotype"),
                                             isPairedEnd = paired_end,
                                             ...)

   message("\t\tAdding feature names and types to count table")
   gene_counts <- cbind(gene_counts$annotation[,c("GeneID","gene_name","gene_biotype")], gene_counts$counts)  %>%
                    magrittr::set_colnames(c("gene_id", "gene_name", "gene_biotype" ,sample_names))


   #create directory for saving feature_count.txt file
   create_dir <- dir.create("~/.autonomics/counts", recursive=TRUE, showWarnings = FALSE)

   write.table(gene_counts,"~/.autonomics/counts/gene_counts.txt", quote=FALSE, sep="\t", row.names = FALSE)

   message("\t\tgene_counts.txt file written under ~/.autonomics/counts/")

}

#' Read rnaseq counts
#' @param file      string: path to rnaseq counts file
#' @param fid_var   string or number: feature id variable
#' @param fname_var string or number: feature name variable
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    file <- 'extdata/stemdiff/rnaseq/gene_counts.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_rnaseq(fid_var = 'gene_id', fname_var = 'gene_name')
#' }
#' @seealso merge_sdata, merge_fdata
#' @importFrom magrittr %>%
#' @export
read_rnaseq <- function(file, fid_var, fname_var = character(0)){

   assertive.files::assert_all_are_existing_files(file)
   dt <- data.table::fread(file, integer64='numeric')

   assertive.sets::assert_is_subset(fid_var, names(dt))
   fid_col <- which(names(dt)==fid_var)
   expr_cols   <- dt %>% vapply(is.integer, logical(1)) %>% unname() %>% which()
   fdata_cols  <- dt[, -fid_col, with = FALSE] %>% vapply(is.integer, logical(1)) %>% magrittr::not() %>% unname() %>% which() %>% magrittr::add(1) %>% c(fid_col, .)

   object <- file %>% read_omics(fid_rows   = 2:nrow(dt),   fid_cols   = fid_col,
                                 sid_rows   = 1,            sid_cols   = expr_cols,
                                 expr_rows  = 2:nrow(dt),   expr_cols  = expr_cols,
                                 fvar_rows  = 1,            fvar_cols  = fdata_cols,
                                 fdata_rows = 2:nrow(dt),   fdata_cols = fdata_cols,
                                 transpose  = FALSE,
                                 verbose    = TRUE)

   sdata(object)$subgroup <- object %>% guess_subgroup_values(verbose = TRUE)

   if (length(fname_var)>0){
      assertive.sets::assert_is_subset(fname_var, fvars(object))
      fdata(object) %<>% (function(x){x$feature_name <- x[[fname_var]];
                                                         x %>% autonomics.support::pull_columns(c('feature_id', 'feature_name'))})
   }

   object

}


#==========================================================
# EXIQON
#==========================================================

#' Read exiqon
#' @param file string: path to exiqon genex file
#' @importFrom magrittr %>%
#' @export
read_exiqon <- function(file){
   assertive.files::assert_all_are_existing_files(file)
   dt <- extract_rectangle(file, sheet=1)
   file %>% read_omics(sheet      = 1,
                       fid_rows   = 1,                       fid_cols   = 2:(ncol(dt)-2),
                       sid_rows   = 2:(nrow(dt)-3),          sid_cols   = 1,
                       expr_rows  = 2:(nrow(dt)-3),          expr_cols  = 2:(ncol(dt)-2),
                       fvar_rows  = (nrow(dt)-2):nrow(dt),   fvar_cols  = 1,
                       svar_rows  = 1,                       svar_cols  = (ncol(dt)-1):ncol(dt),
                       fdata_rows = (nrow(dt)-2):nrow(dt),   fdata_cols = 2:(ncol(dt)-2),
                       sdata_rows = 2:(nrow(dt)-3),          sdata_cols = (ncol(dt)-1):ncol(dt),
                       transpose  = TRUE,
                       verbose    = TRUE)
}


#=======================================================
# SOMASCAN
#=======================================================

#' Read somascan
#'
#' Read data from somascan adat file
#'
#' @param file         string: path to *.adat file
#' @param fid_var      string: feature_id   variable
#' @param sid_var      string: sample_id    variable
#' @param subgroup_var string: subgroup     variable
#' @param fname_var    string: feature_name variable
#' @param ...          provide backward compatibility to deprecated function load_soma
#' @return Summarizedexperiment
#' @seealso prepare_somascan
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    file <- system.file('extdata/stemcomp/soma/stemcomp.adat', package = 'autonomics.data')
#'    file %>% read_somascan()
#' }
#' @importFrom magrittr %>%
#' @export
read_somascan <- function(
   file,
   fid_var      = 'SeqId',
   sid_var      = 'SampleId',
   subgroup_var = 'SampleGroup',
   fname_var    = 'EntrezGeneSymbol'
){
   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.types::assert_is_a_string(fid_var)
   assertive.types::assert_is_a_string(sid_var)

   # Understand file structure
   content <- readLines(file)
   n_row <- length(content)
   n_col <- utils::count.fields(file, quote='', sep='\t') %>% max()

   f_row <- 1 + which(stringi::stri_detect_fixed(content, '^TABLE_BEGIN'))

   s_row <- content                                 %>%
            magrittr::extract(f_row:length(.)) %>%
            purrr::detect_index(function(y) stringi::stri_detect_regex(y, '^\t+', negate=TRUE)) %>%
            magrittr::add(f_row-1)

   f_col <- content                                    %>%
            magrittr::extract(f_row)                   %>%
            stringi::stri_extract_first_regex('^\\t+') %>%
            stringi::stri_count_fixed('\t')            %>%
            magrittr::add(1)

   fid_cols <- (1+f_col):n_col
   fid_rows <- content                             %>%
               magrittr::extract(f_row:(s_row-1))  %>%
               stringi::stri_extract_first_words() %>%
               magrittr::equals(fid_var)           %>%
               which()                             %>%
               magrittr::add(f_row-1)

   sid_rows <- (1+s_row):n_row
   sid_cols <- content                           %>%
               magrittr::extract(s_row)          %>%
               stringi::stri_extract_all_words() %>%
               unlist()                          %>%
               magrittr::equals(sid_var)         %>%
               which()

   # Read
   object <- file %>% read_omics(fid_rows   =  fid_rows,         fid_cols   =  fid_cols,
                                 sid_rows   =  sid_rows,         sid_cols   =  sid_cols,
                                 expr_rows  = (s_row+1):n_row,   expr_cols  = (f_col+1):n_col,
                                 fvar_rows  =  f_row:(s_row-1),  fvar_cols  =  f_col,
                                 fdata_rows =  f_row:(s_row-1),  fdata_cols  = (f_col+1):n_col,
                                 svar_rows  =  s_row,            svar_cols  = 1:(f_col-1),
                                 sdata_rows = (s_row+1):n_row,   sdata_cols = 1:(f_col-1),
                                 transpose  = TRUE,
                                 verbose    = TRUE)

   # sdata
   sdata(object) %<>% (function(y){ y$subgroup <- y[[subgroup_var]]
                                    y %>% autonomics.support::pull_columns(c('sample_id', 'subgroup'))})

   # fdata
   assertive.sets::assert_is_subset(fname_var, fvars(object))
   fdata(object) %<>% (function(y){ y$feature_name <- y[[fname_var]]
                                    y %>% autonomics.support::pull_columns(c('feature_id', 'feature_name'))})

   # Return
   object

}

#' @export
#' @rdname read_somascan
load_soma <- function(file, ...){
   .Deprecated('read_somascan')
   read_somascan(file = file)
}


#======================================================
# METABOLON
#======================================================

#' Read metabolon
#' @param file          string: path to metabolon xlsx file
#' @param sheet         number/string: xls sheet number or name
#' @param fid_var       string: feature_id variable (ideally transcends dataset)
#' @param sid_var       string: sample_id variable
#' @param subgroup_var  string: subgroup variable (human comprehensible)
#' @param fname_var     string: feature_name variable
#' @param ...           enable backward compatibility to deprecated load_metabolon
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    file <- 'extdata/glutaminase/glutaminase.xlsx' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_metabolon()
#' }
#' @importFrom magrittr %>%
#' @export
read_metabolon <- function(
   file,
   sheet        = 2,
   fid_var      = 'COMP_ID',
   sid_var      = 'CLIENT_IDENTIFIER',
   subgroup_var = 'Group',
   fname_var    = 'BIOCHEMICAL'
){

   assertive.files::assert_all_are_existing_files(file)

   d_f <- readxl::read_excel(file, sheet, col_names = FALSE, .name_repair = 'minimal')

   fvar_rows <- which(!is.na(d_f %>% extract_dt_col(1))) %>% magrittr::extract(1)
   svar_cols <- which(!is.na(d_f %>% extract_dt_row(1))) %>% magrittr::extract(1)
   fvar_cols <- fdata_cols <- 1:svar_cols
   svar_rows <- sdata_rows <- 1:fvar_rows

   fid_rows  <- fdata_rows <- expr_rows <- (fvar_rows+1):nrow(d_f)
   sid_cols  <- sdata_cols <- expr_cols <- (svar_cols+1):ncol(d_f)
   fid_cols  <-  d_f %>% extract_dt_row(fvar_rows) %>% magrittr::extract(1:svar_cols) %>% magrittr::equals(fid_var) %>% which()
   sid_rows  <-  d_f %>% extract_dt_col(svar_cols) %>% magrittr::extract(1:fvar_rows) %>% magrittr::equals(sid_var) %>% which()

   object <- file %>% read_omics(sheet      = sheet,
                                 fid_rows   = fid_rows,      fid_cols   = fid_cols,
                                 sid_rows   = sid_rows,      sid_cols   = sid_cols,
                                 expr_rows  = expr_rows,     expr_cols  = expr_cols,
                                 fvar_rows  = fvar_rows,     fvar_cols  = fvar_cols,
                                 svar_rows  = svar_rows,     svar_cols  = svar_cols,
                                 fdata_rows = fdata_rows,    fdata_cols = fdata_cols,
                                 sdata_rows = svar_rows,     sdata_cols = sdata_cols,
                                 transpose  = FALSE,
                                 verbose    = TRUE)

   # sdata
   subgroup_var <- svars(object) %>% magrittr::extract(stringi::stri_detect_fixed(., subgroup_var))
   sdata(object) %<>% (function(y){ y$subgroup <- y[[subgroup_var]]
                                    y %>% autonomics.support::pull_columns(c('sample_id', 'subgroup'))})
   # fdata
   assertive.sets::assert_is_subset(fname_var, fvars(object))
   fdata(object) %<>% (function(y){ y$feature_name <- y[[fname_var]]
                                    y %>% autonomics.support::pull_columns(c('feature_id', 'feature_name'))})

   # return
   object
}

#' @rdname read_metabolon
#' @export
load_metabolon <- function(file, sheet=2, ...){
   .Deprecated('read_metabolon')
   read_metabolon(file = file, sheet=sheet)
}

#====================================
# LCMS PROTEINGROUPS & PHOSPHOSITES
#====================================

#' maxquant patterns
#' @export
maxquant_patterns <- c(`Ratio normalized`             =  '^Ratio ([HM]/[ML]) normalized (.+)$',
                       `Ratio`                        =  '^Ratio ([HM]/[ML]) (?!count|type|variability|iso-count|normalized)(.+)',
                       `LFQ intensity`                =  '^LFQ intensity ([HML])? ?(.+)$',
                       `Reporter intensity corrected` =  '^Reporter intensity corrected ([0-9]+) (.+)$',
                       `Reporter intensity`           =  '^Reporter intensity ([0-9]+) (.+)$',
                       `Intensity labeled`            =  '^Intensity ([HML]) (.+)$',
                       `Intensity`                    =  '^Intensity (.+)$')


#' Guess maxquant quantity from snames
#'
#' charactervector, dataframe, or SummarizedExperiment.
#'
#' @param x charactervector, dataframe, or SummarizedExperiment
#' @param ... used for proper S3 method dispatch
#' @return  character(1): value from names(maxquant_patterns)
#' @examples
#' if (require(autonomics.data)){
#'    require(magrittr)
#'
#'    # file
#'       x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'       guess_maxquant_quantity(x)
#'
#'    # charactervector
#'        guess_maxquant_quantity("Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1")
#'        guess_maxquant_quantity("Ratio M/L STD(L)_EM00(M)_EM01(H)_R1")
#'        guess_maxquant_quantity("LFQ intensity EM00.R1")
#'        guess_maxquant_quantity("Reporter intensity corrected 0 STD(0)EM00(1)EM01(2)_R1")
#'        guess_maxquant_quantity("Reporter intensity 0 STD(0)EM00(1)EM01(2)_R1")
#'        guess_maxquant_quantity("Intensity H STD(L)_EM00(M)_EM01(H)_R1")
#'
#'    # dataframe
#'        x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'              system.file(package = 'autonomics.data') %>%
#'              data.table::fread()
#'        guess_maxquant_quantity(x)
#'
#'    # SummarizedExperiment
#'       x <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data') %>%
#'             read_proteingroups(standardize_snames = FALSE,
#'                                demultiplex_snames = FALSE)
#'       guess_maxquant_quantity(x)
#' }
#' @export
guess_maxquant_quantity <- function(x, ...){
   UseMethod("guess_maxquant_quantity", x)
}

#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.character <- function(x, ...){      # x = character vector of maxquant colnames

   # read if x is filename
   if (assertive.files::is_existing_file(x)){
      x %<>% data.table::fread(header = TRUE, nrows = 1) %>% names()
   }

   # guess from character vector
   for (quantity in names(maxquant_patterns)){
      pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}


#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.data.frame <- function(x, ...){     # x = maxquant dataframe
   x <- names(x)
   for (quantity in names(maxquant_patterns)){
      pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}

#' @rdname guess_maxquant_quantity
#' @export
guess_maxquant_quantity.SummarizedExperiment <- function(x, ...){     # x = SummarizedExperiment
   x <- snames(x)
   for (quantity in names(maxquant_patterns)){
      pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
      if (any(stringi::stri_detect_regex(x, pattern)))   return(quantity)
   }
   stop('quantity could not be infered')
}




#' proteingroups fvars
#' @export
proteingroups_fvars <- c(c('id', 'Majority protein IDs', 'Protein names', 'Gene names', 'Contaminant', 'Potential contaminant', 'Reverse', 'Phospho (STY) site IDs'))

#' Standardize maxquant snames
#'
#' For charactervector or SummarizedExperiment
#'
#' Drop "Ratio normalized", "LFQ intensity" etc from maxquant snames & sample_id values
#'
#' @param x        character(.) or SummarizedExperiment
#' @param quantity maxquant quantity
#' @param verbose  logical(1)
#' @param ...      allow for proper S3 method dispatch
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     standardize_maxquant_snames("Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1")
#'     standardize_maxquant_snames("Ratio M/L STD(L)_EM00(M)_EM01(H)_R1")
#'     standardize_maxquant_snames('LFQ intensity STD_R1')
#'     standardize_maxquant_snames('LFQ intensity L STD(L)_EM00(M)_EM01(H)_R1')
#'     standardize_maxquant_snames('Reporter intensity 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'     standardize_maxquant_snames('Reporter intensity corrected 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'
#' # SummarizedExperiment
#' if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            read_proteingroups(standardize_snames = FALSE,
#'                               demultiplex_snames = FALSE)
#'      x %>% standardize_maxquant_snames(verbose = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
standardize_maxquant_snames <- function (x, ...) {
   UseMethod("standardize_maxquant_snames", x)
}


#' @importFrom magrittr %>%
#' @export
#' @rdname standardize_maxquant_snames
standardize_maxquant_snames.character <- function(
   x,
   quantity = guess_maxquant_quantity(x),
   verbose  = FALSE,
   ...
){
   # x = mix + channel. Return mix if single channel.
   pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$2')
   if (all(channel=='')){ cleanx <- mix
   } else               { cleanx <- sprintf('%s{%s}', mix, channel)
   }
   autonomics.support::cmessage('\t\tStandardize snames: %s  ->  %s', x[1], cleanx[1])
   return(cleanx)
}

#' @importFrom magrittr %>%
#' @export
#' @rdname standardize_maxquant_snames
standardize_maxquant_snames.SummarizedExperiment <- function(
   x,
   quantity = guess_maxquant_quantity(x),
   verbose  = FALSE,
   ...
){
   newsnames <- snames(x) %>% standardize_maxquant_snames(quantity = quantity, verbose=verbose)
   snames(x) <- sdata(x)$sample_id <- newsnames
   x
}


#' Demultiplex snames
#'
#' For charactervector or SummarizedExperiment
#'
#' @param x        character vector or SummarizedExperiment
#' @param verbose  logical
#' @param ...      allow for S3 dispatch
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     # Alternate multiplexing forms supported
#'      demultiplex_snames("STD(L)_EM00(M)_EM01(H)_R1{M/L}")     # Label Ratio
#'      demultiplex_snames('A(0)_B(1)_C(2)_D(3)_R1{0}'     )     # Reporter intensity
#'      demultiplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # Label Intensity
#'
#'    # Alternate separators supported
#'      demultiplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # underscore
#'      demultiplex_snames('STD(L).EM00(M).EM01(H).R1{L}')       # dot
#'      demultiplex_snames('STD(L)EM00(M)EM01(H).R1{L}')         # no separator
#'
#'    # Composite snames supported
#'      demultiplex_snames("WT.t0(L)_WT.t1(M)_WT.t2(H)_R1{H/L}") # composite snames
#'
#'    # Uniqueness ensured by appending labels when necessary
#'      demultiplex_snames(c("STD(L).BM00(M).BM00(H).R10{M/L}",  # implicit uniquification
#'                                          "STD(L).BM00(M).BM00(H).R10{H/L}"))
#'    # Uniplexed snames are returned unchanged
#'      demultiplex_snames(c('STD_R1', 'EM0_R1'))
#'
#' # SummarizedExperiment
#'   if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            read_proteingroups(standardize_snames = FALSE, demultiplex_snames = FALSE) %>%
#'            standardize_maxquant_snames()
#'      x %>% demultiplex_snames(verbose = TRUE)
#'   }
#'
#' @export
demultiplex_snames <- function (x, ...) {
   UseMethod("demultiplex_snames", x)
}

#' @rdname demultiplex_snames
#' @export
demultiplex_snames.character <- function(x, verbose = FALSE, ...){

   # Return unchanged if not multiplexed
   # KD(H)WT(L){H/L}
   pattern <- '(.+)\\{(.+)\\}'
   n_open   <- x %>% stringi::stri_count_fixed('(')
   n_closed <- x %>% stringi::stri_count_fixed(')')
   is_multiplexed <- all(stringi::stri_detect_regex(x, pattern) & (n_open==n_closed) & (n_open>0))
   if (!is_multiplexed) return(x)

   # Separate mix and channel
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$2')

   # Separate labels and samples
   pattern <- '\\(.+?\\)'
   labels  <- mix %>% stringi::stri_extract_all_regex(pattern) %>% lapply(stringi::stri_replace_first_fixed, '(', '') %>%
      lapply(stringi::stri_replace_first_fixed, ')', '')
   samples <- mix %>% stringi::stri_split_regex(pattern) %>%
      # rm sep from samples (but not from replicate - needed to glue back later!)
      lapply(function(y){y[1:length(labels[[1]])] %<>% stringi::stri_replace_first_regex('^[_. ]', ''); y})

   # Return unchanged if mixes differ in no of labels or samples
   are_all_identical <- function(y) if (length(y)==1) TRUE else all(y[-1] == y[1])
   n_samples <- vapply(samples,  length, integer(1))
   n_labels  <- vapply(labels, length, integer(1))
   if (!are_all_identical(n_samples) | !are_all_identical(n_labels)){
      autonomics.support::cmessage('\t\tCannot demultiplexing snames: mixes differ in number of samples or labels')
      return(x)
   }

   # Extract replicate
   n_samples %<>% unique()
   n_labels  %<>% unique()
   if (n_samples > n_labels){ replicate <- mix %>% stringi::stri_split_regex(pattern) %>% vapply((function(y) y %>% magrittr::extract(length(y))), character(1))
                              samples %<>% lapply(magrittr::extract, 1:(n_samples-1))
   } else {                   replicate <- rep('', length(samples))
   }

   # Extract channel samples from mix
   is_ratio <- channel %>% stringi::stri_detect_fixed('/') %>% all()
   samples %<>% mapply(magrittr::set_names, ., labels, SIMPLIFY = FALSE)
   if (is_ratio){
      num_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(magrittr::extract, character(1), 1)
      den_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(magrittr::extract, character(1), 2)
      den_samples <- mapply(magrittr::extract, samples, den_label)
      num_samples <- mapply(magrittr::extract, samples, num_label)
      xdemultiplex <- sprintf('%s_%s%s', num_samples, den_samples, replicate)
   } else {
      samples %<>% mapply(magrittr::extract, ., channel)
      xdemultiplex <- sprintf('%s%s', samples, replicate)
   }
   if (verbose) autonomics.support::cmessage('\t\tDemultiplex snames: %s  ->  %s', x[1], xdemultiplex[1])

   # Ensure uniqueness. Add labels if required.
   idx <- autonomics.support::cduplicated(xdemultiplex) %>% which()
   if (length(idx)>0){
      label_tags <- channel[idx] %>% stringi::stri_replace_first_fixed('/', '')
      if (verbose)   autonomics.support::cmessage('\t\tUniquify snames: %s -> %s%s (for %d/%d snames)',
                                                  xdemultiplex[idx][1], xdemultiplex[idx][1], label_tags[1],
                                                  length(idx), length(xdemultiplex))
      xdemultiplex[idx] %<>% paste0(label_tags)
   }

   # Return
   return(xdemultiplex)
}

#' @rdname demultiplex_snames
#' @importFrom magrittr %>%
#' @export
demultiplex_snames.SummarizedExperiment <- function(x,verbose  = FALSE, ...){
   newsnames <- snames(x) %>% demultiplex_snames(verbose = verbose)
   snames(x) <- sdata(x)$sample_id <- newsnames
   x
}

#' Read proteingroups
#' @param file                character(1). Path to 'proteinGroups.txt'
#' @param quantity            character(1). Expression columns to extract from proteinGroups file: 'Ratio normalized', 'Ratio', 'LFQ intensity',
#'                                          'LFQ intensity', 'Reporter intensity corrected', 'Reporter intensity', 'Intensity labeled', or 'Intensity'.
#' @param fvars               character(n). Annotation columns to extract from proteinGroups file.
#' @param standardize_snames  logical(1).   Standardize maxquant snames: Ratio normalized H/L WT(L).KD(H).R1 -> WT(L).KD(H).R1{H/L} ?
#' @param demultiplex_snames  logical(1).   Demultiplex (standardized) maxquant snames: WT(L).KD(H).R1{H/L} -> KD(H)_WT(L).R1 ?
#' @param verbose             logical(1).   Message progress?
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_proteingroups()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% read_proteingroups()
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_proteingroups <- function(
   file,
   quantity           = guess_maxquant_quantity(file),
   fvars              = proteingroups_fvars,
   standardize_snames = TRUE,
   demultiplex_snames = TRUE,
   verbose            = TRUE
){

   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.sets::assert_is_subset(quantity, names(maxquant_patterns))
   assertive.types::assert_is_a_bool(verbose)

   # Initial Read
   assertive.files::assert_all_are_existing_files(file)
   dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
   fvars %<>% intersect(names(dt))

   # Define components
   fid_rows   <- 2:nrow(dt)
   fid_cols   <- which(names(dt) == 'id')
   sid_rows   <- 1
   sid_cols   <- names(dt) %>% stringi::stri_detect_regex(maxquant_patterns[[quantity]]) %>% which()
   expr_rows  <- 2:nrow(dt)
   expr_cols  <- sid_cols
   fvar_rows  <- 1
   fvar_cols  <- match(fvars, names(dt))
   fdata_rows <- 2:nrow(dt)
   fdata_cols <- fvar_cols

   # Read sumexp
   object <- file %>% read_omics(fid_rows   = fid_rows,     fid_cols   = fid_cols,
                                 sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                 expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                 fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                 fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                 transpose  = FALSE,
                                 verbose    = verbose)
   # Clean sdata
   if (standardize_snames) object %<>% standardize_maxquant_snames(verbose = verbose)
   if (demultiplex_snames) object %<>% demultiplex_snames(verbose = verbose)
   object$subgroup <- object$sample_id %>% guess_subgroup_values(verbose = verbose)
   #object$block    <- object$sample_id %>% guess_subject_values( verbose = TRUE)

   # Clean fdata
   contaminant_var <- c('Contaminant', 'Potential contaminant') %>% intersect(fvars(object))
   fdata(object)[[contaminant_var]] %<>% (function(x){x[is.na(x)] <- ''; x})
   fdata(object)[['Reverse'      ]] %<>% (function(x){x[is.na(x)] <- ''; x})
   fdata(object)$feature_name    <- fdata(object)$`Gene names`
   fdata(object)$feature_uniprot <- fdata(object)$`Majority protein IDs`
   fdata(object) %<>% autonomics.support::pull_columns(c('feature_id', 'feature_name', 'feature_uniprot'))

   # Return
   object

}

#' phosphosites fvars
#' @export
phosphosite_fvars <- c('id', 'Protein group IDs', 'Positions within proteins', 'Localization prob')


#' Read phosphosites
#' @param file                string: phosphosites filepath
#' @param proteingroups_file  string: proteingroups filepath
#' @param quantity            NULL or value in names(maxquant_patterns)
#' @param fvars               string vector
#' @param standardize_snames  logical
#' @param demultiplex_snames  logical
#' @param verbose             logical
#' @examples
#' \dontrun{
#' if (require(autonomics.data)){
#'    require(magrittr)
#'    file <- 'extdata/stemdiff/maxquant/phospho (STY)Sites.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% read_phosphosites()
#' }
#' }
#' @importFrom magrittr %>%
#' @export
read_phosphosites <- function(
   file,
   proteingroups_file = dirname(file) %>% paste0('/phospho (STY)Sites.txt'),
   quantity           = guess_maxquant_quantity(file),
   fvars              = phosphosite_fvars,
   standardize_snames = TRUE,
   demultiplex_snames = TRUE,
   verbose            = FALSE
){

   # Check
   `Protein group IDs` <- `Localization prob` <- NULL

   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.sets::assert_is_subset(quantity, names(maxquant_patterns))
   assertive.types::assert_is_character(fvars)

   # Initial Read
   dt <- data.table::fread(file, integer64 = 'numeric', header = TRUE)
   pattern <- maxquant_patterns %>% magrittr::extract2(quantity)
   value_cols <- names(dt) %>% (function(x) stringi::stri_detect_regex(x,pattern) & !stringi::stri_detect_regex(x, '___[1-3]')) %>% which()
   fvar_cols  <- which(names(dt) %in% fvars)
   #dt %<>% magrittr::extract(, c(fvars, value_cols), with = FALSE)

   # Read phosphosites
   fid_rows   <- 2:nrow(dt)
   fid_cols   <- which(names(dt) == 'id')
   sid_rows   <- 1
   sid_cols   <- value_cols
   expr_rows  <- 2:nrow(dt)
   expr_cols  <- value_cols
   fvar_rows  <- 1
   fvar_cols  <- fvar_cols
   fdata_rows <- 2:nrow(dt)
   fdata_cols <- fvar_cols
   phosphosites  <- file  %>% read_omics(fid_rows   = fid_rows,  fid_cols   = fid_cols,
                                         sid_rows   = sid_rows,     sid_cols   = sid_cols,
                                         expr_rows  = expr_rows,    expr_cols  = expr_cols,
                                         fvar_rows  = fvar_rows,    fvar_cols  = fvar_cols,
                                         fdata_rows = fdata_rows,   fdata_cols = fdata_cols,
                                         transpose  = FALSE,
                                         verbose    = verbose)

   # Calculate occupancies (allows to disentangle phosphorylation and protein expression)
   autonomics.support::cmessage('\t\toccupancies(phosphosites) = exprs(phosphosites) - exprs(proteingroups)')
   proteingroups <- read_proteingroups(proteingroups_file, quantity = quantity, verbose = FALSE) %>%
                    magrittr::extract(phosphosites %>% fvalues("Protein group IDs"), ) %>%
                    magrittr::extract(, phosphosites %>% snames())
   occupancies(phosphosites) <- exprs(phosphosites) - exprs(proteingroups)

   # Simplify snames
   if (standardize_snames) phosphosites %<>% standardize_maxquant_snames(verbose = verbose)
   if (demultiplex_snames) phosphosites %<>% demultiplex_snames(verbose = verbose)

   # Return
   phosphosites

}

