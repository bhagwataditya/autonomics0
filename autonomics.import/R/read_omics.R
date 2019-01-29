#' Read omics data as is from file
#' @param file character(1)
#' @param sheet numeric(1) or character(1)
#' @return data.table
#' @importFrom magrittr %>%
#' @export
read_file <- function(file, sheet){
   assertive.files::assert_all_are_existing_files(file)
   if (tools::file_ext(file) %in% c('xls', 'xlsx')){
      file %>% readxl::read_excel(sheet = sheet, col_names = FALSE) %>% data.table::data.table()
   } else {
      file %>% data.table::fread(na.strings = "", header = FALSE, integer64 = 'numeric')
   }
}

# Extract row
extract_dt_row <- function(dt, i) dt %>% magrittr::extract(i,) %>% as.matrix() %>% magrittr::extract(1, ) %>% unname()
extract_dt_col <- function(dt, i) dt[[i]]


#' Extract rectangle from omics datatable
#' @param dt         datatable with raw omics data
#' @param rows       row selector: logical(n), numeric(n), or character(n)
#' @param cols       col selector: logical(n), numeric(n), or character(n)
#' @param transpose  logical(1)
#' @param drop       logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    dt <- file %>% read_file(sheet = 2)
#'    dt %>% extract_rectangle(rows = 11:401, cols = 15:86)                    %>% extract(1:3, 1:3)           # exprs
#'    dt %>% extract_rectangle(rows = 11:401, cols = 5    )                    %>% extract(1:3, , drop=FALSE)  # fids
#'    dt %>% extract_rectangle(rows = 2,      cols = 15:86)                    %>% extract(, 1:3, drop=FALSE)  # sids
#'    dt %>% extract_rectangle(rows = 10:401, cols = 1:14 )                    %>% extract(1:3, 1:3)           # fdata
#'    dt %>% extract_rectangle(rows = 1:10,   cols = 14:86, transpose = TRUE)  %>% extract(1:3, 1:3)           # sdata
#' }
#' @importFrom magrittr %>% %<>%
#' @export
extract_rectangle <- function(dt, rows, cols, transpose = FALSE, drop = FALSE){
   dt %<>% magrittr::extract(rows, cols, with = FALSE)
   rectangle <- dt %>% as.matrix()
   if (transpose) rectangle %<>% t()
   if (drop) if (nrow(rectangle)==1 | ncol(rectangle)==1) rectangle %<>% as.vector('character')
   rectangle
}

#' Read omics data
#' @param file       character(1): name of text (txt, csv, tsv) or excel (xls, xlsx) file
#' @param sheet      integer(1) or character(1)
#' @param fid_rows   featureid rows: logical(n) or numeric(n)
#' @param fid_cols   featureid cols: logical(n) or numeric(n)
#' @param sid_rows   sampleid  rows: logical(n) or numeric(n)
#' @param sid_cols   sampleid  cols: logical(n) or numeric(n)
#' @param expr_rows  expr      rows: logical(n) or numeric(n)
#' @param expr_cols  expr      rows: logical(n) or numeric(n)
#' @param fvar_rows  fvar      rows: logical(n) or numeric(n)
#' @param fvar_cols  fvar      cols: logical(n) or numeric(n)
#' @param svar_rows  svar      rows: logical(n) or numeric(n)
#' @param svar_cols  svar      cols: logical(n) or numeric(n)
#' @param fdata_rows fdata     rows: logical(n) or numeric(n)
#' @param fdata_cols fdata     cols: logical(n) or numeric(n)
#' @param sdata_rows sdata     rows: logical(n) or numeric(n)
#' @param sdata_cols sdata     cols: logical(n) or numeric(n)
#' @param transpose  logical(1)
#' @param verbose    logical(1)
#' @examples
#' require(magrittr)
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'   file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'   file %>% autonomics.import::read_omics_asis(sheet      = 2,
#'                                          fid_rows   = 11:401,    fid_cols   = 5,
#'                                          sid_rows   = 3,         sid_cols   = 15:86,
#'                                          expr_rows  = 11:401,    expr_cols  = 15:86,
#'                                          fvar_rows  = 10,        fvar_cols  = 1:14,
#'                                          svar_rows  = 1:10,      svar_cols  = 14,
#'                                          fdata_rows = 11:401,    fdata_cols = 1:14,
#'                                          sdata_rows = 1:10,      sdata_cols = 15:86,
#'                                          transpose  = FALSE)
#' }
#'
#' # SOMASCAN
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat', package = 'atkin.2014')
#'    file %>% autonomics.import::read_omics_asis(fid_rows   = 73,       fid_cols   = 27:1343,
#'                                           sid_rows   = 93:520,   sid_cols   = 6:6,
#'                                           expr_rows  = 93:520,   expr_cols  = 27:1343,
#'                                           fvar_rows  = 73:91,    fvar_cols  = 26,
#'                                           svar_rows  = 92,       svar_cols  = 1:25,
#'                                           fdata_rows = 73:91,    fdata_cols = 27:1343,
#'                                           sdata_rows = 93:520,   sdata_cols = 1:25,
#'                                           transpose  = TRUE)
#' }
#'
#' # EXIQON
#' if require(subramanian.2016){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx', package = 'subramanian.2016')
#'    file %>% autonomics.import::read_omics_asis(sheet      = 1,
#'                                           fid_rows   = 1,        fid_cols   = 2:330,
#'                                           sid_rows   = 2:73,     sid_cols   = 1,
#'                                           expr_rows  = 2:73,     expr_cols  = 2:330,
#'                                           fvar_rows  = 74:76,    fvar_cols  = 1,
#'                                           svar_rows  = 1,        svar_cols  = 331:332,
#'                                           fdata_rows = 74:76,    fdata_cols = 2:330,
#'                                           sdata_rows = 2:73,     sdata_cols = 331:332,
#'                                           transpose  = TRUE)
#' }
#'
#' # PROTEINGROUPS LABELFREE
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::read_omics_asis(fid_rows   = 2:5893,   fid_cols   = 275,
#'                                           sid_rows   = 1,        sid_cols   = 227:248,
#'                                           expr_rows  = 2:5893,   expr_cols  = 227:248,
#'                                           fvar_rows  = 1,        fvar_cols  = 1:8,
#'                                           fdata_rows = 2:5893,   fdata_cols = 1:8,
#'                                           transpose  = FALSE)
#' }
#'
#' # PROTEINGROUPS LABELED
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/stemdiff/maxquant/proteinGroups.txt', package = 'autonomics.data')
#'    file %>% autonomics.import::read_omics_asis(fid_rows   = 2:9783,  fid_cols   = 383,
#'                                           sid_rows   = 1,       sid_cols   = seq(124, 316, by = 6),
#'                                           expr_rows  = 2:9783,  expr_cols  = seq(124, 316, by = 6),
#'                                           fvar_rows  = 1,       fvar_cols  = c(2, 6, 7, 383),
#'                                           fdata_rows = 2:9783,  fdata_cols = c(2, 6, 7, 383),
#'                                           transpose  = FALSE)
#' }
#'
#' # RNASEQ
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    file %>% autonomics.import::read_omics_asis(fid_rows   = 2:48821,   fid_cols   = 1,
#'                                           sid_rows   = 1,         sid_cols   = 5:86,
#'                                           expr_rows  = 2:48821,   expr_cols  = 5:86,
#'                                           fvar_rows  = 1,         fvar_cols  = 1:4,
#'                                           fdata_rows = 2:48821,   fdata_cols = 1:4,
#'                                           transpose  = FALSE)
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_omics_asis <- function(
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
   verbose    = FALSE
){
   # Assert
   assertive.files::assert_all_are_existing_files(file)
   assertive.types::assert_is_a_bool(transpose)

   # read
   dt <- autonomics.import::read_file(file, sheet)

   # Extract exprs
   fids1  <- dt %>% autonomics.import::extract_rectangle(fid_rows,  fid_cols,  transpose = transpose, drop = TRUE)
   sids1  <- dt %>% autonomics.import::extract_rectangle(sid_rows,  sid_cols,  transpose = transpose, drop = TRUE)
   exprs1 <- dt %>% autonomics.import::extract_rectangle(expr_rows, expr_cols, transpose = transpose) %>% (function(y){class(y) <- 'numeric'; y})

   # Extract feature annotations
   #    Leave rownames(fdata1) empty: fids1 may contain non-valid values
   #    This happens in MaxQuant files, which sometimes contain missing rows (I think after opening in excell)
   fdata1 <- data.frame(feature_id = fids1, stringsAsFactors = FALSE)
   fdata_available <- !is.null(fvar_rows) & !is.null(fvar_cols)
   if (fdata_available){
      fvars1 <- dt %>% autonomics.import::extract_rectangle(fvar_rows,  fvar_cols,  transpose = transpose, drop = TRUE)
      fdata1 %<>% cbind(autonomics.import::extract_rectangle(dt, fdata_rows, fdata_cols, transpose = transpose) %>%
                        magrittr::set_colnames(fvars1) %>%
                        data.frame(stringsAsFactors = FALSE, check.names = FALSE))
   }

   # Extract sample annotations
   #    Leave rownames(sdata1) empty: sids may contain non-valid values
   #    This happens in SOMA files, where CLIENT_IDENTIFIER is not unique for calibrator and buffer samples
   sdata1 <- data.frame(sample_id = sids1, stringsAsFactors = FALSE)
   sdata_available <- !is.null(svar_rows) & !is.null(svar_cols)
   if (sdata_available){
      svars1 <- dt %>% autonomics.import::extract_rectangle(svar_rows,  svar_cols,  transpose =  transpose, drop = TRUE)
      sdata1 %<>% cbind(autonomics.import::extract_rectangle(dt, sdata_rows, sdata_cols, transpose = !transpose) %>%
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
   autonomics.import::fdata(object) <- fdata1
   autonomics.import::sdata(object) <- sdata1

   # Return
   object
}


#' Add subgroup values to object
#'
#' Three methods: (1) take values from existing svar
#'                (2) guess from sampleids
#'                (3) merge designfile (with one column sample_id) into sdata
#'
#' @param object SummarizedExperiment
#' @param take_from_svar        character(1)
#' @param guess_from_sampleids  logical(1)
#' @param merge_designfile      character(1): designfile path
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    object <- file %>% autonomics.import::read_metabolon_asis()
#' }
#' @importFrom magrittr %>%
#' @export
add_subgroup <- function(
   object,
   take_from_svar        = NULL,
   guess_from_sampleids  = FALSE,
   merge_designfile      = NULL
){

   # Assert
   if (!is.null(take_from_svar)){ assertive.types::assert_is_character(take_from_svar);
                                  assertive.types::assert_is_subset(take_from_svar, autonomics.import::svars(object)) }
   assertive.types::assert_is_logical(guess_from_sampleids)
   if (!is.null(merge_designfile)) assertive.files::assert_all_are_existing_files(merge_design_file)

   # Either take from svar
   if (!is.null(take_from_svar)){
      autonomics.support::cmessage("Add subgroup: take from svar '%s'", take_from_svar)
      autonomics.import::sdata(object) %<>% cbind(subgroup = .[[svar]], .)

   # Or guess from sampleids
   } else if (guess_from_sampleids){
      autonomics.support::cmessage('Add subgroup: guess from sampleids')
      autonomics.import::sdata(object) %<>% cbind(subgroup = autonomics.import::guess_subgroup_values(.$sample_id), .)

   # Or merge in designfile
   } else if (merge_designfile){
      autonomics.support::cmessage("Add subgroup: merge in designfile '%s'", merge_designfile)
      design_df <- autonomics.support::cfread(merge_designfile, data.table = FALSE)
      assertive.sets::assert_is_subset('sample_id', names(design_df))
      autonomics.import::sdata(object) %<>% merge(design_df, by = 'sample_id', sort = FALSE, all.x =TRUE)
   }

   # Return
   return(object)
}

#===============================================================================================================


#' @rdname read_somascan
#' @importFrom magrittr %>%
#' @export
read_somascan_asis <- function(
   file,
   fid_var = 'SeqId',
   sid_var = 'SampleId'
){

   # Assert
   assertive.files::assert_all_are_existing_files(file)

   # Peak
   dt <- data.table::fread(file, header = FALSE)
   fdata_row1 <- 1 + which(dt[[1]] == '^TABLE_BEGIN')[1]
   fvar_col <- which(dt %>% extract_dt_row(fdata_row1) != '')[1]
   svar_row <- which(dt[[1]] != '' &  (1:nrow(dt) > fdata_row1))[1]
   fvars <- dt %>% extract_dt_col(fvar_col) %>% magrittr::extract(fdata_row1:svar_row)
   svars <- dt %>% extract_dt_row(svar_row) %>% magrittr::extract(1:fvar_col-1)
   fid_row <- fdata_row1 + -1 + which(fvars == fid_var)
   sid_col <- which(svars == sid_var)

   # Read
   file %>% read_omics_asis(sheet      =  NULL,
                       fid_rows   =  fid_row,                   fid_cols   = (fvar_col+1):ncol(dt),
                       sid_rows   = (svar_row+1):nrow(dt),      sid_cols   =  sid_col,
                       expr_rows  = (svar_row+1):nrow(dt),      expr_cols  = (fvar_col+1):ncol(dt),
                       fvar_rows  =  fdata_row1:(svar_row-1),   fvar_cols  =  fvar_col,
                       fdata_rows =  fdata_row1:(svar_row-1),   fdata_cols = (fvar_col+1):ncol(dt),
                       svar_rows  =  svar_row,                  svar_cols  = 1:(fvar_col-1),
                       sdata_rows = (svar_row+1):nrow(dt),      sdata_cols = 1:(fvar_col-1),
                       transpose  = TRUE,
                       verbose    = TRUE)
}


#' Read somascan data
#' @param file character(1)
#' @param fid_var  character(1): feature id variable
#' @param sid_var  character(1): sample id variable
#' @param filter_sample_type         character(.): sample  types to be filtered for.     Probably a subset of c('Sample', 'QC', 'Buffer', 'Calibrator').
#' @param filter_feature_type        character(.): feature types to be filtered for.     Probably a subset of c('Protein', 'Hybridization Control Elution', 'Rat Protein').
#' @param filter_sample_quality      character(.): sample  qualities to be filtered for. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param filter_feature_quality     character(.): feature qualities to be filtered for. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param filter_na_svars            logical(1)  : whether to rm NA svars
#' @param filter_single_value_svars  logical(1)  : whether to rm single value svars
#' @param infer_design               logical(1)  : whether to infer design from sampleids
#' @param log2_transform             logical(1)  : whether to log2 transform
#' @return Summarizedexperiment
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::read_somascan_asis()
#'    file %>% autonomics.import::read_somascan()
#' }
#' @importFrom magrittr %>%
#' @export
read_somascan <- function(
   file,
   fid_var                 = 'SeqId',
   sid_var                 = 'SampleId',
   filter_sample_type      = 'Sample',
   filter_feature_type     = 'Protein',
   filter_sample_quality   = c('FLAG', 'PASS'),
   filter_feature_quality  = c('FLAG', 'PASS'),
   rm_na_svars             = TRUE,
   rm_single_value_svars   = TRUE,
   infer_design            = FALSE,
   log2_transform          = TRUE
){

   # Read
   #-----
   autonomics.support::cmessage('\t\tRead %s', file)
   object <- file %>% autonomics.import::read_somascan_asis(fid_var = fid_var, sid_var = sid_var)

   # Filter
   #-------
   if ('SampleType' %in% autonomics.import::svars(object)){ # sample type - older versions don't have it
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Sample types` = autonomics.import::sdata(object)$SampleType))
      object %<>% autonomics.import::filter_samples(SampleType %in% filter_sample_type, verbose = TRUE)
   }
   if ('RowCheck'   %in% autonomics.import::svars(object)){ # sample quality
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Sample qualities ("RowCheck")` = autonomics.import::sdata(object)$RowCheck))
      object %<>% autonomics.import::filter_samples(RowCheck %in% filter_sample_quality, verbose = TRUE)
   }
   if ('Type'       %in% autonomics.import::fvars(object)){ # feature type
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Type` = autonomics.import::fdata(object)$Type))
      object %<>% autonomics.import::filter_features(Type %in% filter_feature_type, verbose = TRUE)
   }
   if ('ColCheck'   %in% autonomics.import::fvars(object)){ # feature quality
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Feature qualities ("ColCheck")` = autonomics.import::fdata(object)$ColCheck))
      object %<>% autonomics.import::filter_features(ColCheck %in% filter_feature_quality, verbose = TRUE)
      message('\t\t=========================================================================')
   }

   # Design
   #-------
   if (autonomics.support::all_svalues_available(object, 'SampleGroup')){
      autonomics.import::sdata(object) %<>% cbind(subgroup = .$SampleGroup, .)
   }

   # Select
   #------
   if (rm_na_svars)            autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns()
   if (rm_single_value_svars)  autonomics.import::sdata(object) %<>% autonomics.support::rm_single_value_columns()

   # Log2 transform
   #---------------
   if (log2_transform){
      autonomics.support::cmessage('\t\tLog2 transform')
      autonomics.import::exprs(object) %<>% log2()
   }


   # Return
   #-------
   object
}



#===========================================================================================================================


#' Read meabolon data
#' @param file     character(1)
#' @param sheet    numeric(1) or character(1)
#' @param fid_var  character(1): feature id variable
#' @param sid_var  character(1): sample id variable
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    file %>% autonomics.import::read_metabolon_asis()
#' }
#' @importFrom magrittr %>%
#' @export
read_metabolon_asis <- function(file, sheet = 2, fid_var = 'COMP_ID', sid_var = 'CLIENT_IDENTIFIER'){

   assertive.files::assert_all_are_existing_files(file)

   d_f <- readxl::read_excel(file, sheet, col_names = FALSE)

   fvar_rows <- which(!is.na(d_f %>% extract_dt_col(1))) %>% magrittr::extract(1)
   svar_cols <- which(!is.na(d_f %>% extract_dt_row(1))) %>% magrittr::extract(1)
   fvar_cols <- fdata_cols <- 1:svar_cols
   svar_rows <- sdata_rows <- 1:fvar_rows

   fid_rows  <- fdata_rows <- expr_rows <- (fvar_rows+1):nrow(d_f)
   sid_cols  <- sdata_cols <- expr_cols <- (svar_cols+1):ncol(d_f)
   fid_cols  <-  d_f %>% extract_dt_row(fvar_rows) %>% magrittr::extract(1:svar_cols) %>% magrittr::equals(fid_var) %>% which()
   sid_rows  <-  d_f %>% extract_dt_col(svar_cols) %>% magrittr::extract(1:fvar_rows) %>% magrittr::equals(sid_var) %>% which()

   file %>% read_omics_asis(sheet      = sheet,
                       fid_rows   = fid_rows,      fid_cols   = fid_cols,
                       sid_rows   = sid_rows,      sid_cols   = sid_cols,
                       expr_rows  = expr_rows,     expr_cols  = expr_cols,
                       fvar_rows  = fvar_rows,     fvar_cols  = fvar_cols,
                       svar_rows  = svar_rows,     svar_cols  = svar_cols,
                       fdata_rows = fdata_rows,    fdata_cols = fdata_cols,
                       sdata_rows = svar_rows,     sdata_cols = sdata_cols,
                       transpose  = FALSE,
                       verbose    = TRUE)
}


#==========================================================================================================================


#' Read rnaseq data
#' @param file     character(1)
#' @param fid_var  character(1): feature id variable
#' @examples
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    file %>% read_rnaseq_asis()
#' }
#' @importFrom magrittr %>%
#' @export
read_rnaseq_asis <- function(file, fid_var = 'gene_id'){

   assertive.files::assert_all_are_existing_files(file)

   dt <- data.table::fread(file, integer64='numeric')
   expr_cols   <- dt %>% vapply(is.integer, logical(1)) %>% unname() %>% which()
   fdata_cols  <- dt %>% vapply(is.integer, logical(1)) %>% magrittr::not() %>% unname() %>% which()
   fid_col     <- which(names(dt)==fid_var)

   file %>% autonomics.import::read_omics_asis(sheet      = NULL,
                                          fid_rows   = 2:nrow(dt),   fid_cols   = fid_col,
                                          sid_rows   = 1,            sid_cols   = expr_cols,
                                          expr_rows  = 2:nrow(dt),   expr_cols  = expr_cols,
                                          fvar_rows  = 1,            fvar_cols  = fdata_cols,
                                          fdata_rows = 2:nrow(dt),   fdata_cols = fdata_cols,
                                          transpose  = FALSE,
                                          verbose    = TRUE)
}


#====================================================================================================================


#' Read exiqon data
#' @param file     character(1)
#' @examples
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx', package = 'subramanian.2016')
#'    file %>% read_exiqon_asis()
#' }
#' @importFrom magrittr %>%
#' @export
read_exiqon_asis <- function(file){
   assertive.files::assert_all_are_existing_files(file)
   dt <- autonomics.import::read_file(file, sheet=1)
   file %>% autonomics.import::read_omics_asis(sheet = 1,
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


