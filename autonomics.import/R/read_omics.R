
#' Read omics data as is from file
#' @param file character(1)
#' @param sheet numeric(1) or character(1)
#' @return data.table
#' @importFrom magrittr %>%
#' @export
read_asis <- function(file, sheet){
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
#'    dt <- file %>% read_asis(sheet = 2)
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
#' @examples
#' require(magrittr)
#'
#' # METABOLON
#' if (require(autonomics.data)){
#'   file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'   file %>% autonomics.import::read_omics(sheet      = 2,
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
#'    file %>% autonomics.import::read_omics(fid_rows   = 73,       fid_cols   = 27:1343,
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
#'    file %>% autonomics.import::read_omics(sheet      = 1,
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
#'    file %>% autonomics.import::read_omics(fid_rows   = 2:5893,   fid_cols   = 275,
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
#'    file %>% autonomics.import::read_omics(fid_rows   = 2:9783,  fid_cols   = 383,
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
#'    file %>% autonomics.import::read_omics(fid_rows   = 2:48821,   fid_cols   = 1,
#'                                           sid_rows   = 1,         sid_cols   = 5:86,
#'                                           expr_rows  = 2:48821,   expr_cols  = 5:86,
#'                                           fvar_rows  = 1,         fvar_cols  = 1:4,
#'                                           fdata_rows = 2:48821,   fdata_cols = 1:4,
#'                                           transpose  = FALSE)
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
   transpose  = FALSE
){
   # Assert
   assertive.types::assert_is_a_bool(transpose)

   # read
   dt <- autonomics.import::read_asis(file, sheet)

   # Extract exprs
   fids1  <- dt %>% autonomics.import::extract_rectangle(fid_rows,  fid_cols,  transpose = transpose, drop = TRUE)
   sids1  <- dt %>% autonomics.import::extract_rectangle(sid_rows,  sid_cols,  transpose = transpose, drop = TRUE)
   exprs1 <- dt %>% autonomics.import::extract_rectangle(expr_rows, expr_cols, transpose = transpose) %>% (function(y){class(y) <- 'numeric'; y})

   # Extract feature and sample annotations
   fdata_available <- !is.null(fvar_rows) & !is.null(fvar_cols)
   sdata_available <- !is.null(svar_rows) & !is.null(svar_cols)
   if (fdata_available){
      fvars1 <- dt %>% autonomics.import::extract_rectangle(fvar_rows,  fvar_cols,  transpose = transpose, drop = TRUE)
      fdata1 <- dt %>% autonomics.import::extract_rectangle(fdata_rows, fdata_cols, transpose = transpose) %>%
                       magrittr::set_colnames(fvars1) %>%
                       data.frame(stringsAsFactors = FALSE, check.names = FALSE)
   }
   if (sdata_available){
      svars1 <- dt %>% autonomics.import::extract_rectangle(svar_rows,  svar_cols,  transpose =  transpose, drop = TRUE)
      sdata1 <- dt %>% autonomics.import::extract_rectangle(sdata_rows, sdata_cols, transpose = !transpose) %>%
                       magrittr::set_colnames(svars1) %>%
                       data.frame(stringsAsFactors = FALSE, check.names = FALSE)
   }

   # Remove features with missing identifiers (happens in MaxQuant files due to empty interspersed lines)
   valid_rows <- !is.na(fids1)
   if (any(!valid_rows)){
      autonomics.support::cmessage('\t\tRm %d features with missing feature_ids', sum(!valid_rows))
      fids1  %<>% magrittr::extract(valid_rows)
      fdata1 %<>% magrittr::extract(valid_rows,)
      exprs1 %<>% magrittr::extract(valid_rows, )
   }

   # Name features and samples
   fids1 %<>% autonomics.support::uniquify('make.unique')
   sids1 %<>% autonomics.support::uniquify('make.unique')
   rownames(exprs1) <- fids1
   colnames(exprs1) <- sids1
   if (fdata_available){   rownames(fdata1) <- fids1;  fdata1 %<>% cbind(feature_id = fids1, .)  }
   if (sdata_available){   rownames(sdata1) <- sids1;  sdata1 %<>% cbind(sample_id  = sids1, .)  }

   # Wrap into Sumexp
   object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = exprs1))
   if (fdata_available) autonomics.import::fdata(object) <- fdata1
   if (sdata_available) autonomics.import::sdata(object) <- sdata1

   # Return
   object
}

#' Read somascan data
#' @param file character(1)
#' @param fid_var  character(1): feature id variable
#' @param sid_var  character(1): sample id variable
#' @return Summarizedexperiment
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::read_somascan()
#' }
#' @importFrom magrittr %>%
#' @export
read_somascan <- function(file, fid_var = 'SeqId', sid_var = 'SampleId'){

   dt <- data.table::fread(file, header = FALSE)
   fdata_row1 <- 1 + which(dt[[1]] == '^TABLE_BEGIN')[1]
   fvar_col <- which(dt %>% extract_dt_row(fdata_row1) != '')[1]
   svar_row <- which(dt[[1]] != '' &  (1:nrow(dt) > fdata_row1))[1]

   fvars <- dt %>% extract_dt_col(fvar_col) %>% magrittr::extract(fdata_row1:svar_row)
   svars <- dt %>% extract_dt_row(svar_row) %>% magrittr::extract(1:fvar_col-1)

   fid_row <- fdata_row1 + -1 + which(fvars == fid_var)
   sid_col <- which(svars == sid_var)

   file %>% read_omics(sheet      =  NULL,
                       fid_rows   =  fid_row,                   fid_cols   = (fvar_col+1):ncol(dt),
                       sid_rows   = (svar_row+1):nrow(dt),      sid_cols   =  sid_col,
                       expr_rows  = (svar_row+1):nrow(dt),      expr_cols  = (fvar_col+1):ncol(dt),
                       fvar_rows  =  fdata_row1:(svar_row-1),   fvar_cols  =  fvar_col,
                       fdata_rows =  fdata_row1:(svar_row-1),   fdata_cols = (fvar_col+1):ncol(dt),
                       svar_rows  =  svar_row,                  svar_cols  = 1:(fvar_col-1),
                       sdata_rows = (svar_row+1):nrow(dt),      sdata_cols = 1:(fvar_col-1),
                       transpose  = TRUE)

}

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
#'    file %>% autonomics.import::read_metabolon()
#' }
#' @importFrom magrittr %>%
#' @export
read_metabolon <- function(file, sheet = 2, fid_var = 'COMP_ID', sid_var = 'CLIENT_IDENTIFIER'){

   d_f <- readxl::read_excel(file, sheet, col_names = FALSE)
   fvar_row <- which(!is.na(d_f %>% extract_dt_col(1))) %>% magrittr::extract(1)
   svar_col <- which(!is.na(d_f %>% extract_dt_row(1))) %>% magrittr::extract(1)
   fid_col   <- d_f %>% extract_dt_row(fvar_row) %>% magrittr::extract(1:svar_col) %>% magrittr::equals(fid_var) %>% which()
   sid_row   <- d_f %>% extract_dt_col(svar_col) %>% magrittr::extract(1:fvar_row) %>% magrittr::equals(sid_var) %>% which()

   file %>% read_omics(sheet      = sheet,
                       expr_rows  = (fvar_row+1):nrow(d_f),    expr_cols  = (svar_col+1):ncol(d_f),
                       fid_rows   = (fvar_row+1):nrow(d_f),    fid_cols   =  fid_col,
                       sid_rows   =  sid_row,                  sid_cols   = (svar_col+1):ncol(d_f),
                       fdata_rows = (fvar_row+1):nrow(d_f),    fdata_cols = 1:svar_col,
                       sdata_rows = 1:fvar_row,                sdata_cols = (svar_col+1):ncol(d_f),
                       transpose  = FALSE)
}

#' Read rnaseq data
#' @param file     character(1)
#' @param fid_var  character(1): feature id variable
#' @examples
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    file %>% read_rnaseq()
#' }
#' @importFrom magrittr %>%
#' @export
read_rnaseq <- function(file, fid_var = 'gene_id'){

   dt <- data.table::fread(file, integer64='numeric')
   expr_cols   <- dt %>% vapply(is.integer, logical(1)) %>% unname() %>% which()
   fdata_cols  <- dt %>% vapply(is.integer, logical(1)) %>% magrittr::not() %>% unname() %>% which()
   fid_col     <- which(names(dt)==fid_var)

   file %>% autonomics.import::read_omics(sheet      = NULL,
                                          fid_rows   = 2:nrow(dt),   fid_cols   = fid_col,
                                          sid_rows   = 1,            sid_cols   = expr_cols,
                                          expr_rows  = 2:nrow(dt),   expr_cols  = expr_cols,
                                          fvar_rows  = 1,            fvar_cols  = fdata_cols,
                                          fdata_rows = 2:nrow(dt),   fdata_cols = fdata_cols,
                                          transpose  = FALSE)
}


#' Read exiqon data
#' @param file     character(1)
#' @examples
#' if require(subramanian.2016){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx', package = 'subramanian.2016')
#'    file %>% read_exiqon()
#' }
#' @importFrom magrittr %>%
#' @export
read_exiqon <- function(file){
   dt <- autonomics.import::read_asis(file, sheet=1)
   file %>% autonomics.import::read_omics(sheet = 1,
                                          fid_rows   = 1,                       fid_cols   = 2:(ncol(dt)-2),
                                          sid_rows   = 2:(nrow(dt)-3),          sid_cols   = 1,
                                          expr_rows  = 2:(nrow(dt)-3),          expr_cols  = 2:(ncol(dt)-2),
                                          fvar_rows  = (nrow(dt)-2):nrow(dt),   fvar_cols  = 1,
                                          svar_rows  = 1,                       svar_cols  = (ncol(dt)-1):ncol(dt),
                                          fdata_rows = (nrow(dt)-2):nrow(dt),   fdata_cols = 2:(ncol(dt)-2),
                                          sdata_rows = 2:(nrow(dt)-3),          sdata_cols = (ncol(dt)-1):ncol(dt),
                                          transpose  = TRUE)
}
