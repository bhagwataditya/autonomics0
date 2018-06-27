
#==========
# MAXQUANT
#==========

#' Create maxquant design df
#'
#' For automated infer of design, use the following sample naming
#' scheme (before running max quant): WT(L).KD(M).OE(H).R1
#'
#' @param proteingroups_file string: full path to protein groups file
#' @param sample_file        string: full path to sample design file
#' @param value_type         string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @param infer_design       logical: should design be infered from sample ids (see details)
#' @param ... backward compatibility to deprecated functions
#' @return sample design dataframe
#' @examples
#' require(magrittr)
#'
#' # LABELED RATIOS and INTENSITIES
#' if (require(autonomics.data)){
#'    file <- system.file(
#'       'extdata/stemcell.differentiation/maxquant/proteinGroups.txt',
#'        package = 'autonomics.data')
#'    file %>% create_maxquant_sample_df()
#'    file %>% create_maxquant_sample_df(infer_design = TRUE)
#'    file %>% create_maxquant_sample_df(infer_design = TRUE, value_type = 'raw.intensity')
#'    file %>% create_maxquant_sample_df(infer_design = TRUE, value_type = 'raw.ratio')
#' }
#'
#' # LFQ INTENSITIES
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% autonomics.import::create_maxquant_sample_df()
#'    file %>% autonomics.import::create_maxquant_sample_df(infer_design = TRUE)
#' }
#'
#' # UNLABELED INTENSITIES
#' if (require(alnoubi.2017)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'alnoubi.2017')
#'    file %>% autonomics.import::create_maxquant_sample_df()
#'    file %>% autonomics.import::create_maxquant_sample_df(infer_design = TRUE)
#' }
#'
#' # REPORTER INTENSITIES
#' if (require(billing.vesicles)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'billing.vesicles')
#'    file %>% autonomics.import::create_maxquant_sample_df()
#'    file %>% autonomics.import::create_maxquant_sample_df(infer_design = TRUE)
#' }
#'
#' @author Aditya Bhagwat
#' @importFrom magrittr %>%
#' @export
create_maxquant_sample_df <- function(
   proteingroups_file,
   value_type   = autonomics.import::infer_maxquant_value_type(proteingroups_file),
   infer_design = FALSE
){

   # Assert
   assertive.files::assert_all_are_dirs(dirname(proteingroups_file))
   #assertive.files::assert_all_are_readable_files(proteingroups_file, warn_about_windows = FALSE)

   # Read
   DT <- autonomics.support::cfread(proteingroups_file)

   # Extract
   sampleids <- DT %>% autonomics.import::get_maxquant_value_columns(value_type) %>% names()

   # Infer design and return
   if (infer_design){
      design <- sampleids %>%
         autonomics.import::designify_maxquant_sampleids()   %>%
         autonomics.import::infer_design_from_sampleids('.')
      design$sample_id <- sampleids
      rownames(design) <- sampleids
      return(design)
   }

   # Return design template if not able to infer
   return(data.frame(
      sample_id = sampleids,
      subgroup  = '',
      replicate = '',
      block     = ''))

}


#'@rdname create_maxquant_sample_df
#'@export
create_maxquant_design_df <- function(...){
   .Deprecated('create_maxquant_sample_df')
   create_maxquant_sample_df(...)
}

#'@rdname create_maxquant_sample_df
#'@export
create_sample_design_df <- function(...){
   .Deprecated('create_maxquant_design_df')
   create_maxquant_design_df(...)
}


#' Create maxquant design file
#'
#' For automated infer of design, use the following sample naming
#' scheme (before running max quant): WT(L).KD(M).OE(H).R1
#'
#' @param proteingroups_file string: full path to protein groups file
#' @param sample_file        string: full path to sample design file
#' @param value_type         string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @param infer_design       logical: should design be infered from sample ids (see details)
#' @param ... backward compatibility to deprecated functions
#' @return sample design dataframe
#' @examples
#' require(magrittr)
#'
#' # LABELED RATIOS and INTENSITIES
#' if (require(autonomics.data)){
#'    file <- system.file(
#'              'extdata/stemcell.differentiation/maxquant/proteinGroups.txt',
#'               package = 'autonomics.data')
#'    file %>% create_maxquant_sample_df()
#'    file %>% create_maxquant_sample_file(sample_file = tempfile())
#' }
#'
#' @importFrom magrittr %>%
#' @export
create_maxquant_sample_file <- function(
   proteingroups_file,
   sample_file  = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   value_type   = autonomics.import::infer_maxquant_value_type(proteingroups_file),
   infer_design = FALSE
){

   # Create
   autonomics.import::create_maxquant_sample_df(proteingroups_file = proteingroups_file,
                                                value_type         = value_type,
                                                infer_design       = infer_design) %>%
   autonomics.import::write_sample_file(sample_file)

   # Return
   return(invisible(sample_file))
}

#' @rdname create_maxquant_sample_file
#' @export
create_maxquant_design_file <- function(...){
   create_maxquant_sample_file(...)
}

#' @rdname create_maxquant_sample_file
#' @export
create_sample_design_file <- function(...){
   .Deprecated('write_maxquant_design')
   write_maxquant_design(...)
}


#========
# EXIQON
#========


#' Create exiqon design dataframe
#' @param exiqon_file  string
#' @param sample_file  string
#' @param value_type   string: any value in autonomics.import::MAXQUANT_VALUE_TYPES
#' @return string: path to design file
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    exiqon_file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                                package = 'subramanian.2016')
#'    exiqon_file %>% autonomics.import::create_exiqon_sample_df()
#'    exiqon_file %>% autonomics.import::create_exiqon_sample_df(infer_design = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
create_exiqon_sample_df <- function(
   exiqon_file,
   infer_design = FALSE
){
   sampleids <- autonomics.import::load_exiqon_sdata(exiqon_file)  %>%
                magrittr::extract2('sample_id')
   if (infer_design){
      sampleids %>% autonomics.import::infer_design_from_sampleids()
   } else {
      return(data.frame(sample_id = sampleids,
                        subgroup  = '',
                        replicate = '',
                        block     = ''))
   }
}


#' Create exiqon design file
#' @param exiqon_file  string: path to exiqon file
#' @param sample_file  string: path to sample file
#' @param infer_design logical: whether to infer design from sampleids
#' @return string: path to design file
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016){
#'    exiqon_file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx',
#'                                package = 'subramanian.2016')
#'    exiqon_file %>% autonomics.import::create_exiqon_sample_file(sample_file = tempfile())
#'
#' }
#' @importFrom magrittr %>%
#' @export
create_exiqon_sample_file <- function(
   exiqon_file,
   sample_file  = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   infer_design = FALSE
){
   autonomics.import::create_exiqon_sample_df(exiqon_file,
                                              infer_design = infer_design) %>%
   autonomics.import::write_sample_file(sample_file)
}


#===========
# METABOLON
#===========


#' Create metabolon design dataframe
#' @param  metabolon_file  string: path to metabolon file
#' @param  sample_file     string: path to sample file
#' @param  infer_design    logical: whether to infer design from CLIENT_IDENTIFIER
#' @return design dataframe
#' @examples
#' if (require(autonomics.data)){
#'    metabolon_file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                                   package = 'autonomics.data')
#'    metabolon_file %>% autonomics.import::create_metabolon_sample_df()
#'    metabolon_file %>% autonomics.import::create_metabolon_sample_df(infer_design = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
create_metabolon_sample_df <- function(
   metabolon_file,
   infer_design = FALSE
){
   sdata1 <- metabolon_file %>% autonomics.import::load_metabolon_sdata()
   design_df <- data.frame(CLIENT_IDENTIFIER = sdata1$CLIENT_IDENTIFIER,
                           row.names         = sdata1$CLIENT_IDENTIFIER)
   if (infer_design){
      design_df %<>% cbind(autonomics.import::infer_design_from_sampleids(.$CLIENT_IDENTIFIER))
   } else {
      design_df$sample_id <- sdata1$CLIENT_IDENTIFIER
      design_df$subgroup  <- sdata1$Group
   }
   design_df
}

#' Create metabolon design file
#' @param metabolon_file  string: path to metabolon file
#' @param sample_file     string: path to sample file
#' @param infer_design    logical: whether to infer design from sampleids
#' @return path to metabolon file
#' @examples
#' if (require(autonomics.data)){
#'    metabolon_file <- system.file('extdata/glutaminase/glutaminase.xlsx',
#'                                   package = 'autonomics.data')
#'    metabolon_file %>% autonomics.import::create_metabolon_sample_file(
#'                          sample_file = tempfile())
#' }
#' @importFrom magrittr %>%
#' @export
create_metabolon_sample_file <- function(
   metabolon_file,
   sample_file     = paste0(dirname(proteingroups_file), '/sample_design.txt'),
   infer_design    = FALSE
){
   df <- metabolon_file %>%
         autonomics.import::create_metabolon_sample_df(infer_design = infer_design)
   df %>% autonomics.import::write_sample_file(sample_file)
   df
}


#======
# SOMA
#======


#' Create soma design dataframe
#' @param soma_file    string: path to soma file
#' @param infer_design logical: whether to infer design from sampleids
#' @return sample design dataframe
#' @examples
#' if (require(autonomics.data)){
#'    soma_file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                              package = 'autonomics.data')
#'    soma_file %>% autonomics.import::create_soma_sample_df()
#'    soma_file %>% autonomics.import::create_soma_sample_df(infer_design = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
create_soma_sample_df <- function(
   soma_file,
   infer_design = FALSE
){

   sdata1 <- autonomics.import::load_soma_sdata(soma_file)
   design_df <- data.frame(SampleId  = sdata1$SampleId,
                           row.names = sdata1$SampleId)
   if (infer_design){
      design_df %<>% cbind(autonomics.import::infer_design_from_sampleids(.$SampleId))
   } else {
      design_df$sample_id <- sdata1$SampleId
      design_df$subgroup  <- sdata1$SampleGroup
   }
   design_df
}


#' Create soma design file
#' @param soma_file     string: path to soma file
#' @param sample_file   string: path to sample file
#' @param infer_design  logical: whether to infer design from sampleids
#' @return path to soma design file
#' @examples
#'    soma_file <- system.file('extdata/stemcell.comparison/stemcell.comparison.adat',
#'                              package = 'autonomics.data')
#'    soma_file %>% autonomics.import::create_soma_sample_file(sample_file = tempfile())
#' @importFrom magrittr %>%
#' @export
create_soma_sample_file <- function(
   soma_file,
   sample_file,
   infer_design = FALSE
){
   df <- soma_file %>%
         autonomics.import::create_soma_sample_df(infer_design = infer_design)
   df %>% autonomics.import::write_sample_file(sample_file)
   df
}

