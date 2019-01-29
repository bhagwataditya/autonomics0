
#====================================================================================
#' Get/set analysis
#' @param object SummarizedExperiment
#' @param value list
#' @return analysis details (get) or updated object (set)
#' @rdname analysis

# Get
#====
#' @rdname analysis
#' @export
setGeneric("analysis",                                   function(object) standardGeneric("analysis"))

#' @rdname analysis
setMethod("analysis", signature("SummarizedExperiment"), function(object) S4Vectors::metadata(object)$analysis)

# Set
#====
#' @rdname analysis
#' @export
setGeneric(      "analysis<-",                                          function(object, value)  standardGeneric("analysis<-"))

#' @rdname analysis
setReplaceMethod("analysis", signature("SummarizedExperiment", "list"), function(object, value){
   S4Vectors::metadata(object)$analysis <- value
   object})


#====================================================================================
#' Get /set annotation
#' @param object eSet
#' @param value  list
#' @return annotation details (get) or updated eset (set)
#' @rdname annotation

# Get
#====
#' @rdname annotation
#' @export
setGeneric("annotation",                                         function(object)  standardGeneric("annotation")        )

#' @rdname annotation
setMethod("annotation", signature("SummarizedExperiment"),       function(object) S4Vectors::metadata(object)$annotation)

#' @rdname annotation
setMethod("annotation", signature("eSet"),                       function(object) Biobase::annotation(object)           )

#' @rdname annotation
setMethod("annotation", signature("EList"),                      function(object) object$annotation                     )

# Set
#====
#' @rdname annotation
#' @export
setGeneric("annotation<-",                                                     function(object, value)  standardGeneric("annotation<-") )

#' @rdname annotation
setReplaceMethod("annotation", signature("SummarizedExperiment", "character"), function(object, value){
   S4Vectors::metadata(object)$annotation <- value
   object })
#' @rdname annotation
setReplaceMethod("annotation", signature("eSet",                 "character"), function(object, value){
   Biobase::annotation(object) <- value
   object })
#' @rdname annotation
setReplaceMethod("annotation", signature("EList",                "character"), function(object, value){
   object$annotation <- value
   object })


#=================================================================================================
#' Get/Set contrast definitions
#' @param object SummarizedExperiment
#' @param value named contrast vector
#' @return updated object
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    object <- subramanian.2016::metabolon
#'    contrastdefs <- subramanian.2016::contrasts.metabolon
#'    autonomics.import::contrastdefs(object)                      # getter
#'    autonomics.import::contrastdefs(object) <- contrastdefs    # setter (conventional)
#'    autonomics.import::set_contrastdefs(object, contrastdefs)  # setter (piping)
#' }
#' @rdname contrastdefs

# Set
#====
#' @rdname contrastdefs
#' @export
setGeneric("contrastdefs<-",                                                     function(object, value)  standardGeneric("contrastdefs<-") )

#' @rdname contrastdefs
setReplaceMethod("contrastdefs", signature("SummarizedExperiment", "character"), function(object, value){ S4Vectors::metadata(object)$contrastdefs <- value; object})

#' @rdname contrastdefs
setReplaceMethod("contrastdefs", signature("SummarizedExperiment", "NULL"),      function(object, value){object})

#' @rdname contrastdefs
#' @export
set_contrastdefs <- function(object, value){
   autonomics.import::contrastdefs(object) <- value
   object
}


# Get
#====
#' @rdname contrastdefs
#' @export
setGeneric("contrastdefs",                                    function(object)   standardGeneric("contrastdefs") )

#' @rdname contrastdefs
setMethod("contrastdefs", signature("SummarizedExperiment"),  function(object) S4Vectors::metadata(object)$contrastdefs )


#=========================================================================

#' @title Get/Set counts
#' @description Get / Set counts matrix
#' @param object SummarizedExperiment
#' @param value count matrix (features x samples)
#' @return count matrix (get) or updated object (set)
#' @rdname counts

# Get
#====
#' @rdname counts
#' @export
setGeneric('counts',                                        function(object)   standardGeneric("counts"))

#' @rdname counts
setMethod("counts",    signature("SummarizedExperiment"),   function(object)   SummarizedExperiment::assays(object)$counts)


# Set
#====
#' @rdname counts
#' @export
setGeneric(      'counts<-',                                                function(object, value) standardGeneric("counts<-"))

#' @rdname counts
setReplaceMethod("counts",    signature("SummarizedExperiment", "matrix"),  function(object, value){
   SummarizedExperiment::assays(object)$counts <- value
   object })
#' @rdname counts
setReplaceMethod("counts",    signature("SummarizedExperiment", "numeric"), function(object, value){
   SummarizedExperiment::assays(object)$counts[] <- value
   object })


#=========================================================================
#' @title Get/Set exprs
#' @description Get / Set exprs matrix
#' @param object SummarizedExperiment
#' @param value ratio matrix (features x samples)
#' @return exprs matrix (get) or updated object (set)
#' @examples
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::exprs(object) <- 0
#'    object
#' }
#' @rdname exprs

# Get
#====
#' @rdname exprs
#' @export
setGeneric('exprs',                                        function(object)   standardGeneric("exprs"))

#' @rdname exprs
setMethod("exprs",    signature("SummarizedExperiment"),   function(object)   SummarizedExperiment::assays(object)$exprs)

#' @rdname exprs
setMethod("exprs",    signature("ExpressionSet"),          function(object)   Biobase::assayDataElement(object, "exprs"))

#' @rdname exprs
setMethod("exprs",    signature("EList"),                  function(object)   object$E)


# Set
#====
#' @rdname exprs
#' @export
setGeneric(      'exprs<-',                                                function(object, value) standardGeneric("exprs<-"))

#' @rdname exprs
setReplaceMethod("exprs",    signature("SummarizedExperiment", "matrix"),  function(object, value){
   SummarizedExperiment::assays(object)$exprs <- value
   object })
#' @rdname exprs
setReplaceMethod("exprs",    signature("SummarizedExperiment", "numeric"), function(object, value){
   SummarizedExperiment::assays(object)$exprs[] <- value
   object })

#' @rdname exprs
setReplaceMethod("exprs",    signature("ExpressionSet", "matrix"),         function(object, value){
   Biobase::assayDataElementReplace(object, "exprs",  value) })
#' @rdname exprs
setReplaceMethod("exprs",    signature("EList", "matrix"),                 function(object, value){   # exprs(obj) <- matrix(...)
   object$E[] <- value
   object })
#' @rdname exprs
setReplaceMethod("exprs",    signature("EList", "numeric"),                function(object, value){   # exprs(obj) <- 0
   object$E[] <- value
   object })


#=============================================================================================================
#' @title Get/Set fdata
#' @description Get/Set feature data
#' @param object SummarizedExperiment, eSet, or EList
#' @param value data.frame
#' @return feature dataframe (get) or updated object (set)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% autonomics.import::fdata() %>% head()
#' }
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::fdata(object) %<>% cbind(z=1)
#'    object
#' }
#' @rdname fdata

# Get
#====
#' @rdname fdata
#' @export
setGeneric('fdata',                                                           function(object)   standardGeneric('fdata'))

#' @rdname fdata
#' @importFrom magrittr %>%
setMethod(        'fdata',  signature('SummarizedExperiment'),                function(object){
                                                                                 object@elementMetadata@listData %>%
                                                                                 as.data.frame(check.names      = FALSE,
                                                                                               row.names        = object@elementMetadata@rownames,
                                                                                               stringsAsFactors = FALSE)})
# NOTES: (1) as.data.frame(object@elementMetadata) doesn't handle check.names correctly!
#        (2) SummarizedExperiment::rowData returns a DataFrame (which is not generic to e.g. EList objects)
#' @rdname fdata
setMethod(        'fdata',  signature('eSet'),                                function(object){
   Biobase::fData(object) })
#' @rdname fdata
setMethod(        'fdata',  signature('EList'),                               function(object){
   object$genes })

# Set
#====
#' @rdname fdata
#' @export
setGeneric(       'fdata<-',                                                  function(object, value)  standardGeneric('fdata<-'))

#' @rdname fdata
setReplaceMethod( 'fdata', signature('SummarizedExperiment', 'data.frame'),   function(object, value){
   SummarizedExperiment::rowData(object) <- S4Vectors::DataFrame(
      value,
      check.names = FALSE)
   object })
#' @rdname fdata
setReplaceMethod( 'fdata', signature('eSet',                 'data.frame'),   function(object, value){
   Biobase::fData(object) <- value
   object })
#' @rdname fdata
setReplaceMethod( 'fdata', signature('EList',                'data.frame'),   function(object, value){
   object$genes <- value
   object })

#==================================================================================
#' Get fvar levels
#' @param  object  SummarizedExperiment, eSet, or EList
#' @param  fvar    feature variable
#' @return fvar values
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% flevels('BIOCHEMICAL')
#' }
#' @importFrom magrittr %>%
#' @export
#' @importFrom magrittr %>%
#' @export
flevels <- function(object, fvar){
   object %>%
   autonomics.import::fvalues(fvar) %>%
  (function(x)if (is.factor(x)) levels(x) else unique(x))
}


#====================================================================================
#' @title Get/Set fnames
#' @description Get/Set feature names
#' @param object SummarizedExperiment
#' @param value character vector with feature names
#' @return feature name vector (get) or updated object (set)
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::fnames(object) %<>% paste0('XXX')
#'    object
#' }
#' @rdname fnames

# Get
#====
#' @rdname fnames
#' @export
setGeneric("fnames",                                                         function(object)          standardGeneric("fnames"))

#' @rdname fnames
setMethod("fnames",         signature("SummarizedExperiment"),               function(object)          rownames(object))

#' @rdname fnames
setMethod("fnames",         signature("eSet"),                               function(object)          Biobase::featureNames(object))

#' @rdname fnames
setMethod("fnames",         signature("EList"),                              function(object)          rownames(object$genes))

# Set
#====
#' @rdname fnames
#' @export
setGeneric("fnames<-",                                                       function(object, value)   standardGeneric("fnames<-"))

#' @rdname fnames
setReplaceMethod("fnames", signature("SummarizedExperiment", "character"),   function(object, value){  rownames(object) <- value
object})

#' @rdname fnames
setReplaceMethod("fnames", signature("eSet",                 "character"),   function(object, value){  Biobase::featureNames(object) <- value
object})

#' @rdname fnames
setReplaceMethod("fnames", signature("EList",                "character"),   function(object, value){  rownames(object) <- rownames(object$genes) <- value
object})


#================================================================================================
#' @title Get fvalues
#' @description Get fvar values
#' @param  object  SummarizedExperiment, eSet, or EList
#' @param  fvar    feature variable
#' @return fvar values
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% autonomics.import::fvalues('BIOCHEMICAL')
#'    subramanian.2016::metabolon %>% autonomics.import::fvalues(NULL)
#' }
#' @importFrom magrittr %>%
#' @export
fvalues <- function(object, fvar){

   # Return NULL output for NULL input
   if (is.null(fvar)) return(NULL)

   # Assert that fvar is present
   assertive.sets::assert_is_subset(fvar, autonomics.import::fvars(object))

   # Extract and return
   object %>%
   autonomics.import::fdata() %>%
   magrittr::extract2(fvar)
}


#================================================================================================
#' @title Get/Set fvars
#' @description Get/Set feature variables
#' @param object SummarizedExperiment
#' @param value character vector with feature variables
#' @return feature variables vector (get) or updated object (set)
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::fvars(object)[1] %<>% paste0('1')
#'    object
#' }
#' @rdname fvars

# Get
#====
#' @rdname fvars
#' @export
setGeneric("fvars",                                                          function(object)        standardGeneric("fvars"))

#' @rdname fvars
setMethod("fvars",          signature("SummarizedExperiment"),               function(object)        names(SummarizedExperiment::rowData(object)))

#' @rdname fvars
setMethod("fvars",          signature("eSet"),                               function(object)        Biobase::fvarLabels(object))

#' @rdname fvars
setMethod("fvars",          signature("EList"),                              function(object)        names(object$genes))

# Set
#====
#' @rdname fvars
#' @export
setGeneric("fvars<-",                                                        function(object, value)  standardGeneric("fvars<-") )

#' @rdname fvars
setReplaceMethod("fvars",  signature("SummarizedExperiment", "character"),   function(object, value){ names(SummarizedExperiment::rowData(object)) <- value
object })
#' @rdname fvars
setReplaceMethod("fvars",  signature("eSet", "character"),                   function(object, value){ Biobase::fvarLabels(object) <- value
object})
#' @rdname fvars
setReplaceMethod("fvars",  signature("EList", "character"),                  function(object, value){ names(object$genes) <- value
object})

#================================================================================================
#'@title Get/set limma results
#'@param object SummarizedExperiment
#'@param value list
#'@return limma results (get) or updated object (set)
#'@rdname limma

# Get
#====
#' @rdname limma
#' @export
setGeneric("limma",                                    function(object)   standardGeneric("limma") )

#' @rdname limma
setMethod("limma", signature("SummarizedExperiment"),  function(object) S4Vectors::metadata(object)$limma )

# Set
#====

#' @rdname limma
#' @export
setGeneric("limma<-",                                                function(object, value)  standardGeneric("limma<-") )

#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "array"), function(object, value){ S4Vectors::metadata(object)$limma <- value; object})

#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "NULL"), function(object, value){object})



#=================================================================================================
#' @title Get/Set prepro
#' @description Get/Set preprocessing details
#' @param object SummarizedExperiment
#' @param value  list
#' @return preprocessing details (get) or updated eset
#' @rdname prepro

# Set
#====
#' @rdname prepro
#' @export
setGeneric("prepro<-",                                                function(object, value)  standardGeneric("prepro<-") )

#' @rdname prepro
setReplaceMethod("prepro", signature("SummarizedExperiment", "list"), function(object, value){ S4Vectors::metadata(object)$prepro <- value; object})

#' @rdname prepro
setReplaceMethod("prepro", signature("SummarizedExperiment", "NULL"), function(object, value){object})

#' @rdname prepro
setReplaceMethod("prepro", signature("eSet", "list"),                 function(object, value){ Biobase::preproc(Biobase::experimentData(object)) <- value; object})

#' @rdname prepro
setReplaceMethod("prepro", signature("EList", "list"),                function(object, value){ object$prepro <- value; object})

# Get
#====
#' @rdname prepro
#' @export
setGeneric("prepro",                                    function(object)   standardGeneric("prepro") )

#' @rdname prepro
setMethod("prepro", signature("SummarizedExperiment"),  function(object) S4Vectors::metadata(object)$prepro )

#' @rdname prepro
setMethod("prepro", signature("eSet"),                  function(object){Biobase::preproc(Biobase::experimentData(object))})

#' @rdname prepro
setMethod("prepro", signature("EList"),                 function(object){object$prepro})


#==================================================================================================================
#' @title Get/Set sdata
#' @description Get/Set sample data
#' @param object SummarizedExperiment, eSet, or EList
#' @param value dataframe
#' @return sample dataframe (get) or updated object (set)
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::sdata(object) %<>% cbind(z=1)
#'    object
#' }
#' @rdname sdata

# Get
#====
#' @rdname sdata
#' @export
setGeneric('sdata',                                    function(object){ standardGeneric('sdata')})

#' @rdname sdata
setMethod('sdata',  signature('eSet'),                 function(object){ Biobase::pData(object) })

#' @rdname sdata
setMethod('sdata',  signature('EList'),                function(object){ object$targets })

#' @rdname sdata
#' @importFrom magrittr %>%
setMethod('sdata',  signature('SummarizedExperiment'), function(object){ object@colData@listData %>%
                                                                         as.data.frame(check.names      = FALSE,
                                                                                       row.names        = object@colData@rownames,
                                                                                       stringsAsFactors = FALSE) })
# NOTES: (1) "as.data.frame(object@colData)" doesn't handle check.names correctly!
#        (2)  SummarizedExperiment::colData returns a DataFrame, which is not generic (to e.g. EList)

# Set
#====

#' @rdname sdata
#' @export
setGeneric(      'sdata<-',                                                   function(object, value)  standardGeneric('sdata<-'))

#' @rdname sdata
setReplaceMethod('sdata',  signature('SummarizedExperiment', 'data.frame'),   function(object, value){
                                                                                 SummarizedExperiment::colData(object) <- S4Vectors::DataFrame(
                                                                                    value,
                                                                                    check.names = FALSE)
                                                                                 object })

#' @rdname sdata
setReplaceMethod('sdata',  signature('eSet',                 'data.frame'),   function(object, value){
                                                                                 Biobase::pData(object) <- value
                                                                                 object })

#' @rdname sdata
setReplaceMethod('sdata',  signature('EList',                'data.frame'),   function(object, value){
                                                                                 object$targets <- value
                                                                                 object })


#=======================================================================
#' Get (exprs) sign matrix
#' @param x SummarizedExperiment
#' @rdname sign-SummarizedExperiment
#' @export
setMethod("sign", signature("SummarizedExperiment"), function(x){
   sign(log(autonomics.import::exprs(x)))
})


#=====================================================================
#' @title Get/Set snames
#' @description Get/Set sample names
#' @param object SummarizedExperiment
#' @param value string vector with sample names
#' @return sample names vector (get) or updated eSet (set)
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::snames(object) %<>% paste0('XXX')
#'    object
#' }
#' @rdname snames

# Get
#====
#' @rdname snames
#' @export
setGeneric("snames",                                       function(object)   standardGeneric("snames"))

#' @rdname snames
setMethod('snames',  signature("SummarizedExperiment"),    function(object)   colnames(object))

#' @rdname snames
setMethod('snames',  signature("eSet"),                    function(object)   Biobase::sampleNames(object))

#' @rdname snames
setMethod('snames',  signature("EList"),                   function(object)   rownames(object$targets))

# Set
#====
#' @rdname snames
#' @export
setGeneric(      "snames<-",                                                function(object, value)  standardGeneric("snames<-"))

#' @rdname snames
setReplaceMethod("snames", signature("SummarizedExperiment", "character"),  function(object, value){
   colnames(object)  <- value
   object })
#' @rdname snames
setReplaceMethod("snames", signature("eSet", "character"),                  function(object, value){
   Biobase::sampleNames(object)  <- value
   object })
#' @rdname snames
setReplaceMethod("snames", signature("EList", "character"),                 function(object, value){
   colnames(object) <- rownames(object$targets) <-value
   object})

#=========================================================
#' @title Get slevels
#' @description Get svar levels
#' @param object SummarizedExperiment, eSet, or eList
#' @param svar sample var (character)
#' @return svar values (character)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% slevels('subgroup')
#'    subramanian.2016::metabolon %>% subgroup_levels()
#' }
#' @rdname slevels

# Get
#====
#
#' @rdname slevels
#' @importFrom magrittr %>%
#' @export
slevels <- function(object, svar){
   object %>%
   autonomics.import::svalues(svar) %>%
   (function(x) if (is.factor(x)) levels(x) else unique(x))
}

#' @rdname slevels
#' @importFrom magrittr %>%
#' @export
subgroup_levels <- function(object){
   object %>% autonomics.import::slevels('subgroup')
}


#=========================================================
#' @title Get/Set svalues
#' @description Get/Set svar values
#' @param object SummarizedExperiment, eSet, or eList
#' @param svar   sample var (character)
#' @param value  value vector
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    subramanian.2016::metabolon %>% svalues('subgroup')
#'    subramanian.2016::metabolon %>% subgroup_values()
#' }
#' @rdname svalues

# Get
#====
#' @rdname svalues
#' @importFrom magrittr %>%
#' @export
svalues <- function(object, svar){
   autonomics.import::sdata(object) %>% magrittr::extract2(svar) %>% (function(x) if (is.factor(x)) as.character(x) else x)
}

#' @rdname svalues
#' @importFrom magrittr %>%
#' @export
subgroup_values <- function(object){
   object %>% autonomics.import::svalues('subgroup')
}

#' @rdname svalues
#' @importFrom magrittr %>%
#' @export
sampleid_values <- function(object){
   object %>% autonomics.import::svalues('sample_id')
}

# Set
#====
#' @rdname svalues
#' @export
setGeneric(      'svalues<-',                                                         function(object, svar, value)  standardGeneric('svalues<-'))

#' @rdname svalues
setReplaceMethod('svalues',  signature('SummarizedExperiment', 'character', "ANY"),   function(object, svar, value){
   SummarizedExperiment::colData(object)[svar] <- value
   object
})


#=========================================================================
#' @title Get/Set svars
#' @description Get/Set sample variables
#' @param object SummarizedExperiment
#' @param value string fector with variable names
#' @return sample variable names (get) or updated eSet (set)
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::svars(object)[1] %<>% paste0('1')
#'    object
#' }
#' @rdname svars

# Get
#====
#' @rdname svars
#' @export
setGeneric(      "svars",                                                  function(object)         standardGeneric("svars") )

#' @rdname svars
setMethod(       "svars",   signature("SummarizedExperiment"),             function(object)         names(SummarizedExperiment::colData((object))))

#' @rdname svars
setMethod(       "svars",   signature("eSet"),                             function(object)         Biobase::varLabels(object))

#' @rdname svars
setMethod(       "svars",    signature("EList"),                           function(object)         names(object$targets))

# Set
#====
#' @rdname svars
#' @export
setGeneric(      "svars<-",                                                function(object, value)  standardGeneric("svars<-") )

#' @rdname svars
setReplaceMethod("svars",  signature("SummarizedExperiment", "character"), function(object, value){ names(SummarizedExperiment::colData(object)) <- value
object })
#' @rdname svars
setReplaceMethod("svars",  signature("eSet", "character"),                 function(object, value){ Biobase::varLabels(object) <- value
object })
#' @rdname svars
setReplaceMethod("svars",  signature("EList", "character"),                function(object, value){ names(object$targets) <- value
object })


#===========================================================================================
#' Get/set tmpplot
#' @param object SummarizedExperiment
#' @param value list
#' @return (temporary) plot (get) or updated object (set)
#' @rdname tmpplot

# Get
#====
#' @rdname tmpplot
#' @export
setGeneric("tmpplot",                                   function(object) standardGeneric("tmpplot"))

#' @rdname tmpplot
setMethod("tmpplot", signature("SummarizedExperiment"), function(object) S4Vectors::metadata(object)$tmpplot)

# Set
#====
#' @rdname tmpplot
#' @export
setGeneric(      "tmpplot<-",                                          function(object, value)  standardGeneric("tmpplot<-"))

#' @rdname tmpplot
setReplaceMethod("tmpplot", signature("SummarizedExperiment", "list"), function(object, value){
   S4Vectors::metadata(object)$tmpplot <- value
   object })


#=========================================================================
#' @title Get/Set weights
#' @description Get/Set weight matrix
#' @param object SummarizedExperiment
#' @param value ratio matrix (features x samples)
#' @param ... addtional params
#' @return weight matrix (get) or updated object (set)
#' @examples
#' if (require(billing.differentiation.data)){
#'    object <- billing.differentiation.data::rna.voomcounts
#'    autonomics.import::weights(object) <- 1
#'    object
#' }
#' @rdname weights

# Get
#====
#' @rdname weights
#' @export
setGeneric('weights',                                        function(object)   standardGeneric("weights"))

#' @rdname weights
setMethod("weights",    signature("SummarizedExperiment"),   function(object)   SummarizedExperiment::assays(object)$weights)

#' @rdname weights
setMethod("weights",    signature("ExpressionSet"),          function(object)   Biobase::assayDataElement(object, "weights"))

#' @rdname weights
setMethod("weights",    signature("EList"),                  function(object)   object$weights)

# Set
#====
#' @rdname weights
#' @export
setGeneric(      'weights<-',                                                function(object, value) standardGeneric("weights<-"))

#' @rdname weights
setReplaceMethod("weights",    signature("SummarizedExperiment", "matrix"),  function(object, value){
                                                                                SummarizedExperiment::assays(object)$weights <- value
                                                                                object })

#' @rdname weights
setReplaceMethod("weights",    signature("SummarizedExperiment", "numeric"), function(object, value){
                                                                                SummarizedExperiment::assays(object)$weights[] <- value
                                                                                object })

#' @rdname weights
setReplaceMethod("weights",    signature("ExpressionSet", "matrix"),         function(object, value){
                                                                                Biobase::assayDataElementReplace(object, "weights",  value) })

#' @rdname weights
setReplaceMethod("weights",    signature("EList", "matrix"),                 function(object, value){   # weights(obj) <- matrix(...)
                                                                                object$E[] <- value
                                                                                object })

#' @rdname weights
setReplaceMethod("weights",    signature("EList", "numeric"),                function(object, value){   # weights(obj) <- 0
                                                                                object$E[] <- value
                                                                                object })
