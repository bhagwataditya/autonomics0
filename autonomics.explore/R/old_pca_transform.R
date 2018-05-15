# pca_transform <- function(object, dim = c(1,2), verbose = TRUE, na.impute = FALSE){
# 
#    # Assert
#    autonomics.import::assert_is_valid_eset(object)
#    assertive.properties::assert_is_non_empty(autonomics.import::exprs(object))
# 
#    # Replace -Inf, arising from log2(0)
#    idx <- autonomics.import::exprs(object)==-Inf
#    if (sum(idx, na.rm=TRUE)>0){
#       if (verbose) autonomics.support::cmessage('\t\tReplace -Inf by NA')
#       autonomics.import::exprs(object)[idx] <- NA
#    }
# 
#    # Impute NAs
#    if (na.impute)   object %<>% autonomics.import::replace_nas_with_zeros(verbose = verbose)
# 
#    # Switch sign if all negative
#    # (required to prevent singularities in PCA transformation)
#    idx <- !is.na(autonomics.import::exprs(object))
#    if (all(sign(autonomics.import::exprs(object)[idx])==-1)){
#       if (verbose){autonomics.support::cmessage('\t\tAll values negative: flip signs to prevent singularities.')}
#       autonomics.import::exprs(object) %<>% magrittr::multiply_by(-1)
#    }
# 
#    # Limit to features with no NA values
#    selector <- !matrixStats::rowAnys(is.na(autonomics.import::exprs(object)))
#    if (verbose){autonomics.support::cmessage('\t\tUse %s/%s features with available values', sum(selector), length(selector))}
#    pcaset <- object[selector, ]
# 
#    # PCA transform
#    prepDF <- data.frame(feature = rownames(pcaset), autonomics.import::exprs(pcaset))
#    mpmRes <- suppressPackageStartupMessages(
#       mpm::mpm(prepDF, logtrans = FALSE, closure = 'none',  center = 'double',
#                normal = 'global', row.weight = 'mean', col.weight = 'constant')
#    )
# 
#    # Convert into plot coordinates
#    plotRes <- suppressPackageStartupMessages(mpm::plot.mpm(mpmRes, do.plot = FALSE, dim = dim))
#    samples <- data.frame(sample = autonomics.import::snames(object), x = plotRes$Columns$X, y = plotRes$Columns$Y)
#    features <- data.frame(feature = autonomics.import::fnames(object), x = NA, y = NA)
#    features$x[selector] <- plotRes$Rows$X
#    features$y[selector] <- plotRes$Rows$Y
#    var <- data.frame(
#       x = round(mpmRes$contrib[dim[1]]*100),
#       y = round(mpmRes$contrib[dim[2]]*100)
#    )
# 
#    # Return as list
#    return(list(samples = samples, features = features, var = var))
# }

