#' Create time contrasts
#' 
#' Create time contrasts from a subgroup matrix (nconc x ntime)
#' 
#' @param subgroup_matrix subgroup matrix: nconc x ntime
#' @return      contrast matrix: nconc x (ntime-1)
#' @noRd
create_time_contrasts <- function(subgroup_matrix, symbol = ' - '){
    contrastmat <- matrix(
        sprintf('%s%s%s', subgroup_matrix[, -1], symbol, subgroup_matrix[, -ncol(subgroup_matrix)]), 
        nrow = nrow(subgroup_matrix),  
        ncol = ncol(subgroup_matrix)-1)
    
    rownames(contrastmat) <- rownames(subgroup_matrix)
    colnames(contrastmat) <- sprintf('%s - %s', 
                    colnames(subgroup_matrix)[-1], colnames(subgroup_matrix)[-ncol(subgroup_matrix)])
    contrastmat
}


#' Create conc contrasts
#' 
#' Create conc contrasts from a subgroup matrix (nconc x ntime)
#' 
#' @param subgroup_matrix subgroup matrix: nconc x ntime
#' @return      contrast matrix: (nconc-1) x ntime
#' @noRd
create_conc_contrasts <- function(subgroup_matrix, symbol = ' - '){
    contrastmat <- matrix(
        sprintf('%s%s%s', subgroup_matrix[-1, ], symbol, subgroup_matrix[-nrow(subgroup_matrix), ]),
        nrow = nrow(subgroup_matrix)-1,  
        ncol = ncol(subgroup_matrix))
    colnames(contrastmat) <- colnames(subgroup_matrix)
    rownames(contrastmat) <- sprintf('%s - %s', 
                     rownames(subgroup_matrix)[-1], rownames(subgroup_matrix)[-nrow(subgroup_matrix)])
    contrastmat
}


#' Aggregate time contrasts across conc (or vice versa)
#' @param contrastmat contrast matrix
aggregate_contrasts <- function(contrastmat, dim){
    apply(contrastmat, 
          dim, 
          function(x){
              paste0(sprintf('(%s)/%d', x, length(x)), collapse = ' + ')})
}

aggregate_time_contrasts_across_conc <- function(time_contrast_matrix){
    aggregate_contrasts(time_contrast_matrix, 2)
}

aggregate_conc_contrasts_across_time <- function(conc_contrast_matrix){
    aggregate_contrasts(conc_contrast_matrix, 1)
}

split_values <- function(x){
    sep <- autonomics.import::guess_sep(x)
    dt <- data.table::data.table(x = x)
    dt[, data.table::tstrsplit(x, sep) ]
}

split_subgroup_values <- function(object){
    subgroupvalues <- autonomics.import::subgroup_values(object)
    cbind(subgroup = subgroupvalues, split_values(subgroupvalues))
}

split_subgroup_levels <- function(object){
    subgrouplevels <- autonomics.import::subgroup_levels(object)
    cbind(subgroup = subgrouplevels, split_values(subgrouplevels))
}

create_subgroup_matrix <- function(object){
    dt <- split_subgroup_levels(object)
    subgroup_matrix <- as.matrix(data.table::dcast(
        dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    subgroup_matrix
}


diff_contrasts <- function(object){
    
    # Subgroup matrix
    subgroup_matrix <- create_subgroup_matrix(object)
    
    # Contrast matrix
    time_contrast_matrix <- create_time_contrasts(subgroup_matrix)
    time_contrast_names  <- create_time_contrasts(subgroup_matrix, '__')
    conc_contrast_matrix <- create_conc_contrasts(subgroup_matrix)
    conc_contrast_names  <- create_conc_contrasts(subgroup_matrix, '__')
    
    # Contrast vector
    time_contrasts <- 
        structure(c(time_contrast_matrix), names = c(time_contrast_names))
    conc_contrasts <- 
        structure(c(conc_contrast_matrix), names = c(conc_contrast_names))
    time_contrasts_across_conc <- 
        aggregate_time_contrasts_across_conc(time_contrast_matrix)
    conc_contrasts_across_time <- 
        aggregate_conc_contrasts_across_time(conc_contrast_matrix)
    
    # Return
    c(  time_contrasts, 
        conc_contrasts, 
        time_contrasts_across_conc, 
        conc_contrasts_across_time)
}

add_limma2 <- function(object, contrastdefs = diff_contrasts(object)){
    S4Vectors::metadata(object)$contrastdefs <- contrastdefs
    object %<>% autonomics.find::add_limma(contrastdefs = contrastdefs)
    object
}

#' Plot contrastogram
#' @param object SummarizedExperiment
#' @param directed TRUE or FALSE: whether to distinguish up and downregulations
#' @param subgroup_colors named color vector (names = subgroups)
#' @examples 
#' if (require(autonomics.data)){
#'     object <- autonomics.data::glutaminase
#'     plot_contrastogram(object)
#'     plot_contrastogram(object, directed = FALSE)
#' }
plot_contrastogram <- function(
    object, directed = TRUE, subgroup_colors = default_color_values2(object), 
    sparse = FALSE
){
    # Perform limma
    object %<>% add_limma2()
    contrastogram_matrices <- compute_connections(
        object, directed = directed, subgroup_colors = subgroup_colors)
    connection_sizes  <- contrastogram_matrices$connection_sizes
    connection_colors <- contrastogram_matrices$connection_colors
    
        widths <- scales::rescale(connection_sizes, c(0.01,30))
        if (sparse) connection_colors[connection_sizes/nrow(object)<0.50] <- "0"
        
    # Plot diagram
    dt <- split_subgroup_levels(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    #dir.create('~/importomicscache/contrastogram')
    #pdf('~/importomicscache/contrastogram/directed_contrastogram.pdf', width = 9, height = 9)
    diagram::plotmat(connection_sizes, nperrow, relsize = 1, box.size = 0.05, 
        name = rownames(connection_sizes), box.col = subgroup_colors, 
        box.type = 'square', arr.lwd = widths, # sqrt(connection_sizes)
        arr.lcol = connection_colors, arr.col = connection_colors) #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
}


default_color_values2 <- function(object){
    autonomics.plot::default_color_values(object)[autonomics.import::subgroup_levels(object)]
}

compute_connections <- function(
    object, 
    directed = FALSE,
    subgroup_colors = default_color_values2(object)
){
    
    # subgroup matrix, difference contrasts, limma
    pvalues <- autonomics.import::limma(object)[, , 'p']
    effects <- autonomics.import::limma(object)[, , 'effect']
    is_significant <-  pvalues < 0.05
    is_up          <- (pvalues < 0.05) & (effects > 0)
    is_down        <- (pvalues < 0.05) & (effects < 0)
    
    n <- apply(if (directed) is_up else is_significant, 2,sum, na.rm = TRUE)
    ndown <- apply(is_down, 2, sum, na.rm = TRUE)
                      
    # Create diagram
    subgrouplevels <- autonomics.import::subgroup_levels(object)
    connection_sizes <- connection_colors <- matrix(
        0,
        nrow = length(subgrouplevels),
        ncol = length(subgrouplevels), 
        dimnames = list(subgrouplevels, subgrouplevels))
    
    # Add time contrasts to diagram
    subgroup_matrix <- create_subgroup_matrix(object)
    for (i in 1:nrow(subgroup_matrix)){
        for (j in 2:ncol(subgroup_matrix)){
            from <- subgroup_matrix[i, j-1]
            to   <- subgroup_matrix[i, j]
            connection_sizes[to, from]         <- n[[paste0(to, '__', from)]]
            connection_colors[to, from] <- subgroup_colors[[to]]
            if (directed){
                connection_sizes[from, to] <- ndown[[paste0(to, '__', from)]]
                connection_colors[from, to] <- subgroup_colors[[to]]
            }
        }
    }
    
    # Add conc contrasts to diagram
    for (j in 1:ncol(subgroup_matrix)){
        for (i in 2:ncol(subgroup_matrix)){
            from <- subgroup_matrix[i-1, j]
            to   <- subgroup_matrix[i, j]
            connection_sizes[to, from]         <- n[[paste0(to, '__', from)]]
            connection_colors[to, from] <- subgroup_colors[[to]]
            if (directed){
                connection_sizes[from, to] <- ndown[[paste0(to, '__', from)]]
                connection_colors[from, to] <- subgroup_colors[[to]]
            }
        }
    }
    
    # Return
    list(connection_sizes = connection_sizes, connection_colors = connection_colors)
}

