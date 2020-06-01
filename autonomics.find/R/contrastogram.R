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

create_subgroup_dt <- function(object){
    subgrouplevels <- autonomics.import::subgroup_levels(object)
    sep <- autonomics.import::guess_sep(subgrouplevels)
    dt <- data.table::data.table(subgroup = subgrouplevels)
    dt %<>% cbind(dt[, data.table::tstrsplit(subgroup, sep) ])
    dt
}

create_subgroup_matrix <- function(object){
    dt <- create_subgroup_dt(object)
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
#' @examples 
#' if (require(autonomics.data)){
#'     object <- autonomics.data::glutaminase
#'     plot_contrastogram(object)
#' }
plot_contrastogram <- function(object, subgroup_colors = default_color_values2(object)){
    
    # Perform limma
    object %<>% add_limma2()
    contrastogram_matrices <- create_contrastogram_matrices(object)
    connections <- contrastogram_matrices$connections
    colors <- contrastogram_matrices$colors
    
        # 
        widths <- scales::rescale(connections, c(0.01,30))
        colors[connections/nrow(object)<0.20] <- "0"
        
    # Plot diagram
    dt <- create_subgroup_dt(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    #dir.create('~/importomicscache/contrastogram')
    #pdf('~/importomicscache/contrastogram/contrastogram.pdf', width = 9, height = 9)
    diagram::plotmat(connections, nperrow, relsize = 1, box.size = 0.05, 
        name = rownames(connections), box.col = subgroup_colors, 
        box.type = 'square', arr.lwd = widths, # sqrt(connections)
        arr.lcol = colors, arr.col = colors) #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
}


default_color_values2 <- function(object){
    autonomics.plot::default_color_values(object)[autonomics.import::subgroup_levels(object)]
}

create_contrastogram_matrices <- function(
    object, 
    subgroup_colors = default_color_values2(object)
){
    
    # subgroup matrix, difference contrasts, limma
    pvalues <- autonomics.import::limma(object)[, , 'p']
    nsignificant <- apply(pvalues, 2, function(x) sum(x<0.05, na.rm=TRUE))
    
    # Create diagram
    subgrouplevels <- autonomics.import::subgroup_levels(object)
    connections <- colors <- matrix(
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
            connections[to, from] <- nsignificant[[paste0(to, '__', from)]]
            colors[to, from] <- subgroup_colors[[to]]}}
    
    # Add conc contrasts to diagram
    for (j in 1:ncol(subgroup_matrix)){
        for (i in 2:ncol(subgroup_matrix)){
            from <- subgroup_matrix[i-1, j]
            to   <- subgroup_matrix[i, j]
            connections[to, from] <- nsignificant[[paste0(to, '__', from)]]
            colors[to, from] <- subgroup_colors[[to]]}}
    
    # Return
    list(connections = connections, colors = colors)
}

