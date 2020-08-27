    diagram::plotmat(connection_sizes, nperrow, relsize = 1, box.size = 0.05, 
        name = rownames(connection_sizes), box.col = subgroup_colors, 
        box.type = 'square', arr.lwd = widths, # sqrt(connection_sizes)
        arr.lcol = connection_colors, arr.col = connection_colors) #, arr.lcol = log2(1+diagram_matrix))

x <- connection_sizes[c(1,2,5,6), c(1,2,5,6)]

sgmat <- create_subgroup_matrix(glutaminase)
coldiffs <- create_time_contrasts(sgmat)
rowdiffs <- create_conc_contrasts(sgmat)


nodes <- matrix(nrow = nrow(sgmat) + nrow(rowdiffs), 
                ncol = ncol(sgmat) + ncol(coldiffs))

columns <- seq(1,ncol(nodes), by=2)
rows    <- seq(1,nrow(nodes),by=2)
nodes[rows, columns] <- sgmat
colnames(nodes)[columns] <- colnames(sgmat)
rownames(nodes)[columns] <- rownames(sgmat)

rows    <- seq(1, nrow(nodes),by=2)
columns <- seq(2, ncol(nodes), by=2)
nodes[rows, columns] <- coldiffs
colnames(nodes)[columns] <- colnames(coldiffs)

rows <- seq(2, nrow(nodes), by=2)
columns <- seq(1, ncol(nodes), by=2)
nodes[rows, columns] <- rowdiffs
rownames(nodes)[rows] <- rownames(rowdiffs)

connections <- matrix(
                0, 
                nrow = nrow(nodes) + ncol(nodes), 
                ncol = nrow(nodes) + ncol(nodes))
rownames(connections) <- c(nodes)

diagram::plotmat(x, c(2,2), box.type='square', box.size=0.05)
