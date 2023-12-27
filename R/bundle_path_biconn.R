#' edge_bundle_path_biconn
#' the function adds column `component`, which indicates which component the edge belongs to.
#' @export
edge_bundle_path_biconn <- function(graph, xy, weightFactor=2, mode="out",
    max_distortion=2, segments=20, bundle_strength=1) {
    el <- as_edgelist(graph, names=FALSE)
    E(graph)$dist <- apply(el, 1, function(row) {
        X1=xy[row[1], 1]
        Y1=xy[row[1], 2]
        X2=xy[row[2], 1]
        Y2=xy[row[2], 2]
        sqrt((X1-X2)**2 + (Y1-Y2)**2)
    })^weightFactor

    bc <- igraph::biconnected_components(graph)

    paths <- do.call(rbind, lapply(seq_along(bc$component_edges), function(ci) {
        xy2 <- xy[sort(as.integer(bc$components[[ci]])),]
        co <- bc$component_edges[[ci]]
        g <- subgraph.edges(graph, co)
        w <- E(g)$dist
        el <- as_edgelist(g, names=FALSE)

        m <- length(co)
        lock <- rep(FALSE, m)
        skip <- rep(FALSE, m)
        sortedEdges <- order(w, decreasing=TRUE)
        cpoints <- vector("list", m)
        

        for (se in sortedEdges) {

            s <- el[se, 1]
            t <- el[se, 2]
            
            cpoints[[se]] <- xy2[c(s, t), ]
            if (lock[se]) {
                next()
            }
            skip[se] <- TRUE
            g1 <- igraph::delete.edges(g, which(skip))
            sp_verts <- suppressWarnings(igraph::shortest_paths(g1, s, t,
                weights = w[!skip], mode=mode)$vpath[[1]])
            
            if (length(sp_verts) < 2) {
                skip[se] <- FALSE
                next
            }
            sp_len <- path_length(sp_verts, xy) ## compatible with the original impl.
            if (sp_len >= max_distortion * w[se]) {
                skip[se] <- FALSE
                next
            }

            lock[igraph::get.edge.ids(g, rep(as.integer(sp_verts), each = 2)[-c(1, 2 * length(sp_verts))])] <- TRUE
            cpoints[[se]] <- xy2[sp_verts, ]

        }
        cpoints <- lapply(cpoints, subdivide, bs = bundle_strength)
        cpoints_bezier <- lapply(cpoints, approximateBezier, n = segments)
        idx <- seq(0, 1, length.out = segments)
        data_bundle <- as.data.frame(cbind(
            do.call("rbind", cpoints_bezier),
            rep(idx, m),
            rep(co[sortedEdges], each = segments)
        ))
        names(data_bundle) <- c("x", "y", "index", "group")
        data_bundle[["component"]] <- ci
        data_bundle
    }))
    paths
}