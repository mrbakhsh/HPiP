    #' plotPPI
    #' @title Plot the Predicted PPI
    #' @param ppi A data.frame containing protein-protein interactions with edge
    #' score.
    #' @param edge.name A character string giving an edge attribute name.
    #' @param node.color The fill color of the node.
    #' @param edge.color The color of the edge.
    #' @param cex.node The size of the node.
    #' @param node.label.dist The distance of the label from the center of the
    #' node.
    #' @return A PPI plot.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom igraph get.adjacency
    #' @importFrom graphics par
    #' @importFrom igraph graph_from_data_frame
    #' @description Plot the predicted PPIs. This function uses the plot
    #' function of the \code{igraph}.
    #' @rdname plotPPI
    #' @export
    #' @examples
    #' df <- data.frame(
    #'     node1 = c("A", "B", "C", "D", "E"),
    #'     node2 = c("C", "E", "E", "E", "A"),
    #'     edge.scores = c(0.5, 0.4, 0.3, 0.2, 0.7)
    #' )
    #' plotPPI(df, edge.name = "edge.scores")
    plotPPI <- function(ppi, edge.name = "ensemble_score",
                        node.color = "grey",
                        edge.color = "orange",
                        cex.node = 4,
                        node.label.dist = 1.5) {
      m <- adjacency_matrix <- # create adjacency matrix
        get.adjacency(
          graph_from_data_frame(ppi, directed = FALSE),
          attr = edge.name,
          sparse = FALSE
        )

      # igraph
      g <-
        igraph::graph.adjacency(m,
          mode = "undirected", weighted = TRUE, diag = TRUE
        )
      graph <-
        igraph::simplify(g)

      par(mar = c(0, 0, 0, 0))
      print(plot(graph,
        edge.color = edge.color, vertex.color = node.color,
        vertex.label.color = "black", vertex.label.dist = node.label.dist,
        vertex.size = cex.node
      ))
    }
