    #' run_clustering
    #' @title Module Detection
    #' @param ppi A data.frame containing pathogen proteins in the
    #' first column,host proteins in the second column, and edge weight in the
    #' third column.
    #' @param method Module detection algorithms including:
    #' \itemize{
    #' \item \code{FC} - fast-greedy algorithm.
    #' \item \code{RW} - walktrap algorithm.
    #' \item \code{ML} - multi-level community algorithm.
    #' \item \code{clp} - label propagation algorithm.
    #' \item \code{MCL} - markov clustering.
    #' }
    #' @param expan Numeric value > 1 for the expansion parameter.
    #' See \code{\link[MCL]{mcl}} for more details.
    #' @param infla Numeric value > 0 for the inflation power coefficient.
    #' See \code{\link[MCL]{mcl}} for more details.
    #' @param iter An integer, the maximum number of iterations for the MCL.
    #' See \code{\link[MCL]{mcl}} for more details.
    #' @return A data.frame with the enrichment analysis results.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom igraph get.adjacency
    #' @importFrom igraph graph.adjacency
    #' @importFrom igraph fastgreedy.community
    #' @importFrom igraph cluster_walktrap
    #' @importFrom igraph cluster_louvain
    #' @importFrom igraph cluster_label_prop
    #' @importFrom igraph simplify
    #' @importFrom igraph membership
    #' @importFrom MCL mcl
    #' @description This function contains five module detection algorithms
    #' including fast-greedy algorithm
    #' (\code{FC}), walktrap algorithm (\code{RW}), multi-level community
    #' algorithm
    #' (\code{ML}), label propagation algorithm (\code{clp}), and markov
    #' clustering (\code{MCL}).
    #' @export
    run_clustering <- function(ppi, method = c("FC",
                                               "RW",
                                               "ML",
                                               "clp",
                                               "MCL"),
                               expan = 2,
                               infla = 5,
                               iter = 50) {

        method <- match.arg(method)

        attr <- colnames(ppi[,3])
        #create adjacency matrix
        m <-
            adjacency_matrix <- #create adjacency matrix
            get.adjacency(
                graph_from_data_frame(ppi, directed = FALSE),
                attr = attr,
                sparse = FALSE
            )
        #igraph
        g <-
            graph.adjacency(m,
                            mode = "undirected", weighted = TRUE, diag = T)
        graph<-
            simplify(g)

        if(method == "FC") {
            fc <-
                fastgreedy.community(graph,merges=TRUE,modularity=TRUE)
            membership<-
                membership(fc)
        }else if(method == "RW") {
        rw <-
            cluster_walktrap(graph,weights = E(graph)$weight,
                             merges=TRUE,modularity=TRUE)
        membership<-
            membership(rw)


        }else if(method == "ML") {
        ml <-
            cluster_louvain(graph,weights = E(graph)$weight)
        membership<-
            membership(ml)

        }else if(method == "clp") {
        clp <-
            cluster_label_prop(graph, weights = E(graph)$weight)
        membership<-
            membership(clp)

        } else if(method == "MCL") {

        mcl.clust <-
            mcl(m, addLoops = T,
                expansion = expan,
                inflation = infla,
                allow1 = FALSE,
                max.iter = iter, ESM = FALSE)
        membership<-
            mcl.clust$Cluster
        names(membership) <- colnames(m)

    }
        dfout <-
            cbind(names(membership),membership)
        colnames(dfout)[1] <-
            "member"
        return(dfout)
    }






