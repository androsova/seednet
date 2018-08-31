#' Load STRING protein-protein interactions
#'
#' This function creates an igraph object from STRING database \code{string_db}.
#' It removes species identifier from the proteins names and returns the
#' largest connected component. Note: you must have an Internet connection to
#' get a STRINGdb reference class.
#'
#' @param string_db STRINGdb reference class.
#'
#' @return igraph object of connected graph from STRING database.
#'
#' @import igraph
#' @import STRINGdb
#'
#' @examples
#' library(STRINGdb)
#' string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
#' graph <- STRING_graph(string_db)
#' print(graph)
#'
#' @export
#'
#' @references Franceschini, A et al. (2013). STRING v9.1: protein-protein
#' interaction networks, with increased coverage and integration. In:’Nucleic Acids Res.
#' 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29
STRING_graph = function(string_db) {
  requireNamespace("STRINGdb", quietly = TRUE)
  string_db_graph = string_db$load() %>%
    igraph::delete_edge_attr("neighborhood") %>%
    igraph::delete_edge_attr("neighborhood_transferred") %>%
    igraph::delete_edge_attr("fusion") %>%
    igraph::delete_edge_attr("cooccurence") %>%
    igraph::delete_edge_attr("homology") %>%
    igraph::delete_edge_attr("coexpression") %>%
    igraph::delete_edge_attr("coexpression_transferred") %>%
    igraph::delete_edge_attr("experiments") %>%
    igraph::delete_edge_attr("experiments_transferred") %>%
    igraph::delete_edge_attr("database") %>%
    igraph::delete_edge_attr("database_transferred") %>%
    igraph::delete_edge_attr("textmining") %>%
    igraph::delete_edge_attr("textmining_transferred")
  igraph::V(string_db_graph)$name = sapply(strsplit(igraph::V(string_db_graph)$name, '\\.'), "[[", 2)
  connected_components = igraph::clusters(string_db_graph)
  string_db_graph = igraph::induced.subgraph(string_db_graph,
                                     igraph::V(string_db_graph)[which(connected_components$membership
                                                              == which.max(connected_components$csize))])
    return(string_db_graph)
}

#' Check for NAs
#'
#' During R update this function might be missing from the base package.
#'
#' @param x R object to be tested.
#'
#' @return A logical value. TRUE in case input has NA, FALSE otherwise.
#'
#' @examples
#' no_NA <- 1:20
#' anyNA(no_NA)
#' present_NA <- rep(c(1:5, NA), 2)
#' anyNA(present_NA)
#'
#' @export
anyNA <- function(x) any(is.na(x))

#' Convert graph to adjacency matrix
#'
#' This function was modified from \code{as_adjacency_matrix()} function from
#' 'igraph' package to yield a faster computation for the large graphs. The graph
#' is represented as symmetric adjacency matrix, where rows and columns correspond
#' to graph vertices.
#'
#' @param graph An igraph graph with weights as edges attribute.
#' @param edge_attrib Name of the edge attribute with edges' weights.
#'
#' @return A symmetric matrix with number of rows and columns corresponding to the
#' number of vertices in the input graph.
#'
#' @note This function results in error if edge attribute with provided
#' \code{edge_attrib} name does not exist.
#'
#' @import igraph
#'
#' @examples
#' graph <- igraph::make_tree(10, 4, mode = "undirected")
#' igraph::E(graph)$combined_score <- 1:igraph::ecount(graph)
#' convert_to_adjacency_matrix(graph, "combined_score")
#'
#' @export
#'
#' @references Gábor Csárdi, Tamás Nepusz: The igraph software package for
#' complex network research. InterJournal Complex Systems, 1695, 2006.
convert_to_adjacency_matrix = function(graph, edge_attrib) {
  # This function was modified from the package 'igraph' to yield a faster computation
  el = igraph::as_edgelist(graph, names = FALSE)
  vc = igraph::vcount(graph)
  value = igraph::edge_attr(graph, name = edge_attrib)
  df = data.frame(Row = pmin(el[, 1], el[, 2]), Col = pmax(el[, 1], el[, 2]), Value = value)
  adj_matrix = matrix(0, nrow = vc, ncol = vc)
  adj_matrix[(df$Col - 1) * nrow(adj_matrix) + df$Row] = df$Value
  diag(adj_matrix) = 0
  adj_matrix = adj_matrix + t(adj_matrix)
  colnames(adj_matrix) <- rownames(adj_matrix) <- igraph::V(graph)$name
  return(adj_matrix)
}

#' Degree-discounting of adjacency matrix
#'
#' This function reduces the weight of interactions (edges) proportionally
#' to the vertex degree. Degree-discounting was adapted from the notion of
#' reduced adjacency matrix in Laplacian matrix theory. This function requires
#' \code{adj_matrix} and \code{alpha} parameter as an input. By default alpha
#' is set to -0.5 as suggested by Satuluri and Parthasarathy (2011).
#'
#' @param adj_matrix A symmetric numeric matrix.
#' @param alpha A numeric value for the alpha parameter.
#'
#' @return A degree-discounted symmetric matrix with number of rows and columns
#' as in input adjacency matrix. The values in the degree-discounted matrix are
#' normalized from 0 to 1.
#'
#' @note This function results in error if the input adjacency matrix is
#' not symmetric.
#'
#' @import igraph
#'
#' @examples
#' adjacency_matrix <- matrix(c(0, 0.85, 0.5, 0.5, 0, 0.9, 0.6,
#' 0, 0, 0.85, rep(0,8), 0.5, rep(0,8), 0.5, 0, 0, 0, 0.5, 0,
#' 0.8, 0.51, rep(0,4), 0.5, 0, 0.9, 0, 0, 0, 0.9, 0, 0, 0,
#' 0.9, rep(0,4), 0.6, 0, 0, 0.8, rep(0,4), 0.6, 0, 0, 0,
#' 0.51, rep(0,11), 0.6, 0, 0), ncol = 9)
#' colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- LETTERS[1:9]
#' degree_discounted_matrix(adjacency_matrix)
#'
#' @export
#' @seealso \url{https://en.wikipedia.org/wiki/Laplacian_matrix} for more
#' information on reduced adjacency matrix.
#' @references Satuluri V. and Parthasarathy S. (2011) Symmetrizations
#' for clustering directed graphs. In Proceedings of the 14th International
#' Conference on Extending Database Technology (EDBT/ICDT '11), ACM, New York, NY,
#' USA, 343-354.
degree_discounted_matrix = function(adj_matrix, alpha = 0.5) {
  diagonal_degree_matrix = diag(colSums(adj_matrix)^alpha)
  penalized_matrix = diagonal_degree_matrix %*% adj_matrix %*% diagonal_degree_matrix
  diag(penalized_matrix) = 0
  #Normalizing data
  normalized_matrix = (penalized_matrix-min(penalized_matrix))/(max(penalized_matrix)-min(penalized_matrix))
  rownames(normalized_matrix) <- colnames(normalized_matrix) <- colnames(adj_matrix)
  return(normalized_matrix)
}

#' Degree-discounting of graph
#'
#' This function reduces the weight of interactions (edges) proportionally
#' to the vertex degree. Degree-discounting was adapted from the notion of
#' reduced adjacency matrix in Laplacian matrix theory. This function requires
#' \code{graph} with weighted edges as an input.
#'
#' @param graph An igraph graph with weights as edges attribute.
#' @param edge_attrib Name of the edge attribute with edges' weights.
#'
#' @return A degree-discounted graph. The edge weights in the degree-discounted
#' graph are normalized from 0 to 1.
#'
#' @note This function results in error if edge attribute with provided
#' \code{edge_attrib} name does not exist.
#'
#' @import igraph
#'
#' @examples
#' graph <- igraph::make_tree(10, 4, mode = "undirected")
#' E(graph)$combined_score <- 1:igraph::ecount(graph)
#' degree_discounted_graph(graph, "combined_score")
#'
#' @export
#' @seealso \url{https://en.wikipedia.org/wiki/Laplacian_matrix} for reduced
#'adjacency matrix in Laplacian matrix theory
#' @references Satuluri V. and Parthasarathy S. (2011) Symmetrizations
#' for clustering directed graphs. In Proceedings of the 14th International
#' Conference on Extending Database Technology (EDBT/ICDT '11), ACM, New York, NY,
#' USA, 343-354.
degree_discounted_graph = function(graph, edge_attrib) {
  # Construction of adjacency matrix
  adj_matrix = convert_to_adjacency_matrix(graph, edge_attrib)

  # Construct reduced adjacency matrix
  penalized_matrix = degree_discounted_matrix(adj_matrix)

  # Degree-discounted graph with filtered edges
  degree_disc_graph = igraph::graph.adjacency(penalized_matrix, mode = "undirected", weighted = TRUE)
  return(degree_disc_graph)
}

#' Table with observed-predicted proteins for permutation setup
#'
#' This function creates a new table with all vertices from the input graph
#' and sets value to 1 for the input \code{seeds} in the "true" column,
#' while second column with be used for permutations.
#'
#' @param graph_vertices A vector with unique values corresponding to the input
#' graph vertex names.
#' @param seeds A vector with unique values corresponding to the input
#' seed names.
#'
#' @return A data.frame with graph vertices in rows and two columns for present
#' seeds and predicted during permutations.
#'
#' @note This function results in error if the input vectors have non-unique entries.
#'
#' @examples
#' library(STRINGdb)
#' string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
#' graph <- STRING_graph(string_db)
#' seeds_df <- data.frame(genes = Batten_disease_seeds)
#' maped_seeds <- string_db$map(seeds_df, "genes")
#' #Remove the species identifier from the protein names
#' seeds <- sapply(strsplit(maped_seeds[!is.na(maped_seeds[,2]),2], '\\.'), "[[", 2)
#' graph_vertices <- igraph::V(graph)$name
#' observed_pred <- observed_pred_table(graph_vertices, seeds)
#' table(observed_pred)
#'
#' @export
observed_pred_table = function(graph_vertices, seeds) {
  observ_pred_df = data.frame(true = rep(0, length(graph_vertices)),
                              predicted = rep(0, length(graph_vertices)))
  rownames(observ_pred_df) = graph_vertices
  observ_pred_df[which(graph_vertices %in% seeds),1] = 1
  return(observ_pred_df)
}

#' Greedy search graph
#'
#' This function creates an igraph graph from the input \code{graph} and
#' set of \code{seeds}. Each step it selects the strongest interacting nodes
#' defined by the edge attribute \code{attrib_name}. It makes number of
#' steps in the graph as defined by \code{steps} and retrieves unique
#' nodes not yet added in the previous step.
#'
#' @param graph An igraph graph with weights as edges attribute.
#' @param seeds A vector with unique values corresponding to the input
#' seed names.
#' @param steps An integer number indicating number of steps of the greedy search
#' @param attrib_name Name of the edge attribute with edges' weights.
#'
#' @return An igraph graph.
#'
#' @note This function results in error if the input vector \code{seeds}
#' has non-unique entries or if edge attribute with provided
#' \code{edge_attrib} name does not exist.
#'
#' @import igraph
#'
#' @examples
#' library(STRINGdb)
#' string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
#' graph <- STRING_graph(string_db)
#' seeds_df <- data.frame(genes = Batten_disease_seeds)
#' maped_seeds <- string_db$map(seeds_df, "genes")
#' #Remove the species identifier from the protein names
#' seeds <- sapply(strsplit(maped_seeds[!is.na(maped_seeds[,2]),2], '\\.'), "[[", 2)
#' greedy_search(graph, seeds, 3, "combined_score")
#'
#' @export
greedy_search = function(graph, seeds, steps, attrib_name) {
  added_nodes = seeds
  for (step in 1:as.numeric(steps)) {
    edges = as.numeric(igraph::E(graph)[from(seeds)])

    adj_edges = data.frame(ends(graph, edges),
                           weignt = igraph::edge_attr(graph, attrib_name, index = edges),
                           stringsAsFactors = FALSE)
    retrieved_nodes = sapply(seeds, function(seed){
      seed_edges = which(adj_edges==seed, arr.ind = TRUE)
      seed_edges_reversed = matrix(c(seed_edges[,1],
                                      ifelse(seed_edges[,2] == 1, 2, 1)),
                                      ncol=2)
      node_df = data.frame(node = adj_edges[seed_edges_reversed],
                           weight = adj_edges[seed_edges[,1],3],
                           stringsAsFactors = FALSE)
      node_df = node_df[which(!(node_df$node %in% added_nodes)),]
      node_df[which.max(node_df$weight),1]
    })
    names(retrieved_nodes) = NULL
    added_nodes = c(added_nodes, retrieved_nodes)

    seeds = retrieved_nodes
  }
  constructed_graph = igraph::induced.subgraph(graph, added_nodes)
  return(constructed_graph)
}

#' First-neighborhood graph
#'
#' This function creates an igraph graph from the input \code{graph} and
#' set of \code{seeds}. It selects all directly connected nodes to the input
#' \code{seeds} disregarding weights of the edges.
#'
#' @param graph An igraph graph with weights as edges attribute.
#' @param seeds A vector with unique values corresponding to the input
#' seed names.
#'
#' @return An igraph graph.
#'
#' @note This function results in error if the input vector \code{seeds}
#' has non-unique entries.
#'
#' @import igraph
#'
#' @examples
#' library(STRINGdb)
#' string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
#' graph <- STRING_graph(string_db)
#' seeds_df <- data.frame(genes = Batten_disease_seeds)
#' maped_seeds <- string_db$map(seeds_df, "genes")
#' #Remove the species identifier from the protein names
#' seeds <- sapply(strsplit(maped_seeds[!is.na(maped_seeds[,2]),2], '\\.'), "[[", 2)
#' first_neighborhood(graph, seeds)
#'
#' @export
first_neighborhood = function(graph, seeds) {
  constructed_graph = igraph::induced.subgraph(graph, igraph::V(subgraph.edges(graph, E(graph)[from(seeds)]))$name)
  return(constructed_graph)
}

#' Shortest path graph
#'
#' This function is a wraper for igraph function \code{get.shortest.paths()}. It creates
#' igraph graph from the input \code{graph} and specified set of vertices \code{seeds}.
#'
#' @param graph An igraph graph with weights as edges attribute.
#' @param seeds A vector with unique values corresponding to the input
#' seed names.
#' @param edge_weights A numeric vector with edge weights (optional, see documentation
#' on \code{igraph::get.shortest.paths()}).
#'
#' @return An igraph graph.
#'
#' @note This function results in error if the input vector \code{seeds}
#' has non-unique entries.
#'
#' @import igraph
#'
#' @examples
#' library(STRINGdb)
#' string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
#' graph <- STRING_graph(string_db)
#' seeds_df <- data.frame(genes = Batten_disease_seeds)
#' maped_seeds <- string_db$map(seeds_df, "genes")
#' # Remove the species identifier from the protein names
#' seeds <- sapply(strsplit(maped_seeds[!is.na(maped_seeds[,2]),2], '\\.'), "[[", 2)
#' #In STRING database, the higher combined score indicates the higher evidence (probability)
#' # of interaction between proteins, thus for shortest path (which selects the minimal
#' # distance between vertices based on the smallest sum of edges weights) we use reversed
#' # weights of the STRING graph.
#' shortest_path(graph, seeds, 1000-E(graph)$combined_score)
#'
#' @export
shortest_path = function(graph, seeds, edge_weights) {
  vertices_between_seeds = igraph::get.shortest.paths(graph, seeds, to = seeds, mode = "all", weights = edge_weights)
  constructed_graph = igraph::induced.subgraph(graph, unique(c(seeds, names(unlist(vertices_between_seeds[[1]])))))
  return(constructed_graph)
}
