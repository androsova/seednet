library(seednet)
library(STRINGdb)

context("Check outputs")

test_that("igraph from STRINGdb reference class is created", {
  string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
  graph <- STRING_graph(string_db)

  expect_is(graph, 'igraph')
})

test_that("NA test returns correct output", {
  expect_equal(anyNA(1:10), FALSE)
  expect_equal(anyNA(c(NA, 1:10, NULL)), TRUE)
})

test_that("adjacency matrix from igraph graph is created", {
  graph <- igraph::make_tree(10, 4, mode = "undirected")
  expect_error(convert_to_adjacency_matrix(graph, "combined_score"))

  igraph::E(graph)$combined_score <- 1:igraph::ecount(graph)
  expect_is(convert_to_adjacency_matrix(graph, "combined_score"), "matrix")
})

test_that("matrix from adjacency matrix is created", {
  graph <- igraph::make_tree(10, 4, mode = "undirected")
  igraph::E(graph)$combined_score <- 1:igraph::ecount(graph)
  adjacency_matrix <- convert_to_adjacency_matrix(graph, "combined_score")

  expect_is(degree_discounted_matrix(adjacency_matrix), "matrix")
})

test_that("degree-discounting creates igraph graph", {
  graph <- igraph::make_tree(10, 4, mode = "undirected")
  E(graph)$combined_score <- 1:ecount(graph)

  expect_is(degree_discounted_graph(graph, "combined_score"), "igraph")
})

test_that("Markov walk creates igraph graph", {
  string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
  graph <- STRING_graph(string_db)
  seeds_df <- data.frame(genes = Batten_disease_seeds)
  maped_seeds <- string_db$map(seeds_df, "genes")
  seeds <- sapply(strsplit(maped_seeds[!is.na(maped_seeds[,2]),2], '\\.'), "[[", 2)

  expect_is(Markov_walk(graph, seeds, 3, "combined_score"), "igraph")
})

test_that("First neighborhood creates igraph graph", {
  string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
  graph <- STRING_graph(string_db)
  seeds_df <- data.frame(genes = Batten_disease_seeds)
  maped_seeds <- string_db$map(seeds_df, "genes")
  seeds <- sapply(strsplit(maped_seeds[!is.na(maped_seeds[,2]),2], '\\.'), "[[", 2)

  expect_is(first_neighborhood(graph, seeds), "igraph")
})

test_that("Shortest path creates igraph graph", {
  string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)
  graph <- STRING_graph(string_db)
  seeds_df <- data.frame(genes = Batten_disease_seeds)
  maped_seeds <- string_db$map(seeds_df, "genes")
  seeds <- sapply(strsplit(maped_seeds[!is.na(maped_seeds[,2]),2], '\\.'), "[[", 2)

  expect_is(shortest_path(graph, seeds, 1000-E(graph)$combined_score), "igraph")
})

test_that("molecule list is non-empty character vector", {
  expect_is(Alzheimer_disease_seeds, "character")
  expect_equal(length(Alzheimer_disease_seeds)>0, TRUE)

  expect_is(Batten_disease_seeds, "character")
  expect_equal(length(Batten_disease_seeds)>0, TRUE)

  expect_is(Diabetes_mellitus_type_2_seeds, "character")
  expect_equal(length(Diabetes_mellitus_type_2_seeds)>0, TRUE)

  expect_is(Epilepsy_seeds, "character")
  expect_equal(length(Epilepsy_seeds)>0, TRUE)

  expect_is(Parkinson_disease_seeds, "character")
  expect_equal(length(Parkinson_disease_seeds)>0, TRUE)

  expect_is(Postoperative_delirium_seeds, "character")
  expect_equal(length(Postoperative_delirium_seeds)>0, TRUE)

  expect_is(Schizophrenia_seeds, "character")
  expect_equal(length(Schizophrenia_seeds)>0, TRUE)
})







