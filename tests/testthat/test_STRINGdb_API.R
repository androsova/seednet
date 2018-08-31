library(seednet)
library(STRINGdb)

context("Internet connection")

string_db = STRINGdb$new(version="10", species=9606, score_threshold=500)

test_that("STRINGdb reference class is created", {
  graph <- STRING_graph(string_db)
  expect_is(graph, 'igraph')
})
