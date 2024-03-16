# statistics for paper - to be sourced and inserted into the text
# Gaspar Jekely 2023

library(tidygraph)
library(dplyr)
library(tibble)
library(igraph)
library(catmaid)
source("code/CATMAID_connection.R")

# statistics from CATMAID ------------

frag_all_annot <- as_tibble(catmaid_get_annotations_for_skeletons(
  "^fragmentum$", pid = 11
))

#all skeletons
with_soma_all_annot <- as_tibble(catmaid_get_annotations_for_skeletons(
  "^with_soma$", pid = 11
))

with_soma_node_count <- catmaid_get_node_count("^with_soma$", pid = 11)
fragment_node_count <- catmaid_get_node_count("^fragmentum$", pid = 11)

with_soma_postsyn <- catmaid_query_connected("^with_soma$", pid = 11)$incoming
with_soma_presyn <- catmaid_query_connected("^with_soma$", pid = 11)$outgoing
fragment_postsyn <- catmaid_query_connected("^fragmentum$", pid = 11)$incoming
fragment_presyn <- catmaid_query_connected("^fragmentum$", pid = 11)$outgoing

with_soma_postsyn_tb <- as_tibble(with_soma_postsyn) %>%
  mutate(connection_type = "postsynapse")
with_soma_presyn_tb <- as_tibble(with_soma_presyn) %>%
  mutate(connection_type = "presynapse")

with_soma_conn <- full_join(with_soma_postsyn_tb, with_soma_presyn_tb)
with_soma_conn <- with_soma_conn %>%
  mutate(skeleton = "with_soma")

fragment_postsyn_tb <- as_tibble(fragment_postsyn) %>%
  mutate(connection_type = "postsynapse")
fragment_presyn_tb <- as_tibble(fragment_presyn) %>%
  mutate(connection_type = "presynapse")

fragment_conn <- full_join(fragment_postsyn_tb, fragment_presyn_tb)
fragment_conn <- fragment_conn %>%
  mutate(skeleton = "fragment")

with_soma_and_fragment_conn <- full_join(with_soma_conn, fragment_conn)

write.table(frag_all_annot, "data/frag_all_annot.txt", sep = "\t")
write.table(with_soma_all_annot, "data/with_soma_all_annot.txt", sep = "\t")
write.table(with_soma_node_count, "data/with_soma_node_count.txt", sep = "\t")
write.table(fragment_node_count, "data/fragment_node_count.txt", sep = "\t")
write.table(with_soma_and_fragment_conn, "data/with_soma_and_fragment_conn.txt", sep = "\t")

frag_all_annot <- read.table("data/frag_all_annot.txt", sep = "\t")
with_soma_all_annot <- read.table("data/with_soma_all_annot.txt", sep = "\t")
with_soma_node_count <- read.table("data/with_soma_node_count.txt", sep = "\t")
fragment_node_count <-  read.table("data/fragment_node_count.txt", sep = "\t")
with_soma_and_fragment_conn <- read.table("data/with_soma_and_fragment_conn.txt", sep = "\t")
