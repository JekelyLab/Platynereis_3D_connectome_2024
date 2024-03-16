# code to generate the connectome graph based on the CATMAID database
# Gaspar Jekely 2023

source("code/Natverse_functions_and_conn.R")

# get all synapses from CATMAID
all_syn_connectors <- catmaid_fetch(
  path = "11/connectors/",
  body = list(
    relation_type = "presynaptic_to",
    relation_type = "postsynaptic_to",
    with_partners = "true"
  )
)
length(all_syn_connectors$connectors)

connector_ids <- unlist(unname(sapply(all_syn_connectors$connectors, "[[", 1)))
connectivity <- catmaid_get_connectors(connector_ids, pid = 11)
length(connector_ids)
dim(connectivity)[1]

# convert connectivity table to tbl graph
conn.tb <- connectivity %>%
  select(pre, post) %>%
  as_tbl_graph()
conn.tb

# retrieve node names, simplify and filter the graph --------------------------------------------------------

# add weighted degree to nodes
conn.tb.str <- conn.tb %>%
  activate(nodes) %>%
  mutate(strength = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE,
    normalized = FALSE
  ))

# test if graph is simple (no multi-edges)
is_simple(as.igraph(conn.tb.str))

# add edge weight of 1
conn.tb.str <- set_edge_attr(conn.tb.str, "weight", value = 1)
conn.tb.str

# merge edges by converting to simple and summing edge weights
conn.tb.str.sum <- conn.tb.str %>%
  activate(edges) %>%
  convert(to_simple) %>%
  mutate(weight = map_dbl(.orig_data, ~ sum(.x$weight))) %>%
  select(-.orig_data, -.tidygraph_edge_index) %>%
  activate(nodes) %>%
  select(-.tidygraph_node_index)

conn.tb.str.sum
is_simple(as.igraph(conn.tb.str.sum))
components(conn.tb.str.sum)

# filter by node strength
conn.tb.str.sum.filt <- conn.tb.str.sum %>%
  activate(nodes) %>%
  filter(strength > 3)

conn.tb.str.sum.filt

# check connected components  ---------------------------

# check connected components, the connectome should only have one component
cl <- components(conn.tb.str.sum.filt)
# size of the largest subnetwork
length(which(cl$membership == 1))

# generate subgraph
conn.tb.str.sum.filt.cl <- subgraph(conn.tb.str.sum.filt, which(cl$membership == 1))

# check if graph is directed
igraph::is_directed(conn.tb.str.sum.filt.cl)

# shorten name, convert to tbl
conn.tb <- conn.tb.str.sum.filt.cl %>% as_tbl_graph()
# export graph for gephi to do force field layout -------------------------

# Serialise the graph to Gephi format
gexf_data <- rgexf::igraph.to.gexf(as.igraph(conn.tb))
# Write the network into a gexf (Gephi) file
write.gexf(gexf_data, output = "data/full_connectome_graph.gexf")

# Force Atlas tool in Gephi 0.1. The inertia was set to 0.1, repulsion strength was 50
# attraction strength was 10, maximum displacement was 5, gravity was 20, speed was 5 and the attraction distribution option was on/off/on.
# The 'auto stabilise function' was off. Nodes were scaled by weighted degree with min 5 max 25 node size.
# Towards the end of the clustering the 'adjust by sizes' option was also selected.
# To prevent node overlap, we then run the 'Noverlap' function

# read gexf file ----------------------------------------------------------

# read the gephi connectome file with nodes positioned by force field clustering
# at gexf export, positions were normalised (0,1)

# import gephi file (will be a directed graph)
conn_gexf <- rgexf::read.gexf("data/full_connectome_force_layout.gexf")

# get coordinates from imported gephi file
coords <- as_tibble(x = conn_gexf$nodesVizAtt$position$x) %>%
  mutate(x = value) %>%
  mutate(y = conn_gexf$nodesVizAtt$position$y) %>%
  mutate(skid = conn_gexf$nodes$label)

# sort by skid
x_coord <- as.list(arrange(coords, desc(skid)) %>%
  select(x))
y_coord <- as.list(arrange(coords, desc(skid)) %>%
  select(y))
conn.tb

# assign coordinates from imported gephi graph to original CATMAID graph
conn.tb.coord <- activate(conn.tb, nodes) %>%
  arrange(desc(name)) %>% # sort by skid to match coordinate list
  mutate(
    x = unlist(x_coord), # add coordinates
    y = unlist(y_coord)
  )

# check names, search for annotations --------------------------------------------------

# list skids
skids <- conn.tb.coord %>%
  activate(nodes) %>%
  select(name) %>%
  pull()

# get neuron names
names <- catmaid_get_neuronnames(skids, pid = 11)

conn.tb.coord.skids.names <- conn.tb.coord %>%
  mutate(names) %>%
  rename(skids = name)

conn.tb.coord.skids.names %>%
  as_tibble() %>%
  filter(names == "INipsidesc_sg1l")
  
#fix dupl name
conn.tb.coord.skids.names <- conn.tb.coord.skids.names %>%
  activate(nodes) %>%
  mutate(names = case_when(
    names == "INipsidesc_sg1l" & skids == "34371" ~ "INipsidesc_sg1l1",
    TRUE ~ names
    )
    )

# list neuron names
cell_names <- conn.tb.coord.skids.names %>%
  activate(nodes) %>%
  select(names) %>%
  pull() %>%
  unname()

# check duplicated names (there should be none)
cell_names[duplicated(cell_names)]
length(cell_names)

# iterate through skid list, get annotations for the first skid in every cell type
# check for the presence of the annotations and add the annotation to the type_of_cell list
type_of_cell <- list()
annot_to_search <- c(
  "Sensory neuron", "interneuron", "motorneuron",
  "effector", "fragmentum"
)

for (i in seq_along(skids)) {
  annot <- catmaid_get_annotations_for_skeletons(skids = skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if (sum(annot$annotation %in% annot_to_search[j]) == 1) {
      type_of_cell[i] <- annot_to_search[j]
      break()
    } else {
      type_of_cell[i] <- "other"
    }
  }
  print(i)
}

# annotations to search for
annot_to_search <- c(
  "episphere", "segment_0", "segment_1",
  "segment_2", "segment_3", "pygidium"
)
segment_of_cell <- list()
for (i in seq_along(skids)) {
  annot <- catmaid_get_annotations_for_skeletons(skids = skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if (sum(annot$annotation %in% annot_to_search[j]) == 1) {
      segment_of_cell[i] <- annot_to_search[j]
      break()
    } else {
      segment_of_cell[i] <- "fragment"
    }
  }
  print(i)
}

# annotations to search for
annot_to_search <- c("left_side", "right_side", "middle")
side_of_cell <- list()
for (i in seq_along(skids)) {
  annot <- catmaid_get_annotations_for_skeletons(skids = skids[i], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if (sum(annot$annotation %in% annot_to_search[j]) == 1) {
      side_of_cell[i] <- annot_to_search[j]
      break()
    } else {
      side_of_cell[i] <- "other"
    }
  }
  print(i)
}

# find all cell type annotations for neuronal and non-neuronal cell types
celltype_annot_skids <- conn.tb.coord.skids.names %>%
  activate(nodes) %>%
  select(skids) %>%
  pull() %>%
  catmaid_get_annotations_for_skeletons(pid = 11) %>%
  filter(grepl("^celltype\\d+$|^celltype_non_neuronal\\d+$", annotation)) %>%
  select(annotation, skid) %>%
  rename(skids = skid) %>%
  mutate(skids = as.character(skids))

# select only nodes
conn.tb.nodes <- as_tibble(conn.tb.coord.skids.names)

# join the cell type annotation table with the nodes tibble
joined_tibble <- left_join(conn.tb.nodes, celltype_annot_skids, by = "skids")
colnames(joined_tibble)

# select the cell type annotations only for all skids (including NA values)
celltype_annotations <- joined_tibble %>%
  select(annotation) %>%
  pull()

# add cell type annotation to the conn.tb
conn.tb.annot <- conn.tb.coord.skids.names %>%
  mutate(celltype_annotation = replace_na(celltype_annotations, "not_celltype"))

# add type of cell to node$group (can be used for visualisation)
conn.tb.annot.type <- conn.tb.annot %>%
  mutate(group = unlist(type_of_cell)) %>%
  mutate(class = unlist(type_of_cell)) %>%
  mutate(segment = unlist(segment_of_cell)) %>%
  mutate(side = unlist(side_of_cell))

# convert to igraph
conn.igraph <- as.igraph(conn.tb.annot.type)

# clustering with Leiden algorithm, run with ModularityVertexPartition and resolution parameter
partition <- leiden(conn.igraph,
  weights = E(conn.igraph)$weight,
  partition_type = "RBConfigurationVertexPartition",
  resolution_parameter = 1.0,
  n_iterations = -1, seed = 31
)
max(partition)

# define colors
col <- (c((brewer.pal(12, "Paired")), Okabe_Ito[1:8])[1:max(partition)])
length(col)
pie(rep(1, max(partition)), col = col)


# assign partition value and color to nodes
# If vertices are in the same cluster/partition, give them the same colour

conn.tb.annot.type.color <- conn.tb.annot.type %>%
  activate(nodes) %>%
  mutate(partition = partition) %>%
  mutate(color = col[partition])

conn.tb.annot.type.color %>%
  select(partition) %>%
  pull() %>%
  max()

# centrality measures --------------------------------------

conn.tb.annot.type.color.cent <- conn.tb.annot.type.color %>%
  activate(nodes) %>%
  mutate("betweenness" = centrality_betweenness(
    weights = E(conn.igraph)$weight,
    directed = TRUE
  )) %>%
  mutate("authority" = centrality_authority(weights = E(conn.igraph)$weight)) %>%
  mutate("pagerank" = centrality_pagerank(
    weights = E(conn.igraph)$weight,
    directed = TRUE
  )) %>%
  mutate("closeness" = centrality_closeness(
    weights = E(conn.igraph)$weight,
    mode = "all"
  )) %>%
  mutate("eigen" = centrality_eigen(directed = TRUE)) %>%
  mutate("hub" = centrality_hub(weights = E(conn.igraph)$weight)) %>%
  mutate("node_eccentricity_out" = node_eccentricity(mode = "out")) %>%
  mutate("node_eccentricity_in" = node_eccentricity(mode = "in")) %>%
  mutate("node_is_sink" = node_is_sink()) %>%
  mutate("node_is_source" = node_is_source()) %>%
  mutate("node_is_cut" = node_is_cut()) %>%
  mutate("local_triangles" = local_triangles())

conn.tb.annot.type.color.cent
# save conn.tb
saveRDS(conn.tb.annot.type.color.cent, "supplements/connectome_graph_tibble.rds")
conn.tb.annot.type.color.cent <- readRDS("supplements/connectome_graph_tibble.rds")

# VisNetwork conversion ---------------------------------------------------

# convert to visNetwork graph
conn_graph.visn <- toVisNetworkData(conn.tb.annot.type.color.cent)

## copy column "weight" to new column "value" in list "edges"
conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$strength

# save vis graph with annotations as R data file and txt file printed with dput() to get an exact copy
saveRDS(conn_graph.visn, "supplements/connectome_graph.rds")
writeLines(capture.output(dput(conn_graph.visn)), "supplements/connectome_graph.txt")
