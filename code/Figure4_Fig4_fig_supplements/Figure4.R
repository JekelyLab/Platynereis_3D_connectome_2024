# code to generate Figure 4 (cell types graph) of the Platynereis connectome paper
# Gaspar Jekely 2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# network plot - read all skids based on cell type annotations ---------------------------

skids_celltypelist <- list()
n_neuron_types <- 202
n_non_neuron_types <- 92

# first we read all skids for Sensory cell types from 1-202
counter <- 0
for (i in c(1:n_neuron_types)) {
  annotation <- paste("annotation:^celltype", i, "$", sep = "")
  # retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "^Sensory neuron$")
  print(i)
  if (length(skids) > 0) {
    counter <- counter + 1
    skids_celltypelist[[counter]] <- skids
    print("Sensory neuron")
  }
}
length(skids_celltypelist)
counter

# we read all skids for interneuron cell types from 1-202
# we do not reset the counter and continue to fill the skids_celltypelist
for (i in c(1:n_neuron_types)) {
  annotation <- paste("annotation:^celltype", i, "$", sep = "")
  # retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "^interneuron$")
  print(i)
  if (length(skids) > 0) {
    counter <- counter + 1
    skids_celltypelist[[counter]] <- skids
    print("interneuron")
  }
}
length(skids_celltypelist)
counter

# we read all skids for motor neuron cell types from 1-202
for (i in c(1:n_neuron_types)) {
  annotation <- paste("annotation:^celltype", i, "$", sep = "")
  # retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "^motorneuron$")
  print(i)
  if (length(skids) > 0) {
    counter <- counter + 1
    skids_celltypelist[[counter]] <- skids
    print("motorneuron")
  }
}
length(skids_celltypelist)
counter
# the skids_celltypelist has the skids ordered by Sensory, inter, motorneuron

# we read the skids for all non-neuronal cell types from 1-92
# in the order 'ciliated cell', 'gland cell', 'pigment cell', 'muscle'
# continue to use the counter
for (i in c(1:n_non_neuron_types)) {
  annotation <- paste("annotation:^celltype_non_neuronal", i, "$", sep = "")
  # retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "^ciliated cell$")
  print(i)
  if (length(skids) > 0) {
    counter <- counter + 1
    skids_celltypelist[[counter]] <- skids
    print("ciliated cell")
  }
}
length(skids_celltypelist)
counter


for (i in c(1:n_non_neuron_types)) {
  annotation <- paste("annotation:^celltype_non_neuronal", i, "$", sep = "")
  # retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "^gland cell$")
  print(i)
  if (length(skids) > 0) {
    counter <- counter + 1
    skids_celltypelist[[counter]] <- skids
    print("gland cell")
  }
}
length(skids_celltypelist)
counter

for (i in c(1:n_non_neuron_types)) {
  annotation <- paste("annotation:^celltype_non_neuronal", i, "$", sep = "")
  # retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "^pigment cell$")
  print(i)
  if (length(skids) > 0) {
    counter <- counter + 1
    skids_celltypelist[[counter]] <- skids
    print("pigment cell")
  }
}
length(skids_celltypelist)
counter

for (i in c(1:n_non_neuron_types)) {
  annotation <- paste("annotation:^celltype_non_neuronal", i, "$", sep = "")
  # retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "^muscle$")
  print(i)
  if (length(skids) > 0) {
    counter <- counter + 1
    skids_celltypelist[[counter]] <- skids
    print("muscle")
  }
}
length(skids_celltypelist)
counter

# we ignore other non-neuronal cell types

# retrieve connectivity from CATMAID --------------------------------------

# define empty synapse matrix
synapse_matrix <- matrix(NA,
  nrow = (length(skids_celltypelist)),
  ncol = (length(skids_celltypelist))
)
dim(synapse_matrix)

# iterate through cell type list, retrieve and count connections and populate connectivity matrix
for (i in 1:(length(skids_celltypelist))) {
  for (j in 1:(length(skids_celltypelist))) {
    cellgroup_conn <- catmaid_get_connectors_between(
      pre = unlist(skids_celltypelist[i]),
      post = unlist(skids_celltypelist[j]), pid = 11, conn = conn_http1
    )
    N_synapses <- nrow(cellgroup_conn)
    if (length(N_synapses) == 0) N_synapses <- 0 # if there are no synapses, need to add 0
    synapse_matrix[i, j] <- N_synapses
  }
  print(i)
}

# get neuron names and rename matrix col and row ---------------------------------

# use the Catmaid neuron name before the "_" to parse the generic cell type name

# get the neuron name of the first skid in the skids list
celltype_names <- list()
skids_celltypelist
for (i in 1:length(skids_celltypelist)) {
  name <- catmaid_get_neuronnames(skids_celltypelist[[i]][1], pid = 11)
  name <- sub("_.*$", "", name)
  celltype_names[i] <- name
}
celltype_names
# check duplicated names
celltype_names[duplicated(celltype_names)]
# check length
length(celltype_names)

# assign column and row names
rownames(synapse_matrix) <- celltype_names
colnames(synapse_matrix) <- celltype_names

# search for annotations --------------------------------------------------

# these are the annotations to search for
annot_to_search <- c("Sensory neuron", "interneuron", "motorneuron", "effector")

# iterate through skid list, get annotations for the first skid in every cell type
# check for the presence of the annotations and add the annotation to the type_of_cell list
type_of_cell <- list()
for (i in seq_along(skids_celltypelist)) {
  annot <- catmaid_get_annotations_for_skeletons(skids = skids_celltypelist[[i]][1], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if (sum(annot[, "annotation"] %in% annot_to_search[j]) == 1) {
      type_of_cell[i] <- annot_to_search[j]
      break()
    } else {
      type_of_cell[i] <- "NULL"
    }
  }
}

# convert to tibble
syn_tb <- as_tibble(synapse_matrix) %>%
  mutate(presyn_cell = unlist(celltype_names)) %>%
  mutate(type_presyn = unlist(type_of_cell)) %>%
  pivot_longer(-c(presyn_cell, type_presyn),
    names_to = "postsyn_cell",
    values_to = "synapses"
  ) %>%
  filter(synapses > 0) %>%
  group_by(postsyn_cell) %>%
  mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))

# annotation and saving -----------------------------------------------------

# convert to igraph form
syn.igraph <- syn_tb %>%
  select(presyn_cell, postsyn_cell, synapses) %>%
  filter(synapses > 0) %>%
  as_tbl_graph()

# node names
node_names <- V(syn.igraph)$name

# create df of names vs type of cell
name_to_type <- data.frame(name = unlist(celltype_names), type = unlist(type_of_cell))
node_types <- lapply(node_names, function(x) name_to_type[which(name_to_type$name == x), ]$type)

# assign type to nodes in igraph node list
V(syn.igraph)$group <- unlist(node_types)

# calculate node weighted degree
weights <- strength(syn.igraph, vids = V(syn.igraph), mode = c("all"), loops = TRUE)
# assign weighted degree to nodes
V(syn.igraph)$value <- weights

# delete two unconnected vertices (MBmouth and INface)
syn.igraph <- delete_vertices(syn.igraph, "MBmouth")
syn.igraph <- delete_vertices(syn.igraph, "INface")


# convert back to tibble for saving
syn_tb <- syn.igraph %>%
  as_tbl_graph()

# determine source, sink and cut nodes

syn_tb <- syn_tb %>%
  activate(nodes) %>%
  mutate("node_eccentricity_out" = node_eccentricity(mode = "out")) %>%
  mutate("node_eccentricity_in" = node_eccentricity(mode = "in")) %>%
  mutate("node_is_sink" = node_is_sink()) %>%
  mutate("node_is_source" = node_is_source()) %>%
  mutate("node_is_cut" = node_is_cut())

syn_tb

syn_tb <- syn_tb %>%
  mutate("betweenness" = centrality_betweenness(
    weights = E(syn.igraph)$synapses,
    directed = TRUE, normalized = TRUE
  )) %>%
  mutate("closeness" = centrality_closeness(
    weights = E(syn.igraph)$synapses,
    mode = "out"
  )) %>%
  mutate("PageRank" = centrality_pagerank(
    weights = E(syn.igraph)$synapses,
    directed = TRUE
  ))

syn_tb <- syn_tb %>%
  mutate(group = case_when(
    group == "Sensory neuron" ~ "sensory neuron",
    TRUE ~ group
  )) %>%
  mutate(group = case_when(
    group == "motorneuron" ~ "motoneuron",
    TRUE ~ group
  ))

# save cell type connectivity data as source data
saveRDS(syn_tb, "source_data/Figure4_source_data1.rds")

# continue here if you want to load the network from data/ ----------------
syn_tb <- readRDS("source_data/Figure4_source_data1.rds")
syn.igraph <- as.igraph(syn_tb)

# Serialise the graph to Gephi format
gexf_data <- rgexf::igraph.to.gexf(syn.igraph)
# Write the network into a gexf (Gephi) file
write.gexf(gexf_data, output = "data/celltype_connectome_graph.gexf")


# Force Atlas tool in Gephi 0.9.6. The inertia was set to 0.1, repulsion strength was 1000
# attraction strength was 5, maximum displacement was 5, gravity was 1, speed was 5 and the attraction distribution option was off.
# The 'auto stabilise function' was off. Nodes were scaled by weighted degree with min 5 max 25 node size.
# Towards the end of the clustering the 'adjust by sizes' option was also selected.
# To prevent node overlap, we then run the 'Noverlap' and 'Label Adjust' functions

# read gexf file ----------------------------------------------------------

# read the gephi connectome file with nodes positioned by force field clustering
# at gexf export, positions were normalised (0,1)

# import gephi file (will be a directed graph)
conn_gexf <- rgexf::read.gexf("data/celltype_connectome_force_layout.gexf")

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

# assign coordinates from imported gephi graph to original CATMAID graph
syn_tb <- activate(syn_tb, nodes) %>%
  arrange(desc(name)) %>% # sort by skid to match coordinate list
  mutate(
    x = unlist(x_coord), # add coordinates
    y = unlist(y_coord)
  )

# graph visualisation -----------------------------------------------------

# convert to visNet form
syn.vis <- syn_tb %>%
  toVisNetworkData()
syn.vis

# assign sqrt of number of synapses to edge 'value'
syn.vis$edges$value <- sqrt(syn.vis$edges$synapses)
syn.vis$nodes$value

#coordinates as matrix
coords <- matrix(c(syn.vis$nodes$x, syn.vis$nodes$y), ncol = 2)
coords

visNet <- visNetwork(syn.vis$nodes, syn.vis$edges) %>%
  visIgraphLayout(type = "full", layoutMatrix = coords) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0.1),
    scaling = list(min = 0.1, max = 20),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.8, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = syn.vis$nodes$color, border = "black"),
    scaling = list(min = 5, max = 50),
    font = list(color = "black", size = 40),
  ) %>%
  visGroups(
    groupname = "sensory neuron", color = "#E69F00", shape = "dot",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "interneuron", shape = "square",
    opacity = 1, color = "#CC79A7"
  ) %>%
  visGroups(
    groupname = "motoneuron", shape = "diamond",
    opacity = 1, color = "#0072B2"
  ) %>%
  visGroups(
    groupname = "effector", shape = "triangle",
    opacity = 1, color = "#CCCCCC"
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 3000, height = 1500) %>%
  addFontAwesome()

visNet

# save as html
saveNetwork(visNet, "pictures/Figure4_celltype_network.html",
  selfcontained = TRUE
)
# save from web browser
webshot::webshot(
  url = "pictures/Figure4_celltype_network.html",
  file = "pictures/Figure4_celltype_network.png",
  vwidth = 3000, vheight = 1500, # define the size of the browser window
  cliprect = c(120, 130, 2800, 1360), zoom = 2, delay = 2
)


# plot node source/sink  --------------------------------------------------

all_cells <- syn_tb %>%
  activate(nodes) %>%
  as_tibble() %>%
  filter(group %in% c(
    "sensory neuron", "interneuron",
    "motoneuron", "effector"
  )) %>%
  mutate(with_soma = TRUE) %>%
  ggplot() +
  geom_histogram(aes(x = with_soma, fill = group), stat = "count") +
  labs(x = "all", y = "# of nodes", fill = "cell class", title = "") +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.background = element_blank(),
    axis.ticks = element_line(size = 0.2)
  ) +
  scale_fill_manual(
    values = list(
      `sensory neuron` = "#E69F00",
      interneuron = "#CC79A7",
      motoneuron = "#0072B2",
      effector = "grey40"
    )
  )
all_cells

cut_nodes <- syn_tb %>%
  activate(nodes) %>%
  as_tibble() %>%
  filter(group %in% c(
    "sensory neuron", "interneuron",
    "motoneuron", "effector"
  )) %>%
  mutate(with_soma = TRUE) %>%
  filter(node_is_cut == TRUE) %>%
  ggplot() +
  geom_histogram(aes(x = node_is_cut, fill = group), stat = "count") +
  labs(x = "cut", y = "", fill = "cell class", title = "") +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_line(size = 0.2),
    panel.background = element_blank(),
    axis.ticks = element_line(size = 0.2)
  ) +
  scale_fill_manual(
    values = list(
      interneuron = "#CC79A7",
      motoneuron = "#0072B2",
      `sensory neuron` = "#E69F00",
      effector = "grey40"
    )
  ) +
  guides(fill = "none")
cut_nodes

source <- syn_tb %>%
  activate(nodes) %>%
  as_tibble() %>%
  filter(group %in% c(
    "sensory neuron", "interneuron",
    "motoneuron", "effector"
  )) %>%
  mutate(with_soma = TRUE) %>%
  filter(node_is_source == TRUE) %>%
  ggplot() +
  geom_histogram(aes(x = node_is_source, fill = group), stat = "count") +
  labs(x = "source", y = "", fill = "cell class", title = "") +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_line(size = 0.2),
    panel.background = element_blank(),
    axis.ticks = element_line(size = 0.2)
  ) +
  scale_fill_manual(
    values = list(
      interneuron = "#CC79A7",
      motoneuron = "#0072B2",
      `sensory neuron` = "#E69F00",
      effector = "grey40"
    )
  ) +
  guides(fill = "none")
source

sink <- syn_tb %>%
  activate(nodes) %>%
  as_tibble() %>%
  filter(group %in% c(
    "sensory neuron", "interneuron",
    "motoneuron", "effector"
  )) %>%
  mutate(with_soma = TRUE) %>%
  filter(node_is_sink == TRUE) %>%
  ggplot() +
  geom_histogram(aes(x = node_is_sink, fill = group), stat = "count") +
  labs(x = "sink", y = "", fill = "cell class", title = "") +
  theme_minimal_hgrid() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_line(size = 0.2),
    panel.background = element_blank(),
    axis.ticks = element_line(size = 0.2)
  ) +
  scale_fill_manual(
    values = list(
      `sensory neuron` = "#E69F00",
      interneuron = "#CC79A7",
      motoneuron = "#0072B2",
      effector = "grey40"
    )
  ) +
  guides(fill = "none")
sink

sink_cut_source <- all_cells + source + sink + cut_nodes + plot_layout(
  guides = "collect",
  nrow = 1
) +
  plot_annotation(tag_levels = NULL) &
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9),
    legend.key.size = unit(3, "mm"),
    panel.grid.minor = element_blank()
  )

sink_cut_source

# plot graph coloured by partition -----------------------------------------


# clustering with Leiden algorithm, run with ModularityVertexPartition and resolution parameter
partition <- leiden(syn.igraph,
  partition_type = "ModularityVertexPartition",
  resolution_parameter = 1, weights = E(syn.igraph)$synapses,
  n_iterations = 1000
)
max(partition)

# color by partition
# define some colors
col <- c(brewer.pal(12, "Paired"), sample(blues, (max(partition)), replace = TRUE))

# If vertices are in the same cluster, give them the same color
syn.vis$nodes$color <- col[partition]
syn.vis$nodes$partition <- partition

visNet_mod <- visNetwork(syn.vis$nodes, syn.vis$edges) %>%
  visIgraphLayout(type = "full", layoutMatrix = coords) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0.1),
    scaling = list(min = 0.1, max = 20),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.8, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = syn.vis$nodes$color, border = "black"),
    scaling = list(min = 5, max = 50),
    font = list(color = "black", size = 35),
  ) %>%
  visGroups(
    groupname = "sensory neuron", shape = "dot",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "interneuron", shape = "square",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "motoneuron", shape = "diamond",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "effector", shape = "triangle",
    opacity = 1
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 3000,
    height = 1500) %>%
  addFontAwesome()

visNet_mod

# save as html
saveNetwork(visNet_mod, "pictures/visNetwork_celltype_connectome2.html")
# save from web browser
webshot::webshot(
  url = "pictures/visNetwork_celltype_connectome2.html",
  file = "pictures/visNetwork_celltype_connectome2.png",
  vwidth = 3000, vheight = 1500, # define the size of the browser window
  cliprect = c(120, 130, 2800, 1360), zoom = 2, delay = 2
)


# create labels ------------------------------------------------------------


label_celltype_distr <- tibble(
  data = c(1, 2, 3, 4),
  type = c(
    "sensory neuron", "interneuron",
    "motor neuron", "muscle"
  )
)

label1 <- label_celltype_distr %>%
  ggplot(aes(type, data, fill = type)) +
  geom_tile() +
  scale_fill_manual(
    values = c("#E69F00", "#CC79A7", "#0072B2", "#BD0026"),
    breaks = c(
      "sensory neuron", "interneuron",
      "motor neuron", "muscle"
    )
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

# crop and save
label1 <- label1 &
  theme(plot.margin = margin(t = -5, r = -55, b = -190, l = -80))

ggsave("pictures/Fig4_label_celltype_distr.png",
  width = 1400, height = 76,
  units = c("px"), label1
)

label_graph <- tibble(
  data = c(1, 2, 3, 4),
  type = c(
    "sensory neuron", "interneuron",
    "motor neuron", "effector"
  )
)
label2 <- label_graph %>%
  ggplot(aes(type, data, fill = type)) +
  geom_tile() +
  scale_fill_manual(
    values = c("#E69F00", "#CC79A7", "#0072B2", "grey40"),
    breaks = c(
      "sensory neuron", "interneuron",
      "motor neuron", "effector"
    )
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

# crop and save
label2 <- label2 &
  theme(plot.margin = margin(t = -5, r = -191, b = -190, l = -222))

ggsave("pictures/Fig4_label_graph.png",
  width = 1225, height = 76,
  units = c("px"), label2
)


# graph statistics --------------------------------------------------

# calculate node centrality measures

graph_size <- list()

my_theme <- theme_minimal() +
  theme(axis.title.x = element_text(size = 11)) +
  theme(axis.title.y = element_text(size = 11)) +
  theme(axis.text = element_text(size = 8)) +
  theme(legend.text = element_text(size = 7)) +
  theme(legend.title = element_text(size = 9)) +
  theme(
    legend.key.size = unit(3, "mm"),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.background = element_blank(),
    axis.ticks = element_line(linewidth = 0.2)
  )


# network size with different edge thresholds
for (i in 1:100) {
  # create subgraph with edges above a threshold
  conn_graph_filt <- delete_edges(syn_tb, which(E(syn.igraph)$synapses < i))
  # check size of largest subgraph
  cl <- components(conn_graph_filt)
  print(max(cl$csize))
  graph_size[i] <- max(cl$csize)
}

# plot synapse cut-off graph
syn_cutoff <- as.data.frame(graph_size) %>%
  pivot_longer(everything()) %>%
  mutate(syn_cutoff = c(1:100)) %>%
  rename("nodes" = value) %>%
  select(syn_cutoff, nodes) %>%
  ggplot(aes(x = syn_cutoff, y = nodes)) +
  geom_col() +
  labs(x = "synapse cutoff", y = "# of nodes", title = " ") +
  theme_minimal_hgrid() +
  my_theme
syn_cutoff


# plot edge weight distribution
edge_weight <- as.data.frame(E(syn.igraph)$synapses) %>%
  ggplot(aes(x = E(syn.igraph)$synapses)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10(breaks = c(1, 2, 10, 100)) +
  labs(x = "edge weight", y = "# of edges", title = " ") +
  theme_minimal_hgrid() +
  my_theme
edge_weight


# plot node weighted degree distribution
syn_tb
w_degree <- syn_tb %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "effector"),
    aes(value),
    fill = "grey40", alpha = 0.6, bins = 48,
    color = "grey", linewidth = 0.1
  ) +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "sensory neuron"),
    aes(value),
    fill = "#E69F00", alpha = 1, bins = 54,
    color = "grey", linewidth = 0.1
  )  +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "interneuron"),
    aes(value),
    fill = "#CC79A7", alpha = 0.5, bins = 52,
    color = "grey", linewidth = 0.1
  )  +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "motoneuron"),
    aes(value),
    fill = "#0072B2", alpha = 0.6, bins = 50,
    color = "grey", linewidth = 0.1
  ) +
  labs(x = "weighted degree", y = "# of nodes", title = " ") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 25, 80)) +
  theme_minimal_hgrid() +
  my_theme
w_degree


# plot node betweenness centrality
between <- syn_tb %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "effector"),
    aes(betweenness),
    fill = "grey40", alpha = 0.6, bins = 48,
    color = "grey", linewidth = 0.1
  ) +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "sensory neuron"),
    aes(betweenness),
    fill = "#E69F00", alpha = 1, bins = 54,
    color = "grey", linewidth = 0.1
  )  +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "interneuron"),
    aes(betweenness),
    fill = "#CC79A7", alpha = 0.5, bins = 52,
    color = "grey", linewidth = 0.1
  )  +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "motoneuron"),
    aes(betweenness),
    fill = "#0072B2", alpha = 0.6, bins = 50,
    color = "grey", linewidth = 0.1
  ) +
  labs(x = "betweenness", y = "# of nodes", title = " ") +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) +
  theme_minimal_hgrid() +
  my_theme
between


# plot node PageRank centrality

page <- syn_tb %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "effector"),
    aes(PageRank),
    fill = "grey40", alpha = 0.6, bins = 48,
    color = "grey", linewidth = 0.1
  ) +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "sensory neuron"),
    aes(PageRank),
    fill = "#E69F00", alpha = 1, bins = 54,
    color = "grey", linewidth = 0.1
  )  +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "interneuron"),
    aes(PageRank),
    fill = "#CC79A7", alpha = 0.5, bins = 52,
    color = "grey", linewidth = 0.1
  )  +
  geom_histogram(
    data = syn_tb %>% as_tibble() %>% filter(group == "motoneuron"),
    aes(PageRank),
    fill = "#0072B2", alpha = 0.6, bins = 50,
    color = "grey", linewidth = 0.1
  ) +
  labs(x = "PageRank", y = "# of nodes", title = " ") +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.02)) +
  theme_minimal_hgrid() +
  my_theme
page

# plot node weighted degree against PageRank

degree_rank <- syn_tb %>%
  activate(nodes) %>%
  as_tibble() %>%
  ggplot(aes(
    x = value, y = PageRank, color = group, shape = group, alpha = group,
    size = betweenness
  )) +
  geom_point() +
  scale_color_manual(
    values = list(
      interneuron = "#CC79A7",
      motoneuron = "#0072B2",
      `sensory neuron` = "#E69F00",
      effector = "grey40"
    )
  ) +
  scale_shape_manual(
    values = list(
      interneuron = 16,
      motoneuron = 15,
      `sensory neuron` = 17,
      effector = 18
    )
  ) +
  scale_alpha_manual(
    values = c(0.6, 0.6, 0.8, 1)
  ) +
  scale_size_area(max_size = 5) +
  scale_x_continuous(limits = c(-6, 85)) +
  labs(
    x = "weighted degee", color = "cell class",
    alpha = "cell class", shape = "cell class"
  ) +
  geom_text(
    aes(label = name),
    size = 2.5, alpha = 1, nudge_x = 0.003, nudge_y = 0.0005,
    check_overlap = T, col = "black"
  ) +
  my_theme

degree_rank

# distances plot ---------------------


syn_tb_filt <- syn_tb %>%
  activate(edges) %>%
  filter(synapses>3) %>%
  activate(nodes)

# distances from SN to effectors ------------
d_eff <- distances(syn_tb_filt, v = syn_tb_filt %>%
                     filter(group == "sensory neuron") %>% pull(name),
                   to = syn_tb_filt %>%
                     filter(group == "effector") %>% pull(name),
                   mode = "out")

SN_dist_to_eff_plot <- as.data.frame(d_eff) %>% rownames_to_column("input") %>%
  pivot_longer(
    -input, names_to = "effector", values_to = "distance"
  ) %>%
  group_by(input) %>%
  mutate(shortest_path_to_effector = min(distance)) %>%
  select(input, shortest_path_to_effector) %>%
  unique() %>%
  ggplot(aes(shortest_path_to_effector)) +
  geom_histogram(bins = 10) +
  my_theme +
  labs(x = "distance to effector", y = "# sensory neurons")
SN_dist_to_eff_plot


# table on cell type numbers -------------------------------------

nNodes <- igraph::gorder(as.igraph(syn_tb))
nEdges <- igraph::gsize(as.igraph(syn_tb))
nSyn <- syn_tb %>%
  activate(edges) %>%
  select(synapses) %>%
  pull() %>%
  sum()

trans <- transitivity(as.igraph(syn_tb), type = "localaverageundirected")
mean_dist <- igraph::mean_distance(as.igraph(syn_tb))
eDens <- edge_density(as.igraph(syn_tb), loops = FALSE)

# the numbers come from running the Suppl_Figure_celltype_conn_matrix.R that collects cell types based on annotations

table <- plot_ly(
  type = "table",
  columnwidth = c(10, 2.5),
  columnorder = c(0, 1),
  header = list(
    values = c("cell type or network property", "#"),
    align = c("center", "center"),
    line = list(width = 1, color = "black"),
    fill = list(color = c("#CCCCCC", "#AAAAAA")),
    font = list(family = "Arial", size = 14, color = "black")
  ),
  cells = list(
    values = rbind(
      c(
        "total neuronal cell types", "sensory neuron types",
        "interneuron types", "motor neuron types",
        "total non-neuronal cell types", "muscle types",
        "locomotor ciliated-cell types", "gland-cell types",
        "pigment-cell types", "other non-neuronal types",
        "nodes", "edges",
        "in-graph synapse",
        "density", "mean distance", "local avg. transitivity"
      ),
      c(
        "202", "84", "84", "34", "92", "53",
        "6", "7", "9", "17",
        nNodes, nEdges, nSyn, format(round(eDens, 4)),
        format(round(mean_dist, 2)), format(round(trans, 2))
      )
    ),
    align = c("center", "center"),
    line = list(color = "black", width = 0.3),
    font = list(family = "Arial", size = 12, color = c("black"))
  )
)

table
saveNetwork(table, "pictures/celltype_stats_table.html")
webshot::webshot(
  url = "pictures/celltype_stats_table.html",
  file = "pictures/celltype_stats_table.png",
  vwidth = 300, vheight = 450, # define the size of the browser window
  cliprect = c(20, 50, 250, 360), zoom = 10
)


# assemble figure ---------------------------------------------------------

#check max of network
network_temp <- read_rds("source_data/Figure4_source_data1.rds")
network_temp %>%
  activate(edges) %>%
  as_tibble %>%
  select(synapses) %>%
  max()

  

panel_glands_network <- ggdraw() + 
  draw_image(readPNG("pictures/visNetwork_glands.png")) +
  draw_label("MNgland-head circuit", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11) + 
  draw_label("# of synapses", x = 0.94, y = 0.9, size = 8, hjust = 1) +
  draw_label("1", x = 0.87, y = 0.82, size = 8, hjust = 1) + 
  draw_label("31", x = 0.87, y = 0.76, size = 8, hjust = 1) +
  draw_line(x = c(0.88, 0.93), y = c(0.82, 0.82), size = 0.3, color = 'grey') +
  draw_line(x = c(0.88, 0.93), y = c(0.76, 0.76), size = 2, color = 'grey')


img_network <- readPNG("pictures/Figure4_celltype_network.png")
img_table <- readPNG("pictures/celltype_stats_table.png")
panel_network <- ggdraw() + draw_image(img_network, scale = 1.05) +
    draw_label("# of synapses", x = 1, y = 0.9, size = 8, hjust = 1) +
  draw_label("1", x = 0.94, y = 0.87, size = 8, hjust = 1) + 
  draw_label("521", x = 0.94, y = 0.84, size = 8, hjust = 1) +
  draw_line(x = c(0.95, 1), y = c(0.87, 0.87), size = 0.2, color = 'grey') +
  draw_line(x = c(0.95, 1), y = c(0.84, 0.84), size = 1.6, color = 'grey')

panel_table <- ggdraw() + draw_image(img_table, scale = 1.3)

layout <- "
AAAAAAAAAAAAAAAAAAAAAA
BBBBBB################
BBBBBB#CCDDDEEEFFFFFFF
GGGGGGGGHHHHHHHIIIIIII
"

Figure4 <- panel_network + 
  panel_table  +  edge_weight + syn_cutoff + SN_dist_to_eff_plot + sink_cut_source +
  w_degree + page + degree_rank  +  
  plot_layout(
    design = layout, 
    heights = c(2, 0.25, 0.75, 1)
    )+
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "", "", "", "G", "H", "I"))) &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure4.png",
  limitsize = FALSE,
  units = c("px"), Figure4, width = 4000, height = 4000, bg = "white"
)

ggsave("Figures/Figure4.pdf",
  limitsize = FALSE,
  units = c("px"), Figure4, width = 4000, height = 4000
)
