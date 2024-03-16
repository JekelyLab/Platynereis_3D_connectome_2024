# R code to generate the connectome graph images in the 3d Platynereis connectome paper
# Gaspar Jekely March 2022

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# plot graph with coordinates from gephi ----------------------------------

# read graph
connect.tb <- readRDS("supplements/connectome_graph_tibble.rds")
conn_graph.visn <- readRDS("supplements/connectome_graph.rds")
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol = 2)

partitions <- connect.tb %>%
  select(partition) %>%
  pull() %>%
  max()

coords_rotated <- autoimage::rotate(
  coords,
  -pi / 3,
  pivot = c(0, 0)
)


#flip along x axis
coords_rotated[,1] = coords_rotated[,1] * -1

visNet <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 0.05, max = 35),
    color = list(inherit = TRUE, opacity = 0.5),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.1,
    color = list(border = "black"),
    opacity = 1,
    scaling = list(min = 10, max = 50),
    font = list(size = 20)
  ) %>%
  visOptions(
    highlightNearest = list(
      enabled = TRUE, degree = 1,
      algorithm = "hierarchical", labelOnly = FALSE
    ),
    width = 2500, height = 2500, autoResize = FALSE
  ) %>%
  visInteraction(
    dragNodes = TRUE, dragView = TRUE,
    zoomView = TRUE, hover = TRUE,
    multiselect = TRUE
  ) %>%
  addFontAwesome() %>%
  visEvents(selectNode = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['label', 'label_long'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, label : label_info[0].label_long, label_long : label_info[0].label});
            }") %>%
  visEvents(blurNode = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['label', 'label_long'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, label : label_info[0].label_long, label_long : label_info[0].label});
  }")

# save as html
saveNetwork(visNet, "pictures/Full_connectome_modules.html", selfcontained = TRUE)

webshot::webshot(
  url = "pictures/Full_connectome_modules.html",
  file = "pictures/Full_connectome_modules_webshot.png",
  vwidth = 2550, vheight = 2550, # define the size of the browser window
  cliprect = c(90, 105, 2420, 2400), zoom = 2
)


# plot neurons by colours matching network pic -----------------------------

#background plot
plot_neuron_modules <- function(neurons, color, alpha){
  plot3d(neurons, soma=TRUE, lwd=2,
         add=T, alpha=alpha,
         col=color)
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
}

plot_background_ventral()
plot3d(scalebar_50um_ventral, lwd = 3, color = "black")
open_device <- cur3d()

# plot modules  in the colour that was used in the network plot
for (i in 1:partitions) {
  print(i)
  skids <- connect.tb %>%
  filter(partition == i) %>%
  pull(skids)
  
  color <- connect.tb %>%
  filter(partition == i) %>%
  pull(color) %>%
  unique()
  
  neurons = nlapply(read.neurons.catmaid(skids, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  
  set3d(open_device)
  plot_neuron_modules(
    neurons, 
    color, 
    runif(length(skids), min = 0.7, max = 0.9)
  )
  
  plot_background_ventral()
  plot_neuron_modules(
    neurons, 
    color, 
    runif(length(skids), min = 0.7, max = 0.9)
  )
  
  text <- ifelse(sum(as.numeric(skids %in% "1333316")) == 1, "mech_l", 
          ifelse(sum(as.numeric(skids %in% "1333515")) == 1, "mech_r",
          ifelse(sum(as.numeric(skids %in% "12115")) == 1, "eye",
          ifelse(sum(as.numeric(skids %in% "1298426")) == 1, "ciliomotor",
          ifelse(sum(as.numeric(skids %in% "228289")) == 1, "pigment",
          ifelse(sum(as.numeric(skids %in% "251790")) == 1, "NS",
          ifelse(sum(as.numeric(skids %in% "77507")) == 1, "MB",
          ifelse(sum(as.numeric(skids %in% "1421839")) == 1, "MUSsp_l",      
          ifelse(sum(as.numeric(skids %in% "955074")) == 1, "MUSsp_r",      
          ifelse(sum(as.numeric(skids %in% "1732111")) == 1, "postural",
          ifelse(sum(as.numeric(skids %in% "1307524")) == 1, "MUSbow_r",
          ifelse(sum(as.numeric(skids %in% "893804")) == 1, "MUS2_r",
          ifelse(sum(as.numeric(skids %in% "1398955")) == 1, "MUS3_r",
                 i)))))))))))))
      
  rgl.snapshot(paste("pictures/conn_module_", text, ".png", sep = ""))
  close3d()
  
}

set3d(open_device)
rgl.snapshot("pictures/conn_module_all.png")

close3d()
# add partition name based on occurrence of key cells --------------

key_skids <- c(
  "251790", "77507", 
  "1333316", "1333515", "12115",
  "1298426", "228289", "1732111", "1421839",
  "955074", "1307524", "893804",
  "1398955")
partition_names <- c(
  "NS", "MB", 
  "mech (l)", "mech (r)", "visual",
  "ciliomotor", "pigment", 
  "postural", "muscle 1",
  "muscle 2","muscle 3","muscle 4",
  "muscle 5")

connect.tb <- connect.tb %>%
  mutate(partition_name = "")

for(i in 1:13){
part <- connect.tb %>%
  filter(skids  == key_skids[i]) %>%
  pull(partition)
connect.tb <- connect.tb %>%
  mutate(partition_name = ifelse(
    partition == part, partition_names[i], partition_name)
    )
}

connect.tb %>%
  select(partition_name) %>%
  pull()

# cell class name adjustment --------------
connect.tb <- connect.tb %>%
  mutate(class = ifelse(class == "Sensory neuron", "sensory neuron", 
                        ifelse(class == "motorneuron", "motoneuron", class)))
  
connect.tb

# plot showing which cell types are in which module ----------

segmental_colors <- brewer.pal(6, "Paired")

connect.tb %>%
  activate(nodes) %>%
  as_tibble() %>%
  group_by(partition_name) %>% #
  filter(
    class %in% c(
      "sensory neuron", "interneuron",
      "motoneuron", "effector"
      )
  ) %>%
  ggplot() +
  geom_histogram(aes(x = partition_name, fill = class), stat = "count") +
  labs(x = "module", y = "# of cells", fill = "cell class", title = "") +
  scale_fill_manual(
    values = list(
      interneuron = "#CC79A7",
      motoneuron = "#0072B2",
      `sensory neuron` = "#E69F00",
      effector = "grey40"
    )
  ) +
  scale_x_discrete(limits = partition_names) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(3, "mm"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave("pictures/connectome_partitions_plot.png", width = 1600, height = 1200, units = "px", bg = "white")

# print graph with cut nodes highlighted ----------------------------------

# overwrite group value (partition) with node_is_cut value (for colouring)
conn_graph.visn$nodes$group <- as.character(conn_graph.visn$nodes$node_is_cut)
# remove colour info (which takes precedence over group colour)
conn_graph.visn$nodes$color <- c()

visNet_cut <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 1, max = 15),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = conn_graph.visn$nodes$color, border = "black"),
    opacity = 1,
    scaling = list(min = 15, max = 50),
    font = list(size = 20)
  ) %>%
  visOptions(
    highlightNearest = list(
      enabled = TRUE, degree = 1,
      algorithm = "hierarchical", labelOnly = FALSE
    ),
    width = 2500, height = 2500, autoResize = FALSE
  ) %>%
  visGroups(
    groupname = "TRUE", color = "#D55E00", shape = "dot",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "FALSE", shape = "diamond",
    opacity = 0.5, color = "#56B4E9"
  )

# save as html
saveNetwork(visNet_cut, "pictures/Full_connectome_cut_nodes.html", selfcontained = TRUE)
webshot::webshot(
  url = "pictures/Full_connectome_cut_nodes.html",
  file = "pictures/Full_connectome_cut_nodes.png",
  vwidth = 2550, vheight = 2550, # define the size of the browser window
  cliprect = c(70, 75, 2420, 2400), zoom = 1
)


# print graph with sink nodes highlighted ----------------------------------

# overwrite group value (partition) with node_is_sink value (for colouring)
conn_graph.visn$nodes$group <- as.character(conn_graph.visn$nodes$node_is_sink)
# remove colour info (which takes precedence over group colour)
conn_graph.visn$nodes$color <- c()

visNet_sink <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 1, max = 15),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = conn_graph.visn$nodes$color, border = "black"),
    opacity = 1,
    scaling = list(min = 15, max = 50),
    font = list(size = 20)
  ) %>%
  visOptions(
    highlightNearest = list(
      enabled = TRUE, degree = 1,
      algorithm = "hierarchical", labelOnly = FALSE
    ),
    width = 2500, height = 2500, autoResize = FALSE
  ) %>%
  visGroups(
    groupname = "TRUE", color = "#D55E00", shape = "dot",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "FALSE", shape = "diamond",
    opacity = 0.5, color = "#56B4E9"
  )

# save as html
saveNetwork(visNet_sink, "pictures/Full_connectome_sink_nodes.html", selfcontained = TRUE)
webshot::webshot(
  url = "pictures/Full_connectome_sink_nodes.html",
  file = "pictures/Full_connectome_sink_nodes.png",
  vwidth = 2550, vheight = 2550, # define the size of the browser window
  cliprect = c(70, 75, 2420, 2400), zoom = 1
)


# print graph with source nodes highlighted ----------------------------------

# overwrite group value (partition) with node_is_source value (for colouring)
conn_graph.visn$nodes$group <- as.character(conn_graph.visn$nodes$node_is_source)
# remove colour info (which takes precedence over group colour)
conn_graph.visn$nodes$color <- c()

visNet_source <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 1, max = 15),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = conn_graph.visn$nodes$color, border = "black"),
    opacity = 1,
    scaling = list(min = 15, max = 50),
    font = list(size = 20)
  ) %>%
  visOptions(
    highlightNearest = list(
      enabled = TRUE, degree = 1,
      algorithm = "hierarchical", labelOnly = FALSE
    ),
    width = 2500, height = 2500, autoResize = FALSE
  ) %>%
  visGroups(
    groupname = "TRUE", color = "#D55E00", shape = "dot",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "FALSE", shape = "diamond",
    opacity = 0.5, color = "#56B4E9"
  )

# save as html
saveNetwork(visNet_source, "pictures/Full_connectome_source_nodes.html", selfcontained = TRUE)
webshot::webshot(
  url = "pictures/Full_connectome_source_nodes.html",
  file = "pictures/Full_connectome_source_nodes.png",
  vwidth = 2550, vheight = 2550, # define the size of the browser window
  cliprect = c(70, 75, 2420, 2400), zoom = 1
)

# print graph with cells coloured by cell class ----------------------------------

{
  # overwrite group value (partition) with segment value (for colouring)
  conn_graph.visn$nodes$group <- as.character(conn_graph.visn$nodes$class)
  # remove colour info (which takes precedence over group colour)
  conn_graph.visn$nodes$color <- c()
  
  conn_graph.visn$nodes$group <- gsub("Sensory neuron", "SN", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group <- gsub("interneuron", "IN", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group <- gsub("motorneuron", "MN", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group <- gsub("fragmentum", "fragment", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group
  
  visNet_class <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
    visEdges(
      smooth = list(type = "curvedCW", roundness = 0),
      scaling = list(min = 1, max = 15),
      color = list(inherit = FALSE, opacity = 0.1),
      arrows = list(to = list(
        enabled = TRUE,
        scaleFactor = 0.5, type = "arrow"
      ))
    ) %>%
    visNodes(
      borderWidth = 0.3,
      color = list(background = conn_graph.visn$nodes$group, border = "black"),
      opacity = 1,
      font = list(size = 0)
    ) %>%
    visOptions(
      highlightNearest = list(
        enabled = TRUE, degree = 1,
        algorithm = "hierarchical", labelOnly = FALSE
      ),
      width = 700, height = 700, autoResize = FALSE
    ) %>%
    visGroups(
      groupname = "SN", color = "#E69F00", shape = "dot",
      size = 25, opacity = 1
    ) %>%
    visGroups(
      groupname = "IN", shape = "dot",
      opacity = 1, size = 25, color = "#CC79A7"
    ) %>%
    visGroups(
      groupname = "MN", shape = "dot",
      opacity = 1, size = 25, color = "#0072B2"
    ) %>%
    visGroups(
      groupname = "effector", shape = "triangle",
      opacity = 0.4, size = 25, color = "#000000"
    ) %>%
    visGroups(
      groupname = "fragment", shape = "dot",
      opacity = 1, size = 5, color = "#aaaaaa"
    ) %>%
    visGroups(
      groupname = "other", shape = "dot",
      opacity = 1, size = 10, color = "#aaaaaa"
    ) %>%
    visLegend(
      useGroups = TRUE,
      width = 0.1,
      ncol = 1,
      position = "left"
    )
  
  # save as html
  saveNetwork(visNet_class, "pictures/Full_connectome_cell_classes.html", selfcontained = TRUE)
  webshot::webshot(
    url = "pictures/Full_connectome_cell_classes.html",
    file = "pictures/Full_connectome_cell_classes.png",
    vwidth = 1000, vheight = 1000, # define the size of the browser window
    cliprect = c(60, 80, 740, 680), zoom = 10, delay = 5
  )
}


# table on connectome stats -----------------------------------------------

# overwrite group value with class value
conn_graph.visn$nodes$group <- as.character(conn_graph.visn$nodes$class)

SN <- as_tibble(conn_graph.visn$nodes) %>%
  filter(group == "Sensory neuron") %>%
  select(group)
nSN <- length(unlist(SN))

IN <- as_tibble(conn_graph.visn$nodes) %>%
  filter(group == "interneuron") %>%
  select(group)
nIN <- length(unlist(IN))

MN <- as_tibble(conn_graph.visn$nodes) %>%
  filter(group == "motorneuron") %>%
  select(group)
nMN <- length(unlist(MN))

effector <- as_tibble(conn_graph.visn$nodes) %>%
  filter(group == "effector") %>%
  select(group)
nEff <- length(unlist(effector))

fragment <- as_tibble(conn_graph.visn$nodes) %>%
  filter(group == "fragmentum") %>%
  select(group)
nFrag <- length(unlist(fragment))

nSyn <- sum(conn_graph.visn$edges$weight)


nNodes <- igraph::gorder(as.igraph(connect.tb))
nEdges <- igraph::gsize(as.igraph(connect.tb))

trans <- transitivity(as.igraph(connect.tb), type = "localaverageundirected")
mean_dist <- igraph::mean_distance(as.igraph(connect.tb))
eDens <- edge_density(as.igraph(connect.tb), loops = FALSE)

table <- plot_ly(
  type = "table",
  columnwidth = c(6, 2.5),
  columnorder = c(0, 1),
  header = list(
    values = c("graph statistics", "#"),
    align = c("center", "center"),
    line = list(width = 1, color = "black"),
    fill = list(color = c("#EEEEEE", "#DDDDDD")),
    font = list(family = "Arial", size = 14, color = "black")
  ),
  cells = list(
    values = rbind(
      c(
        "nodes", "edges",
        "in-graph synapse",
        "density", "mean distance", "local avg. transitivity"
      ),
      c(
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
saveNetwork(table, "pictures/connectome_stats_cells_table.html")
webshot::webshot(
  url = "pictures/connectome_stats_cells_table.html",
  file = "pictures/connectome_stats_cells_table.png",
  vwidth = 240, vheight = 240, # define the size of the browser window
  cliprect = c(23, 48, 230, 154), zoom = 10
)

# graph statistics --------------------------------------------------


# create subgraph with edges above a threshold
conn.igraph <- as.igraph(connect.tb)
graph_size <- list()
for (i in 1:20) {
  conn_graph_filt <- delete.edges(conn.igraph, which(E(conn.igraph)$weight < i))
  # check size of largest subgraph
  cl <- components(conn_graph_filt)
  print(max(cl$csize))
  graph_size[i] <- max(cl$csize)
}

# plot edge cutoff
as.data.frame(graph_size) %>%
  pivot_longer(everything()) %>%
  mutate(syn_cutoff = c(1:20)) %>%
  rename("nodes" = value) %>%
  select(syn_cutoff, nodes) %>%
  ggplot(aes(x = syn_cutoff, y = nodes)) +
  geom_col() +
  labs(x = "synapse cutoff", y = "# of nodes", title = " ") +
  theme_minimal_hgrid() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("pictures/connectome_size_syn_cutoff.png",
  width = 1200, height = 1200, units = "px", bg = "white"
)



# plot edge weight distribution

as_tibble(E(conn.igraph)$weight) %>%
  ggplot(aes(x = E(conn.igraph)$weight)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10(breaks = c(1, 2, 5, 30)) +
  labs(x = "edge weight", y = "# of edges", title = " ") +
  theme_minimal_hgrid() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
ggsave("pictures/connectome_weight_distr.png",
  width = 1200, height = 1200, units = "px", bg = "white"
)


# plot node weighted degree distribution

as.data.frame(V(conn.igraph)$strength) %>%
  ggplot(aes(x = V(conn.igraph)$strength)) +
  geom_histogram(binwidth = 0.05) +
  scale_x_log10(breaks = c(1, 2, 10, 100, 500)) +
  labs(x = "weighted degree", y = "# of nodes", title = " ") +
  theme_minimal_hgrid() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("pictures/connectome_degree_distr.png",
  width = 1200, height = 1200, units = "px", bg = "white"
)


# plot node weighted degree against eccentricity

connect.tb %>%
  as_tibble() %>%
  filter(class %in% c(
    "sensory neuron", "interneuron",
    "motoneuron", "effector"
  )) %>%
  ggplot(aes(
    x = strength,
    y = unlist((node_eccentricity_in - node_eccentricity_out) /
      (node_eccentricity_in + node_eccentricity_out)),
    color = class,
    shape = class,
    size = pagerank
  )) +
  geom_point() +
  scale_x_log10() +
  labs(
    x = "weighted degree", y = "eccentricity",
    title = " ", color = "cell class",
    alpha = "cell class", shape = "cell class"
  ) +
  theme_minimal_hgrid() +
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
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(3, "mm"),
    panel.grid.minor = element_blank()
  )

ggsave("pictures/connectome_weight_vs_eccentr.png",
  width = 1600, height = 1200, units = "px", bg = "white"
)


# plot node source/sink

all_nodes <- as_tibble(connect.tb) %>%
  filter(side == "left_side" | side == "right_side" | side == "middle") %>%
  filter(class %in% c(
    "sensory neuron", "interneuron",
    "motoneuron", "effector"
  )) %>%
  mutate(with_soma = TRUE) %>%
  ggplot() +
  geom_histogram(aes(x = with_soma, fill = class), stat = "count") +
  labs(x = "all", y = "# of nodes", fill = "cell class", title = "") +
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
  )
all_nodes

cut <- as_tibble(connect.tb) %>%
  group_by(partition) %>%   
  filter(
    class %in% c(
      "sensory neuron", "interneuron",
      "motoneuron", "effector"
    )
    ) %>%
  filter(node_is_cut == TRUE)  %>%
  ggplot() +
  geom_histogram(aes(x = node_is_cut, fill = class), stat = "count") +
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
  )
cut

source <- as_tibble(connect.tb) %>%
  group_by(partition) %>% 
  filter(
    class %in% c(
      "sensory neuron", "interneuron",
      "motoneuron", "effector"
    )) %>%
  filter(node_is_source == TRUE) %>%
  ggplot() +
  geom_histogram(aes(x = node_is_source, fill = class), stat = "count") +
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
    ), guide = "none"
  )
source

sink <- as_tibble(connect.tb) %>%
  group_by(partition) %>% # Specify group indicator
  filter(
    class %in% c(
      "sensory neuron", "interneuron",
      "motoneuron", "effector"
    )) %>%
  filter(node_is_sink == TRUE) %>%
  ggplot() +
  geom_histogram(aes(x = node_is_sink, fill = class), stat = "count") +
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
      interneuron = "#CC79A7",
      motoneuron = "#0072B2",
      `sensory neuron` = "#E69F00",
      effector = "grey40"
    )
  )
sink

sink_cut_source <- all_nodes + source + sink + cut + plot_layout(
  guides = "collect",
  nrow = 1
) &
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(3, "mm"),
    panel.grid.minor = element_blank()
  )
sink_cut_source

ggsave(
  "pictures/connectome_sink_cut_source.png",
  limitsize = FALSE,
  units = c("px"), sink_cut_source, width = 1600, height = 1200
)


# assemble figure ---------------------------------------------------------

# read with image_read (magick package) and rotate
img_conn <- readPNG("pictures/Full_connectome_modules_webshot.png")

panel_mod1 <- ggdraw() + draw_image(readPNG("pictures/conn_module_mech_l.png")) +
  draw_label("mechanosensory (l)", x = 0.5, y = 0.95, size = 9)
panel_mod2 <- ggdraw() + draw_image(readPNG("pictures/conn_module_MUSsp_r.png")) +
  draw_label("muscle motor", x = 0.5, y = 0.05, size = 9)
panel_mod3 <- ggdraw() + draw_image(readPNG("pictures/conn_module_MB.png")) +
  draw_label("MB and central brain", x = 0.5, y = 0.95, size = 9)
panel_mod4 <- ggdraw() + draw_image(readPNG("pictures/conn_module_NS.png")) +
  draw_label("anterior NS", x = 0.5, y = 0.95, size = 9)
panel_mod5 <- ggdraw() + draw_image(readPNG("pictures/conn_module_mech_r.png")) +
  draw_label("mechanosensory (r)", x = 0.5, y = 0.95, size = 9)
panel_mod6 <- ggdraw() + draw_image(readPNG("pictures/conn_module_MUSsp_l.png")) +
  draw_label("muscle motor", x = 0.5, y = 0.05, size = 9)
panel_mod7 <- ggdraw() + draw_image(readPNG("pictures/conn_module_ciliomotor.png")) +
  draw_label("ciliomotor", x = 0.6, y = 0.95, size = 9)
panel_mod_post <- ggdraw() + draw_image(readPNG("pictures/conn_module_postural.png")) +
  draw_label("postural control", x = 0.5, y = 0.95, size = 9)

panel_mod8 <- ggdraw() + draw_image(readPNG("pictures/conn_module_MUSbow_r.png")) +
  draw_label("muscle motor", x = 0.5, y = 0.05, size = 9)
panel_mod9 <- ggdraw() + draw_image(readPNG("pictures/conn_module_eye.png")) +
  draw_label("visual", x = 0.5, y = 0.95, size = 9)
panel_mod11 <- ggdraw() + draw_image(readPNG("pictures/conn_module_pigment.png")) +
  draw_label("pigment motor", x = 0.5, y = 0.95, size = 9)
panel_mod10 <- ggdraw() + draw_image(readPNG("pictures/conn_module_MUS2_r.png")) +
  draw_label("muscle motor", x = 0.5, y = 0.05, size = 9)
panel_mod12 <- ggdraw() + draw_image(readPNG("pictures/conn_module_MUS3_r.png")) +
  draw_label("muscle motor", x = 0.5, y = 0.05, size = 9)

panel_all <- ggdraw() + draw_image(readPNG("pictures/conn_module_all.png")) +
  draw_label("all modules", x = 0.5, y = 0.95, size = 9) +
  draw_label(expression(paste("50 ", mu, " m")), x = 0.84, y = 0.08, size = 10)
  

panel_conn <- ggdraw() + draw_image(img_conn) +
  draw_label("anterior NS", x = 0.63, y = 0.93, size = 12) +
  draw_label("MB and central brain", x = 0.25, y = 0.65, size = 12) +
  draw_label("visual", x = 0.18, y = 0.78, size = 12) +
  draw_label("ciliomotor", x = 0.52, y = 0.56, size = 12) +
  draw_label("postural control", x = 0.25, y = 0.5, size = 12) +
  draw_label("pigment motor", x = 0.75, y = 0.6, size = 12) +
  draw_label("mechanosensory (l)", x = 0.65, y = 0.35, size = 12) +
  draw_label("mechanosensory (r)", x = 0.8, y = 0.45, size = 12) +
  draw_label("muscle motor", x = 0.2, y = 0.15, size = 12) +
  draw_label("muscle motor", x = 0.3, y = 0.32, size = 12) +
  draw_label("muscle motor", x = 0.5, y = 0.15, size = 12)

# define layout with textual representation
layout_A <- "
aABcC
#EEE#
DEEEF
#EEE#
GEEEH
#EEE#
IJKLM
"

Figure2 <- panel_all + panel_mod9 + panel_mod4 + panel_mod3 + panel_mod_post +
  panel_mod7 + panel_conn +
  panel_mod11 +
  panel_mod5 + panel_mod1 +
  panel_mod6 + panel_mod8 + panel_mod2 + panel_mod12 + panel_mod10 +
  plot_layout(design = layout_A,  heights = c(8, 0.7, 8, 0.7, 8, 0.7, 8))

ggsave("Figures/Figure2.png",
  limitsize = FALSE,
  units = c("px"), Figure2, width = 3000, height = 3400, bg = "white"
)

ggsave("Figures/Figure2.pdf",
       limitsize = FALSE,
       units = c("px"), Figure2, width = 3000, height = 3400
)

# suppl fig --------------

img_syn_cut <- readPNG("pictures/connectome_size_syn_cutoff.png")
img_weight <- readPNG("pictures/connectome_weight_distr.png")
img_degree <- readPNG("pictures/connectome_degree_distr.png")
img_between <- readPNG("pictures/connectome_betweenness.png")

panel_degree <- ggdraw() + draw_image(img_degree)
panel_weight <- ggdraw() + draw_image(img_weight)

panel_partitions <- ggdraw() + draw_image(readPNG("pictures/connectome_partitions_plot.png"))
panel_syn_cut <- ggdraw() + draw_image(img_syn_cut)
# read with image_read (magick package) and rotate
panel_source <- ggdraw() + draw_image(readPNG("pictures/Full_connectome_source_nodes.png")) +
  draw_label("source nodes", x = 0.6, y = 0.97, size = 9)
panel_sink <- ggdraw() + draw_image(readPNG("pictures/Full_connectome_sink_nodes.png")) +
  draw_label("sink nodes", x = 0.67, y = 0.97, size = 9)
panel_cut <- ggdraw() + draw_image(readPNG("pictures/Full_connectome_cut_nodes.png")) +
  draw_label("cut nodes", x = 0.6, y = 0.97, size = 9)

panel_cell_classes <- ggdraw() +
  draw_image(readPNG("pictures/Full_connectome_cell_classes.png")) +
  draw_label("cell class", x = 0.6, y = 0.97, size = 9)

panel_sink_cut_source <- ggdraw() +
  draw_image(readPNG("pictures/connectome_sink_cut_source.png"))
panel_weight_vs_ecc <- ggdraw() +
  draw_image(readPNG("pictures/connectome_weight_vs_eccentr.png"))

layout <- "
AABBEE
AABBEE
AABBFF
CCDDFF
CCDDGG
CCDDGG
######
HHIIJJ
"

Fig_conn <- panel_cell_classes + panel_source + panel_sink + panel_cut +
  panel_syn_cut + panel_weight + panel_degree +
  panel_weight_vs_ecc + panel_sink_cut_source + panel_partitions +
  plot_layout(
    design = layout, 
    heights = c(1, 1, 1, 1, 1, 1, 0.05, 2.4)
    ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))


ggsave("Figures/Figure2_fig_suppl1.pdf",
  limitsize = FALSE,
  units = c("px"), Fig_conn, width = 2400, height = 2200
)

ggsave("Figures/Figure2_fig_suppl1.png",
  limitsize = FALSE,
  units = c("px"), Fig_conn, width = 2400, height = 2200, bg = "white"
)

# save source data ------------

#save annotated connectome tibble
saveRDS(connect.tb, "source_data/Figure2_source_data1.rds")

