# R code to generate Fig10 fig suppl3 of the 3d Platynereis connectome paper
# Gaspar Jekely 2024

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# plot graph with coordinates from gephi ----------------------------------

# read graph
connect.tb <- readRDS("supplements/connectome_graph_tibble.rds")

filter_region_class <- function(region_to_filter, class_to_filter){
 connect.tb %>%
    as_tibble() %>%
    filter(segment %in% region_to_filter) %>%
    filter(class %in% class_to_filter) %>%
    select(skids) %>%
    pull()
}

# head skids --------------------
head_SN_skids <- filter_region_class("episphere", "Sensory neuron")
head_IN_skids <- filter_region_class("episphere", "interneuron")
head_MN_skids <- filter_region_class("episphere", "motorneuron")
head_effector_skids <- filter_region_class("episphere", "effector")

# right side skids --------------------
torso_SN_skids <- filter_region_class(c("segment_0", "segment_1", "segment_2", "segment_3", "pygidium"), "Sensory neuron")
torso_IN_skids <- filter_region_class(c("segment_0", "segment_1", "segment_2", "segment_3", "pygidium"), "interneuron")
torso_MN_skids <- filter_region_class(c("segment_0", "segment_1", "segment_2", "segment_3", "pygidium"), "motorneuron")
torso_effector_skids <- filter_region_class(c("segment_0", "segment_1", "segment_2", "segment_3", "pygidium"), "effector")

cell_groups <- list(
  head_SN_skids, head_IN_skids, head_MN_skids, head_effector_skids,
  torso_SN_skids, torso_IN_skids, torso_MN_skids, torso_effector_skids
)

# connectivity function
number_of_syns <- function(pre, post){
  syn_number <- catmaid_get_connectors_between(
  unlist(pre), 
  unlist(post),
  pid = 11) 
  ifelse(length(syn_number) == 0, 0, syn_number %>%
           as_tibble %>%
           select(connector_id) %>%
           pull() %>% length())
}

head_trunk_matrix <- sapply(cell_groups, function(x)
  sapply(cell_groups, function(y)
  number_of_syns(x, y)
    )
  )

names <- c(
  "SN head", "IN head", "MN head", "effector head",
  "SN trunk", "IN trunk", "MN trunk", "effector trunk"
  )

colnames(head_trunk_matrix) <- names
rownames(head_trunk_matrix) <- names

# plot matrix -------------
as.data.frame(head_trunk_matrix) %>%
  rownames_to_column(var = "presynaptic") %>%
  pivot_longer(-presynaptic, names_to = "postsynaptic",
               values_to = "synapses") %>%
  ggplot(aes(presynaptic, postsynaptic, fill = sqrt(synapses))) +
  geom_tile() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 0, vjust = 0.5, size = 14
    ),
    axis.text.y = element_text(
      angle = 0, hjust = 1, vjust = 0.5, size = 14
    ),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.ticks = element_line(linewidth = 0.1)
  ) +
  scale_fill_gradientn(
    colours = c("white", "#0072B2",  "#D55E00")
  ) +
  labs(
    x = "postsynaptic cell groups",
    y = "presynaptic cell groups",
    title = " ", fill = "sqrt\nsynapses"
  ) +
  scale_x_discrete(
    limits = names,
    position = "top"
  ) +
  scale_y_discrete(
    limits = rev(names)
  ) +
  geom_text(aes(label = synapses, size = sqrt(synapses))) +
  guides(size = "none")

# Saving R ggplot with R ggsave Function
ggsave("pictures/head_trunk_matrix.png",
       width = 2600,
       height = 2600, limitsize = TRUE,
       units = c("px")
)  

write.csv(
  t(head_trunk_matrix), 
  "source_data/Figure10_fig_suppl3_source_data1.txt"
)


# Sankey plot ---------------------

head_trunk_matrix_filt <- head_trunk_matrix

# edge weight filtering on the matrix
head_trunk_matrix_filt[head_trunk_matrix_filt < 20] <- 0

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  sqrt(t(head_trunk_matrix_filt)),
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

# Convert to object suitable for networkD3
Connectome_graph_d3 <- igraph_to_networkD3(Conn_graph)
# weights of the edges
E(Conn_graph)$weight

# calculate node weighted degree
degree <- degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)

# Convert to object suitable for networkD3
Connectome_graph_d3 <- igraph_to_networkD3(Conn_graph)

# assign weights to the nodes
Connectome_graph_d3$nodes$weight <- degree
Connectome_graph_d3$links$group <- Connectome_graph_d3$links$source
Connectome_graph_d3$links$group
Connectome_graph_d3$nodes$group <- c("SN", "IN", "MN",  "effector", "SN", "IN", "MN",  "effector")

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["SN", "IN", "MN",  "effector"]) .range(["#E69F00", "#CC79A7", "#0072B2", "#cccccc"])'

# Plot sankeyNetwork
head_trunk_modules <- networkD3::sankeyNetwork(
  Links = Connectome_graph_d3$links, 
  Nodes = Connectome_graph_d3$nodes, Source = "source",
  Target = "target", NodeID = "name", Value = "value",
  NodeGroup = "group",
  colourScale = my_color, LinkGroup = NULL, units = "", 
  fontSize = 56,
  fontFamily = "Arial", nodeWidth = 30, 
  nodePadding = 750, margin = NULL,
  height = 1000, width = 1000, 
  iterations = 1000, sinksRight = TRUE
)
head_trunk_modules

saveNetwork(head_trunk_modules, "pictures/Sankey_head_trunk.html")
webshot2::webshot(url="pictures/Sankey_head_trunk.html",
                  file="pictures/Sankey_head_trunk.png",
                  vwidth=1000, vheight=300, #define the size of the browser window
                  cliprect = c(0, 0, 1000, 1000), zoom=2)

# figure assembly ----------------

panelA <- ggdraw() + 
  draw_image(
    readPNG("pictures/head_trunk_matrix.png"), 
    scale = 1
  )

panelB <- ggdraw() + 
  draw_image(
    readPNG("pictures/Sankey_head_trunk.png"), 
    scale = 1
  )

layout = "
AAAAABBB
"

Figure10_fig_suppl3 <- panelA + panelB + 
  plot_layout(
    design = layout, 
    heights = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "plain")
  )

ggsave("Figures/Figure10_fig_suppl3.png",
       limitsize = FALSE,
       units = c("px"), Figure10_fig_suppl3,
       width = 4400, height = 2600, bg = "white"
)

ggsave("Figures/Figure10_fig_suppl3.pdf",
       limitsize = FALSE,
       units = c("px"), Figure10_fig_suppl3, 
       width = 4400, height = 2600
)

