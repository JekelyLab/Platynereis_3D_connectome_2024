# R code to generate Fig10 fig suppl4 of the 3d Platynereis connectome paper
# Gaspar Jekely 2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# plot graph with coordinates from gephi ----------------------------------

# read graph
connect.tb <- readRDS("supplements/connectome_graph_tibble.rds")

filter_side_class <- function(side_to_filter, class_to_filter){
 connect.tb %>%
    filter(side %in% side_to_filter) %>%
    filter(class %in% class_to_filter) %>%
    select(skids) %>%
    pull()
}

# left side skids --------------------
left_SN_skids <- filter_side_class("left_side", "Sensory neuron")
left_IN_skids <- filter_side_class("left_side", "interneuron")
left_MN_skids <- filter_side_class("left_side", "motorneuron")
left_effector_skids <- filter_side_class("left_side", "effector")

# right side skids --------------------
right_SN_skids <- filter_side_class("right_side", "Sensory neuron")
right_IN_skids <- filter_side_class("right_side", "interneuron")
right_MN_skids <- filter_side_class("right_side", "motorneuron")
right_effector_skids <- filter_side_class("right_side", "effector")

cell_groups <- list(
  left_SN_skids, left_IN_skids, left_MN_skids, left_effector_skids,
  right_SN_skids, right_IN_skids, right_MN_skids, right_effector_skids
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

left_right_matrix <- sapply(cell_groups, function(x)
  sapply(cell_groups, function(y)
  number_of_syns(x, y)
    )
  )

names <- c(
  "SN left", "IN left", "MN left", "effector left",
  "SN right", "IN right", "MN right", "effector right"
  )

colnames(left_right_matrix) <- names
rownames(left_right_matrix) <- names

# plot matrix -------------
as.data.frame(left_right_matrix) %>%
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
ggsave("pictures/left_right_Matrix.png",
       width = 2600,
       height = 2600, limitsize = TRUE,
       units = c("px")
)  

write.csv(
  t(left_right_matrix), 
  "source_data/Figure10_fig_suppl4_source_data1.txt"
)


# Sankey plot ---------------------

left_right_matrix_filt <- left_right_matrix

# edge weight filtering on the matrix
left_right_matrix_filt[left_right_matrix_filt < 20] <- 0

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  sqrt(t(left_right_matrix_filt)),
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

Connectome_graph_d3$nodes$group <- c("SN", "IN", "MN",  "effector", "SN", "IN", "MN",  "effector")
# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["SN", "IN", "MN",  "effector"]) .range(["#E69F00", "#CC79A7", "#0072B2", "#cccccc"])'

# Plot sankeyNetwork
left_right_modules <- networkD3::sankeyNetwork(
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
left_right_modules

saveNetwork(left_right_modules, "pictures/Sankey_left_right.html")
#save from browser after node position adjustment

# figure assembly ----------------

panelA <- ggdraw() + 
  draw_image(
    readPNG("pictures/left_right_Matrix.png"), 
    scale = 1
  )

panelB <- ggdraw() + 
  draw_image(
    readPNG("pictures/Sankey_left_right.png"), 
    scale = 1
  )

layout = "
AAAAABBB
"

Figure10_fig_suppl4 <- panelA + panelB + 
  plot_layout(
    design = layout, 
    heights = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "plain")
  )

ggsave("Figures/Figure10_fig_suppl4.png",
       limitsize = FALSE,
       units = c("px"), Figure10_fig_suppl4,
       width = 4400, height = 2600, bg = "white"
)

ggsave("Figures/Figure10_fig_suppl4.pdf",
       limitsize = FALSE,
       units = c("px"), Figure10_fig_suppl4, 
       width = 4400, height = 2600
)

