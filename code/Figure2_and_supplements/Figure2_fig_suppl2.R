# R code to generate Fig2 fig suppl2 of the 3d Platynereis connectome paper
# Gaspar Jekely 2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# plot graph with coordinates from gephi ----------------------------------

# read graph
connect.tb <- readRDS("supplements/connectome_graph_tibble.rds")

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

# function to filter by partition --------------
filter_partition <- function(name){
 connect.tb %>%
    filter(partition_name == name) %>%
    select(skids) %>%
    pull()
}

# skids of partitions -----------
NS_skids <- filter_partition("NS")
MB_skids <- filter_partition("MB")
Mech_l_skids <- filter_partition("mech (l)")
Mech_r_skids <- filter_partition("mech (r)")
Visual_skids <- filter_partition("visual")
cMN_skids <- filter_partition("ciliomotor")
pigment_skids <- filter_partition("pigment")
postural_skids <- filter_partition("postural")
MUS1_skids <- filter_partition("muscle 1")
MUS2_skids <- filter_partition("muscle 2")
MUS3_skids <- filter_partition("muscle 3")
MUS4_skids <- filter_partition("muscle 4")
MUS5_skids <- filter_partition("muscle 5")

cell_groups <- list(
  NS_skids, MB_skids, Mech_l_skids, Mech_r_skids,
  Visual_skids, cMN_skids, pigment_skids, postural_skids,
  MUS1_skids, MUS2_skids, MUS3_skids, MUS4_skids, MUS5_skids
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

partition_conn_matrix <- sapply(cell_groups, function(x)
  sapply(cell_groups, function(y)
  number_of_syns(x, y)
    )
  )

colnames(partition_conn_matrix) <- partition_names
rownames(partition_conn_matrix) <- partition_names
max(partition_conn_matrix)

# plot matrix -------------
as.data.frame(partition_conn_matrix) %>%
  rownames_to_column(var = "presynaptic") %>%
  pivot_longer(-presynaptic, names_to = "postsynaptic",
               values_to = "synapses") %>%
  ggplot(aes(presynaptic, postsynaptic, fill = sqrt(synapses))) +
  geom_tile() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 0, vjust = 0.5, size = 12
    ),
    axis.text.y = element_text(
      angle = 0, hjust = 1, vjust = 0.5, size = 12
    ),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.ticks = element_line(linewidth = 0.1)
  ) +
  scale_fill_gradientn(
    colours = c("white", "#0072B2",  "#D55E00")
  ) +
  labs(
    x = "postsynaptic module",
    y = "presynaptic module",
    title = " ", fill = "sqrt\nsynapses"
  ) +
  scale_x_discrete(
    limits = partition_names,
    position = "top"
  ) +
  scale_y_discrete(
    limits = rev(partition_names)
  ) +
  geom_text(aes(label = synapses, size = sqrt(synapses))) +
  guides(size = "none")

# Saving R ggplot with R ggsave Function
ggsave("pictures/Module_Matrix.png",
       width = 3000,
       height = 3000, limitsize = TRUE,
       units = c("px")
)  


write.csv(t(partition_conn_matrix), 
  "source_data/Figure2_fig_suppl2_source_data1.txt"
  )

# Sankey plot ---------------------

partition_conn_matrix_filt <- partition_conn_matrix
partition_conn_matrix_filt
# edge weight filtering on the matrix
partition_conn_matrix_filt[partition_conn_matrix_filt < 20] <- 0

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  sqrt(t(partition_conn_matrix_filt)),
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

colors <- connect.tb %>%
  as_tibble() %>%
  select(partition, partition_name, color) %>%
  arrange(partition) %>%
  unique() %>%
  slice(c(1:14)) %>% pull(color)
colors

partition_name <- connect.tb %>%
  as_tibble() %>%
  select(partition, partition_name, color) %>%
  arrange(partition) %>%
  unique() %>%
  slice(c(1:14)) %>% pull(partition_name)
partition_name

Connectome_graph_d3$nodes$name

# Define 'group' based on partition:
Connectome_graph_d3$nodes$group <- as.factor(gsub("\\s", "\\_", Connectome_graph_d3$nodes$name))
Connectome_graph_d3$nodes$group

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["mech_(l)", "NS", "muscle_2", "mech_(r)", "MB", "postural", "muscle_1", "muscle_4", "ciliomotor", "visual", "muscle_5", "muscle_3", "pigment"]) .range(["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#E69F00", "#56B4E9"])'

# Plot sankeyNetwork
SN_modules <- networkD3::sankeyNetwork(
  Links = Connectome_graph_d3$links, Nodes = Connectome_graph_d3$nodes, Source = "source",
  Target = "target", NodeID = "name", Value = "value",
  LinkGroup = NULL, units = "", NodeGroup = "group",
  colourScale = my_color, fontSize = 56,
  fontFamily = "Arial", nodeWidth = 30, nodePadding = 250, margin = NULL,
  height = 1000, width = 1000, 
  iterations = 1000, sinksRight = TRUE
)
SN_modules

saveNetwork(SN_modules, "pictures/Sankey_modules.html")

#save png from browser

# module members -------------------

connect.tb %>%
  as_tibble() %>%
  filter(grepl("chaeMe", names)) %>%
  select(names, partition) 

connect.tb %>%
  as_tibble() %>%
  filter(partition==4) %>%
  select(names) %>%
  pull()

# module sensitivity to resolution  ------------------

conn.igraph <- as.igraph(connect.tb)

num_of_partitions <- list()
for(resol_param in 1:100){
  # clustering with Leiden algorithm
  partition <- leiden(conn.igraph,
                    weights = E(conn.igraph)$weight,
                    partition_type = "RBConfigurationVertexPartition",
                    resolution_parameter = resol_param/50,
                    n_iterations = -1, seed = 31
                    )
  num_of_partitions[resol_param] <- max(partition)
  print(max(partition))
}

partition_data <- tibble(
  num_modules = unlist(num_of_partitions),
  resol_param = c(1:100/50)
)

plot_theme <-  theme_minimal() + 
  theme(
  text = element_text(size = 10),
  axis.text = element_text(angle = 90, hjust = 1, size = 10),
  axis.title = element_text(size = 12),
  axis.ticks.length.y = unit(1, "mm"),
  legend.key.size = unit(3, "mm"),
  legend.position = "top",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 8),
  axis.line = element_blank()
)

partition_plot <- partition_data %>%
  ggplot(aes(resol_param, num_modules)) +
  geom_point() +
  plot_theme +
  labs(x = "Leiden resolution parameter", y = "# of modules")

partition_plot
ggsave("pictures/partition_plot.png", width = 1200, height = 1200, units = "px", bg = "white")

# figure assembly -----names()# figure assembly ----------------

panelA <- ggdraw() + 
  draw_image(
    readPNG("pictures/Module_Matrix.png"), 
    scale = 1
    )

panelB <- ggdraw() + 
  draw_image(
    readPNG("pictures/Sankey_modules.png"), 
    scale = 1
  )

panel_partition_plot <- ggdraw() +
  draw_image(readPNG("pictures/partition_plot.png"))

layout = "
AB
AC
"

Figure2_fig_suppl2 <- panelA + panelB + panel_partition_plot +
  plot_layout(
    design = layout, 
    heights = c(1,1), 
    widths = c(3,1.5)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "plain")
    )

ggsave("Figures/Figure2_fig_suppl2.png",
       limitsize = FALSE,
       units = c("px"), Figure2_fig_suppl2,
       width = 4400, height = 3000, bg = "white"
)

ggsave("Figures/Figure2_fig_suppl2.pdf",
       limitsize = FALSE,
       units = c("px"), Figure2_fig_suppl2, 
       width = 4400, height = 3000
)


