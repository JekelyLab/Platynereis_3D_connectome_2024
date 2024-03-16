# R code to generate the head-trunk connectivity figure 10 and Fig 10 fig suppl 2of the Platynereis connectome paper
# Gaspar Jekely 2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

#color scheme for segments
segmental_colors <- brewer.pal(6, 'Paired')
pie(rep(1,6),col=segmental_colors,segmental_colors)

head_connectome_skids <- skids_by_2annotations("episphere", "connectome")
trunk_connectome_skids <- skids_by_2annotations("torso", "connectome")

head_trunk_desc_skids <- skids_by_2annotations("episphere", "head_trunk")
head_trunk_asc_skids <- skids_by_2annotations("torso", "head_trunk")

length(head_trunk_desc_skids)
length(head_trunk_asc_skids)

scalebar_50um_v = read.neurons.catmaid("^scalebar_50um_ventral$", pid=11)
outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                              invertFaces = T, conn = NULL, pid = 11)


head_connectome_cells <- nlapply(
  read.neurons.catmaid(head_connectome_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
  )

trunk_connectome_cells <- nlapply(
  read.neurons.catmaid(trunk_connectome_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
)

head_trunk_desc_cells <- nlapply(
  read.neurons.catmaid(head_trunk_desc_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
)

head_trunk_asc_cells <- nlapply(
  read.neurons.catmaid(head_trunk_asc_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
)

# plotting ----------------------------------------------------------------

#plot all cells
plot_background_ventral()
nview3d("ventral", extramat=rotationMatrix(0.4, 1, 0.048, 0.5))
plot3d(
  scalebar_50um_v, alpha = 1, col = "black",
  lwd=3, add = T
)
       
par3d(zoom=0.52)
#z-axis clip
clipplanes3d(0, 0, 1, -10000)

plot3d(
  head_connectome_cells, alpha = 0.6, col = segmental_colors[1], 
  WithConnectors = F, WithNodes = F, soma=T,
  lwd=1, add = T
)

plot3d(
  trunk_connectome_cells, alpha = 0.5, col = segmental_colors[6], 
  WithConnectors = F, WithNodes = F, soma=T,
  lwd=1, add = T
)

plot3d(outline, add=T, alpha=0.02, col="#E2E2E2") 

# make snapshot
rgl.snapshot("pictures/Head_trunk_cells.png")
close3d()

#plot head-trunk projecting cells only

plot_background_ventral()
nview3d("ventral", extramat=rotationMatrix(0.4, 1, 0.048, 0.5))
par3d(zoom=0.52)
#z-axis clip
clipplanes3d(0, 0, 1, -10000)

plot3d(
  head_trunk_desc_cells, alpha = 1, col = segmental_colors[1], 
  WithConnectors = F, WithNodes = F, soma=T,
  lwd=1, add = T
)

plot3d(outline, add=T, alpha=0.02, col="#E2E2E2") 

# make snapshot
rgl.snapshot("pictures/Head_trunk_desc_cells.png")
close3d()

plot_background_ventral()
nview3d("ventral", extramat=rotationMatrix(0.4, 1, 0.048, 0.5))

par3d(zoom=0.52)
#z-axis clip
clipplanes3d(0, 0, 1, -10000)

plot3d(
  head_trunk_asc_cells, alpha = 1, col = segmental_colors[6], 
  WithConnectors = F, WithNodes = F, soma=T,
  lwd=1, add = T
)

plot3d(outline, add=T, alpha=0.02, col="#E2E2E2") 

# make snapshot
rgl.snapshot("pictures/Head_trunk_asc_cells.png")
close3d()


# connectome coloured by head trunk ---------------------------------------

#read the connectome visNetwork R data file
#this file was generated with /Figure2-connectome/Figure_connectome.R

conn_graph.visn <- readRDS("supplements/connectome_graph.rds", refhook = NULL)
#shifting coordinates
#conn.tb <-  as_tbl_graph(conn_graph.visn)
#conn.tb.shifted <- conn.tb %>% 
#  mutate(y = ifelse(
#    segment == 'episphere', y+1, y))
#conn_graph.visn <- toVisNetworkData(as_tbl_graph(conn.tb.shifted))

coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol = 2)
library(autoimage)

conn_graph.visn$nodes$segment

coords_rotated <- autoimage::rotate(
  coords, 
  pi/3.5, 
  pivot = c(0, 0)
)

{
  #overwrite group value (partition) with segment value (for colouring)
  conn_graph.visn$nodes$group <-  as.character(conn_graph.visn$nodes$segment)
  
  #rename segments
  conn_graph.visn$nodes$group <- gsub("episphere", "head", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group <- gsub("segment_", "sg", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group
  
  #remove colour info (which takes precedence over group colour)
  conn_graph.visn$nodes$color <- c()
  
  visNet_head_trunk <- visNetwork(conn_graph.visn$nodes,conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0),
             scaling=list(min=0.1, max=1),
             color = list(inherit=TRUE, opacity=0.2),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 0.5, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(background=conn_graph.visn$nodes$group, border='black'),
             opacity = 1, 
             font = list(size = 0)) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                       algorithm = 'hierarchical',labelOnly=FALSE), 
               width = 2000, height = 2000, autoResize = FALSE) %>%
    visGroups(groupname = "head", color=segmental_colors[1], shape = "dot", 
              opacity=1, size=20) %>%
    visGroups(groupname = "sg0", shape = "dot", 
              opacity=0.6, size=10, color=segmental_colors[6]) %>%
    visGroups(groupname = "sg1", shape = "dot", 
              opacity=0.6, size=10, color=segmental_colors[6]) %>%
    visGroups(groupname = "sg2", shape = "dot", 
              opacity=0.6, size=10, color=segmental_colors[6]) %>%
    visGroups(groupname = "sg3", shape = "dot", 
              opacity=0.6, size=10, color=segmental_colors[6]) %>%
    visGroups(groupname = "pygidium", shape = "dot", 
              opacity=0.6, size=10, color=segmental_colors[6]) %>%
    visGroups(groupname = "fragment", shape = "dot", 
              opacity=0, size=2, color='#aaaaaa')  %>%
    addFontAwesome()
  
  visNet_head_trunk
 
  
  #save as html
  saveNetwork(visNet_head_trunk, "pictures/Full_connectome_head_trunk.html", selfcontained = TRUE)
  webshot::webshot(url="pictures/Full_connectome_head_trunk.html",
                    file="pictures/Full_connectome_head_trunk.png",
                    vwidth = 2000, vheight = 2000, #define the size of the browser window
                    cliprect = c(130, 140, 1860, 1880), zoom=2, delay = 1)
}

# retrieve partners --------------------

head_trunk_connections <- catmaid_get_connectors_between(
  pre_skids = head_connectome_skids, 
  post_skids = trunk_connectome_skids, pid = 11
)


#select trunk targets of head neurons
head_targets_in_trunk_skids <- as_tibble(head_trunk_connections) %>%
  select(post_skid) %>%
  unique() %>%
  pull()
#trunk targets of head neurons
length(head_targets_in_trunk_skids)

head_targets_in_trunk_skids_with_annot <- 
  as_tibble(
    catmaid_get_annotations_for_skeletons(
      head_targets_in_trunk_skids, pid = 11
      )
    )

SN_trunk_skids <- head_targets_in_trunk_skids_with_annot %>%
  filter(annotation == "Sensory neuron") %>%
  select(skid) %>%
  pull()
IN_trunk_skids <- head_targets_in_trunk_skids_with_annot %>%
  filter(annotation == "interneuron") %>%
  select(skid) %>%
  pull()
MN_trunk_skids <- head_targets_in_trunk_skids_with_annot %>%
  filter(annotation == "motorneuron") %>%
  select(skid) %>%
  pull()
effector_trunk_skids <- head_targets_in_trunk_skids_with_annot %>%
  filter(annotation == "effector") %>%
  select(skid) %>%
  pull()

#plot all cells
plot_background_ventral()
nview3d("ventral", extramat=rotationMatrix(0.4, 1, 0.048, 0.5))
par3d(zoom=0.52)
#z-axis clip
clipplanes3d(0, 0, 1, -10000)

plot3d(nlapply(
  read.neurons.catmaid(SN_trunk_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
), alpha = 1, col = "#E69F00", 
WithConnectors = F, soma=T,
lwd=1, add = T
)

plot3d(nlapply(
  read.neurons.catmaid(IN_trunk_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
  ), alpha = 0.5, col = "#CC79A7", 
  WithConnectors = F, soma=T, lwd=1, add = T
  )

plot3d(nlapply(
  read.neurons.catmaid(MN_trunk_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
  ), alpha = 0.8, col = "#0072B2", 
  WithConnectors = F, soma=T, lwd=1, add = T
  )

plot3d(nlapply(
  read.neurons.catmaid(effector_trunk_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
  ), alpha = 1, col = "grey80", 
  WithConnectors = F, soma=T, lwd=1, add = T
  )

plot3d(outline, add=T, alpha=0.02, col="#E2E2E2") 

# make snapshot
rgl.snapshot("pictures/head_targets_in_trunk.png")
close3d()

#select trunk targets of head neurons

trunk_head_connections <- catmaid_get_connectors_between(
  pre_skids = trunk_connectome_skids, 
  post_skids = head_connectome_skids, pid = 11
)

#select trunk targets of head neurons
trunk_targets_in_head_skids <- as_tibble(trunk_head_connections) %>%
  select(post_skid) %>%
  unique() %>%
  pull()
length(trunk_targets_in_head_skids)

trunk_targets_in_head_skids_with_annot <- 
  as_tibble(
    catmaid_get_annotations_for_skeletons(
      trunk_targets_in_head_skids, pid = 11
    )
  )

SN_head_skids <- trunk_targets_in_head_skids_with_annot %>%
  filter(annotation == "Sensory neuron") %>%
  select(skid) %>%
  pull()
IN_head_skids <- trunk_targets_in_head_skids_with_annot %>%
  filter(annotation == "interneuron") %>%
  select(skid) %>%
  pull()
MN_head_skids <- trunk_targets_in_head_skids_with_annot %>%
  filter(annotation == "motorneuron") %>%
  select(skid) %>%
  pull()
effector_head_skids <- trunk_targets_in_head_skids_with_annot %>%
  filter(annotation == "effector") %>%
  select(skid) %>%
  pull()

#plot all cells
plot_background_ventral()
nview3d("ventral", extramat=rotationMatrix(0.4, 1, 0.048, 0.5))
par3d(zoom=0.52)
#z-axis clip
clipplanes3d(0, 0, 1, -10000)

plot3d(nlapply(
  read.neurons.catmaid(SN_head_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
), alpha = 1, col = "#E69F00", 
WithConnectors = F, soma=T,
lwd=1, add = T
)

plot3d(nlapply(
  read.neurons.catmaid(IN_head_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
), alpha = 0.5, col = "#CC79A7", 
WithConnectors = F, soma=T, lwd=1, add = T
)

plot3d(nlapply(
  read.neurons.catmaid(MN_head_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
), alpha = 0.8, col = "#0072B2", 
WithConnectors = F, soma=T, lwd=1, add = T
)

plot3d(nlapply(
  read.neurons.catmaid(effector_head_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
), alpha = 1, col = "grey80", 
WithConnectors = F, soma=T, lwd=1, add = T
)

plot3d(outline, add=T, alpha=0.02, col="#E2E2E2") 

# make snapshot
rgl.snapshot("pictures/trunk_targets_in_head.png")
close3d()


# rank by syn number ------------------------------------------------------

#get neuron names
head_neurons <- as_tibble(head_trunk_connections) %>%
  select(pre_skid) %>%
  pull() %>%
  catmaid_get_neuronnames(pid = 11)

head_cells_annot <- 
  as_tibble(
    catmaid_get_annotations_for_skeletons(
      names(head_connectome_cells), pid = 11
    )
  )

SN_all_head_skids <- head_cells_annot %>%
  filter(annotation == "Sensory neuron") %>%
  select(skid) %>%
  pull()
IN_all_head_skids <- head_cells_annot %>%
  filter(annotation == "interneuron") %>%
  select(skid) %>%
  pull()
MN_all_head_skids <- head_cells_annot %>%
  filter(annotation == "motorneuron") %>%
  select(skid) %>%
  pull()
effector_all_head_skids <- head_cells_annot %>%
  filter(annotation == "effector") %>%
  select(skid) %>%
  pull()

head_trunk_connections_annot <- head_trunk_connections %>%
  mutate(postsyn_class = ifelse(
    post_skid %in% SN_all_trunk_skids, 'SN', 
    ifelse(post_skid %in% IN_all_trunk_skids, 'IN',
           ifelse(post_skid %in% MN_all_trunk_skids, 'MN', 
                  ifelse(post_skid %in% effector_all_trunk_skids, 'effector', 'n.a.')))
  )) %>%
  mutate(presyn_class = ifelse(
    pre_skid %in% SN_all_head_skids, 'SN', 
    ifelse(pre_skid %in% IN_all_head_skids, 'IN',
           ifelse(pre_skid %in% MN_all_head_skids, 'MN', 
                  ifelse(pre_skid %in% effector_all_head_skids, 'effector', 'n.a.')))
  )) %>%
  mutate(presyn_color = ifelse(
    pre_skid %in% SN_all_head_skids, '#E69F00', 
    ifelse(pre_skid %in% IN_all_head_skids, '#CC79A7',
           ifelse(pre_skid %in% MN_all_head_skids, '#0072B2', 
                  ifelse(pre_skid %in% effector_all_head_skids, 'grey50', 'n.a.')))
  )) %>%
  mutate(names = head_neurons) %>%
  group_by(names) %>%
  mutate(num_occurrences = n()) %>%
  arrange(desc(num_occurrences)) %>%
  select(names, postsyn_class, num_occurrences, presyn_color)

head_trunk_connections_annot_reduced <- head_trunk_connections_annot %>%
  group_by(names, postsyn_class) %>%
  mutate(num_class_target = n())

theme_head_trunk <- theme(
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  axis.text.x = element_text(
    angle = 90, hjust = 1, vjust = 0.5,
    size = 10, color = "black"
  ),
  axis.text.y = element_text(
    angle = 0, hjust = 1, vjust = 0.5,
    size = 10, color = "black"
  ),
  axis.title = element_text(
    size = 12, color = "black"
  ),
  axis.ticks = element_line(linewidth = 0.2),
  axis.ticks.length = unit(0.4, "mm"),
  axis.line = element_line(linewidth = 0.2),
  legend.text = element_text(size = 8, hjust = 0, margin = margin(r = 2)),
  legend.key.height = unit(2.5, "mm"),
  legend.key.width = unit(2, "mm"),
  legend.box.spacing = unit(0, "mm"),
  legend.title = element_text(size = 10)
)

#top 50 names
top_50_names <- as_tibble(head_trunk_connections_annot) %>%
  group_by(names) %>%
  summarize(num_occurrences = n()) %>%
  arrange(desc(num_occurrences)) %>%
  top_n(50, num_occurrences) %>%
  select(names) %>%
  pull()

#color list
top_50_colors <- as_tibble(head_trunk_connections_annot) %>%
  distinct(names, presyn_color, .keep_all = TRUE) %>%
  filter(names %in% top_50_names) %>%
  pull(presyn_color)

top_50_colors

head_trunk_connections_annot_reduced %>%
  distinct(names, postsyn_class, .keep_all = TRUE) %>%
  group_by(names) %>%
  ggplot(aes(
    x = names, y = postsyn_class,
    fill = sqrt(num_class_target)
    ))+
  geom_tile() +
#  scale_fill_manual(values = c("SN" = "#E69F00", "IN" = "#CC79A7",
#                               "MN" = "#0072B2", "effector" = "grey50")) +
  scale_fill_gradientn(colours = c("grey97", "#0072B2")) +
  scale_x_discrete(limits = top_50_names) +
  scale_y_discrete(limits = c('effector', 'MN', 'IN', 'SN')) +
  theme_minimal() +
  theme_head_trunk +
  labs(x = "head neurons", y = "target class (trunk only)", fill = "sqrt\n(synapses)") +
  theme(axis.text.x = element_text(color = top_50_colors))

# Save source data
write.table(head_trunk_connections_annot_reduced, "source_data/Figure10_source_data1.txt", sep = "\t")
head_trunk_connections_annot_reduced <- read.table("source_data/Figure10_source_data1.txt", sep = "\t")

# Saving plot
ggsave("pictures/Head_cells_trunk_targets_by_class.png",
       width = 3600, height = 1400, limitsize = TRUE,
       units = c("px"), bg = "white")


# rank for trunk neurons ------------

#get neuron names
trunk_neurons <- as_tibble(trunk_head_connections) %>%
  select(pre_skid) %>%
  pull() %>%
  catmaid_get_neuronnames(pid = 11)

trunk_cells_annot <- 
  as_tibble(
    catmaid_get_annotations_for_skeletons(
      names(trunk_connectome_cells), pid = 11
    )
  )

SN_all_trunk_skids <- trunk_cells_annot %>%
  filter(annotation == "Sensory neuron") %>%
  select(skid) %>%
  pull()
IN_all_trunk_skids <- trunk_cells_annot %>%
  filter(annotation == "interneuron") %>%
  select(skid) %>%
  pull()
MN_all_trunk_skids <- trunk_cells_annot %>%
  filter(annotation == "motorneuron") %>%
  select(skid) %>%
  pull()
effector_all_trunk_skids <- trunk_cells_annot %>%
  filter(annotation == "effector") %>%
  select(skid) %>%
  pull()

trunk_head_connections_annot <- trunk_head_connections %>%
  mutate(postsyn_class = ifelse(
    post_skid %in% MN_all_head_skids, 'SN', 
    ifelse(post_skid %in% IN_head_skids, 'IN',
           ifelse(post_skid %in% MN_all_head_skids, 'MN', 
                  ifelse(post_skid %in% effector_head_skids, 'effector', 'n.a.')))
  )) %>%
  mutate(presyn_color = ifelse(
    pre_skid %in% SN_all_trunk_skids, '#E69F00', 
    ifelse(pre_skid %in% IN_all_trunk_skids, '#CC79A7',
           ifelse(pre_skid %in% MN_all_trunk_skids, '#0072B2', 
                  ifelse(pre_skid %in% effector_all_trunk_skids, 'grey50', 'n.a.')))
  )) %>%
  mutate(presyn_class = ifelse(
    pre_skid %in% SN_all_trunk_skids, 'SN', 
    ifelse(pre_skid %in% IN_all_trunk_skids, 'IN',
           ifelse(pre_skid %in% MN_all_trunk_skids, 'MN', 
                  ifelse(pre_skid %in% effector_all_trunk_skids, 'effector', 'n.a.')))
  )) %>%
  mutate(names = trunk_neurons) %>%
  group_by(names) %>%
  mutate(num_occurrences = n()) %>%
  arrange(desc(num_occurrences)) %>%
  filter(presyn_class != "n.a") %>%
  select(pre_skid, names, presyn_class, postsyn_class, num_occurrences, presyn_color)

trunk_head_connections_annot_reduced <- trunk_head_connections_annot %>%
  group_by(names, postsyn_class) %>%
  mutate(num_class_target = n())


MN_all_trunk_skids %in% '1277759'

#top 50 names
top_50_names_trunk <- as_tibble(trunk_head_connections_annot) %>%
  group_by(names) %>%
  summarize(num_occurrences = n()) %>%
  arrange(desc(num_occurrences)) %>%
  top_n(50, num_occurrences) %>%
  select(names) %>%
  pull()
top_50_names_trunk
#color list
top_50_colors_trunk <- as_tibble(trunk_head_connections_annot) %>%
  distinct(names, presyn_color, .keep_all = TRUE) %>%
  filter(names %in% top_50_names_trunk) %>%
  pull(presyn_color)

top_50_colors_trunk

as_tibble(trunk_head_connections_annot_reduced) %>%
  group_by(names) %>%
  ggplot(aes(
    x = names, y = postsyn_class,
    fill = sqrt(num_class_target)
  ))+
  geom_tile() +
  #  scale_fill_manual(values = c("SN" = "#E69F00", "IN" = "#CC79A7",
  #                               "MN" = "#0072B2", "effector" = "grey50")) +
  scale_fill_gradientn(colours = c("white", "#0072B2")) +
  scale_x_discrete(limits = top_50_names_trunk) +
  scale_y_discrete(limits = c('effector', 'MN', 'IN', 'SN')) +
  theme_minimal() +
  theme_head_trunk +
  labs(x = "trunk neurons", y = "target class (head only)", fill = "sqrt\n(synapses)") +
  theme(axis.text.x = element_text(color = top_50_colors_trunk))


# Save source data
write.table(trunk_head_connections_annot_reduced, "source_data/Figure10_source_data2.txt", sep = "\t")
trunk_head_connections_annot_reduced <- read.table("source_data/Figure10_source_data2.txt", sep = "\t")

# Saving plot
ggsave("pictures/trunk_cells_head_targets_by_class.png",
       width = 3600, height = 1400, limitsize = TRUE,
       units = c("px"), bg = "white")


# Make multi-panel figure -------------------------------------


# read png and convert to image panel

panel_full <- ggdraw() + draw_image(
  readPNG("pictures/Head_trunk_cells.png")
) + 
  draw_label("head cells", x = 0.3, y = 0.99, size = 10, color = segmental_colors[1]) + 
  draw_label("trunk cells", x = 0.7, y = 0.99, size = 10, color = segmental_colors[6]) +
  draw_label(expression(paste("50 ", mu, " m")), x = 0.86, y = 0.12, size = 10)
  

panel_conn <- ggdraw() + draw_image(
  readPNG("pictures/Full_connectome_head_trunk.png")
) + 
  draw_label("head cells", x = 0.3, y = 0.99, size = 10, color = segmental_colors[1]) + 
  draw_label("trunk cells", x = 0.7, y = 0.99, size = 10, color = segmental_colors[6])

panel_desc <- ggdraw() + draw_image(
  readPNG("pictures/Head_trunk_desc_cells.png")
) +
  draw_label("head -> trunk", x = 0.5, y = 0.99, size = 10)

panel_asc <- ggdraw() + draw_image(
  readPNG("pictures/Head_trunk_asc_cells.png")
) +
  draw_label("trunk -> head", x = 0.5, y = 0.99, size = 10)

panel_head_trunk <- ggdraw() + draw_image(
  readPNG("pictures/head_targets_in_trunk.png")
) +
  draw_label("head targets in trunk", x = 0.5, y = 0.99, size = 10)

panel_trunk_head <- ggdraw() + draw_image(
  readPNG("pictures/trunk_targets_in_head.png")
) +
  draw_label("trunk targets in head", x = 0.5, y = 0.99, size = 10) +
  draw_label("sensory neuron", x = 0.23, y = 0.12, size = 9, color = "#E69F00") +
  draw_label("interneuron", x = 0.2, y = 0.07, size = 9, color = "#CC79A7") +
  draw_label("motor neuron", x = 0.76, y = 0.12, size = 9, color = "#0072B2") +
  draw_label("effector", x = 0.75, y = 0.07, size = 9, color = "grey50")

panel_head_trunk_targets <- ggdraw() + draw_image(
  readPNG("pictures/Head_cells_trunk_targets_by_class.png")
)

panel_trunk_head_targets <- ggdraw() + draw_image(
  readPNG("pictures/trunk_cells_head_targets_by_class.png")
) 

# define layout with textual representation for pathchwork assembly of figure
layout <- "
ABCD
EFFF
GHHH
"

Fig10 <- panel_full +panel_conn + panel_desc + panel_asc +
  panel_head_trunk + panel_head_trunk_targets +
  panel_trunk_head + panel_trunk_head_targets +
  plot_layout(design = layout, widths = c(0.6, 0.8, 0.6, 0.6)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12, face = "plain"))
  
ggsave("Figures/Figure10.png",
         limitsize = FALSE,
         units = c("px"), Fig10, width = 2400, height = 2400, bg = "white"
)
  

ggsave("Figures/Figure10.pdf",
       limitsize = FALSE,
       units = c("px"), Fig10, width = 2400, height = 2400
)


# Figure supplement -------------------------------------------------------

trunk_sg0_skids <- trunk_cells_annot %>%
  filter(annotation == "segment_0") %>%
  select(skid) %>%
  pull()
trunk_sg1_skids <- trunk_cells_annot %>%
  filter(annotation == "segment_1") %>%
  select(skid) %>%
  pull()
trunk_sg2_skids <- trunk_cells_annot %>%
  filter(annotation == "segment_2") %>%
  select(skid) %>%
  pull()
trunk_sg3_skids <- trunk_cells_annot %>%
  filter(annotation == "segment_3") %>%
  select(skid) %>%
  pull()
trunk_pyg_skids <- trunk_cells_annot %>%
  filter(annotation == "pygidium") %>%
  select(skid) %>%
  pull()


head_trunk_connections_annot_seg <- head_trunk_connections %>%
  mutate(postsyn_seg = ifelse(
    post_skid %in% trunk_sg0_skids, 'segment_0', 
    ifelse(post_skid %in% trunk_sg1_skids, 'segment_1',
           ifelse(post_skid %in% trunk_sg2_skids, 'segment_2',
           ifelse(post_skid %in% trunk_sg3_skids, 'segment_3', 
                  ifelse(post_skid %in% trunk_pyg_skids, 'pygidium', 'n.a.'))))
  )) %>%
  mutate(presyn_class = ifelse(
    pre_skid %in% SN_all_head_skids, 'SN', 
    ifelse(pre_skid %in% IN_all_head_skids, 'IN',
           ifelse(pre_skid %in% MN_all_head_skids, 'MN', 
                  ifelse(pre_skid %in% effector_all_head_skids, 'effector', 'n.a.')))
  )) %>%
  mutate(presyn_color = ifelse(
    pre_skid %in% SN_all_head_skids, '#E69F00', 
    ifelse(pre_skid %in% IN_all_head_skids, '#CC79A7',
           ifelse(pre_skid %in% MN_all_head_skids, '#0072B2', 
                  ifelse(pre_skid %in% effector_all_head_skids, 'grey50', 'n.a.')))
  )) %>%
  mutate(names = head_neurons) %>%
  group_by(names) %>%
  mutate(num_occurrences = n()) %>%
  arrange(desc(num_occurrences)) %>%
  select(names, postsyn_seg, num_occurrences, presyn_color)


head_trunk_connections_annot_seg_reduced <- head_trunk_connections_annot_seg %>%
  group_by(names, postsyn_seg) %>%
  mutate(num_class_target = n())


head_trunk_connections_annot_seg_reduced %>%
  distinct(names, postsyn_seg, .keep_all = TRUE) %>%
  group_by(names) %>%
  ggplot(aes(
    x = names, y = postsyn_seg,
    fill = sqrt(num_class_target)
  ))+
  geom_tile() +
  #  scale_fill_manual(values = c("SN" = "#E69F00", "IN" = "#CC79A7",
  #                               "MN" = "#0072B2", "effector" = "grey50")) +
  scale_fill_gradientn(colours = c("grey97", "#0072B2")) +
  scale_x_discrete(limits = top_50_names) +
  scale_y_discrete(limits = rev(c(
    'segment_0', 'segment_1', 'segment_2', 
    'segment_3', 'pygidium'))) +
  theme_minimal() +
  theme_head_trunk +
  labs(x = "head neurons", y = "target segment (trunk only)", fill = "sqrt\n(synapses)") +
  theme(axis.text.x = element_text(color = top_50_colors))

# Save source data
write.table(head_trunk_connections_annot_seg_reduced, "source_data/Figure10_fig_suppl2_source_data1.txt", sep = "\t")
head_trunk_connections_annot_seg_reduced <- read.table("source_data/Figure10_fig_suppl2_source_data1.txt", sep = "\t")

# Saving plot
ggsave("Figures/Figure10_fig_suppl2.png",
       width = 3600, height = 1400, limitsize = TRUE,
       units = c("px"), bg = "white")

ggsave("Figures/Figure10_fig_suppl2.pdf",
       width = 3600, height = 1400, limitsize = TRUE,
       units = c("px"), bg = "white")


