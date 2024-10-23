#Platynereis 3d larva connectome paper Figure on gland circuits
#Gaspar Jekely

source("code/Natverse_functions_and_conn.R")

# load cells ---------------------------------------

#Glands
spinGland = nlapply(read.neurons.catmaid("^celltype_non_neuronal7$", pid=11),
               function(x) smooth_neuron(x, sigma=6000))
ciliatedGland = nlapply(read.neurons.catmaid("^celltype_non_neuronal9$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
MVGland = nlapply(read.neurons.catmaid("^celltype_non_neuronal17$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
HeadGland = nlapply(read.neurons.catmaid("^celltype_non_neuronal29$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
outline_HeadGland = nlapply(read.neurons.catmaid(891030, pid=11),
                    function(x) smooth_neuron(x, sigma=200))
InterparaGland = nlapply(read.neurons.catmaid("^celltype_non_neuronal30$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
spinMicroGland = nlapply(read.neurons.catmaid("^celltype_non_neuronal31$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))

#MNs
MNgland_head = nlapply(read.neurons.catmaid("^celltype166$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
MNspinning = nlapply(read.neurons.catmaid("^celltype47$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
INDLSO = nlapply(read.neurons.catmaid("^celltype186$", pid=11),
          function(x) smooth_neuron(x, sigma=6000))
INpreSer = nlapply(read.neurons.catmaid("^celltype4$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
Ser_h1 = nlapply(read.neurons.catmaid("^celltype8$", pid=11),
          function(x) smooth_neuron(x, sigma=6000))

# plot cells --------------------------

plot_background_ventral_no_ac()
par3d(windowRect = c(0, 0, 800, 1066))
plot3d(spinGland, soma=T, lwd=2, add=T, alpha=0.4, col=Okabe_Ito[1])
plot3d(ciliatedGland, soma=T, lwd=2, add=T, alpha=0.5, col=Okabe_Ito[2])
plot3d(MVGland, soma=T, lwd=2, add=T, alpha=1, col=Okabe_Ito[6])
plot3d(HeadGland, soma=T, lwd=2, add=T, alpha=0.6, col=Okabe_Ito[8])
plot3d(outline_HeadGland, soma=F, lwd=1, add=T, alpha=0.3, col=Okabe_Ito[8])
plot3d(InterparaGland, soma=T, lwd=2, add=T, alpha=1, col=Okabe_Ito[7])
plot3d(spinMicroGland, soma=T, lwd=2, add=T, alpha=0.5, col=Okabe_Ito[3])
plot3d(scalebar_50um_ventral2, lwd=3, add=T, alpha=2, col=Okabe_Ito[8])
texts3d(98000, 110000, 182500, text = "50 μm", col = "black", cex = 2)
texts3d(48000, 20000, 130000, text = "spinGland", col = Okabe_Ito[1], cex = 2)
texts3d(48000, 20000, 25000, text = "HeadGland", col = "black", cex = 2)
texts3d(113000, 20000, 160000, text = "InterparaGland", col = Okabe_Ito[7], cex = 2)
texts3d(110000, 20000, 120000, text = "spinMicroGland", col = Okabe_Ito[3], cex = 2)
texts3d(110000, 20000, 74000, text = "ciliatedGland", col=Okabe_Ito[2], cex = 2)
texts3d(42000, 50000, 84000, text = "MVGland", col=Okabe_Ito[6], cex = 2)

rgl.snapshot("pictures/all_glands.png")
close3d()

# plot MNgland-head --------------------------

plot_background_ventral_no_ac()
par3d(windowRect = c(0, 0, 800, 1066))
plot3d(ciliatedGland, soma=T, lwd=2, add=T, alpha=0.5, col=Okabe_Ito[2])
plot3d(MVGland, soma=T, lwd=2, add=T, alpha=1, col=Okabe_Ito[6])
plot3d(MNgland_head, soma=T, lwd=c(3,5), add=T, alpha=1, col=bluepurple[c(9,8)])

texts3d(110000, 20000, 74000, text = "ciliatedGland", col=Okabe_Ito[2], cex = 2)
texts3d(42000, 50000, 84000, text = "MVGland", col=Okabe_Ito[6], cex = 2)
texts3d(55000, 50000, 12000, text = "MNgland-head", col=bluepurple[9], cex = 2)

rgl.snapshot("pictures/MNgland_head_glands_v.png")
close3d()

# anterior view

plot_background()
plot3d(ciliatedGland, soma=T, lwd=2, add=T, alpha=0.5, col=Okabe_Ito[2])
plot3d(MVGland, soma=T, lwd=2, add=T, alpha=1, col=Okabe_Ito[6])
plot3d(MNgland_head, soma=T, lwd=c(3,5), add=T, alpha=1, col=bluepurple[c(9,8)])

plot3d(INDLSO, soma=T, lwd=c(3,5), add=T, alpha=1, col=bluepurple[c(3,4)])
plot3d(INpreSer, soma=T, lwd=c(4,5), add=T, alpha=1, col=bluepurple[c(5,6)])
plot3d(Ser_h1, soma=T, lwd=c(2,3), add=T, alpha=0.6, col=oranges[c(6,8)])
plot3d(scalebar_50um_anterior, soma=T, lwd=3, add=T, col="black")
par3d(zoom = 0.55)

texts3d(120000, 115000, 44000, text = "ciliatedGland", col=Okabe_Ito[2], cex = 2)
texts3d(105000, 123500, 44000, text = "50 μm", col=Okabe_Ito[8], cex = 1.8)
texts3d(72000, 110000, 44000, text = "MVGland", col=Okabe_Ito[6], cex = 2)
texts3d(68000, 37000, 8000, text = "MNgland-head", col=bluepurple[9], cex = 2)
texts3d(97000, 52000, 8000, text = "INDLSO", col=bluepurple[5], cex = 2)
texts3d(45000, 32000, 8000, text = "INpreSer", col=bluepurple[5], cex = 2)
texts3d(33000, 46000, 8000, text = "Ser-h1", col=oranges[8], cex = 2)

rgl.snapshot("pictures/MNgland_head_glands_a.png")
close3d()

# retrieve connectivity  -------------

cell_groups <-  list(
  ciliatedGland,
  MVGland,
  MNgland_head,
  INDLSO, INpreSer, Ser_h1
  )

N_cell_groups <- length(cell_groups)

#get connectivity all against all
synapse_list <- vector("list", N_cell_groups*N_cell_groups)
for (i in 1:N_cell_groups) {
  for (j in 1:N_cell_groups){
    #get connectors between two cell groups
    presyn_skids <- attr(cell_groups[i][[1]],"df")$skid
    postsyn_skids <- attr(cell_groups[j][[1]],"df")$skid
    connectivity <- catmaid_get_connectors_between(pre=presyn_skids, 
                                                   post=postsyn_skids, pid=11)
    #check the number of synapses from group1 -> group2
    N_synapses <-  dim(connectivity)[1]
    #change "NULL" to 0
    if(is.null(N_synapses)){N_synapses = 0}
    
    print ((i*N_cell_groups-N_cell_groups)+j)
    print (N_synapses)
    #add value to synapse list
    synapse_list [[(i*N_cell_groups-N_cell_groups)+j]] <- N_synapses
  }
}

#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=N_cell_groups)

cell_group_names <-  list(
  "ciliatedGland",
  "MVGland", 
  "MNgland_head",
  "INDLSO", "INpreSer",
  "Ser_h1"
  )
rownames(synapse_matrix) <- as.character(cell_group_names)
colnames(synapse_matrix) <- as.character(cell_group_names)
synapse_matrix

write.csv(synapse_matrix, 
          "source_data/Figure7_fig_suppl1_source_data1.txt"
)
# graph conversion --------------------------------------------------------

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

#calculate node weighted degree
degree=degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)

# use visNetwork to plot the network --------------------------------------

## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(Conn_graph)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight

Conn_graph.visn$nodes$value = degree
Conn_graph.visn$nodes
node_colors <- c("#cccccc","#cccccc",
                 "#56B4E9",
                "#CC79A7","#CC79A7", "#BD0026")

Conn_graph.visn$nodes$label
Conn_graph.visn$nodes$group <- list("gland", "gland",
                                    "MN", 
                                    "IN", "IN", "MN")
Conn_graph.visn$nodes$color <- node_colors

#hierarchical layout
Conn_graph.visn$nodes$level <- c("3", "3",
                                 "2",
                                 "1","1", "0")
visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout_nicely", type = "full", physics = FALSE) %>%
  visHierarchicalLayout(levelSeparation=200, 
           nodeSpacing=150,
           treeSpacing = 10,
           direction='LR',
           sortMethod=NULL,
           shakeTowards=NULL) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0.1),
           scaling=list(min=2, max=12),
           color = list(inherit=TRUE, opacity=0.7),
           arrows = list(
             to = list(enabled = TRUE, 
                       scaleFactor = 1.2, type = 'arrow'))
           ) %>%
  visNodes(borderWidth=0.3, 
           color = list(background=Conn_graph.visn$nodes$color, border='black'),
           opacity=0.9,
           shape='dot', 
           font=list(color='black', size=24),
           scaling = list(label=list(enabled=TRUE, min=23, max=28)),
           level= Conn_graph.visn$nodes$level) %>%
  visInteraction(navigationButtons = FALSE,
           dragNodes = TRUE, dragView = FALSE,
           zoomView = TRUE) %>%
  visGroups(groupname = "gland", shape = "square", opacity=1)
visNet

saveNetwork(visNet, "pictures/visNetwork_glands.html")
webshot2::webshot(url = "pictures/visNetwork_glands.html",
                  file = "pictures/visNetwork_glands.png",
                  vwidth = 1000, vheight = 500, #define the size of the browser window
                  cliprect = c(150, 140, 750, 310), zoom = 5, delay = 1)


# table ---------------------

table <- plot_ly(
  type = "table",
  columnwidth = c(8, 6, 6),
  columnorder = c(0, 1, 2),
  header = list(
    values = c(
      "gland type", 
      "segment",
      "# of cells"),
    align = c("center", "center", "center"),
    line = list(width = 1, color = "black"),
    fill = list(color = c("#CCCCCC", "#AAAAAA", "#AAAAAA")),
    font = list(family = "Arial", size = 14, color = "black")
  ),
  cells = list(
    values = rbind(
      c( "spinGland", "ciliatedGland",
         "MVGland", "HeadGland",
         "InterparaGland", "spinMicroGland"
      ),
      c(
        "sg2, sg3", "sg1",
        "sg1", "head",
        "sg2, sg3", "sg2, sg3"
      ),
      c(
        "4", "24",
        "6", "5",
        "4", "19"
      )
    ),
    align = c("center", "center","center","center","center","center"),
    line = list(color = "black", width = 0.3),
    font = list(family = "Arial", size = 12, color = c("black"))
  )
)

table
saveNetwork(table, "pictures/glands_stats_table.html")
webshot::webshot(
  url = "pictures/glands_stats_table.html",
  file = "pictures/glands_stats_table.png",
  vwidth = 300, vheight = 200, # define the size of the browser window
  cliprect = c(23, 58, 235, 152), zoom = 10
)


# multi-panel figure -------------------------------------------------

panel_glands <- ggdraw() + 
  draw_image(readPNG("pictures/all_glands.png")) +
  draw_label("all gland cell types", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)  +
    geom_segment(aes(x = 0.1,
                     y = 0.9,
                     xend = 0.1,
                     yend = 0.82),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
    geom_segment(aes(x = 0.1,
                     y = 0.82,
                     xend = 0.1,
                     yend = 0.9),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
    draw_label("a", x = 0.1, y = 0.93, size = 8) +
    draw_label("p", x = 0.1, y = 0.79, size = 8) 

panel_MNgland_head_v <- ggdraw() + 
  draw_image(readPNG("pictures/MNgland_head_glands_v.png")) +
  draw_label("ventral view", x = 0.6, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_MNgland_head_a <- ggdraw() + 
  draw_image(readPNG("pictures/MNgland_head_glands_a.png")) +
  draw_label("anterior view", x = 0.6, y = 0.99, 
             fontfamily = "sans", size = 11) +
    geom_segment(aes(x = 0.1,
                     y = 0.9,
                     xend = 0.1,
                     yend = 0.82),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
    geom_segment(aes(x = 0.1,
                     y = 0.82,
                     xend = 0.1,
                     yend = 0.9),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
    draw_label("d", x = 0.1, y = 0.93, size = 8) +
    draw_label("v", x = 0.1, y = 0.79, size = 8) 

panel_table <- ggdraw() + 
  draw_image(readPNG("pictures/glands_stats_table.png"))


#check max of network
network_temp <- read_csv("source_data/Figure7_fig_suppl1_source_data1.txt")
network_temp %>%
  select(-...1) %>% min()

panel_glands_network <- ggdraw() + 
  draw_image(readPNG("pictures/visNetwork_glands.png")) +
  draw_label("MNgland-head circuit", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11) + 
  draw_label("# of synapses", x = 0.94, y = 0.9, size = 8, hjust = 1) +
  draw_label("1", x = 0.87, y = 0.82, size = 8, hjust = 1) + 
  draw_label("31", x = 0.87, y = 0.76, size = 8, hjust = 1) +
  draw_line(x = c(0.88, 0.93), y = c(0.82, 0.82), size = 0.3, color = 'grey') +
  draw_line(x = c(0.88, 0.93), y = c(0.76, 0.76), size = 2, color = 'grey')


panel_headGland1 <- ggdraw() + 
  draw_image(readPNG("images_notR/headGland_2um_1.png")) +
  draw_label("headGland", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_headGland2 <- ggdraw() + 
  draw_image(readPNG("images_notR/headGland_2um_2.png")) +
  draw_label("HeadGland", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_ciliatedGland <- ggdraw() + 
  draw_image(readPNG("images_notR/ciliatedGland.png")) +
  draw_label("ciliatedGland", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_MVGland <- ggdraw() + 
  draw_image(readPNG("images_notR/MVGland.png")) +
  draw_label("MVGland", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_spinmicroGland <- ggdraw() + 
  draw_image(readPNG("images_notR/spinmicroGland.png")) +
  draw_label("spinMicroGland", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_interparaGland <- ggdraw() + 
  draw_image(readPNG("images_notR/interparaGland.png")) +
  draw_label("InterparaGland", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)


layout <- "
AAAABBBBCCC
###########
DDDDEEEEEEE
###########
F#G#H#I#J#K
"
Fig_glands  <- panel_glands + panel_MNgland_head_v + panel_MNgland_head_a +
  panel_table + panel_glands_network +
  panel_headGland1 + panel_headGland2 +
  panel_ciliatedGland + panel_MVGland + panel_spinmicroGland +
  panel_interparaGland +
  plot_layout(design = layout, 
              widths = c(1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1), 
              heights = c(1, 0.05, 0.5, 0.05, 0.5)
              ) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))
  
ggsave("Figures/Figure7_fig_suppl1.png", limitsize = FALSE, 
       units = c("px"), Fig_glands, width = 2800, height = 2150, bg='white')

ggsave("Figures/Figure7_fig_suppl1.pdf", limitsize = FALSE, 
       units = c("px"), Fig_glands, width = 2800, height = 2150)
