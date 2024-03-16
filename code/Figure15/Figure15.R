# R/natverse code to generate Figure 15 for the Platynereis 3d connectome paper
# Gaspar Jekely Feb-Dec 2022

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")


# read neurons from CATMAID -----------------------------------------------

# read chaeMech and partner neuron groups
{
  chaeMech <- nlapply(
    read.neurons.catmaid("^celltype71$", pid = 11),
    function(x) smooth_neuron(x, sigma = 8000)
  )
  SNbronto <- nlapply(
    read.neurons.catmaid("^celltype168$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  MNhose <- nlapply(
    read.neurons.catmaid("^celltype66$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MNring <- nlapply(
    read.neurons.catmaid("^celltype64$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MNspinning <- nlapply(
    read.neurons.catmaid("^celltype47$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MNcrab <- nlapply(
    read.neurons.catmaid("^celltype65$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MNbow <- nlapply(
    read.neurons.catmaid("^celltype67$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MNacicX <- nlapply(
    read.neurons.catmaid("^celltype69$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MNspider <- nlapply(
    read.neurons.catmaid("^celltype61$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  INchaeMech <- nlapply(
    read.neurons.catmaid("^celltype70$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INsplitCR <- nlapply(
    read.neurons.catmaid("^celltype73$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INsplitBronto <- nlapply(
    read.neurons.catmaid("^celltype149$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INsplitCRATO <- nlapply(
    read.neurons.catmaid("^celltype74$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INasc_pyg <- nlapply(
    read.neurons.catmaid("^celltype147$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INCM <- nlapply(
    read.neurons.catmaid("^celltype60$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  INFVa_pyg <- nlapply(
    read.neurons.catmaid("^celltype72$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  INrope <- nlapply(
    read.neurons.catmaid("^celltype58$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )


  CR <- nlapply(
    read.neurons.catmaid("^CRneurons$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  PU <- nlapply(
    read.neurons.catmaid("^PUneurons$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
}

# load other cell clusters
{
  stomodeum <- nlapply(
    read.neurons.catmaid("^stomodeum$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MUSring_pyg <- nlapply(
    read.neurons.catmaid("^celltype_non_neuronal78$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  gut <- nlapply(
    read.neurons.catmaid("^gut$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  girdle <- nlapply(
    read.neurons.catmaid("^mechanosensory_girdle$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  chaeta <- nlapply(
    read.neurons.catmaid("^chaeta$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
}

# background plot function ------------------------------------------------

# adjusted background plotting function
plot_background_mech <- function() {
  plot_background_ventral()
  clear3d()
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  plot3d(yolk,
    WithConnectors = F, WithNodes = F, soma = TRUE, lwd = 0,
    rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.05,
    col = "grey80"
  )
  plot3d(stomodeum,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 3,
    rev = FALSE, fixup = F, add = T, alpha = 0.1, col = "grey50"
  )
}

# plot chaeMech cells -----------------------------------------------------

{
  plot_background_mech()
  plot3d(chaeMech,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 3,
    rev = FALSE, fixup = F, add = T, alpha = 0.8,
    col = sample(oranges[4:9], 20, replace = TRUE)
  )
  par3d(zoom = 0.55)
  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(chaeta,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    rev = FALSE, fixup = F, add = T, alpha = 0.3, col = "grey50"
  )
  plot3d(scalebar_50um_ventral, add = T, alpha = 1,
         lwd = 2, col = "black"
  )

  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1) %*% rotationMatrix(-pi / 2, 0, 0, 1)))
  filename <- paste("pictures/Figure_mec_chaeMech_lat.png")
  rgl.snapshot(filename)
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  texts3d(78000, 139000, 25000, text = "stomodeum", col = "black", cex = 2.2)
  texts3d(140000, 139000, 45000, text = "chaetae", col = "black", cex = 2.2)
  plot3d(scalebar_50um_ventral,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 5,
    rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 1,
    col = "black"
  )
  filename <- paste("pictures/Figure_mec_chaeMech_ventr.png")
  rgl.snapshot(filename)

  clear3d()
  plot3d(chaeMech,
    WithConnectors = T, WithNodes = F, soma = T, lwd = 2,
    rev = FALSE, fixup = F, add = T, alpha = 0.5,
    col = "grey80"
  )
  par3d(zoom = 0.55)
  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(yolk,
    WithConnectors = F, WithNodes = F, soma = TRUE, lwd = 0,
    rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.05,
    col = "grey80"
  )
  plot3d(stomodeum,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 3,
    rev = FALSE, fixup = F, add = T, alpha = 0.1, col = "grey50"
  )
  filename <- paste("pictures/Figure_mec_chaeMech_syn.png")
  rgl.snapshot(filename)
  close3d()
}

# circuit analysis and plotting -------------------------------------------

# get connectivity between mech cell clusters and their output clusters
cell_groups <- list(
  chaeMech, CR, SNbronto, INchaeMech,
  INsplitCR, INsplitBronto, INsplitCRATO, INasc_pyg,
  INCM, INFVa_pyg, INrope, MNhose,
  MNring, MNspinning, MNcrab, MNacicX,
  MNspider
)

N_cell_groups <- length(cell_groups)
N_cell_groups

cell_group_attr <- data.frame(
  cell_group_names = c(
    "chaeMech", "CR", "SNbronto", "INchaeMech",
    "INsplitCR", "INsplitBronto", "INsplitCRATO", "INasc-pyg",
    "INCM", "INFVa-pyg", "INrope", "MNhose",
    "MNring", "MNspinning", "MNcrab", "MNacicX",
    "MNspider"
  ),
  type = c(
    "SN", "SN", "SN", "IN",
    "INsplit", "INsplit", "INsplit", "IN",
    "INsplit", "IN", "IN", "MN",
    "MN", "MN", "MN", "MN",
    "MN"
  ),
  level = c(
    "1", "1", "1", "3",
    "2", "2", "2", "3",
    "2", "3", "3", "4",
    "4", "4", "4", "4",
    "4"
  )
)

dim(cell_group_attr)



# iterate through cell group neuron lists and get connectivity for all against all
{
  # define empty synapse list with the right dimensions
  synapse_list <- vector("list", N_cell_groups * N_cell_groups)
  for (i in 1:N_cell_groups) {
    for (j in 1:N_cell_groups) {
      # get connectors between two cell groups
      presyn_skids <- attr(cell_groups[i][[1]], "df")$skid
      postsyn_skids <- attr(cell_groups[j][[1]], "df")$skid
      connectivity <- catmaid_get_connectors_between(
        pre = presyn_skids,
        post = postsyn_skids, pid = 11
      )
      # check the number of synapses from group1 -> group2
      N_synapses <- dim(connectivity)[1]
      # change "NULL" to 0
      if (is.null(N_synapses)) {
        N_synapses <- 0
      }
      print((i * N_cell_groups - N_cell_groups) + j)
      print(N_synapses)
      # add value to synapse list
      synapse_list[[(i * N_cell_groups - N_cell_groups) + j]] <- N_synapses
    }
  }
}

# convert synapse list into a matrix of appropriate dimensions
synapse_matrix <- matrix(unlist(synapse_list), byrow = TRUE, nrow = N_cell_groups)
rownames(synapse_matrix) <- cell_group_attr$cell_group_names
colnames(synapse_matrix) <- cell_group_attr$cell_group_names
synapse_matrix

write.csv(as.data.frame(synapse_matrix), "source_data/Figure15_source_data1.txt")


# graph conversion --------------------------------------------------------

# edge weight filtering on the matrix to remove weak edges
synapse_matrix[synapse_matrix < 5] <- 0
max(synapse_matrix)

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

# calculate node weighted degree
degree <- degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)


# visNetwork plotting -----------------------------------------------------

## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(Conn_graph)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight
Conn_graph.visn$nodes$value <- degree
Conn_graph.visn$nodes$group <- cell_group_attr$type

# hierarchical layout

# level	: Number. Default to undefined. When using the hierarchical layout, the level determines where the node is going to be positioned.
Conn_graph.visn$nodes$level <- cell_group_attr$level
# hierarchical layout
{
  visNet <- visNetwork(Conn_graph.visn$nodes, Conn_graph.visn$edges) %>%
    visHierarchicalLayout(
      levelSeparation = 250,
      direction = "LR",
      sortMethod = "hubsize",
      shakeTowards = "roots"
    ) %>%
    visEdges(
      smooth = list(type = "curvedCW", roundness = 0.2),
      scaling = list(min = 2, max = 12),
      color = list(inherit = TRUE, opacity = 0.5),
      arrows = list(to = list(
        enabled = TRUE,
        scaleFactor = 1.2, type = "arrow"
      ))
    ) %>%
    visNodes(
      borderWidth = 0.3,
      color = list(background = Conn_graph.visn$nodes$color, border = "black"),
      opacity = 0.9,
      shape = "dot",
      font = list(color = "black", size = 44),
      scaling = list(label = list(enabled = TRUE, min = 34, max = 44)),
      level = Conn_graph.visn$nodes$level
    ) %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1, algorithm = "hierarchical", labelOnly = FALSE)) %>%
    visInteraction(
      dragNodes = TRUE, dragView = TRUE,
      zoomView = TRUE, hover = TRUE,
      multiselect = TRUE
    ) %>%
    visGroups(
      groupname = "SN", shape = "ellipse",
      opacity = 1, color = "#E69F00"
    ) %>%
    visGroups(
      groupname = "INsplit", shape = "square",
      opacity = 1, color = "#56B4E9"
    ) %>%
    visGroups(
      groupname = "IN", shape = "square",
      opacity = 1, color = "#CC79A7"
    ) %>%
    visGroups(
      groupname = "MN", shape = "dot",
      opacity = 1, color = "#cccccc"
    ) %>%
    addFontAwesome() %>%
    visLegend(
      addNodes = list(
        list(
          label = "148 synapses      ", shape = "icon",
          icon = list(code = "f2d1", size = 80, color = "#E69F00")
        )
      ),
      useGroups = TRUE, width = 0.08, ncol = 1,
      position = "right", stepX = 50, stepY = 75, zoom = TRUE
    )


  visNet
}

# save as html
saveNetwork(visNet, "pictures/visNetwork_mech_chaeMech_grouped_circuit.html")
webshot2::webshot(
  url = "pictures/visNetwork_mech_chaeMech_grouped_circuit.html",
  file = "pictures/visNetwork_mech_chaeMech_grouped_circuit.png",
  cliprect = c(100, 50, 1000, 480), zoom = 5
)

# plot with ggplot
{
  as.data.frame((synapse_matrix)) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses") %>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE)) %>%
    ggplot(aes(x = postsyn_cell_group, y = presyn_cell_group)) +
    geom_raster(aes(fill = sqrt(synapses))) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 10),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    labs(x = "postsynaptic cell groups", y = "presynaptic cell groups", title = " ") +
    scale_x_discrete(limits = cell_group_attr$cell_group_names) +
    scale_y_discrete(limits = cell_group_attr$cell_group_names) +
    scale_fill_gradientn(colours = c("white", "#0072B2")) +
    geom_text(aes(label = synapses, size = synapses / (synapses + 0.1))) +
    scale_radius(range = c(0, 2)) +
    guides(size = "none")


  # Saving R ggplot with R ggsave Function
  ggsave("pictures/mech_girdle_chaeMech_syn_matrix.png",
    width = 1700,
    height = 1300, limitsize = TRUE,
    units = c("px")
  )

  write.csv2(synapse_matrix, file = "supplements/mech_girdle_chaeMech_synapse_matrix.csv")
}


# assemble figure ---------------------------------------------------------
{
  imgA <- readPNG("pictures/Figure_mec_chaeMech_ventr.png")
  imgB <- readPNG("pictures/Figure_mec_chaeMech_syn.png")
  imgC <- readPNG("pictures/Figure_mec_chaeMech_lat.png")
  imgD <- readPNG("pictures/chaeMech_EM_2.5um.png")
  imgE <- readPNG("pictures/mech_girdle_chaeMech_syn_matrix.png")
  imgF <- readPNG("pictures/visNetwork_mech_chaeMech_grouped_circuit.png")

  # convert png to image panel
  panelA <- ggdraw() + draw_image(imgA, scale = 1) +
    draw_label("chaeMech, ventral",
      x = 0.4, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label(expression(paste("50 ", mu, " m")), x = 0.8, y = 0.05, size = 10)
  
  panelB <- ggdraw() + draw_image(imgB, scale = 1) +
    draw_label("chaeMech, synapses",
      x = 0.45, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("outgoing",
      x = 0.15, y = 0.9, fontfamily = "sans", fontface = "plain",
      color = oranges[8], size = 9, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("incoming",
      x = 0.15, y = 0.85, fontfamily = "sans", fontface = "plain",
      color = blues[4], size = 9, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelC <- ggdraw() + draw_image(imgC, scale = 1) +
    draw_label("chaeMech, lateral",
      x = 0.4, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelD <- ggdraw() + draw_image(imgD, scale = 1) +
    draw_label(expression(paste("2.5 ", mu, "m", sep = "")),
      x = 0.21, y = 0.18, fontfamily = "sans", fontface = "bold",
      color = "white", size = 11)
  panelE <- ggdraw() + draw_image(imgE, scale = 1)
  panelF <- ggdraw() + draw_image(imgF, scale = 1)
}

{
# define layout with textual representation
layout <- "
AAAABBBBCCCCDDDD
EEEEEEEEFFFFFFFF
EEEEEEEEFFFFFFFF
"

Fig_mech_circuits <- panelA + panelB + panelC + panelD +
    panelE + panelF +
    plot_layout(design = layout, heights = c(1, 0.6, 0.6)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12, face = "plain"))


ggsave("Figures/Figure15.pdf",
    limitsize = FALSE,
    units = c("px"), Fig_mech_circuits, width = 2600, height = 1900
  )


ggsave("Figures/Figure15.png",
    limitsize = FALSE,
    units = c("px"), Fig_mech_circuits, width = 2600, height = 1900, bg = "white"
  )
}
