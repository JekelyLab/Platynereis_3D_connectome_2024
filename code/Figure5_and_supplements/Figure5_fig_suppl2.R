#R/natverse code to generate Figure 5 - figure supplement 2 of the Platynereis 3d connectome paper
#Gaspar Jekely 2024

#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

load_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(
    annotation, pid=11), 
    function(x) smooth_neuron(x, sigma=6000))
}

# load PDF neurons ---------------
cMNPDF <- load_neuron("^celltype10$")
SNnuch <- load_neuron("^celltype13$")
INMBPDF <- load_neuron("^celltype184$")
SNantlerPDF <- load_neuron("^celltype26$")
SNPDF_pyg <- load_neuron("^celltype115$")
SNPDF_dc <- load_neuron("^celltype27$")
SN_DLSO3_PDF <- load_neuron("^celltype40$")
SN_DLSO1.2_4 <- load_neuron("^celltype28$")
cMNdc <- load_neuron("^celltype12$")
INdescLuqinPDF <- load_neuron("^celltype175$")
hPU2l_asymPDF <- load_neuron("^celltype99$")

# load PDF partners ------------------

INarc1 <- load_neuron("^celltype55$")
INarc2 <- load_neuron("^celltype56$")
INdecussM <- load_neuron("^celltype199$")
MN3 <- load_neuron("^celltype84$")
MNant <- load_neuron("^celltype19$")
SNpygM <- load_neuron("^celltype170$")
SN_DLSO1.3 <- load_neuron("^celltype163$")
cMNATO <- load_neuron("^celltype11$")
MNspinning <- load_neuron("^celltype47$")
INrope <- load_neuron("^celltype58$")

# load ATO-Leuco neurons -------------

INATOpyg <- load_neuron("^celltype145$")
MNspider_ant <- load_neuron("^celltype61$")
INsplitCRATO <- load_neuron("^celltype74$")
INpreLadderATO <- load_neuron("^celltype154$")

INleucoPU <- load_neuron("^celltype150$")

# load ATO-leuco partners -------------

Ser_tr1 <- load_neuron("^celltype142$")
INrope <- load_neuron("^celltype58$")
MUSac_notP <- load_neuron("^celltype_non_neuronal38$")
MNladder <- load_neuron("^celltype151$")

trochPU <- load_neuron("^celltype91$")
INsplitPB <- load_neuron("^celltype79$")

INtorii <- load_neuron("^celltype191$")
SNtorii <- load_neuron("^celltype187$")
MNwave <- load_neuron("^celltype68$")
MNcommUpL <- load_neuron("^celltype106$")
troch  <- load_neuron("^ciliated cell$")
MUS <- load_neuron("^muscle$")

# retrieve connectivity  for PDF -------------

cell_groups <-  list(
  cMNPDF, SNnuch, INMBPDF, 
  SNantlerPDF, SNPDF_pyg, SNPDF_dc, 
  SN_DLSO3_PDF, cMNdc, 
  INdescLuqinPDF, hPU2l_asymPDF, INarc1, 
  INarc2, INdecussM, MN3, 
  MNant, SNpygM, SN_DLSO1.3, 
  cMNATO, MNspinning, INrope
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
  "cMNPDF", "SNnuch", "INMBPDF", 
  "SNantlerPDF", "SNPDF-pyg", "SNPDF-dc", 
  "SN-DLSO3-PDF", "cMNdc", 
  "INdescLuqinPDF", "hPU2l-asymPDF", "INarc1", 
  "INarc2", "INdecussM", "MN3", 
  "MNant", "SNpygM", "SN-DLSO1.3", 
  "cMNATO", "MNspinning", "INrope"
)
rownames(synapse_matrix) <- as.character(cell_group_names)
colnames(synapse_matrix) <- as.character(cell_group_names)
synapse_matrix

write.csv(synapse_matrix, 
          "source_data/Figure5_fig_suppl2_source_data1.txt"
)

# graph conversion --------------------------------------------------------

synapse_matrix_filt <- synapse_matrix
synapse_matrix_filt[synapse_matrix_filt < 3] <- 0

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix_filt,
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

node_colors <- c("#4D004B", "#4D004B", "#4D004B", 
                 "#4D004B", "#4D004B", "#4D004B", 
                 "#4D004B", "#4D004B", 
                 "#4D004B", "#4D004B", "#CC79A7", 
                 "#CC79A7", "#CC79A7", "#56B4E9", 
                 "#56B4E9", "#E69F00", "#E69F00", 
                 "#56B4E9", "#56B4E9", "#CC79A7"
                 )

Conn_graph.visn$nodes$label

Conn_graph.visn$nodes$group <- list("PDF", "PDF", "PDF", 
                                    "PDF", "PDF", "PDF", 
                                    "PDF", "PDF", 
                                    "PDF", "PDF", "IN", 
                                    "IN", "IN", "MN", 
                                    "MN", "SN", "SN", 
                                    "MN", "MN", "IN")
Conn_graph.visn$nodes$color <- node_colors

#hierarchical layout

Conn_graph.visn$nodes$level <- c("2", "1", "2", 
                                 "1", "1", "1", 
                                 "1",  "2", 
                                 "2", "1", 
                                 "3", 
                                 "3", "3", "4", 
                                 "4", "0", "0", 
                                 "4", "4", "3")
visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout_nicely", type = "full", physics = TRUE) %>%
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
  visOptions(width=1000, height=1000) %>%
  visGroups(groupname = "PDF", shape = "square", opacity=1)
visNet

saveNetwork(visNet, "pictures/visNetwork_PDF.html")
webshot2::webshot(url = "pictures/visNetwork_PDF.html",
                  file = "pictures/visNetwork_PDF.png",
                  vwidth = 950, vheight = 900, #define the size of the browser window
                  cliprect = c(50, 100, 950, 900), zoom = 5, delay = 1)



# retrieve connectivity  for ATO -------------

cell_groups <-  list(
  INATOpyg, MNspider_ant, INsplitCRATO,
  INpreLadderATO, Ser_tr1, INrope, MUS, 
  MNladder,
  INleucoPU, trochPU,
  INsplitPB, INtorii, SNtorii,
  MNwave, MNcommUpL, troch
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
  "INATOpyg", "MNspider-ant", 
  "INsplitCRATO", "INpreLadderATO", 
  "Ser-tr1", "INrope", 
  "muscle", "MNladder",
  "INleucoPU", "trochPU",
  "INsplitPB", "INtorii", "SNtorii",
  "MNwave", "MNcommUpL", "ciliary band"
)
rownames(synapse_matrix) <- as.character(cell_group_names)
colnames(synapse_matrix) <- as.character(cell_group_names)
synapse_matrix

write.csv(synapse_matrix, 
          "source_data/Figure5_fig_suppl2_source_data2.txt"
)

# graph conversion --------------------------------------------------------

synapse_matrix_filt <- synapse_matrix
synapse_matrix_filt[synapse_matrix_filt < 3] <- 0

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix_filt,
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


node_colors <- c("#4D004B", "#4D004B", "#4D004B", 
                 "#4D004B", "#56B4E9", 
                 "#CC79A7", "#cccccc", "#56B4E9", 
                 "#CC79A7", "#E69F00", "#56B4E9",
                 "#56B4E9", "#E69F00", 
                 "#56B4E9", "#56B4E9", "#cccccc"
)

Conn_graph.visn$nodes$label

Conn_graph.visn$nodes$group <- list(
  "ATO", "ATO", "ATO", 
  "ATO", "MN", "IN", 
  "effector", "MN",
  "Leuco", "SN", "IN",  "IN",  "SN",
  "MN", "MN", "effector")

Conn_graph.visn$nodes$color <- node_colors

"INATOpyg", "MNspider-ant", 
"INsplitCRATO", "INpreLadderATO", 
"Ser-tr1", "INrope", 
"muscle", "MNladder",
"INleucoPU", "trochPU",
"INsplitPB", "INtorii", "SNtorii",
"MNwave", "MNcommUpL", "Troch"

#hierarchical layout
Conn_graph.visn$nodes$level <- c(
  "3", "4", "2", 
  "3", "4", "3", 
  "5",  "4",
  "2", "0", "1",
  "1", "0",
  "4", "4", "5")

visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout_nicely", type = "full", physics = TRUE) %>%
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
  visOptions(width=1000, height=1000) %>%
  visGroups(groupname = "ATO", shape = "square", opacity=1) %>%
  visGroups(groupname = "Leuco", shape = "square", opacity=1)
visNet

saveNetwork(visNet, "pictures/visNetwork_ATO_Leuco.html")
webshot2::webshot(url = "pictures/visNetwork_ATO_Leuco.html",
                  file = "pictures/visNetwork_ATO_Leuco.png",
                  vwidth = 950, vheight = 500, #define the size of the browser window
                  cliprect = c(50, 300, 950, 500), zoom = 5, delay = 1)


# assemble figure ----------------------

panel_PDF <- ggdraw() + 
  draw_image(readPNG("pictures/visNetwork_PDF.png")) + 
  draw_label(
    "circuit of PDF+ neurons", x = 0.3, y = 0.99, 
    color = "black", size = 11
  )

panel_ATO <- ggdraw() + 
  draw_image(readPNG("pictures/visNetwork_ATO_Leuco.png")) + 
  draw_label(
    "circuit of ATO+ and Leuco+ neurons", x = 0.4, y = 0.99, 
    color = "black", size = 11
  )


layout = "
AAAABBBBB
"

Figure5_fig_suppl2 <- panel_PDF + panel_ATO+
  plot_layout(
    design = layout, 
    heights = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure5_fig_suppl2.png",
       limitsize = FALSE,
       units = c("px"), Figure5_fig_suppl2,
       width = 3000, height = 1200, bg = "white"
)

ggsave("Figures/Figure5_fig_suppl2.pdf",
       limitsize = FALSE,
       units = c("px"), Figure5_fig_suppl2, 
       width = 3000, height = 1200
)
