#R/natverse code to generate Figure 16 and Fig 16 suppl 1 for the Platynereis 3d connectome paper
#Gaspar Jekely March 2022


#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")


# load cell clusters ------------------------------------------------------

{
  SNbronto = nlapply(read.neurons.catmaid("^celltype168$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  
  INsplitBronto = nlapply(read.neurons.catmaid("^celltype149$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  
  hPU = nlapply(read.neurons.catmaid("^celltype96$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  INsplitPUh = nlapply(read.neurons.catmaid("^celltype83$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
 
  interparaPM = nlapply(read.neurons.catmaid("^celltype81$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  INsplitPB = nlapply(read.neurons.catmaid("^celltype79$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  INsplitCR = nlapply(read.neurons.catmaid("^celltype73$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  SNblunt = nlapply(read.neurons.catmaid("^celltype148$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INsplitVent = nlapply(read.neurons.catmaid("^celltype157$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  
  CR = nlapply(read.neurons.catmaid("^CRneurons$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  PU = nlapply(read.neurons.catmaid("^PUneurons$", pid=11),
               function(x) smooth_neuron(x, sigma=6000)) 
  PB = nlapply(read.neurons.catmaid("^Biciliated_penetrating_cell$", pid=11),
               function(x) smooth_neuron(x, sigma=6000)) 
  INsplit = nlapply(read.neurons.catmaid("^INsplit$", pid=11),
               function(x) smooth_neuron(x, sigma=6000)) 
  MUSlongV = nlapply(read.neurons.catmaid("^celltype_non_neuronal76$", pid=11),
               function(x) smooth_neuron(x, sigma=6000)) 
  MUStrans = nlapply(read.neurons.catmaid("^celltype_non_neuronal74$", pid=11),
               function(x) smooth_neuron(x, sigma=6000)) 
  
  Loop = nlapply(read.neurons.catmaid("^celltype59$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  MC3cover = nlapply(read.neurons.catmaid("^celltype87$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  MNant = nlapply(read.neurons.catmaid("^celltype19$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  MNbiramous = nlapply(read.neurons.catmaid("^celltype63$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  MNspinning = nlapply(read.neurons.catmaid("^celltype47$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  MNspider_ant = nlapply(read.neurons.catmaid("^celltype61$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  MNspider_post = nlapply(read.neurons.catmaid("^celltype62$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  MNcrab = nlapply(read.neurons.catmaid("^celltype65$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  spinGland = nlapply(read.neurons.catmaid("^celltype_non_neuronal7$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  ciliary_band = nlapply(read.neurons.catmaid("^ciliated cell$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  
  cover_cells = nlapply(read.neurons.catmaid("^celltype_non_neuronal8$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  INCM = nlapply(read.neurons.catmaid("^celltype60$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  MUSac = nlapply(read.neurons.catmaid("^acicular muscle$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  MUSchae = nlapply(read.neurons.catmaid("^chaetal sac muscle$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  MUSlong = nlapply(read.neurons.catmaid("^longitudinal muscle$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  MUSob = nlapply(read.neurons.catmaid("^oblique muscle$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  
  stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))  
  girdle = nlapply(read.neurons.catmaid("^mechanosensory_girdle$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
}


#adjusted background plotting function
plot_background_mech <- function(){plot_background_ventral()
  clear3d()
  nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=0,
         add=T, forceClipregion = F, alpha=0.05,
         col='grey80')
  plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=F, lwd=3,
         add=T, alpha=0.1, col="grey50")
  par3d(zoom=0.58)
}


# plot SN - INsplit cells  ----------

plot_background_mech()
plot3d(CR, soma=T, lwd=2,
       add=T, alpha=0.7, col=Okabe_Ito[1])
plot3d(INsplitCR, soma=T, lwd=3,
       add=T, alpha=0.7, col=Okabe_Ito[5])
texts3d(121000,38000, 24000, text = "CR", col=Okabe_Ito[1], cex = 2.6)
texts3d(33000,140000, 99000, text = "INsplitCR", col=Okabe_Ito[5], cex = 2.6)

plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=0.15, col="grey50")
plot3d(scalebar_50um_ventral, WithConnectors = F, WithNodes = F, soma=F, lwd=5,
       rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
       col="black")
par3d(zoom=0.55)
rgl.snapshot("pictures/Figure_mec_INsplit_CR.png")
close3d()

plot_background_mech()
plot3d(PB, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=0.7, col=Okabe_Ito[6])
plot3d(INsplitPB, WithConnectors = F, WithNodes = F, soma=T, lwd=4,
       rev = FALSE, fixup = F, add=T, alpha=0.9, col=Okabe_Ito[2])
texts3d(44000,140000, 66000, text = "INsplitPB", col=Okabe_Ito[2], cex = 2.6)
texts3d(110000,18000, 24000, text = "PB", col=Okabe_Ito[6], cex = 2.6)

plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=0.15, col="grey50")
par3d(zoom=0.55)
rgl.snapshot("pictures/Figure_mec_INsplit_PB.png")
close3d()

plot_background_mech()
plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=0.15, col="grey50")
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       add=T, alpha=0.3, col="grey50")
plot3d(SNbronto, soma=T, lwd=5, add=T, alpha=1, 
       col=Okabe_Ito[1])
plot3d(INsplitBronto, soma=T, lwd=10, add=T, alpha=1, 
       col=Okabe_Ito[2])
texts3d(120000, 105000, 2000, text = "SNbronto", col=Okabe_Ito[1], cex = 2.6)
texts3d(55000, 145000, 65000, text = "INsplitBronto", col=Okabe_Ito[2], cex = 2.6)
rgl.snapshot("pictures/Figure_mec_INsplit_Bronto.png")
close3d()


plot_background_mech()
plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=0.15, col="grey50")
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       add=T, alpha=0.3, col="grey50")
plot3d(hPU, soma=T, lwd=5, add=T, alpha=1, 
       col=Okabe_Ito[6])
plot3d(INsplitPUh, soma=T, lwd=10, add=T, alpha=1, 
       col=Okabe_Ito[5])
texts3d(105000, 95000, 3000, text = "hPU", col=Okabe_Ito[6], cex = 2.6)
texts3d(45000, 145000, 35000, text = "INsplitPUh", col=Okabe_Ito[5], cex = 2.6)
par3d(zoom=0.55)
rgl.snapshot("pictures/Figure_mec_INsplit_PUh.png")
close3d()


plot_background_mech()
plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=0.15, col="grey50")
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       add=T, alpha=0.3, col="grey50")
plot3d(interparaPM, soma=T, lwd=5, add=T, alpha=1, 
       col=Okabe_Ito[1])
plot3d(INsplitPB, soma=T, lwd=3, add=T, alpha=1, 
       col=Okabe_Ito[5])
texts3d(42000, 125000, 55000, text = "interparaPM", col=Okabe_Ito[1], cex = 2.6)
texts3d(42000, 145000, 95000, text = "INsplitPB", col=Okabe_Ito[5], cex = 2.6)
par3d(zoom=0.55)
rgl.snapshot("pictures/Figure_mec_interparaPM.png")
close3d()


plot_background_mech()
plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       add=T, alpha=0.3, col="grey50")
plot3d(INsplitCR, soma=T, lwd=3, add=T, alpha=0.5, 
       col=Okabe_Ito[1])
plot3d(INsplitPB, soma=T, lwd=2, add=T, alpha=1, 
       col=Okabe_Ito[5])
plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=0.15, col="grey50")
texts3d(36000, 125000, 63000, text = "INsplitCR", col=Okabe_Ito[1], cex = 2.6)
texts3d(42000, 145000, 115000, text = "INsplitPB", col=Okabe_Ito[5], cex = 2.6)
par3d(zoom=0.55)
rgl.snapshot("pictures/Figure_mec_INsplit_CR_PB.png")
close3d()

plot_background_mech()
plot3d(SNblunt, soma=T, lwd=4,
       add=T, alpha=0.7, col=Okabe_Ito[6])
plot3d(INsplitPUh, soma=T, lwd=3,
       add=T, alpha=0.9, col=Okabe_Ito[5])
plot3d(INsplitVent, soma=T, lwd=4,
       add=T, alpha=0.9, col=Okabe_Ito[2])
texts3d(105000,18000, 25000, text = "SNblunt", col=Okabe_Ito[6], cex = 2.6)
texts3d(28000,110000, 42000, text = "INsplitPUh", col=Okabe_Ito[5], cex = 2.6)
texts3d(38000,110000, 122000, text = "INsplitVent", col=Okabe_Ito[2], cex = 2.6)
plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=0.15, col="grey50")
par3d(zoom=0.55)
rgl.snapshot("pictures/Figure_mec_INsplit_SNblunt.png")
close3d()

# individual INsplit cells ----------

par3d(windowRect = c(0, 0, 300, 400))
plot3d(INsplitPUh[2], soma=T, lwd=4, add=T, alpha=1, 
       col=Okabe_Ito[5])
nview3d("ventral")
par3d(zoom=0.74)
rgl.snapshot("pictures/Figure_mec_INsplitPUh2.png")
close3d()

par3d(windowRect = c(0, 0, 300, 400))
plot3d(INsplitBronto[1], soma=T, lwd=4, add=T, alpha=1, 
       col=Okabe_Ito[5])
nview3d("ventral")
par3d(zoom=0.74)
rgl.snapshot("pictures/Figure_mec_INsplitBronto1.png")
close3d()

par3d(windowRect = c(0, 0, 200, 400))
plot3d(INsplitVent[2], soma=T, lwd=4, add=T, alpha=1, 
       col=Okabe_Ito[5])
nview3d("ventral")
par3d(zoom=0.5)
rgl.snapshot("pictures/Figure_mec_INsplitVent2.png")
close3d()

par3d(windowRect = c(0, 0, 200, 400))
plot3d(INsplitPB[1], soma=T, lwd=4, add=T, alpha=1, 
       col=Okabe_Ito[5])
nview3d("ventral")
par3d(zoom=0.5)
rgl.snapshot("pictures/Figure_mec_INsplitPB1.png")
close3d()

par3d(windowRect = c(0, 0, 200, 400))
plot3d(INsplitCR[9], soma=T, lwd=4, add=T, alpha=1, 
       col=Okabe_Ito[5])
nview3d("ventral")
par3d(zoom=0.5)
rgl.snapshot("pictures/Figure_mec_INsplitCR9.png")
close3d()

# get connectivity between mech cell clusters -----------------

{
cell_groups <- list(CR, PB, SNbronto, hPU, interparaPM, SNblunt,
                    INsplitPUh, INsplitCR, INsplitBronto, INsplitPB, INCM, 
                    Loop, MC3cover, MNant, MNbiramous, MNspinning,
                    MNspider_ant, MNspider_post, MNcrab,
                    spinGland, ciliary_band, MUSac, MUSchae, MUSlong, 
                    MUSob, cover_cells)
N_cell_groups <- length(cell_groups)
N_cell_groups

cell_group_attr <- data.frame(
  cell_group_names  = c('CR', 'PB', 'SNbronto', 'hPU', 'interparaPM', 'SNblunt',
                        'INsplitPUh', 'INsplitCR', 'INsplitBronto', 'INsplitPB', 'INCM', 
                        'Loop', 'MC3cover', 'MNant', 'MNbiramous', 'MNspinning',
                        'MNspider_ant', 'MNspider_post', 'MNcrab',
                        'spinGland', 'ciliary_band', 'MUSac', 'MUSchae', 'MUSlong', 
                        'MUSob', 'cover_cells'),
  type = c('SN', 'SN', 'SN', 'SN', 'SN', 'SN',
           'IN', 'IN', 'IN', 'IN', 'IN',
           'MN', 'MN', 'MN', 'MN', 'MN',
           'MN', 'MN', 'MN',
           'effector', 'effector', 'effector', 'effector', 'effector', 
           'effector', 'effector'),
  level = c('1', '1', '1', '1', '1', '1',
            '2', '2', '2', '2', '2',
            '3', '3', '3', '3', '3',
            '3', '3', '3',
            '4', '4', '4', '4','4','4','4')
)

dim(cell_group_attr)

#iterate through cell group neuron lists and get connectivity for all against all
{
  #define empty synapse list with the right dimensions
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
}

#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=N_cell_groups)
rownames(synapse_matrix) <- cell_group_attr$cell_group_names
colnames(synapse_matrix) <- cell_group_attr$cell_group_names
synapse_matrix
write.csv(as.data.frame(synapse_matrix), "source_data/Figure16_source_data1.txt")

#plot with ggplot
{
  as.data.frame((synapse_matrix[1:19, ])) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
    ggplot(aes(x = postsyn_cell_group, y = presyn_cell_group)) +
    geom_raster(aes(fill=sqrt(synapses)))+ 
    theme(
      axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
      axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
      axis.title.x = element_text (size=15),
      axis.title.y = element_text (size=15)
    )+
    labs(x="postsynaptic cell groups",y="presynaptic cell groups",title=" ") +
    scale_x_discrete(limits = cell_group_attr$cell_group_names) +
    scale_y_discrete(limits = rev(cell_group_attr$cell_group_names[1:19]))+
    scale_fill_gradientn(colours=c("white","#0072B2")) +
    geom_text(aes(label=synapses, size=synapses/(synapses+0.1))) +
    scale_radius(range = c(0,2)) +
    guides(size = 'none') 
  
  
  # Saving R ggplot with R ggsave Function
  ggsave("Figures/Figure16_fig_suppl1.pdf", 
         width = nrow(synapse_matrix)/1.3, 
         height = ncol(synapse_matrix)/1.6, limitsize = TRUE, 
         units = c("cm"))
  
  # Saving R ggplot with R ggsave Function
  ggsave("Figures/Figure16_fig_suppl1.png", 
         width = 2200, 
         height = 1400, limitsize = TRUE, 
         units = c("px"))
  
  synapse_matrix[1:19,]
  
  # Save source data
  write.table(synapse_matrix[1:19,], "source_data/Figure16_fig_suppl1_source_data1.txt", sep = "\t")
  synapse_matrix <- read.table("source_data/Figure16_fig_suppl1_source_data1.txt", sep = "\t")
  
}



# graph conversion -----------

#edge weight filtering on the matrix to remove weak edges
synapse_matrix[synapse_matrix < 5] <- 0
max(synapse_matrix)

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

#calculate node weighted degree
degree=degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)

# use visNetwork to plot the network ------------

## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(Conn_graph)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight
Conn_graph.visn$nodes$value = degree
Conn_graph.visn$nodes$group <- cell_group_attr$type

#hierarchical layout
#level	: Number. Default to undefined. When using the hierarchical layout, the level determines where the node is going to be positioned.
Conn_graph.visn$nodes$level <- cell_group_attr$level
#hierarchical layout
{
  visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
    visHierarchicalLayout(levelSeparation=350, 
                          direction='LR',
                          sortMethod='hubsize',
                          shakeTowards='roots') %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0.1),
             scaling=list(min=2, max=18),
             color = list(inherit=TRUE, opacity=0.5),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 1.2, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(background=Conn_graph.visn$nodes$color, border='black'),
             opacity=0.9,
             shape='dot', 
             font=list(color='black', size=44),
             scaling = list(label=list(enabled=TRUE, min=34, max=44)),
             level= Conn_graph.visn$nodes$level) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, algorithm='hierarchical',labelOnly=FALSE)) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    visGroups(groupname = "SN", shape = "ellipse", 
              opacity=0.7, color=Okabe_Ito[1]) %>%
    visGroups(groupname = "IN", shape = "square", 
              opacity=1, color="#56B4E9") %>%
    visGroups(groupname = "MN", shape = "dot", 
              opacity=1, color="#cccccc")  %>%
    visGroups(groupname = "effector", shape = "dot", 
              opacity=1, color="grey60") %>%
    addFontAwesome()
  
  
  visNet
}


#save as html
saveNetwork(visNet, "pictures/visNetwork_mech_INsplit_circuit.html")
webshot2::webshot(url="pictures/visNetwork_mech_INsplit_circuit.html",
                  file="pictures/visNetwork_mech_INsplit_circuit.png",
                  cliprect = c(200,50,700, 480), zoom=5)


}

# assemble figure ----------------------------------------

{
  imgCR <- readPNG("pictures/Figure_mec_INsplit_CR.png")
  imgPB <- readPNG("pictures/Figure_mec_INsplit_PB.png")
  imgM <- readPNG("pictures/Figure_mec_INsplit_Bronto.png")
  imgN <- readPNG("pictures/Figure_mec_INsplit_PUh.png")
  imgO <- readPNG("pictures/Figure_mec_interparaPM.png")
  imgP <- readPNG("pictures/Figure_mec_INsplit_CR_PB.png")
  img_blunt <- readPNG("pictures/Figure_mec_INsplit_SNblunt.png")
  imgQ <- readPNG("pictures/visNetwork_mech_INsplit_circuit.png")
  
  
  imgPUh2 <- readPNG("pictures/Figure_mec_INsplitPUh2.png")
  imgBronto1 <- readPNG("pictures/Figure_mec_INsplitBronto1.png")
  imgVent2 <- readPNG("pictures/Figure_mec_INsplitVent2.png")
  imgPB1 <- readPNG("pictures/Figure_mec_INsplitPB1.png")
  imgCR9 <- readPNG("pictures/Figure_mec_INsplitCR9.png")
  
  
  
  panelCR <- ggdraw() + draw_image(imgCR, scale=1) +
    draw_label('CR', x=0.2, y=0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1) +
    draw_label(expression(paste("50 ", mu, " m")), x = 0.79, y = 0.04, size = 10)
  
  panelPB <- ggdraw() + draw_image(imgPB, scale=1) +
    draw_label('PB', x=0.2, y=0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
    
  panelM <- ggdraw() + draw_image(imgM, scale=1) +
    draw_label('SNbronto', x=0.27, y=0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1) 
  
  panelN <- ggdraw() + draw_image(imgN, scale=1) +
    draw_label('hPU', x=0.18, y=0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  panelO <- ggdraw() + draw_image(imgO, scale=1) +
    draw_label('interparaPM', x=0.33, y=0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  panelP <- ggdraw() + draw_image(imgP, scale=1) +
    draw_label('INsplitCR/PB', x=0.35, y=0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  panel_blunt <- ggdraw() + draw_image(img_blunt, scale=1) +
    draw_label('SNblunt', x=0.27, y=0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  

  panel_PUh2 <- ggdraw() + draw_image(imgPUh2, scale=1) +
    draw_label('INsplitPU', x=0.5, y=0.6, size = 9,
               color=Okabe_Ito[5])
  panel_Bronto1 <- ggdraw() + draw_image(imgBronto1, scale=1) +
    draw_label('INsplitBronto', x=0.5, y=0.7, size = 9,
               color=Okabe_Ito[5])
  panel_Vent2 <- ggdraw() + draw_image(imgVent2, scale=1) +
    draw_label('INsplitVent', x=0.5, y=0.5, size = 9,
               color=Okabe_Ito[5])
  panel_PB1 <- ggdraw() + draw_image(imgPB1, scale=1) +
    draw_label('INsplitPB', x=0.5, y=0.55, size = 9,
               color=Okabe_Ito[5])
  panel_CR9 <- ggdraw() + draw_image(imgCR9, scale=1) +
    draw_label('INsplitCR', x=0.5, y=0.6, size = 9,
               color=Okabe_Ito[5])

#check network max  
network_temp <- read_csv("source_data/Figure16_source_data1.txt")
network_temp %>%
  select(-...1) %>%
  max()

  panelQ <- ggdraw() + draw_image(imgQ, scale=1) +
    draw_label("# of synapses", x = 0.3, y = 0.94, size = 8, hjust = 1) +
    draw_label("6", x = 0.24, y = 0.9, size = 8, hjust = 1) + 
    draw_label("471", x = 0.24, y = 0.86, size = 8, hjust = 1) +
    draw_line(x = c(0.25, 0.3), y = c(0.9, 0.9), size = 0.2, color = 'grey') +
    draw_line(x = c(0.25, 0.3), y = c(0.86, 0.86), size = 1.6, color = 'grey')

  
rm(imgCR, imgPB, imgM, imgN, imgO, imgP, imgPUh2,
   imgBronto1, imgVent2, imgPB1, imgCR9, imgQ)

layout <- "
ABCCCCCCDE
##########
FGHHHIIIMM
FGJJKKLLMM

"
  
Figure16 <-  panelCR + panelPB + panelM + panelN + panelO + 
  panel_blunt  + panelP + panel_PUh2 + panel_Bronto1 +
  panel_Vent2 + panel_CR9 + panel_PB1 + panelQ +
    plot_layout(
      design=layout, heights=c(1, 0.05, 0.5, 0.5),
      widths = c(0.6, 0.6, 0.1,0.1,0.1,0.1,0.1,0.1, 0.6, 0.6)) +
    plot_annotation(tag_levels = list(c('A','B','C','D','E','F','G','H',
                                   '','','','','I'))) & 
    theme(plot.tag = element_text(size = 12, face='plain'))
  
ggsave("Figures/Figure16.png", limitsize = FALSE, 
       units = c("px"), Figure16, width = 3000, height = 1600, bg='white')

ggsave("Figures/Figure16.pdf", limitsize = FALSE, 
         units = c("px"), Figure16, width = 3000, height = 1600)
  
}
