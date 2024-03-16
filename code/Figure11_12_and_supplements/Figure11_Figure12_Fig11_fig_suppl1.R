#Platynereis 3d larva connectome paper Figure on segment 0 and segment 1
#Gaspar Jekely

source("code/Natverse_functions_and_conn.R")

# load all cell clusters by segment ---------------------------------------

{
  head = nlapply(read.neurons.catmaid("^episphere$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  sg0 = nlapply(read.neurons.catmaid("^segment_0$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  sg1 = nlapply(read.neurons.catmaid("^segment_1$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  sg2 = nlapply(read.neurons.catmaid("^segment_2$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  sg3 = nlapply(read.neurons.catmaid("^segment_3$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  pyg = nlapply(read.neurons.catmaid("^pygidium$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
}

#color scheme for segments
segmental_colors <- brewer.pal(6, 'Paired')
pie(rep(1,6),col=segmental_colors,segmental_colors)

# plot all cells in body with segmental coloring --------------------------

{
plot_background_ventral()
  
plot3d(head, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         rev = FALSE, fixup = F, add=T, alpha=0.3, col=segmental_colors[1])
plot3d(sg0, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         rev = FALSE, fixup = F, add=T, alpha=0.3, col=segmental_colors[2])
plot3d(sg1, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         rev = FALSE, fixup = F, add=T, alpha=0.3, col=segmental_colors[3])
plot3d(sg2, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         rev = FALSE, fixup = F, add=T, alpha=0.3, col=segmental_colors[4])
plot3d(sg3, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         rev = FALSE, fixup = F, add=T, alpha=0.3, col=segmental_colors[5])
plot3d(pyg, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         rev = FALSE, fixup = F, add=T, alpha=0.3, col=segmental_colors[6])
  
rgl.snapshot("pictures/Figure11_seg_seg_color.png")
close3d()
}


# print graph with cells coloured by segment ----------------------------------

# read the connectome visNetwork R data file
conn_graph.visn <- readRDS("supplements/connectome_graph.rds", refhook = NULL)

coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol = 2)

coords_rotated <- autoimage::rotate(
  coords,
  -pi / 3,
  pivot = c(0, 0)
)

#flip along x axis
coords_rotated[,1] = coords_rotated[,1] * -1

{
  # overwrite group value (partition) with segment value (for colouring)
  conn_graph.visn$nodes$group <- as.character(conn_graph.visn$nodes$segment)
  
  # rename segments
  conn_graph.visn$nodes$group <- gsub("episphere", "head", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group <- gsub("segment_", "sg", conn_graph.visn$nodes$group)
  conn_graph.visn$nodes$group
  
  # remove colour info (which takes precedence over group colour)
  conn_graph.visn$nodes$color <- c()
  segmental_colors <- rev(brewer.pal(6, "Paired"))
  
  visNet_segment <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
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
      groupname = "head", color = segmental_colors[1], shape = "dot",
      opacity = 1, size = 25
    ) %>%
    visGroups(
      groupname = "sg0", shape = "dot",
      opacity = 1, size = 25, color = segmental_colors[2]
    ) %>%
    visGroups(
      groupname = "sg1", shape = "dot",
      opacity = 1, size = 25, color = segmental_colors[3]
    ) %>%
    visGroups(
      groupname = "sg2", shape = "dot",
      opacity = 1, size = 25, color = segmental_colors[4]
    ) %>%
    visGroups(
      groupname = "sg3", shape = "dot",
      opacity = 1, size = 25, color = segmental_colors[5]
    ) %>%
    visGroups(
      groupname = "pygidium", shape = "dot",
      opacity = 1, size = 30, color = segmental_colors[6]
    ) %>%
    visGroups(
      groupname = "fragment", shape = "dot",
      opacity = 1, size = 5, color = "#aaaaaa"
    ) %>%
    addFontAwesome() %>%
    visLegend(
      useGroups = TRUE,
      width = 0.2,
      ncol = 1,
      position = "left"
    )
  
  # save as html
  saveNetwork(visNet_segment, "pictures/Full_connectome_segments.html", selfcontained = TRUE)
  webshot::webshot(
    url = "pictures/Full_connectome_segments.html",
    file = "pictures/Full_connectome_segments.png",
    vwidth = 1000, vheight = 1000, # define the size of the browser window
    cliprect = c(60, 140, 800, 680), zoom = 10, delay = 1
  )
}


#read sg0 cell clusters
{
SNcirri <-  nlapply(read.neurons.catmaid("anterior pair of tentacular cirri", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
SNstomod <-nlapply(read.neurons.catmaid("SN stomodeum", pid=11),
                                function(x) smooth_neuron(x, sigma=6000)) 
stomodeum <-  nlapply(read.neurons.catmaid("stomodeum", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
SNstomod <-nlapply(read.neurons.catmaid("SN stomodeum", pid=11),
                   function(x) smooth_neuron(x, sigma=6000)) 
SNsg0_to_head <-nlapply(read.neurons.catmaid("SN_sg0_to_head", pid=11),
                   function(x) smooth_neuron(x, sigma=6000)) 
metatroch <-nlapply(read.neurons.catmaid("celltype_non_neuronal5", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
}

#plot sg0 cell clusters
{
plot_background()
plot3d(scalebar_50um_anterior, WithConnectors = F, WithNodes = F, soma=F, lwd=5,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="black")
#plot sensory cells without soma
plot3d(SNcirri, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
           rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[1])
plot3d(SNstomod, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[2])
plot3d(SNsg0_to_head, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[7])
plot3d(metatroch, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="grey80")
plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="grey80")
}
rgl.snapshot("pictures/Figure_sg0.png")
close3d()

# Segment 1 ---------------------------------------------------------------

#read cell groups in sg1
{
doCRunp <-nlapply(read.neurons.catmaid("celltype102", pid=11),
                        function(x) smooth_neuron(x, sigma=6000)) 
MS5 <-nlapply(read.neurons.catmaid("celltype38", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 

INCM <-nlapply(read.neurons.catmaid("celltype60", pid=11),
                  function(x) smooth_neuron(x, sigma=6000)) 
INsplitPB_RF_Ya <-nlapply(read.neurons.catmaid("celltype77", pid=11),
              function(x) smooth_neuron(x, sigma=6000)) 
INsplitPBant <-nlapply(read.neurons.catmaid("celltype78", pid=11),
               function(x) smooth_neuron(x, sigma=6000)) 
INsplitBronto <-nlapply(read.neurons.catmaid("celltype149", pid=11),
                          function(x) smooth_neuron(x, sigma=6000)) 
INleucoPU <-nlapply(read.neurons.catmaid("celltype150", pid=11),
                       function(x) smooth_neuron(x, sigma=6000)) 
INcomm_DownL <-nlapply(read.neurons.catmaid("celltype158", pid=11),
                        function(x) smooth_neuron(x, sigma=6000)) 


MC3cover <-nlapply(read.neurons.catmaid("celltype87", pid=11),
                          function(x) smooth_neuron(x, sigma=6000)) 
Ser_tr1 <-nlapply(read.neurons.catmaid("celltype142", pid=11),
                       function(x) smooth_neuron(x, sigma=6000)) 
Loop <-nlapply(read.neurons.catmaid("celltype59", pid=11),
                        function(x) smooth_neuron(x, sigma=6000)) 
MNsmile <-nlapply(read.neurons.catmaid("celltype182", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
MNladder <-nlapply(read.neurons.catmaid("celltype151", pid=11),
                       function(x) smooth_neuron(x, sigma=6000)) 
}

#read all cell types per segment
{
sg0_celltypes <-nlapply(read.neurons.catmaid(skids_by_2annotations(
  "segment_0", "celltype"), 
  pid=11), function(x) smooth_neuron(x, sigma=6000)) 
sg1_celltypes <-nlapply(read.neurons.catmaid(skids_by_2annotations(
  "segment_1", "celltype"), 
  pid=11), function(x) smooth_neuron(x, sigma=6000)) 
sg2_celltypes <-nlapply(read.neurons.catmaid(skids_by_2annotations(
  "segment_2", "celltype"), 
  pid=11), function(x) smooth_neuron(x, sigma=6000)) 
sg3_celltypes <-nlapply(read.neurons.catmaid(skids_by_2annotations(
  "segment_3", "celltype"), 
  pid=11), function(x) smooth_neuron(x, sigma=6000)) 
pyg_celltypes <-nlapply(read.neurons.catmaid(skids_by_2annotations(
  "pygidium", "celltype"), 
  pid=11), function(x) smooth_neuron(x, sigma=6000)) 
}

#read MNspinning and partners
{
MNspinning <- nlapply(read.neurons.catmaid("celltype47", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
INsplitSpin <- nlapply(read.neurons.catmaid("celltype202", pid=11),
                        function(x) smooth_neuron(x, sigma=6000)) 
chaeMech <- nlapply(read.neurons.catmaid("celltype71", pid=11),
                      function(x) smooth_neuron(x, sigma=6000)) 
INsplitBronto <- nlapply(read.neurons.catmaid("celltype149", pid=11),
                      function(x) smooth_neuron(x, sigma=6000)) 
INsplitPB <- nlapply(read.neurons.catmaid("celltype79", pid=11),
                         function(x) smooth_neuron(x, sigma=6000)) 
PB <- nlapply(read.neurons.catmaid("celltype80", pid=11),
                     function(x) smooth_neuron(x, sigma=6000)) 
spinGland <- nlapply(read.neurons.catmaid("celltype_non_neuronal7", pid=11),
              function(x) smooth_neuron(x, sigma=6000)) 

}

#read pyg cellgroups and partners
{
SNpygM <- nlapply(read.neurons.catmaid("celltype170", pid=11),
                     function(x) smooth_neuron(x, sigma=6000)) 
SNPDFpyg <- nlapply(read.neurons.catmaid("celltype115", pid=11),
                   function(x) smooth_neuron(x, sigma=6000)) 
cioMNcover <- nlapply(read.neurons.catmaid("celltype82", pid=11),
                   function(x) smooth_neuron(x, sigma=6000)) 
INascpyg <- nlapply(read.neurons.catmaid("celltype147", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
INdescLuqinPDF <- nlapply(read.neurons.catmaid("celltype175", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
INATOpyg <- nlapply(read.neurons.catmaid("celltype145", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
INFVapyg <- nlapply(read.neurons.catmaid("celltype72", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
MNladder <- nlapply(read.neurons.catmaid("celltype151", pid=11),
                   function(x) smooth_neuron(x, sigma=6000)) 
chaeMech <- nlapply(read.neurons.catmaid("celltype71", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
covercell <-nlapply(read.neurons.catmaid("celltype_non_neuronal8", pid=11),
                   function(x) smooth_neuron(x, sigma=6000)) 
pyg_cirrus_CR <- nlapply(read.neurons.catmaid("celltype101", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
pygPBunp <- nlapply(read.neurons.catmaid("celltype54", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
SerTr1 <- nlapply(read.neurons.catmaid("celltype142", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
MNring <- nlapply(read.neurons.catmaid("celltype64", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
INsplitPB <- nlapply(read.neurons.catmaid("celltype79", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
}

#read synaptic inputs and targets
{
MNsmile_mus <-nlapply(read.neurons.catmaid("MNsmile_mus", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
MNladder_mus <-nlapply(read.neurons.catmaid("MNladder_mus", pid=11),
                function(x) smooth_neuron(x, sigma=6000)) 

Loop_Ser_tr_ciliary_targets <-nlapply(read.neurons.catmaid("Loop_Ser_tr_ciliary_targets", pid=11),
                                      function(x) smooth_neuron(x, sigma=6000)) 
covercell_MC3target <- nlapply(read.neurons.catmaid("covercell_MC3target", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))


SNblunt <-nlapply(read.neurons.catmaid("celltype148", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
SNbronto <-nlapply(read.neurons.catmaid("celltype168", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
INMC3 <-nlapply(read.neurons.catmaid("celltype75", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
headCR <-nlapply(read.neurons.catmaid("celltype100", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
MNant <-nlapply(read.neurons.catmaid("celltype19", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
INUturn <-nlapply(read.neurons.catmaid("celltype144", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
SN47Ach <-nlapply(read.neurons.catmaid("celltype52", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
INrope <-nlapply(read.neurons.catmaid("celltype58", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
INpreLadderATO <-nlapply(read.neurons.catmaid("celltype154", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
INasc_pyg <-nlapply(read.neurons.catmaid("celltype147", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
chaeMech <-nlapply(read.neurons.catmaid("celltype71", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
MNspinning <-nlapply(read.neurons.catmaid("celltype47", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
MNspider <-nlapply(read.neurons.catmaid("celltype61", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
VentralParaPU <-nlapply(read.neurons.catmaid("celltype91", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))

INsplitPB <-nlapply(read.neurons.catmaid("celltype79", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
INtorii <-nlapply(read.neurons.catmaid("celltype191", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
MNwave <-nlapply(read.neurons.catmaid("celltype68", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
MNmouth <-nlapply(read.neurons.catmaid("celltype106", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))

}

#plot sg1 cMN and ciliated targets
{
plot_background_ventral2()
plot3d(Ser_tr1, WithConnectors = F, WithNodes = F, soma=T, lwd=4,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[c(2,5)])
plot3d(Loop, WithConnectors = F, WithNodes = F, soma=T, lwd=c(3,5),
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[c(6,1)])
plot3d(Loop_Ser_tr_ciliary_targets, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=0.3, col="grey80") 
texts3d(63000, 145000, 48000, text = "Loop", col='black', cex = 3)
texts3d(82000, 145000, 76000, text = "Ser-tr", col='black', cex = 3)
plot3d(sg1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=0.1, col='grey90') 
plot3d(scalebar_50um_ventral, lwd=3,
       add=T, alpha=1, col='black') 
par3d(zoom=0.53)
rgl.snapshot("pictures/Figure_sg1_Loop_Ser_trochs.png")
close3d()
}

#plot sg1 MN and targets
{
plot_background_ventral2()
plot3d(MNsmile, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[8])
plot3d(MNladder, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[7])
plot3d(MNladder_mus, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=0.5, col=Okabe_Ito[1])
plot3d(MC3cover, WithConnectors = F, WithNodes = F, soma=T, lwd=c(5, 4),
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[c(1, 6)])
plot3d(covercell_MC3target, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T,  alpha=0.5, col="grey80") 
texts3d(81000, 145000, 49000, text = "MNladder", col='black', cex = 3)
texts3d(42000, 145000, 61000, text = "MC3cover", col='black', cex = 3)
texts3d(49000, 145000, 71000, text = "MNsmile", col='black', cex = 3)
plot3d(sg1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=0.1, col='grey90') 
rgl.snapshot("pictures/Figure_sg1_MN_targets.png")
close3d()
}

#plot sg1 INs
{
plot_background_ventral2()
plot3d(INCM, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[8])
plot3d(INsplitPB_RF_Ya, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[2])
plot3d(INsplitPBant, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[7])
plot3d(INsplitBronto, WithConnectors = F, WithNodes = F, soma=T, lwd=7,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[4])
plot3d(INleucoPU, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[5])
plot3d(INcomm_DownL, WithConnectors = F, WithNodes = F, soma=T, lwd=4,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[6])
plot3d(sg1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=0.1, col='grey90') 
rgl.snapshot("pictures/Figure_sg1_IN.png")
close3d()
}

#plot segmental cell types
{
plot_background_ventral2()
plot3d(sg1_celltypes, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=1, col=segmental_colors[3])
rgl.snapshot("pictures/Figure11_sg1_celltypes.png")
close3d()

plot_background_ventral2()
plot3d(sg2_celltypes, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=1, col=segmental_colors[4])
rgl.snapshot("pictures/Figure11_sg2_celltypes.png")
close3d()

plot_background_ventral2()
plot3d(sg3_celltypes, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=1, col=segmental_colors[5])
rgl.snapshot("pictures/Figure11_sg3_celltypes.png")
close3d()

plot_background_ventral2()
plot3d(pyg_celltypes, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
       rev = FALSE, fixup = F, add=T, alpha=1, col=segmental_colors[6])
par3d(zoom=0.45)
rgl.snapshot("pictures/Figure11_pyg_celltypes.png")
close3d()
}

# get connectivity between sg1-specific cell types  --------

{
#these are the sg1 cell groups and their targets
cell_groups <-  list(MC3cover, covercell_MC3target, 
                     MNladder, MNladder_mus, MNsmile, MNsmile_mus, INsplitBronto, INsplitPB_RF_Ya,
                     INsplitPBant, INCM, INcomm_DownL, INleucoPU,
                     SNblunt, SNbronto, INMC3, headCR,
                     MNant, INrope, INpreLadderATO, INasc_pyg,
                     chaeMech, MNspinning, MNspider, VentralParaPU, INsplitPB,
                     INtorii, MNwave)
N_cell_groups <- length(cell_groups)

#iterate through cell group neuron lists and get connectivity for all agains all
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

#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=N_cell_groups)

cell_group_names <-  list("MC3cover (sg1)", "covercell_MC3target", "MNladder (sg1)","MNladder_mus",
                          "MNsmile (sg1)","MNsmile_mus","INsplitBronto (sg1)","INsplitPB_RF_Ya (sg1)",
                          "INsplitPBant (sg1)","INCM (sg1)","INcomm_DownL (sg1)","INleucoPU (sg1)",
                          "SNblunt","SNbronto","INMC3","headCR", "MNant","INrope","INpreLadderATO",
                          "INasc_pyg", "chaeMech","MNspinning","MNspider","VentralParaPU",
                          "INsplitPB","INtorii","MNwave")
rownames(synapse_matrix) <- as.character(cell_group_names)
colnames(synapse_matrix) <- as.character(cell_group_names)
synapse_matrix
write.csv(as.data.frame(synapse_matrix), "source_data/Figure12_source_data1.txt")

#plot with ggplot
as.data.frame((synapse_matrix)) %>%
  rownames_to_column(var = "presyn_cell_group") %>%
  pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
  group_by(postsyn_cell_group) %>%
  mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
  ggplot() +
  geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction), 
             stroke = 0) + 
  scale_x_discrete(limits = (c("MC3cover (sg1)", "covercell_MC3target", "MNladder (sg1)","MNladder_mus",
                               "MNsmile (sg1)","MNsmile_mus","INsplitBronto (sg1)","INsplitPB_RF_Ya (sg1)",
                               "INsplitPBant (sg1)","INCM (sg1)","INcomm_DownL (sg1)","INleucoPU (sg1)",
                               "SNblunt","SNbronto","INMC3","headCR", "MNant","INrope","INpreLadderATO",
                               "INasc_pyg", "chaeMech","MNspinning","MNspider","VentralParaPU",
                               "INsplitPB","INtorii","MNwave"))) +
  scale_y_discrete(limits = rev(c("MC3cover (sg1)", "covercell_MC3target", "MNladder (sg1)","MNladder_mus",
                                  "MNsmile (sg1)","MNsmile_mus","INsplitBronto (sg1)","INsplitPB_RF_Ya (sg1)",
                                  "INsplitPBant (sg1)","INCM (sg1)","INcomm_DownL (sg1)","INleucoPU (sg1)",
                                  "SNblunt","SNbronto","INMC3","headCR", "MNant","INrope","INpreLadderATO",
                                  "INasc_pyg", "chaeMech","MNspinning","MNspider","VentralParaPU",
                                  "INsplitPB","INtorii","MNwave"))) +
  #  coord_flip() +
  theme(
    axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
    axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
    axis.title.x = element_text (size=15),
    axis.title.y = element_text (size=15)
  ) +
  labs(x="postsynaptic cell groups",y="presynaptic cell groups",title=" ") +
  scale_size_area(max_size=5) +
  guides(color = 'legend') +
  scale_colour_gradient2(
    low = "#0072B2",
    mid = "#D55E00",
    high ="#D55E00",
    midpoint = 0.5,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  )+
  #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
  theme(panel.background = element_rect(fill = "grey98", color = "white"))

# Saving R ggplot with R ggsave Function
ggsave("pictures/Head_cellgroups_syn_matrix.pdf", 
       width = nrow(synapse_matrix)/1.3, 
       height = ncol(synapse_matrix)/1.6, limitsize = TRUE, 
       units = c("cm"))

# Saving R ggplot with R ggsave Function
ggsave("pictures/Head_cellgroups_syn_matrix.png", 
       width = 1700, 
       height = 1300, limitsize = TRUE, 
       units = c("px"))
}

# graph conversion --------------------------------------------------------
{
#edge weight filtering on the matrix to remove weak edges
synapse_matrix[synapse_matrix < 5] <- 0

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

#calculate node weighted degree
degree=degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)
}

# use visNetwork to plot the network --------------------------------------
{
## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(Conn_graph)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight

Conn_graph.visn$nodes$value = degree
Conn_graph.visn$nodes
node_colors <- c("#56B4E9", "#cccccc","#56B4E9","#cccccc","#56B4E9","#cccccc","#CC79A7","#CC79A7",
                "#CC79A7","#CC79A7","#CC79A7","#CC79A7","#E69F00","#E69F00","#CC79A7","#E69F00",
                "#56B4E9","#CC79A7","#CC79A7","#CC79A7","#E69F00","#56B4E9","#56B4E9","#E69F00",
                "#CC79A7","#CC79A7","#56B4E9")
length(node_colors)
length(Conn_graph.visn$nodes$id)
Conn_graph.visn$nodes$label
Conn_graph.visn$nodes$group <- list("sg1", "effector", 
                                    "sg1","effector","sg1","effector","sg1","sg1",
                                    "sg1","sg1","sg1","sg1",
                                    "head","head","head","head",
                                    "head","head","INpreLadderATO","INasc_pyg",
                                    "chaeMech","MNspinning","MNspider","VentralParaPU",
                                    "INsplitPB","head","MNwave")
Conn_graph.visn$nodes$color <- node_colors

("MC3cover (sg1)", "covercell_MC3target", "MNladder (sg1)","MNladder_mus",
  "MNsmile (sg1)","MNsmile_mus","INsplitBronto (sg1)","INsplitPB_RF_Ya (sg1)",
  "INsplitPBant (sg1)","INCM (sg1)","INcomm_DownL (sg1)","INleucoPU (sg1)",
  "SNblunt","SNbronto","INMC3","headCR", 
  "MNant","INrope","INpreLadderATO","INasc_pyg", 
  "chaeMech","MNspinning","MNspider","VentralParaPU",
  "INsplitPB","INtorii","MNwave")

#hierarchical layout
#level	: Number. Default to undefined. When using the hierarchical layout, the level determines where the node is going to be positioned.
Conn_graph.visn$nodes$level <- c("2", "6", "2","3", 
                                 "2","3","1","1",
                                 "0","0","0","1",
                                 "6","6","5","6",
                                 "5","5","3","4",
                                 "4","3","3","4",
                                 "4","5","4")

visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 42, type="square") %>%
  visHierarchicalLayout(levelSeparation=220, 
           nodeSpacing=120,
           direction='LR',
           sortMethod='hubsize',
           shakeTowards='roots') %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0.1),
           scaling=list(min=2, max=12),
           color = list(inherit=TRUE, opacity=0.7),
           arrows = list(to = list(enabled = TRUE, 
                                   scaleFactor = 1.2, type = 'arrow'))) %>%
  visNodes(borderWidth=0.3, 
           color = list(background=Conn_graph.visn$nodes$color, border='black'),
           opacity=0.9,
           shape='dot', 
           font=list(color='black', size=24),
           scaling = list(label=list(enabled=TRUE, min=23, max=28)),
           level= Conn_graph.visn$nodes$level) %>%
  visOptions(highlightNearest = TRUE, width = 1800, height = 800) %>%
  visInteraction(navigationButtons = FALSE,
           dragNodes = TRUE, dragView = FALSE,
           zoomView = TRUE) %>%
  visGroups(groupname = "sg1", shape = "square", opacity=1)
visNet

saveNetwork(visNet, "pictures/visNetwork_sg1_specific_input_output.html")
webshot2::webshot(url = "pictures/visNetwork_sg1_specific_input_output.html",
                  file = "pictures/visNetwork_sg1_specific_input_output.png",
                  vwidth = 1800, vheight = 800, #define the size of the browser window
                  cliprect = c(160, 40, 1560, 740), zoom = 5, delay = 3)
}


# retrieve segment-to-segment (SN IN MN effector) connectivity --------
{
  #retrieve all annotations for the cells in the connectome
  connectome_cells <- catmaid_get_annotations_for_skeletons("annotation:^connectome_with_soma$", pid = 11)
  
  #define the six body regions, matching the catmaid annotations
  body_regions <- c('episphere','segment_0', 'segment_1', 'segment_2', 'segment_3', 'pygidium')
  #define the cell categories, matching the catmaid annotations
  cell_category <- c('Sensory neuron','interneuron', 'motorneuron', 'effector')
  
  list_position=0; synapse_list_seg=list(); matrix_lists=list()
  #these iterated loops will query the connectome_cells data based on annotations and retrieve the connectivity between sets of cells (defined by their skids)
  for (i in c(1:6)){  #iterate through the six body regions
    body_region_i = body_regions[i]
    cycle=0
    for (j in c(1:4)){  #iterate through the cell categories
      cell_category_j = cell_category[j]
      print (body_regions[i]); print (cell_category[j])
      presyn_skids_list1 <- list(); presyn_skids_list2 <- list(); counter_k1=0; counter_k2=0
      
      for (k in c(1:length(connectome_cells$annotation))){
        if (connectome_cells$annotation[k] == body_region_i)
        {#here we collect the presyn skids that match the first annotation in the nested loops
          counter_k1 <- counter_k1+1; presyn_skids_list1[[counter_k1]] <- connectome_cells$skid[k]}
        if (connectome_cells$annotation[k] == cell_category_j)
        {#here we collect the presyn skids that match the second annotation in the nested loops
          counter_k2 <- counter_k2+1; presyn_skids_list2[[counter_k2]] <- connectome_cells$skid[k]}
        #find skids that are shared between the two lists
      }
      presyn_skids <- intersect(presyn_skids_list1,presyn_skids_list2)
      
      #we do these nested iterations to get the postsyn skids per category
      for (l in c(1:6)){  #iterate through the six body regions
        body_region_l = body_regions[l]
        for (m in c(1:4)){  #iterate through the cell categories
          cell_category_m = cell_category[m]
          postsyn_skids_list1 <- list(); postsyn_skids_list2 <- list(); counter_n1=0; counter_n2=0  
          
          for (n in c(1:length(connectome_cells$annotation))){
            if (connectome_cells$annotation[n] == body_region_l)
            {#here we collect the presyn skids that match the first annotation in the nested loops
              counter_n1 <- counter_n1+1; postsyn_skids_list1[[counter_n1]] <- connectome_cells$skid[n]}
            if (connectome_cells$annotation[n] == cell_category_m)
            {#here we collect the presyn skids that match the second annotation in the nested loops
              counter_n2 <- counter_n2+1; postsyn_skids_list2[[counter_n2]] <- connectome_cells$skid[n]}
            #find skids that are shared between the two lists
          }
          postsyn_skids <- intersect(postsyn_skids_list1,postsyn_skids_list2)
          cellgroup_conn=list()
          # get connectors between neurons of interest
          if (length(presyn_skids)==0 | length(postsyn_skids)==0) {
            N_synapses <- 0
          } else {
            cellgroup_conn = catmaid_get_connectors_between(pre=presyn_skids, post=postsyn_skids, pid=11)
            N_synapses=nrow(cellgroup_conn)
            if(length(cellgroup_conn) == 0) {N_synapses <- 0}
          }
          list_position <- list_position + 1; print (list_position); print(N_synapses)
          synapse_list_seg [[list_position]] <- N_synapses
          
        }
      }
    }
  }
  
  #convert synapse list into a matrix of appropriate dimensions
  synapse_matrix_seg = matrix(unlist(synapse_list_seg), byrow=TRUE, nrow=24 )
  
  #define the names for plotting of the six body regions and cell categories
  body_regions_name <- c('head','sg0', 'sg1', 'sg2', 'sg3', 'pyg')
  #define the cell categories, matching the catmaid annotations
  cell_category_name <- c('SN','IN', 'MN', 'eff')
  
  #we make a cell group name list 
  cell_group_names=list();counter=0
  for (i in c(1:6)){  #iterate through the six body regions
    body_region_i = body_regions_name[i]
    for (j in c(1:4)){  #iterate through the cell categories
      counter <- counter+1
      cell_group_names[[counter]] <- paste(cell_category_name[j], body_regions_name[i], sep="-")
    }}  
  
  synapse_matrix_seg = as.data.frame(synapse_matrix_seg)
  
  #assign column names to matrix
  synapse_matrix_seg=setNames(synapse_matrix_seg, as.character(cell_group_names))
  
  #assign row names to matrix
  rownames(synapse_matrix_seg) <- as.character(cell_group_names)
  synapse_matrix_seg = as.matrix(synapse_matrix_seg)
}


# segment-to-segment matrix -----------------------------------------------
{
#plot with ggplot
synapse_matrix_seg_plot <- as.data.frame((synapse_matrix_seg)) %>%
  rownames_to_column(var = "presyn_cell_group") %>%
  pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
  group_by(postsyn_cell_group) %>%
  mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
  ggplot() +
  geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = synapses, color = synapse_fraction), 
             stroke = 0) + 
  scale_x_discrete(limits = as.character(cell_group_names),
                   position = "top") +
  scale_y_discrete(limits = rev(as.character(cell_group_names)) ) +
  #  coord_flip() +
  theme(
    axis.text.x = element_text (
      angle = 90, hjust = 0, vjust = 0.5, size = 8.5,
      color = c(rep(segmental_colors[1],4),
                rep(segmental_colors[2],4),
                rep(segmental_colors[3],4),
                rep(segmental_colors[4],4),
                rep(segmental_colors[5],4),
                rep(segmental_colors[6],4) )
      ), 
    axis.text.y = element_text (
      angle = 0, hjust = 1, vjust = 0.5, size = 8.5,
      color = rev(c(rep(segmental_colors[1],4),
                rep(segmental_colors[2],4),
                rep(segmental_colors[3],4),
                rep(segmental_colors[4],4),
                rep(segmental_colors[5],4),
                rep(segmental_colors[6],4) ))
      ),
    axis.title = element_text (size=11)
  ) +
  labs(
    x = "postsynaptic cell groups", 
    y = "presynaptic cell groups", 
    title= " ",
    color = "fraction \nsynapses",
    size = "# synapses") +
  scale_size_area(max_size=4.2) +
  guides(color = 'legend') +
  scale_colour_gradient2(
    low = "#0072B2",
    mid = "#D55E00",
    high ="#D55E00",
    midpoint = 0.5,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
  theme(panel.background = element_rect(fill = "grey98", color = "white"))

synapse_matrix_seg_plot

write.table(
  synapse_matrix_seg, 
  "source_data/Figure11_source_data1.txt", sep = "\t"
  )

}

#seg to seg matrix - combined by type of cell
{
synapses <- c(); c <- 0
for(i in c(1, 5, 9, 13, 17, 21)){
  for(j in c(1, 5, 9, 13, 17, 21)){
    c <- c+1
    synapses[c] <- sum(synapse_matrix_seg[i:(i+3), j:(j+3)])
    }
}

combined_seg_matrix <- tibble(
  synapses,
  presyn = c(rep("head", 6), rep("sg0", 6),
                                rep("sg1", 6), rep("sg2", 6),
                                rep("sg3", 6), rep("pyg", 6)
                                ),
  postsyn = c(rep(c("head", "sg0", "sg1", "sg2", "sg3", "pyg"), 6))
)

combined_seg_matrix_plot <- combined_seg_matrix %>%
  ggplot() +
  geom_text(aes(
    x = presyn, 
    y = postsyn, 
    label = synapses, 
    size = log2(synapses),
    color = synapses)
    ) + 
  scale_x_discrete(limits = c("head", "sg0", "sg1", "sg2", "sg3", "pyg"),
                   position = "top") +
  scale_y_discrete(limits = rev(c("head", "sg0", "sg1", "sg2", "sg3", "pyg")) ) +
  theme(
    axis.text.x = element_text (
      angle = 90, hjust = 0, vjust = 0.5, size=8.5,
      color = segmental_colors[1:6]
    ), 
    axis.text.y = element_text (
      angle = 0, hjust = 1, vjust = 0.5, size=8.5,
      color = segmental_colors[6:1]
    ),
    axis.title.x = element_text (size=11)
  ) +
  labs(
    x = "postsynaptic segment", 
    y = "presynaptic segment") +
  scale_size_area(max_size=3, guide = NULL) +
  scale_colour_gradient2(
    low = "black",
    mid = segmental_colors[5],
    high = segmental_colors[6],
    midpoint = 2500,
    space = "Lab",
    na.value = "grey50",
    guide = NULL,
    aesthetics = "colour"
  ) +
  theme(panel.background = element_rect(fill = "grey98", color = "white"))

combined_seg_matrix_plot

write.table(combined_seg_matrix, "source_data/Figure11_source_data2.txt", sep = "\t")


}

# graph visualisation -----------------------------------------------------
{
  #with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
  celltype_conn_graph <- graph_from_adjacency_matrix(synapse_matrix_seg,
                                                     mode = c("directed"),
                                                     weighted = TRUE,  diag = TRUE, add.colnames = NULL, add.rownames = NA)
  
  # Convert to object suitable for networkD3
  celltype_conn_graph_d3 <- igraph_to_networkD3(celltype_conn_graph)
  
  #filter values by synapse cutoff
  synapse_cutoff <- 0
  celltype_conn_graph_d3$links$value[celltype_conn_graph_d3$links$value < synapse_cutoff] <- 0
  celltype_conn_graph_d3$links$value
  
  # Define 'group' based on SN/IN type:
  celltype_conn_graph_d3$nodes$group <- as.factor(c("SN","IN","MN","eff",
                                                    "SN","IN","MN","eff",
                                                    "SN","IN","MN","eff",
                                                    "SN","IN","MN","eff",
                                                    "SN","IN","MN","eff",
                                                    "SN","IN","MN","eff"))
  
  celltype_conn_graph_d3$nodes$group
# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["SN", "IN","MN", "eff"]) .range(["#E69F00", "#CC79A7", "#0072B2", "grey40"])'

#add segment position as factor to nodes
celltype_conn_graph_d3$nodes$seg  <- as.factor(c("head","head","head","head",
                                                 "sg0","sg0","sg0","sg0",
                                                 "sg1","sg1","sg1","sg1",
                                                 "sg2","sg2","sg2","sg2",
                                                 "sg3","sg3","sg3","sg3",
                                                 "pygidium","pygidium","pygidium","pygidium"))

#define segment colors to nodes, can be use as an alternative color scheme
my_color_seg <- 'd3.scaleOrdinal() .domain(["head", "sg0","sg1", "sg2","sg3","pygidium"]) .range(["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C"])'

# Plot
sankey_seg <-   sankeyNetwork(Links = celltype_conn_graph_d3$links, Nodes = celltype_conn_graph_d3$nodes, Source = "source",
                Target = "target",  NodeID = "name",Value = "value",
                LinkGroup = NULL, units = "", NodeGroup = "group",
                colourScale = my_color, fontSize = 24,
                fontFamily = "Arial", nodeWidth = 50, nodePadding = 30, margin = NULL,
                height = 400, width = 1400, iterations = 1000, sinksRight = TRUE
  )
sankey_seg
saveNetwork(sankey_seg, "pictures/Sankey_seg_conn.html")
webshot2::webshot(url="pictures/Sankey_seg_conn.html",
                    file="pictures/Sankey_seg_conn.png",
                    vwidth=2000, vheight=750, #define the size of the browser window
                    cliprect = c(150,0,1100, 425), zoom=5, delay = 2)
  
}

# calculate connectivity statistics and plot ------------------------------
{
celltype_conn_graph_d3$links <- as_tibble(celltype_conn_graph_d3$links) %>%
    mutate(source_class = 'empty') %>%
    mutate(source_class = ifelse(source == '0' | source == '4' | source == '8' | source == '12' | source == '16' | source == '20', "SN", source_class)) %>%
    mutate(source_class = ifelse(source == '1' | source == '5' | source == '9' | source == '13'| source == '17' | source == '21', "IN", source_class)) %>%
    mutate(source_class = ifelse(source == '2' | source == '6' | source == '10' | source == '14'| source == '18' | source == '22', "MN", source_class)) %>%
    mutate(source_class = ifelse(source == '3' | source == '7' | source == '11' | source == '15'| source == '19' | source == '23', "effector", source_class)) %>%
    mutate(target_class = 'empty') %>%
    mutate(target_class = ifelse(target == '0' | target == '4' | target == '8' | target == '12' | target == '16' | target == '20', "SN", target_class)) %>%
    mutate(target_class = ifelse(target == '1' | target == '5' | target == '9' | target == '13'| target == '17' | target == '21', "IN", target_class)) %>%
    mutate(target_class = ifelse(target == '2' | target == '6' | target == '10' | target == '14'| target == '18' | target == '22', "MN", target_class)) %>%
    mutate(target_class = ifelse(target == '3' | target == '7' | target == '11' | target == '15'| target == '19' | target == '23', "effector", target_class)) %>%
    mutate(sourcesegment = 'empty') %>%
    mutate(sourcesegment = ifelse(source == '0' | source == '1' | source == '2' | source == '3', "head", sourcesegment)) %>%
    mutate(sourcesegment = ifelse(source == '4' | source == '5' | source == '6' | source == '7', "sg0", sourcesegment)) %>%
    mutate(sourcesegment = ifelse(source == '8' | source == '9' | source == '10' | source == '11', "sg1", sourcesegment)) %>%
    mutate(sourcesegment = ifelse(source == '12' | source == '13' | source == '14' | source == '15', "sg2", sourcesegment)) %>%
    mutate(sourcesegment = ifelse(source == '16' | source == '17' | source == '18' | source == '19', "sg3", sourcesegment)) %>%
    mutate(sourcesegment = ifelse(source == '20' | source == '21' | source == '22' | source == '23', "pygidium", sourcesegment)) %>%
    mutate(targetsegment = 'empty') %>%
    mutate(targetsegment = ifelse(target == '0' | target == '1' | target == '2' | target == '3', "head", targetsegment)) %>%
    mutate(targetsegment = ifelse(target == '4' | target == '5' | target == '6' | target == '7', "sg0", targetsegment)) %>%
    mutate(targetsegment = ifelse(target == '8' | target == '9' | target == '10' | target == '11', "sg1", targetsegment)) %>%
    mutate(targetsegment = ifelse(target == '12' | target == '13' | target == '14' | target == '15', "sg2", targetsegment)) %>%
    mutate(targetsegment = ifelse(target == '16' | target == '17' | target == '18' | target == '19', "sg3", targetsegment)) %>%
    mutate(targetsegment = ifelse(target == '20' | target == '21' | target == '22' | target == '23', "pygidium", targetsegment))
  
#add column with source name
celltype_conn_graph_d3$links$source_name <- lapply(celltype_conn_graph_d3$links$source, 
                                                     function(x) print(celltype_conn_graph_d3$nodes$name[x+1]))
  
#plot connections by target class
celltype_conn_graph_target_plot <- as_tibble(celltype_conn_graph_d3$links) %>%
    ggplot(aes(x=unlist(source_name), y=value, fill = target_class)) +
    geom_col() +
    scale_fill_manual(values = list(SN = "#E69F00", 
                                    IN = "#CC79A7", 
                                    MN = "#0072B2", 
                                    effector = 'grey40') ) +
    scale_x_discrete(limits = c(celltype_conn_graph_d3$nodes$name)) +
    labs(x="source cell class", y="# of synapses",
         title="", fill = 'target \ncell class') +
    theme_minimal() +
    theme(legend.position='top', 
          axis.title.x = element_text(size=11),
          axis.text.y = element_text(size=8.5),
          axis.text.x = element_text(size=8.5, angle = 90),
          legend.text = element_text(size=9),
          legend.title = element_text(size=10),
          legend.key.size = unit(4, "mm"),
          panel.grid.minor = element_blank() )

celltype_conn_graph_target_plot

#plot connections by target segment
celltype_conn_graph_target_seg_plot <- as_tibble(celltype_conn_graph_d3$links) %>%
    ggplot(aes(x=unlist(source_name), y = value,  fill = targetsegment)) +
    geom_col() +
    scale_fill_manual(values = list(head = brewer.pal(6, 'Paired')[1], 
                                    sg0 = brewer.pal(6, 'Paired')[2], 
                                    sg1 = brewer.pal(6, 'Paired')[3], 
                                    sg2 = brewer.pal(6, 'Paired')[4],
                                    sg3 = brewer.pal(6, 'Paired')[5],
                                    pygidium = brewer.pal(6, 'Paired')[6])) +
    scale_x_discrete(limits = c(celltype_conn_graph_d3$nodes$name)) +
    labs(x="source cell class", y="# of synapses",
         title="", fill = 'target \nsegment') +
    theme_minimal() +
    theme(legend.position='top', 
          axis.title.x = element_text(size=11),
          axis.text.y = element_text(size=8.5),
          axis.text.x = element_text(size=8.5, angle = 90, vjust = 0.5, hjust = 1),
          legend.text = element_text(size=9),
          legend.title = element_text(size=10),
          legend.key.size = unit(4, "mm"),
          panel.grid.minor = element_blank() )
  
celltype_conn_graph_target_seg_plot

write_csv(
  as_tibble(celltype_conn_graph_d3$links), 
  "source_data/Figure11_source_data4.txt"
  )

}

# segment 2 cells and circuits ------------------------------

#plot sg2 MNs
{
plot_background_ventral_no_ac()
read_plot_neuron_add('celltype47', blues[c(5,8)], 1)
read_plot_neuron_add('celltype67', bluepurple[c(7,8)], 1)
texts3d(40000,141000, 108000, text = "MNspinning", col='black', cex = 3)
texts3d(112000,141000, 97000, text = "MNbow", col='black', cex = 3)
par3d(zoom=0.38)
plot3d(sg2, WithConnectors = F, soma=TRUE, lwd=1, add=T, alpha=0.1,
       col='grey90')
rgl.snapshot("pictures/Figure12_sg2_MNs.png")
close3d()
}

#plot sg2 INs
{
plot_background_ventral_no_ac()
read_plot_neuron_add('celltype108', blues[c(5,6,8)], 1)
read_plot_neuron_add('celltype157', bluepurple[c(7,8)], 1)
read_plot_neuron_add('celltype169', oranges[c(5,6,7)], 1)
read_plot_neuron_add('celltype180', oranges[c(3,5)], 1)  
plot3d(sg2, WithConnectors = F, soma=TRUE, lwd=1, add=T, alpha=0.1,
       col='grey90')
par3d(zoom=0.38)
texts3d(40000,161000, 102000, text = "INbiaxLeg", col='black', cex = 3)
texts3d(104000,161000, 127000, text = "INsplitVent", col='black', cex = 3)
texts3d(95000,151000, 100000, text = "INcommascFV", col=oranges[8], cex = 3)
texts3d(80000,145000, 88000, text = "INbiaxHmid", col=oranges[6], cex = 3) 
rgl.snapshot("pictures/Figure11_sg2_INs.png")
close3d()
}


# spinGland  cells and circuits ------------------------------

#plot spinGland and circuits
{
plot_background_ventral_no_ac()
read_plot_neuron_add('celltype47', blues[c(5,8)], 1)
read_plot_neuron_add('celltype202', oranges[c(5,8)], 0.5)
read_plot_neuron_add('celltype71', blues[c(3,7)], 0.5)
read_plot_neuron_add('celltype149', 'black', 0.8)
read_plot_neuron_add('celltype79', oranges[6], 0.8)
read_plot_neuron_add('celltype_non_neuronal7', oranges[7], 0.5)
plot3d(sg2, WithConnectors = F, soma=TRUE, lwd=1, add=T, alpha=0.1,
         col='grey90')
par3d(zoom=0.38)
texts3d(72000,141000, 127000, text = "MNspinning", col=blues[8], cex = 3)
texts3d(65000,141000, 70000, text = "INsplitSpin", col= oranges[6], cex = 3)
texts3d(33000,141000, 62000, text = "chaeMech", col=blues[6], cex = 3)
texts3d(100000,141000, 54000, text = "INsplitBronto", col='black', cex = 3)
texts3d(68000, 151000, 104000, text = "INsplitPB", col=oranges[6], cex = 3)
texts3d(32000, 141000, 145000, text = "spinGland", col=oranges[7], cex = 3)
rgl.snapshot("pictures/Figure12_MNspin_in_out.png")
close3d()
}

#spinGland circuit graph

# get connectivity between pyg-specific cell types  --------

{
#these are the sg1 cell groups and their targets
cell_groups <-  list(
  MNspinning, INsplitSpin, chaeMech, 
  INsplitBronto, INsplitPB, spinGland,
  PB, INdescLuqinPDF
)
N_cell_groups <- length(cell_groups)
  
#iterate through cell group neuron lists and get connectivity for all agains all
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
  
#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=N_cell_groups)
  
cell_group_names <-  list(
  "MNspinning", "INsplitSpin", "chaeMech", 
  "INsplitBronto", "INsplitPB", "spinGland",
  "PB", "INdescLuqinPDF"
)
rownames(synapse_matrix) <- as.character(cell_group_names)
colnames(synapse_matrix) <- as.character(cell_group_names)
synapse_matrix
write.csv(as.data.frame(synapse_matrix), "source_data/Figure12_source_data2.txt")

#plot with ggplot
as.data.frame((synapse_matrix)) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
    ggplot() +
    geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction), 
               stroke = 0) + 
    scale_x_discrete(limits = (c(
      "MNspinning", "INsplitSpin", "chaeMech", 
      "INsplitBronto", "INsplitPB", "spinGland",
      "PB", "INdescLuqinPDF"
    ))) +
    scale_y_discrete(limits = rev(c(
      "MNspinning", "INsplitSpin", "chaeMech", 
      "INsplitBronto", "INsplitPB", "spinGland",
      "PB", "INdescLuqinPDF"
    ))) +
    #  coord_flip() +
    theme(
      axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
      axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
      axis.title.x = element_text (size=15),
      axis.title.y = element_text (size=15)
    ) +
    labs(x="postsynaptic cell groups",y="presynaptic cell groups",title=" ") +
    scale_size_area(max_size=5) +
    guides(color = 'legend') +
    scale_colour_gradient2(
      low = "#0072B2",
      mid = "#D55E00",
      high ="#D55E00",
      midpoint = 0.5,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    )+
    #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
    theme(panel.background = element_rect(fill = "grey98", color = "white"))
  
# Saving R ggplot with R ggsave Function
ggsave("pictures/MNspin_cellgroups_syn_matrix.pdf", 
         width = nrow(synapse_matrix)/1.1, 
         height = ncol(synapse_matrix)/1.45, limitsize = TRUE, 
         units = c("cm"))
  
# Saving R ggplot with R ggsave Function
ggsave("pictures/MNspin_cellgroups_syn_matrix.png", 
         width = 1700, 
         height = 1300, limitsize = TRUE, 
         units = c("px"))
}

# graph conversion --------------------------------------------------------
{
#edge weight filtering on the matrix to remove weak edges
synapse_matrix[synapse_matrix < 1] <- 0

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

#calculate node weighted degree
degree=degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)
}

# use visNetwork to plot the network --------------------------------------
{
## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(Conn_graph)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight
  
Conn_graph.visn$nodes$value = degree
Conn_graph.visn$nodes
  
node_colors <- c(
    "#56B4E9", "#CC79A7", "#E69F00",
    "#CC79A7", "#CC79A7", "#cccccc",
    "#E69F00", "#CC79A7"
  )
length(node_colors)
length(Conn_graph.visn$nodes$id)
Conn_graph.visn$nodes$label
Conn_graph.visn$nodes$group <- list(
    "sg2", "sg1", "sg1_3",
    "sg1", "sg1_3", "sg2_3",
    "sg1_3", "head"
)
Conn_graph.visn$nodes$color <- node_colors
  
#hierarchical layout
  
#level	: Number. Default to undefined. When using the hierarchical layout, the level determines where the node is going to be positioned.
Conn_graph.visn$nodes$level <- c(
    "3", "4", "5",
    "4", "4", "3",
    "5", "2"
)

visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 42, type="square") %>%
    visHierarchicalLayout(levelSeparation=250, 
                          nodeSpacing=150,
                          direction='RL',
                          sortMethod='directed',
                          shakeTowards='roots') %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0.1),
             scaling=list(min=2, max=12),
             color = list(inherit=TRUE, opacity=0.7),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 1.2, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(background=Conn_graph.visn$nodes$color, border='black'),
             opacity=0.9,
             shape='box', 
             font=list(color='black', size=24),
             scaling = list(label=list(enabled=TRUE, min=23, max=28)),
             level= Conn_graph.visn$nodes$level) %>%
    visOptions(highlightNearest = TRUE) %>%
    visInteraction(navigationButtons = FALSE,
                   dragNodes = TRUE, dragView = FALSE,
                   zoomView = TRUE) %>%
    visGroups(groupname = "sg2", shape = "dot", opacity=1)
visNet 
  
saveNetwork(visNet, "pictures/visNetwork_MNspinning_input_output.html")
webshot2::webshot(url="pictures/visNetwork_MNspinning_input_output.html",
                    file="pictures/visNetwork_MNspinning_input_output.png",
                    vwidth=1200, vheight=500, #define the size of the browser window
                    cliprect = c(40,80,930,400), zoom=5, delay = 3)
  
}

# pygidial cells and circuits ------------------------------

#plot cioMNcover and covercell
{
plot_background_ventral_no_ac()
read_plot_neuron_add('celltype82', blues[c(5, 6, 8)], 1)
read_plot_neuron_add('celltype_non_neuronal8', oranges[8], 0.2)
texts3d(62000, 141000, 182000, text = "cioMNcover", col=blues[8], cex = 3)
texts3d(28000,141000, 17000, text = "cover cells", col=oranges[8], cex = 3)
plot3d(pyg, WithConnectors = F, soma=TRUE, lwd=1, add=T, alpha=0.1,
         col='grey90')
par3d(zoom=0.45)
rgl.snapshot("pictures/Figure12_pyg_cioMNcover.png")
close3d()
}

#plot SNpygM and other pyg cells
{
plot_background_ventral_no_ac()
read_plot_neuron_add('celltype170', blues[c(5,8)], 1)
read_plot_neuron_add('celltype115', oranges[7], 0.5)
read_plot_neuron_add('celltype147', oranges[8], 0.5)
read_plot_neuron_add('celltype145', oranges[8], 0.5)
read_plot_neuron_add('celltype72', oranges[5], 0.5)
texts3d(55000,148000, 183000, text = "SNpygM", col=blues[7], cex = 3)
texts3d(30000,141000, 169000, text = "SNPDF-pyg", col=oranges[7], cex = 3)
texts3d(92000,141000, 170000, text = "INasc-pyg", col=oranges[8], cex = 3)
texts3d(75000,141000, 160000, text = "INATOpyg", col=oranges[9], cex = 3)
texts3d(72000,141000, 147000, text = "INFVa-pyg", col=oranges[6], cex = 3)
plot3d(pyg, WithConnectors = F, soma=TRUE, lwd=1, add=T, alpha=0.1,
         col='grey90')
par3d(zoom=0.45)
rgl.snapshot("pictures/Figure12_pyg_SN.png")
close3d()
}

#plot pygPBunp

{
plot_background_ventral_no_ac()
read_plot_neuron_add('celltype54', blues[8], 1)
read_plot_neuron_add('celltype8', oranges[8], 1)
read_plot_neuron_add('celltype9', Okabe_Ito[8], 0.5)
read_plot_neuron_add('celltype_non_neuronal33', blues[c(4, 6, 9)], 0.5)
read_plot_neuron_add('celltype_non_neuronal8', oranges[8], 0.2)
texts3d(79000,148000, 180000, text = "pygPBunp", col=blues[7], cex = 3)
texts3d(53000,146000, 40600, text = "Ser-h1", col=oranges[7], cex = 3)
texts3d(81000,141000, 33200, text = "MC", col=Okabe_Ito[8], cex = 3)
texts3d(100000,141000, 8200, text = "vacuolar cells", col=blues[9], cex = 3)
texts3d(28000,141000, 17400, text = "cover cells", col=oranges[8], cex = 3)
plot3d(pyg, WithConnectors = F, soma=TRUE, lwd=1, add=T, alpha=0.1,
         col='grey90')
par3d(zoom=0.45)
rgl.snapshot("pictures/Figure12_pygPBunp.png")
close3d()
}

# get connectivity between pyg-specific cell types  --------
{
#these are the sg1 cell groups and their targets
cell_groups <-  list(
  SNpygM, SNPDFpyg, INascpyg, 
  INdescLuqinPDF, INATOpyg, INFVapyg,
  MNladder, chaeMech, covercell,
  cioMNcover, pyg_cirrus_CR, pygPBunp,
  SerTr1, MNring, INsplitPB
)
N_cell_groups <- length(cell_groups)
  
#iterate through cell group neuron lists and get connectivity for all agains all
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
  
#convert synapse list into a matrix of appropriate dimensions
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=N_cell_groups)
  
cell_group_names <-  list(
  "SNpygM", "SNPDFpyg", "INasc-pyg", 
  "INdescLuqinPDF", "INATOpyg", "INFVa-pyg",
  "MNladder", "chaeMech", "covercell",
  "cioMNcover", "pyg_cirrus_CR", "pygPBunp",
  "Ser-tr1", "MNring", "INsplitPB"
)
rownames(synapse_matrix) <- as.character(cell_group_names)
colnames(synapse_matrix) <- as.character(cell_group_names)
synapse_matrix
write.csv(as.data.frame(synapse_matrix), "source_data/Figure12_source_data3.txt")

#plot with ggplot
as.data.frame((synapse_matrix)) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
    ggplot() +
    geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction), 
               stroke = 0) + 
    scale_x_discrete(limits = (c(
      "SNpygM", "SNPDFpyg", "INasc-pyg", 
      "INdescLuqinPDF", "INATOpyg", "INFVa-pyg",
      "MNladder", "chaeMech", "covercell",
      "cioMNcover", "pyg_cirrus_CR", "pygPBunp",
      "Ser-tr1", "MNring", "INsplitPB"
    ))) +
    scale_y_discrete(limits = rev(c(
      "SNpygM", "SNPDFpyg", "INasc-pyg", 
      "INdescLuqinPDF", "INATOpyg", "INFVa-pyg",
      "MNladder", "chaeMech", "covercell",
      "cioMNcover", "pyg_cirrus_CR", "pygPBunp",
      "Ser-tr1", "MNring", "INsplitPB"
    ))) +
    #  coord_flip() +
    theme(
      axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
      axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
      axis.title.x = element_text (size=15),
      axis.title.y = element_text (size=15)
    ) +
    labs(x="postsynaptic cell groups",y="presynaptic cell groups",title=" ") +
    scale_size_area(max_size=5) +
    guides(color = 'legend') +
    scale_colour_gradient2(
      low = "#0072B2",
      mid = "#D55E00",
      high ="#D55E00",
      midpoint = 0.5,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    )+
    #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
    theme(panel.background = element_rect(fill = "grey98", color = "white"))
  
# Saving R ggplot with R ggsave Function
ggsave("pictures/pyg_cellgroups_syn_matrix.pdf", 
         width = nrow(synapse_matrix)/1.1, 
         height = ncol(synapse_matrix)/1.45, limitsize = TRUE, 
         units = c("cm"))
  
# Saving R ggplot with R ggsave Function
ggsave("pictures/pyg_cellgroups_syn_matrix.png", 
         width = 1700, 
         height = 1300, limitsize = TRUE, 
         units = c("px"))
}

# graph conversion --------------------------------------------------------
{
#edge weight filtering on the matrix to remove weak edges
synapse_matrix[synapse_matrix < 1] <- 0

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

#calculate node weighted degree
degree=degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)
}

# use visNetwork to plot the network --------------------------------------
{
## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(Conn_graph)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight
  
Conn_graph.visn$nodes$value = degree
Conn_graph.visn$nodes

node_colors <- c(
  "#E69F00", "#E69F00", "#CC79A7",
  "#CC79A7", "#CC79A7", "#CC79A7",
  "#56B4E9", "#E69F00", "#cccccc",
  "#56B4E9", "#E69F00", "#E69F00", 
  "#56B4E9","#56B4E9", "#CC79A7"
)
length(node_colors)
length(Conn_graph.visn$nodes$id)
Conn_graph.visn$nodes$label
Conn_graph.visn$nodes$group <- list(
  "pyg", "pyg", "pyg",
  "head", "pyg", "pyg",
  "sg1", "sg1_3", "head",
  "pyg", "trunk", "pyg",
  "sg1", "trunk", "sg1_3"
)
Conn_graph.visn$nodes$color <- node_colors

#hierarchical layout
  
#level	: Number. Default to undefined. When using the hierarchical layout, the level determines where the node is going to be positioned.
Conn_graph.visn$nodes$level <- c(
  "5", "5", "4",
  "1", "4", "4",
  "2", "3", "1",
  "4", "3", "5",
  "2", "3", "3"
)

visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout_nicely", physics = TRUE, randomSeed = 42, type="square") %>%
    visHierarchicalLayout(levelSeparation=250, 
                          nodeSpacing=150,
                          direction='RL',
                          sortMethod='hubsize',
                          shakeTowards='roots') %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0.1),
             scaling=list(min=2, max=12),
             color = list(inherit=TRUE, opacity=0.7),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 1.2, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(background=Conn_graph.visn$nodes$color, border='black'),
             opacity=0.9,
             shape='box', 
             font=list(color='black', size=24),
             scaling = list(label=list(enabled=TRUE, min=23, max=28)),
             level= Conn_graph.visn$nodes$level) %>%
    visOptions(highlightNearest = TRUE) %>%
    visInteraction(navigationButtons = FALSE,
                   dragNodes = TRUE, dragView = FALSE,
                   zoomView = TRUE) %>%
    visGroups(groupname = "pyg", shape = "dot", opacity=1)
visNet 

saveNetwork(visNet, "pictures/visNetwork_pyg_specific_input_output.html")
webshot2::webshot(url="pictures/visNetwork_pyg_specific_input_output.html",
                    file="pictures/visNetwork_pyg_specific_input_output.png",
                    vwidth=1200, vheight=500, #define the size of the browser window
                    cliprect = c(30,80,930,400), zoom=5, delay = 3)
  
}



# quantify effector inputs from different segments ---------------

head <- skids_by_2annotations("episphere", "connectome")
sg0 <- skids_by_2annotations("segment_0", "connectome")
sg1 <- skids_by_2annotations("segment_1", "connectome")
sg2 <- skids_by_2annotations("segment_2", "connectome")
sg3 <- skids_by_2annotations("segment_3", "connectome")
pyg <- skids_by_2annotations("pygidium", "connectome")

pigment <- skids_by_2annotations("^pigment cell$", "connectome")
gland <- skids_by_2annotations("^gland cell$", "connectome")
ciliated <- skids_by_2annotations("^ciliated cell$", "connectome")
muscle <- skids_by_2annotations("^muscle$", "connectome")


seg_conn <- tibble(presyn = c("head", "sg0", "sg1", "sg2", "sg3", "pyg"),
                   'pigment cell' = c(
                     length(catmaid_get_connectors_between(head, pigment, pid = 11)$pre_skid),
                     length(catmaid_get_connectors_between(sg0, pigment, pid = 11)$pre_skid),
                     length(catmaid_get_connectors_between(sg1, pigment, pid = 11)$pre_skid),
                     length(catmaid_get_connectors_between(sg2, pigment, pid = 11)$pre_skid),
                     length(catmaid_get_connectors_between(sg3, pigment, pid = 11)$pre_skid),
                     length(catmaid_get_connectors_between(pyg, pigment, pid = 11)$pre_skid)
                   )
)

seg_conn <- seg_conn %>%
  mutate('gland cell' =  c(
    length(catmaid_get_connectors_between(head, gland, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg0, gland, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg1, gland, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg2, gland, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg3, gland, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(pyg, gland, pid = 11)$pre_skid)
  )
  )

seg_conn <- seg_conn %>%
  mutate(muscle =  c(
    length(catmaid_get_connectors_between(head, muscle, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg0, muscle, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg1, muscle, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg2, muscle, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg3, muscle, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(pyg, muscle, pid = 11)$pre_skid)
  )
  )

seg_conn <- seg_conn %>%
  mutate('ciliated cell' =  c(
    length(catmaid_get_connectors_between(head, ciliated, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg0, ciliated, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg1, ciliated, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg2, ciliated, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(sg3, ciliated, pid = 11)$pre_skid),
    length(catmaid_get_connectors_between(pyg, ciliated, pid = 11)$pre_skid)
  )
  )

seg_con_plot <- seg_conn %>%
  pivot_longer(-presyn, names_to = "postsyn_type", values_to = "synapses") %>%
  ggplot(aes(presyn, synapses, fill = presyn)) +
  geom_col() +
  facet_wrap(~postsyn_type) +
  theme_minimal() +
  coord_flip() +
  scale_x_discrete(limits = rev(c("head", "sg0", "sg1", "sg2", "sg3", "pyg"))) +
  scale_fill_manual(values = segmental_colors) +
  theme(
    axis.text.x = element_text (
      angle = 0, hjust = 0, vjust = 0.5, size = 8.5
    ), 
    axis.text.y = element_text (
      angle = 0, hjust = 1, vjust = 0.5, size = 8.5
    ),
    axis.title = element_text (size=11),
    strip.text = element_text(size = 11)
  ) +
  labs(
    x = "presynaptic segment", 
    y = "synapses",
    facet = "a"
  ) +
  guides(fill = "none") 

seg_con_plot

write_csv(seg_conn, "source_data/Figure11_source_data3.txt")

# plot global reach neurons -------------

load_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(
    annotation, pid=11), 
    function(x) smooth_neuron(x, sigma=6000))
}

head_global_skids <- skids_by_2annotations("episphere", "global_reach")
sg1_global_skids <- skids_by_2annotations("segment_1", "global_reach")
pyg_global_skids <- skids_by_2annotations("pygidium", "global_reach")

# load PDF neurons ---------------
head_global <- load_neuron(head_global_skids)
sg1_global <- load_neuron(sg1_global_skids)
pyg_global <- load_neuron(pyg_global_skids)

plot_background_ventral_no_ac()
plot3d(head_global, soma = TRUE, lwd = c(2:4), color = blues[c(3:9)])
rgl.snapshot("pictures/head_global.png")
close3d()

plot_background_ventral_no_ac()
plot3d(sg1_global, soma = TRUE, lwd = c(2:4), color = bluepurple[c(3:9)])
rgl.snapshot("pictures/sg1_global.png")
close3d()

plot_background_ventral_no_ac()
plot3d(pyg_global, soma = TRUE, lwd = c(2:4), color = oranges[c(3:9)])
rgl.snapshot("pictures/pyg_global.png")
close3d()


# multi-panel figure -------------------------------------------------

panel_all_seg <- ggdraw() + 
  draw_image(readPNG  ("pictures/Figure11_seg_seg_color.png")) +
  draw_label("head", x = 0.55, y = 0.94, size = 11) +
  draw_label("sg0", x = 0.54, y = 0.79, size = 11) +
  draw_label("sg1", x = 0.5, y = 0.65, size = 11) +
  draw_label("sg2", x = 0.5, y = 0.47, size = 11) +
  draw_label("sg3", x = 0.48, y = 0.3, size = 11) +
  draw_label("pygidium", x = 0.44, y = 0.13, size = 11) +
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


panel_segments <- ggdraw() +
  draw_image(readPNG("pictures/Full_connectome_segments.png")) 
 
panel_sankey_seg <- ggdraw() + 
    draw_image(readPNG("pictures/Sankey_seg_conn.png"))

#segmental cell types panels
panel_sg1_celltypes <- ggdraw() + 
  draw_image(readPNG("pictures/Figure11_sg1_celltypes.png")) +
  draw_label("sg1 all cell types", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_sg2_celltypes <- ggdraw() + 
  draw_image(readPNG("pictures/Figure11_sg2_celltypes.png")) +
  draw_label("sg2 all cell types", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_sg3_celltypes <- ggdraw() + 
  draw_image(readPNG("pictures/Figure11_sg3_celltypes.png")) +
  draw_label("sg3 all cell types", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_pyg_celltypes <- ggdraw() + 
  draw_image(readPNG("pictures/Figure11_pyg_celltypes.png")) +
  draw_label("pyg all cell types", x = 0.5, y = 0.99, 
             fontfamily = "sans", size = 11)

panel_head_global <-ggdraw() + draw_image(readPNG("pictures/head_global.png")) + 
  draw_label(
    "head global reach", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )
panel_sg1_global <-ggdraw() + draw_image(readPNG("pictures/sg1_global.png")) + 
  draw_label(
    "sg1 global reach", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_pyg_global <-ggdraw() + draw_image(readPNG("pictures/pyg_global.png")) + 
  draw_label(
    "pyg global reach", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_overview <- plot_grid(
  panel_all_seg, panel_segments, panel_sankey_seg,
  ncol=3, rel_widths = c(1, 2.2, 3.8), rel_heights = c(1),
  labels = c("A", "B", "C"),
  label_size = 12, label_y = 1.03, label_x = 0,
  label_fontfamily = "sans", label_fontface = "plain"
)

panel_mid <- plot_grid(
  combined_seg_matrix_plot, seg_con_plot,
  celltype_conn_graph_target_plot, celltype_conn_graph_target_seg_plot,
  ncol=4, rel_widths = c(2, 3, 3, 3), rel_heights = c(1),
  labels = c("D", "E", "F", "G"),
  label_size = 12, label_y = 1.03, label_x = 0,
  label_fontfamily = "sans", label_fontface = "plain"
)

#assemble anatomy panels
panels_anatomy <- plot_grid(
  panel_sg1_celltypes, panel_sg2_celltypes, panel_sg3_celltypes, panel_pyg_celltypes,
  panel_head_global, panel_sg1_global, panel_pyg_global,
  ncol=7, rel_widths = c(1), rel_heights = c(1),
  labels = c("H", "I", "J", "K", "L", "M", "N"),
  label_size = 12, label_y = 1.03, label_x = 0,
  label_fontfamily = "sans", label_fontface = "plain"
)

layout <- "
A
#
B
#
C
"

Fig11  <- panel_overview + 
  panel_mid +
  panels_anatomy + 
  plot_layout(design = layout, heights = c(1, 0.05, 0.8, 0.05, 1)) & 
  theme(plot.tag = element_text(size = 12, face='plain'))
  
ggsave("Figures/Figure11.png", limitsize = FALSE, 
       units = c("px"), Fig11, width = 4600, height = 3200, bg='white')

ggsave("Figures/Figure11.pdf", limitsize = FALSE, 
       units = c("px"), Fig11, width = 4600, height = 3200)

 


panel_sg0 <- ggdraw() + draw_image(readPNG("pictures/Figure_sg0.png")) +
  draw_label("sg0 and mouth neurons", x = 0.45, y = 0.96, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1) +
  #use expression to define mu with its latex code
  draw_label(expression(paste("50 ", mu, "m")), x = 0.73, y = 0.14, fontfamily = "sans", fontface = "plain",
             color = "black", size = 10)

panel_cil <- ggdraw() + draw_image(readPNG("pictures/Figure_sg1_Loop_Ser_trochs.png")) +
  draw_label("sg1 ciliomotors", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11) 
panel_MC3 <- ggdraw() + draw_image(readPNG("pictures/Figure_sg1_MN_targets.png")) +
  draw_label("sg1 other motor", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11)
panel_sg1IN <- ggdraw() + draw_image(readPNG("pictures/Figure_sg1_IN.png")) +
  draw_label("sg1 INs", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11)
panel_sg1_network <- ggdraw() + draw_image(readPNG("pictures/visNetwork_sg1_specific_input_output.png")) +
  draw_label("sg1 networks", x = 0.2, y = 0.97, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11) +
  draw_line(c(0.44, 0.44), c(1, 0.03), size=0.3) +
  draw_line(c(0.69, 0.69), c(1, 0.03), size=0.3) +
  draw_label("sg1", x = 0.37, y = 0.99, size = 11) +
  draw_label("trunk", x = 0.53, y = 0.99, size = 11) +
  draw_label("head", x = 0.77, y = 0.99, size = 11)

#sg2 panels
panel_sg2_MNs <- ggdraw() + 
  draw_image(readPNG("pictures/Figure12_sg2_MNs.png")) 
panel_sg2_INs <- ggdraw() + 
  draw_image(readPNG("pictures/Figure11_sg2_INs.png"))


panel_MNspin_partners <- ggdraw() + 
  draw_image(readPNG("pictures/Figure12_MNspin_in_out.png"))
panel_MNspin_network <- ggdraw() + 
  draw_image(readPNG("pictures/visNetwork_MNspinning_input_output.png")) +
  draw_label("MNspinning network", x = 0.3, y = 0.97,
             color = "black", size = 11)

panel_spinGland <- ggdraw() + 
  draw_image(readPNG("images_notR/Spinning_Gland_TEM.png")) +
  draw_label("spinGland", 
             x = 0.63, y = 0.96, size = 12, color = "white") +
  draw_label(expression(paste("2 ", mu, " m")), 
             x = 0.23, y = 0.54, size = 9, color = "white") +
  draw_label(expression(paste("5 ", mu, " m")), 
             x = 0.23, y = 0.04, size = 9, color = "white")


#pyg panels
panel_pyg1 <- ggdraw() + 
  draw_image(readPNG("pictures/Figure12_pyg_cioMNcover.png")) 
panel_pyg2 <- ggdraw() + 
  draw_image(readPNG("pictures/Figure12_pyg_SN.png"))
panel_pyg3 <- ggdraw() + 
  draw_image(readPNG("pictures/Figure12_pygPBunp.png"))


panel_pyg_cover_TEM <- ggdraw() + 
  draw_image(readPNG("images_notR/pygidium_covercell_TEM.png")) +
  draw_label("covercell", 
             x = 0.67, y = 0.96, size = 12, color = "white") +
  draw_label("pygidium", 
             x = 0.65, y = 0.45, size = 12, color = "white") +
  draw_label(expression(paste("2 ", mu, " m")), 
             x = 0.25, y = 0.55, size = 9, color = "white") +
  draw_label(expression(paste("20 ", mu, " m")), 
             x = 0.24, y = 0.05, size = 9, color = "white")

panel_pyg_network <- ggdraw() + 
  draw_image(readPNG("pictures/visNetwork_pyg_specific_input_output.png")) +
  draw_line(c(0.38, 0.38), c(1, 0.07), size=0.3) +
  draw_line(c(0.60, 0.60), c(1, 0.07), size=0.3) +
  draw_line(c(0.75, 0.75), c(1, 0.07), size=0.3) +
  draw_label("pygidium", x = 0.31, y = 0.99, size = 11) +
  draw_label("sg1-3", x = 0.48, y = 0.99, size = 11) +
  draw_label("sg1", x = 0.67, y = 0.99, size = 11) +
  draw_label("head", x = 0.8, y = 0.99, size = 11) +
  draw_label("pygidium\nnetwork", x = 0.13, y = 0.95, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11)



layout <- "
JJKKLLMMMMMMMM
NNOOPPQQ#RRRRR
##############
SSTTUUVVWWWWWW
"

Fig12  <- panel_cil + panel_MC3 + panel_sg1IN + panel_sg1_network +
  panel_spinGland + panel_MNspin_partners + panel_sg2_MNs + panel_sg2_INs + 
  panel_MNspin_network + 
  panel_pyg_cover_TEM + panel_pyg1 + panel_pyg2 + panel_pyg3 + panel_pyg_network +
  plot_layout(design = layout, heights = c(1, 1, 0.02, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("Figures/Figure12.png", limitsize = FALSE, 
       units = c("px"), Fig12, width = 3800, height = 2500, bg='white')

ggsave("Figures/Figure12.pdf", limitsize = FALSE, 
       units = c("px"), Fig12, width = 3800, height = 2500)

# Fig 11 fig suppl 1 ------------

ggsave(
  "Figures/Figure11_fig_suppl1.png", 
  limitsize = FALSE, 
  units = c("px"), 
  synapse_matrix_seg_plot, 
  width = 1500, height = 1500)

ggsave(
  "Figures/Figure11_fig_suppl1.pdf", 
  limitsize = FALSE, 
  units = c("px"), 
  synapse_matrix_seg_plot, 
  width = 1500, height = 1500)




