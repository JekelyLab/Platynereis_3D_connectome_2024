#R/natverse code to generate Figure mech girdle anatomy overview for the Platynereis 3d connectome paper
#Gaspar Jekely March-Dec 2022


#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# functions ---------------------------------------------------------------

#adjusted background plotting function
plot_background_mech <- function(){
  plot_background_ventral()
  clear3d()
  nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=0,
         add=T, alpha=0.05,
         col='grey80')
  plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=F, lwd=3,
         add=T, alpha=0.1, col="grey50")
}

read_neurons <- function(annotation){
  nlapply(read.neurons.catmaid(annotation, pid=11),
  function(x) smooth_neuron(x, sigma=8000))
}


# read neurons ------------------------------------------------------------

#read mech neuron types
{
doCRunp <- read_neurons("^celltype102$")
h_antCR <- read_neurons("^celltype100$")
cirrusCR <- read_neurons("^celltype101$")
pygPBunp <- read_neurons("^celltype54$")
PB <- read_neurons("^celltype80$")
antPU <- read_neurons("^celltype88$")
hCirrus_interparaPU <- read_neurons("^celltype94$")
pygCirrusPU <- read_neurons("^celltype98$")
hPU <- read_neurons("^celltype96$")
neuro_notoPU <- read_neurons("^celltype95$")
dsoPU <- read_neurons("^celltype97$")
paraPU <- read_neurons("^celltype91$")
spinPU <- read_neurons("^celltype89$")
hPU2l_asymPDF <- read_neurons("^celltype99$")
pygCirrusPUSM <- read_neurons("^celltype93$")
SNblunt <- read_neurons("^celltype148$")
SNpygM <- read_neurons("^celltype170$")
SNantler <- read_neurons("^celltype26$")
SNPDF_pyg <- read_neurons("^celltype115$")
SNFVa <- read_neurons("^celltype143$")
SNbronto <- read_neurons("^celltype168$")
interparaPM <- read_neurons("^celltype81$")
chaeMech <- read_neurons("^celltype71$")
palp <- read_neurons("^palp sensory neuron$")
antenna <- read_neurons("^antennae_cell$")
INsplitPUh <- read_neurons("^celltype83$")
INMC3 <- read_neurons("^celltype75$")
INrope <- read_neurons("^celltype58$")
INsplitCR <- read_neurons("^celltype73$")
INsplitCRATO <- read_neurons("^celltype74$")
INsplitBronto <- read_neurons("^celltype149$")
INsplitVent <- read_neurons("^celltype157$")
INasc_pyg <- read_neurons("^celltype147$")
INsplitPBant <- read_neurons("^celltype78$")
INsplitPB <- read_neurons("^celltype79$")
INsplitPB_RFYa <- read_neurons("^celltype77$")
INCM <- read_neurons("^celltype60$")
MNant <- read_neurons("^celltype19$")
Loop <- read_neurons("^celltype59$")
MC3cover <- read_neurons("^celltype87$")

CR <- read_neurons("^CRneurons$")
PU <- read_neurons("^PUneurons$") 
PB <- read_neurons("^Biciliated_penetrating_cell$") 
INsplit <- read_neurons("^INsplit$") 
MUSlongV <- read_neurons("^celltype_non_neuronal76$") 
MUStrans <- read_neurons("^celltype_non_neuronal74$") 
}

#load other cell clusters
{
stomodeum <- nlapply(
  read.neurons.catmaid("^stomodeum$", pid=11),
  function(x) smooth_neuron(x, sigma=8000)
)
eye <- read_neurons("^celltype1$")
MUSring_pyg <- read_neurons("^celltype_non_neuronal78$")
gut<- read_neurons("^gut$")
girdle <- read_neurons("^mechanosensory_girdle$")
acicula <- read_neurons("^acicula$")
}


# plot neurons ------------------------------------------------------------

#plot mechanosensory neurons in girdle
{
plot_background_ventral()
clear3d()
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))

plot3d(yolk, soma=TRUE, lwd = 0,
         add = TRUE, alpha = 0.05,
         col='grey70')
  
plot3d(doCRunp, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(h_antCR, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[4])
plot3d(cirrusCR, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[5])
plot3d(pygPBunp, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.6, col = oranges[7])
plot3d(PB, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.4, col = oranges[8])
plot3d(antPU, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.8, col = oranges[4])
plot3d(hCirrus_interparaPU, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[2])
plot3d(hPU, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(neuro_notoPU, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = oranges[4])
plot3d(dsoPU, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.4, col = Okabe_Ito[8])
plot3d(paraPU, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.5, col = oranges[7])
plot3d(spinPU, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.4, col = oranges[8])
plot3d(hPU2l_asymPDF, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.7, col = oranges[5])
plot3d(pygCirrusPUSM, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[2])
plot3d(SNblunt, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(SNpygM, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = oranges[9])
plot3d(SNantler, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.2, col = Okabe_Ito[8])
plot3d(SNPDF_pyg, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.6, col = oranges[7])
plot3d(SNFVa, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.5, col = oranges[8])
plot3d(SNbronto, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.8, col = oranges[4])
plot3d(interparaPM, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[2])
plot3d(chaeMech, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(antenna, soma=TRUE, lwd = 2,
         add = TRUE, alpha = 0.3,
         col=oranges[5])
plot3d(palp, soma=TRUE, lwd = 2,
         add = TRUE, alpha = 0.3,
         col="#D55E00")
plot3d(stomodeum, soma=F, lwd = 2,
         add = TRUE, alpha = 0.2, col = "grey50")
plot3d(MUSring_pyg, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.3, col = "grey40")
plot3d(gut, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.2, col = "grey50")
  
#  texts3d(78000,132000, 24000, text = "stomodeum", col = 'black', cex = 2.2)
par3d(zoom=0.55)
filename <- paste("pictures/Figure_mech_overview_SN.png")
rgl.snapshot(filename)
close3d()
}

#plot interneurons in girdle
{
plot_background_ventral()
clear3d()
plot3d(yolk, soma=TRUE, lwd = 0,
         add = TRUE, alpha = 0.05,
         col='grey70')
plot3d(INsplitPUh, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[8])
plot3d(INMC3, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[7])
plot3d(INrope, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[5])
plot3d(INsplitCR, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[9])
plot3d(INsplitCRATO, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = blues[4])
plot3d(INCM, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = blues[3])
plot3d(INsplitPB_RFYa, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[7])
plot3d(INsplitBronto, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[7])
plot3d(INsplitVent, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[5])
plot3d(INasc_pyg, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[2])
plot3d(INsplitPBant, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.5, col = bluepurple[9])
plot3d(INsplitPB, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.8, col = bluepurple[9])
plot3d(stomodeum, soma=F, lwd = 2,
         add = TRUE, alpha = 0.2, col = "grey50")
plot3d(MUSring_pyg, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.3, col = "grey40")
plot3d(gut, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.2, col = "grey50")
  
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
par3d(zoom=0.55)
texts3d(78000,119000, 30000, text = "stomodeum", col = 'black', cex = 3)
filename <- paste("pictures/Figure_mec_overview_girdle_vis_ventr.png")
rgl.snapshot(filename)

nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)%*%rotationMatrix(-pi/2, 0, 0, 1)))
par3d(zoom=0.55)
texts3d(78000,95000, 178000, text = "hindgut", col = 'black', cex = 4)
filename <- paste("pictures/Figure_mec_overview_girdle_vis_lat.png")
rgl.snapshot(filename)  
  
#define a z clipping plane for the frontal view
nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5)%*%rotationMatrix(0.14, 0, 1, 0)%*%rotationMatrix(-0.05, 1, 0, 0))
#z-axis clip
#clipplanes3d(0, 0, -1, 140000)
#y-axis clip
clipplanes3d(1, 0, 0.16, 7000)
#x-axis clip
clipplanes3d(0, -1, 0.16, 130000)
par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
par3d(zoom=0.52)

clear3d()
plot3d(yolk, soma=TRUE, lwd = 0,
         add = TRUE, alpha = 0.05,
         col='grey70')
plot3d(INsplitPUh, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[8])
plot3d(INMC3, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[7])
plot3d(INrope, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[5])
plot3d(INsplitCR, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[9])
plot3d(INsplitCRATO, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = blues[4])
plot3d(INCM, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = blues[3])
plot3d(INsplitPB_RFYa, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[7])
plot3d(INsplitBronto, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[7])
plot3d(INsplitVent, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[5])
plot3d(INasc_pyg, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[2])
plot3d(INsplitPBant, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.5, col = bluepurple[9])
plot3d(INsplitPB, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.8, col = bluepurple[9])
plot3d(stomodeum, soma=F, lwd = 2,
         add = TRUE, alpha = 0.2, col = "grey50")
plot3d(MUSring_pyg, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.3, col = "grey40")
plot3d(gut, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.2, col = "grey50")
plot3d(eye, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.2, col = Okabe_Ito[8])
texts3d(36000,42000, 24000, text = "eye", col = 'black', cex = 4)
texts3d(96000,18000, 25000, text = "yolk", col = 'grey30', cex = 4)
  
filename <- paste("pictures/Figure_mec_overview_girdle_vis_ant.png")
rgl.snapshot(filename)
close3d()
}  

#plot MN in girdle
{
plot_background_ventral()
clear3d()
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))

plot3d(yolk, soma=TRUE, lwd = 0,
         add = TRUE, alpha = 0.05,
         col='grey70')
plot3d(Loop, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[8])
plot3d(MC3cover, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.6, col = blues[9])
plot3d(pygCirrusPUSM, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 1, col = Okabe_Ito[7])
plot3d(MNant, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.6, col = 'black')
plot3d(pygPBunp, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.7, col = blues[9])
plot3d(stomodeum, soma=F, lwd = 2,
         add = TRUE, alpha = 0.2, col = "grey50")
plot3d(MUSring_pyg, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.3, col = "grey40")
plot3d(gut, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.2, col = "grey50")
par3d(zoom=0.55)
  
filename <- paste("pictures/Figure_mech_overviewMN.png")
rgl.snapshot(filename)
close3d()
}  

#plot all mechanosensory girdle neurons
{
plot_background_ventral()
clear3d()
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
par3d(zoom=0.6)
plot3d(yolk, soma=TRUE, lwd = 0,
         add = TRUE, alpha = 0.05,
         col='grey70')
plot3d(scalebar_50um_ventral, soma=F, lwd = 5,
       add = TRUE, alpha = 1,
       col="black")
  
plot3d(Loop, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[8])
plot3d(MC3cover, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.6, col = blues[9])
plot3d(pygCirrusPUSM, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 1, col = Okabe_Ito[7])
plot3d(MNant, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.6, col = 'black')
plot3d(pygPBunp, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.7, col = blues[9])
plot3d(doCRunp, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(h_antCR, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[4])
plot3d(cirrusCR, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[5])
plot3d(pygPBunp, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.6, col = oranges[7])
plot3d(PB, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.4, col = oranges[8])
plot3d(antPU, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.8, col = oranges[4])
plot3d(hCirrus_interparaPU, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[2])
plot3d(hPU, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(neuro_notoPU, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = oranges[4])
plot3d(dsoPU, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.4, col = Okabe_Ito[8])
plot3d(paraPU, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.5, col = oranges[7])
plot3d(spinPU, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.4, col = oranges[8])
plot3d(hPU2l_asymPDF, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.7, col = oranges[5])
plot3d(pygCirrusPUSM, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[2])
plot3d(SNblunt, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(SNpygM, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = oranges[9])
plot3d(SNantler, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.2, col = Okabe_Ito[8])
plot3d(SNPDF_pyg, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.6, col = oranges[7])
plot3d(SNFVa, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.5, col = oranges[8])
plot3d(SNbronto, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.8, col = oranges[4])
plot3d(interparaPM, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = oranges[2])
plot3d(chaeMech, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = oranges[3])
plot3d(antenna, soma=TRUE, lwd = 2,
         add = TRUE, alpha = 0.3,
         col=oranges[5])
plot3d(palp, soma=TRUE, lwd = 2,
         add = TRUE, alpha = 0.3,
         col="#D55E00")
plot3d(INsplitPUh, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[8])
plot3d(MC3cover, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.6, col = blues[7])
plot3d(INMC3, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[7])
plot3d(INrope, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[5])
plot3d(INsplitCR, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.4, col = blues[9])
plot3d(INsplitCRATO, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = blues[4])
plot3d(INCM, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = blues[3])
plot3d(INsplitPB_RFYa, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[7])
plot3d(INsplitBronto, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = blues[7])
plot3d(INsplitVent, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[5])
plot3d(INasc_pyg, soma = TRUE, lwd = 4,
         add = TRUE, alpha = 0.8, col = Okabe_Ito[2])
plot3d(INsplitPBant, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.5, col = bluepurple[9])
plot3d(INsplitPB, soma = TRUE, lwd = 1,
         add = TRUE, alpha = 0.8, col = bluepurple[9])
plot3d(stomodeum, soma=F, lwd = 2,
         add = TRUE, alpha = 0.2, col = "grey50")
plot3d(MUSring_pyg, soma = TRUE, lwd = 2,
         add = TRUE, alpha = 0.3, col = "grey40")
plot3d(gut, soma = TRUE, lwd = 3,
         add = TRUE, alpha = 0.2, col = "grey50")
par3d(zoom=0.55)
filename <- paste("pictures/Figure_mec_overview_all_ventr.png")
rgl.snapshot(filename)
  
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)%*%rotationMatrix(-pi/2, 0, 0, 1)))
filename <- paste("pictures/Figure_mec_overview_all_lat.png")
rgl.snapshot(filename)
close3d()
}


# connectome graph coloured by girdle cell type ---------------------------


#read the connectome visNetwork R data file
#this file was generated with /Figure2-connectome/Figure_connectome.R
conn_graph.visn <- readRDS("supplements/connectome_graph.rds", refhook = NULL)

#mechanosensory girdle skids
girdle_skids <- catmaid_query_by_annotation("mechanosensory_girdle", pid = 11)$skid
girdle_SN_skids <- skids_by_2annotations("mechanosensory_girdle", "Sensory neuron")
girdle_IN_skids <- skids_by_2annotations("mechanosensory_girdle", "interneuron")
girdle_MN_skids <- skids_by_2annotations("mechanosensory_girdle", "motorneuron")

# print graph with cells coloured by girdle ----------------------------------

# create a list showing if a skid is in the girdle or not
if_girdle <- ifelse(as.numeric(conn_graph.visn$nodes$skids) %in% girdle_SN_skids, "SN girdle",
                    ifelse(as.numeric(conn_graph.visn$nodes$skids) %in% girdle_IN_skids, "IN girdle",
                           ifelse(as.numeric(conn_graph.visn$nodes$skids) %in% girdle_MN_skids, "MN girdle", "not girdle")))

#overwrite group value (partition) with segment value (for colouring)
conn_graph.visn$nodes$group <-  if_girdle

#rotate coordinates
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol = 2)
library(autoimage)
coords_rotated <- autoimage::rotate(
  coords, 
  -pi / 3, 
  pivot = c(0, 0)
)

#remove colour info (which takes precedence over group colour)
conn_graph.visn$nodes$color <- c()

visNet_girdle <- visNetwork(conn_graph.visn$nodes,conn_graph.visn$edges) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(smooth = list(type = 'curvedCW', roundness=0),
           scaling=list(min=1, max=15),
           color = list(inherit=FALSE, opacity=0.1),
           arrows = list(to = list(enabled = TRUE, 
                                   scaleFactor = 0.5, type = 'arrow'))) %>%
  visNodes(borderWidth=0.3, 
           color = list(background=conn_graph.visn$nodes$group, border='black'),
           opacity = 1, 
           font = list(size = 0)) %>%
  visOptions(highlightNearest = list(enabled=TRUE, degree=1, 
                                     algorithm = 'hierarchical',labelOnly=FALSE), 
             width = 800, height = 800, autoResize = FALSE) %>%
  visGroups(groupname = "SN girdle", color=Okabe_Ito[1], shape = "dot", 
            opacity=1, size=30) %>%
  visGroups(groupname = "IN girdle", color=Okabe_Ito[7], shape = "square", 
            opacity=1, size=30) %>%
  visGroups(groupname = "MN girdle", color=Okabe_Ito[5], shape = "triangle", 
            opacity=1, size=45) %>%
  visGroups(groupname = "not girdle", shape = "dot", 
            opacity=1, size=15, color="#DDDDDD") %>%
  addFontAwesome() %>%
  visLegend(
    useGroups = TRUE, 
    width = 0.1,
    ncol = 1,
    position = "left"
  )
visNet_girdle 

#save as html
saveNetwork(visNet_girdle, "pictures/Full_connectome_girdle.html", selfcontained = TRUE)
webshot2::webshot(url="pictures/Full_connectome_girdle.html",
                  file="pictures/Full_connectome_girdle.png",
                  vwidth = 800, vheight = 800, #define the size of the browser window
                  cliprect = c(80, 60, 850, 800), zoom=1, delay = 1)


# assemble figure ---------------------------------------------------------


imgSN <- readPNG("pictures/Figure_mech_overview_SN.png")
imgMN <- readPNG("pictures/Figure_mech_overviewMN.png")
imgAll <- readPNG("pictures/Figure_mec_overview_all_ventr.png")
imgAll_lat <- readPNG("pictures/Figure_mec_overview_all_lat.png")
  
imgIN_ventr <- readPNG("pictures/Figure_mec_overview_girdle_vis_ventr.png")
imgIN_lat <- readPNG("pictures/Figure_mec_overview_girdle_vis_lat.png")
imgIN_ant <- readPNG("pictures/Figure_mec_overview_girdle_vis_ant.png")


panelSN <- ggdraw() + draw_image(imgSN, scale = 1) + 
    draw_label("mechanosensory", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11)

panelMN <- ggdraw() + draw_image(imgMN, scale = 1) + 
    draw_label("girdle MN", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11)
panelAll <- ggdraw() + draw_image(imgAll, scale = 1) + 
    draw_label("girdle all, ventral", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11) +
  draw_label(expression(paste("50 ", mu, " m")), x = 0.75, y = 0.04, size = 10)  +
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

panelAll_lat <- ggdraw() + draw_image(imgAll_lat, scale = 1) + 
    draw_label("girdle all, lateral", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11)  +
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
    draw_label("p", x = 0.1, y = 0.79, size = 8)   +
    geom_segment(aes(x = 0.78,
                     y = 0.93,
                     xend = 0.88,
                     yend = 0.93),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
    geom_segment(aes(x = 0.88,
                     y = 0.93,
                     xend = 0.78,
                     yend = 0.93),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
    draw_label("d", x = 0.76, y = 0.93, size = 8) +
    draw_label("v", x = 0.9, y = 0.93, size = 8) 

  
panelIN_ventr <- ggdraw() + draw_image(imgIN_ventr, scale = 1) + 
    draw_label("girdle IN, ventral", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11)  +
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
panelIN_ant <- ggdraw() + draw_image(imgIN_ant, scale = 1) + 
    draw_label("girdle IN, anterior", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11) +
    geom_segment(aes(x = 0.05,
                     y = 0.9,
                     xend = 0.05,
                     yend = 0.82),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
    geom_segment(aes(x = 0.05,
                     y = 0.82,
                     xend = 0.05,
                     yend = 0.9),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
    draw_label("d", x = 0.05, y = 0.93, size = 8) +
    draw_label("v", x = 0.05, y = 0.79, size = 8) 
panelIN_lat <- ggdraw() + draw_image(imgIN_lat, scale = 1) + 
    draw_label("girdle IN, lateral", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11)  +
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
    draw_label("p", x = 0.1, y = 0.79, size = 8)   +
    geom_segment(aes(x = 0.78,
                     y = 0.93,
                     xend = 0.88,
                     yend = 0.93),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
    geom_segment(aes(x = 0.88,
                     y = 0.93,
                     xend = 0.78,
                     yend = 0.93),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
    draw_label("d", x = 0.76, y = 0.93, size = 8) +
    draw_label("v", x = 0.9, y = 0.93, size = 8) 

panel_girdle <- ggdraw() + 
  draw_image(image_read("pictures/Full_connectome_girdle.png")) +
  draw_label("girdle cells", x = 0.5, y = 0.99, size = 11)

layout <- "
ABCD
####
EFGH
"
  
Fig14 <- panel_girdle + panelAll + panelAll_lat + panelSN + 
  panelIN_ant + panelIN_ventr + panelIN_lat + panelMN + 
  plot_layout(design = layout, heights = c(1, 0.05, 1),
              widths = c(1.1, 1, 1, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("Figures/Figure14.png", limitsize = FALSE, 
       units = c("px"), Fig14, width = 3000, height = 1640, bg='white')

  
ggsave("Figures/Figure14.pdf", limitsize = FALSE, 
         units = c("px"), Fig14, width = 3000, height = 1640)
  
  
  
