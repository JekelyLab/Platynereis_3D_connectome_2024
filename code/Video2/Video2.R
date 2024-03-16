#Generate Video of head neuropils for the Platynereis 3d connectome paper
#Gaspar Jekely 2021-2023

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)

# load nat and all associated packages, incl catmaid
{
  library(catmaid)
  library(magick)
  library(nat)
  options(nat.plotengine = 'rgl')
  require("graphics")
  library(RColorBrewer)
  library(png)
  library(cowplot)
  library(ggplot2)
}

mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)

# catmaid connection, needs username, password AND token - weird!
{
  # can run this separate file using source function
  conn <- source("~/R/conn.R")
  #for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
  #for this we configure to http/1.1
  conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
}

#see the available volumes
catmaid_get_volumelist(conn = NULL, pid = 11)

#read volumes
{
  outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                                invertFaces = T, conn = NULL, pid = 11)
  yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                             invertFaces = T, conn = NULL, pid = 11)
  #these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
  bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  scalebar_50um_anterior = read.neurons.catmaid("^scalebar_50um_anterior$", pid=11)
}

#read sensory cell clusters
{
  Dorsal_SO = nlapply(read.neurons.catmaid("^Dorsal_sensory_organ$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  Adult_eye = nlapply(read.neurons.catmaid("^celltype1$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  nuchal = nlapply(read.neurons.catmaid("^nuchal organ$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  DLSO = nlapply(read.neurons.catmaid("Dorsolateral_sense_organs", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  palp = nlapply(read.neurons.catmaid("^palp sensory neuron$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  antenna = nlapply(read.neurons.catmaid("^antenna$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  MBSN = nlapply(read.neurons.catmaid("^MBSN$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  eyespot = nlapply(read.neurons.catmaid("^eyespot_PRC$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  ANS_SN = nlapply(read.neurons.catmaid("^apical sensory organ$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  HSG3 = nlapply(read.neurons.catmaid("head sensory group 3", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  cMS_cell = nlapply(read.neurons.catmaid("central_MS_cell", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  headPU = nlapply(read.neurons.catmaid("headPU", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  antler = nlapply(read.neurons.catmaid("celltype26", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  vMIP = nlapply(read.neurons.catmaid("celltype24", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  SNbicil = nlapply(read.neurons.catmaid("celltype50", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SNcirri = nlapply(read.neurons.catmaid("anterior pair of tentacular cirri", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
}

#read head interneuron clusters
{
  visualIN = nlapply(read.neurons.catmaid("visual_interneuron", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  ANS_IN = nlapply(read.neurons.catmaid("ns_plexusIN", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INdesc_decus = nlapply(read.neurons.catmaid("INdesc-decuss-head", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  INcentrHead = nlapply(read.neurons.catmaid("central head IN cluster", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  MBprojIN = nlapply(read.neurons.catmaid("MBprojIN", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  MBintrIN = nlapply(read.neurons.catmaid("MBintrIN", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  MBprojMouth = nlapply(read.neurons.catmaid("MBprojMouth", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
}

#read head MN clusters
{
  vMN = nlapply(read.neurons.catmaid("vMN", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  ciliomotor = nlapply(read.neurons.catmaid("ciliomotor_head", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  MNgland = nlapply(read.neurons.catmaid("celltype166", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  MNheadV = nlapply(read.neurons.catmaid("celltype178", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
}

#read MB cells and minicircuit cell groups
{
  SNhorn = nlapply(read.neurons.catmaid("^celltype17$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SNhook = nlapply(read.neurons.catmaid("^celltype18$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SNlasso = nlapply(read.neurons.catmaid("^celltype20$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SNtorii = nlapply(read.neurons.catmaid("^celltype187$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  INlasso = nlapply(read.neurons.catmaid("^celltype21$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INfoot = nlapply(read.neurons.catmaid("^celltype116$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INdecussfoot = nlapply(read.neurons.catmaid("^celltype152$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  INPDF_lc = nlapply(read.neurons.catmaid("^celltype184$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))  
  INhorn = nlapply(read.neurons.catmaid("^celltype183$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INdecusshook = nlapply(read.neurons.catmaid("^celltype153$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  INlat2 = nlapply(read.neurons.catmaid("^celltype121$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INRGWa = nlapply(read.neurons.catmaid("^celltype6$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INW = nlapply(read.neurons.catmaid("^celltype117$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  INZ = nlapply(read.neurons.catmaid("^celltype127$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  INproT2 = nlapply(read.neurons.catmaid("^celltype125$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INsn = nlapply(read.neurons.catmaid("^celltype57$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  INUturn = nlapply(read.neurons.catmaid("^celltype144$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INUturnMB = nlapply(read.neurons.catmaid("^celltype190$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INbigloop = nlapply(read.neurons.catmaid("^celltype140$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMB = nlapply(read.neurons.catmaid("^mushroom body neuron type$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  INSturn = nlapply(read.neurons.catmaid("^celltype122$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  MNant = nlapply(read.neurons.catmaid("^celltype19$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  vMN1_2 = nlapply(read.neurons.catmaid("^vMN1-2$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  MBprojIN = nlapply(read.neurons.catmaid("^MBprojIN$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  MBmouth = nlapply(read.neurons.catmaid("^celltype189$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  MBintrIN = nlapply(read.neurons.catmaid("^MBintrIN$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
}

#read all episphere cells
episphere = nlapply(read.neurons.catmaid("^episphere$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))

#define background plotting function
plot_background <- function(x){
  nopen3d() # opens a pannable 3d window
  #plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
  #      rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.07,
  #     col="#E2E2E2") 
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  #  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  #  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
  #we define a z clipping plane for the frontal view
  par3d(zoom=0.52)
  nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
  #z-axis clip
  clipplanes3d(0, 0, -1, 65000)
  #y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
}

#same without view orientation  -- to be used in movie snapshot generation
plot_background_no_view_change <- function(x){
  clear3d()  
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
  #z-axis clip
  clipplanes3d(0, 0, -1, 65000)
  #y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
}

plot_background_no_view_change_no_yolk <- function(x){
  clear3d()  
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  #z-axis clip
  clipplanes3d(0, 0, -1, 65000)
  #y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
}

#cb friendly color palette
{
  #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                 "#CC79A7", "#000000")
  Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933',
                 '#CC6677', '#882255', '#AA4499', '#DDDDDD')
  library(RColorBrewer)
  display.brewer.all(colorblindFriendly = TRUE)
  brewer12 <- brewer.pal(12, 'Paired')
}

# plot sensory neurons in head ------------------

#plot sensory cells without soma
plot_SN_No_soma <- function(x){
  plot3d(Dorsal_SO, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[2])
  plot3d(MBSN, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[1])
  plot3d(Adult_eye, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[3])
  plot3d(nuchal, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[4])
  plot3d(DLSO, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[7])
  plot3d(palp, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[6])
  plot3d(antenna, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[5])
  plot3d(eyespot, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, alpha=0.5, col="black")
  plot3d(ANS_SN, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Tol_muted[3])
  plot3d(SNcirri, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col=Tol_muted[4])
}

#add a text label
add_text_SN <- function(x){
  texts3d(69000,16000, 30000, text = "DSO", col='black', cex = 1.5)
  texts3d(35000,49000, 27000, text = "eye", col='black', cex = 1.5)
  texts3d(23000,48000, 52000, text = "nuchal", col='black', cex = 1.5)
  texts3d(46000,20000, 30000, text = "DLSO", col='black', cex = 1.5)
  texts3d(42000,115000, 7000, text = "palp", col='black', cex = 1.5)
  texts3d(87500,81000, 4000, text = "antenna", col='black', cex = 1.5)
  texts3d(20500,78500, 15000, text = "MBSN", col='black', cex = 1.5)
  texts3d(72500,42500, 5000, text = "ANS_SN", col='black', cex = 1.5)
  texts3d(24000,67000, 20000, text = "eyespot", col='black', cex = 1.5)
  texts3d(27300,114500, 47000, text = "cirrus", col='black', cex = 1.5)
}

#plot sensory cells with soma
plot_SN_soma <- function(x, alpha){ alpha=alpha
plot3d(Dorsal_SO, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[2])
plot3d(MBSN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[1])
plot3d(Adult_eye, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[3])
plot3d(nuchal, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[4])
plot3d(DLSO, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[7])
plot3d(palp, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[6])
plot3d(antenna, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[5])
plot3d(eyespot, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col="black")
plot3d(ANS_SN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[3])
plot3d(SNcirri, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[4])
}

plot_background()
plot_SN_No_soma()
add_text_SN()
um <- par3d()$userMatrix
rotation <- 101

#export rotation by frame for video
for (l in 1:30){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin1_SN_1_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  #reopen image and add label (this will not rotate in 3d)
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain sensory neurons", x = 0.3, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}

for (l in 30:1){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin1_SN_2_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain sensory neurons", x = 0.3, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}


plot_SN_soma(x, alpha=1)

for (l in 1:30){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, -1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin1_SN_3_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)

  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain sensory neurons", x = 0.3, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}



for (l in 30:1){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, -1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin1_SN_4_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain sensory neurons", x = 0.3, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}


close3d()

# plot interneurons in head ----------

#plot IN without soma
plot_IN_no_soma <- function(x){ alpha=1 
plot3d(visualIN, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[2])
plot3d(ANS_IN, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[8])
plot3d(INdesc_decus, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[4])
plot3d(MBintrIN, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[1])
plot3d(MBprojIN, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[5])
plot3d(INcentrHead, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[3])
}

#plot cells with soma
plot_IN_soma <- function(x,alpha){ alpha=alpha 
plot3d(visualIN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha+0.2, col=Okabe_Ito[2])
plot3d(ANS_IN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha+0.2, col=Okabe_Ito[8])
plot3d(INdesc_decus, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[4])
plot3d(MBintrIN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[1])
plot3d(MBprojIN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[5])
plot3d(INcentrHead, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[3])
}


add_text_IN <- function(x){
  texts3d(59000,32000, 20000, text = "INcentr", col='black', cex = 1.5)
  texts3d(35500,44000, 42000, text = "INeye", col='black', cex = 1.5)
  texts3d(87500,78000, 4000, text = "INdesc-decuss", col='black', cex = 1.5)
  texts3d(22500,83500, 30000, text = "MBintrIN", col='black', cex = 1.5)
  texts3d(25500,68500, 35000, text = "MBprojIN", col='black', cex = 1.5)
  texts3d(74500,36000, 20000, text = "ANS_IN", col='black', cex = 1.5)
}

plot_background()
plot_IN_no_soma()
add_text_IN()



#export rotation by frame for video
for (l in 1:30){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin2_IN_1_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  #reopen image and add label (this will not rotate in 3d)
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain interneurons", x = 0.25, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}

for (l in 30:1){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin2_IN_2_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain interneurons", x = 0.25, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}


plot_IN_soma(x, alpha=1)

for (l in 1:30){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, -1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin2_IN_3_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain interneurons", x = 0.25, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}

for (l in 30:1){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, -1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin2_IN_4_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain interneurons", x = 0.25, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}

close3d()


# plot MNs in head --------------------


add_text_MN <- function(x){
  texts3d(52500,53800, 30000, text = "MNant", col='black', cex = 1.5)
  texts3d(38500,44700, 33000, text = "Ser-h1", col='black', cex = 1.5)
  texts3d(62000,95300, 7000, text = "vMN", col='black', cex = 1.5)
  texts3d(75000,85000, 9000, text = "MNheadV", col='black', cex = 1.5)
  texts3d(38000,120500, 30000, text = "cMN", col='black', cex = 1.5)
  texts3d(85500,38000, 20000, text = "MNgland", col='black', cex = 1.5)
  texts3d(78500,51400, 20000, text = "MC", col='black', cex = 1.5)
  texts3d(68000,100300, 7000, text = "MNakro", col='black', cex = 1.5)
}

#plot cells with soma
plot_MN_soma <- function(x, alpha){ alpha=alpha
plot3d(vMN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[1])
plot3d(ciliomotor, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[5])
plot3d(MNheadV, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[7])
plot3d(MNgland, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[8])
}


plot_background()
plot_MN_no_soma()
add_text_MN()

plot_MN_soma(x, alpha=1)

#export rotation by frame for video
for (l in 1:30){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin3_MN_1_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain motoneurons", x = 0.24, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}

for (l in 30:1){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin3_MN_2_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("brain motoneurons", x = 0.24, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}

close3d()


# plot MB neurons in head --------------


#plot cells with soma
plot_MB_soma <- function(x, alpha){ alpha=alpha
plot3d(MBintrIN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[1])
plot3d(MBprojIN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha+0.3, col=Okabe_Ito[2])
plot3d(MBmouth, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col=Okabe_Ito[3])
plot3d(MBSN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Okabe_Ito[7])
}


add_text_MB <- function(x){
  texts3d(52500,72000, 17000, text = "MBintrIN", col='black', cex = 1.5)
  texts3d(49500,64500, 35000, text = "MBprojIN", col='black', cex = 1.5)
  texts3d(32500,89000, 18000, text = "MBSN", col='black', cex = 1.5)
}

plot_background()
plot_MB_no_soma()
add_text_MB()

plot_MB_soma(x, alpha=1)

#export rotation by frame for video
for (l in 1:30){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, -1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin4_MB_1_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("mushroom body", x = 0.23, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}

for (l in 30:1){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, -1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin4_MB_2_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  
  img <- readPNG(filename)
  panel <- ggdraw() + draw_image(img) +
    draw_label("mushroom body", x = 0.23, y = 0.98, size = 8)
  ggsave(filename, units = c("px"), panel, width = 800, height = 800)
  
  rotation = rotation + 1 
}


close3d()


# unique head SNs --------------

plot_SN_unique_soma <- function(x, alpha){ alpha=alpha
plot3d(HSG3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[4], depth_mask = FALSE)
plot3d(cMS_cell, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col='red', depth_mask = FALSE)
plot3d(headPU, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[6], depth_mask = FALSE)
plot3d(antler, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[7], depth_mask = FALSE)
plot3d(vMIP, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[8], depth_mask = FALSE)
plot3d(SNbicil, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=alpha, col='#44BB99', depth_mask = FALSE)
}



# plot everything and rotate ---------------
plot_background()
#export rotation by frame for video

plot3d(episphere, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=0.2, col="grey")
plot_MB_soma(x,alpha=1) 
plot_SN_soma(x, alpha=1)
plot_IN_soma(x, alpha=1)
plot_MN_soma(x,alpha=1)
plot_SN_unique_soma(x,alpha=1)

for (l in 1:180){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 1, 0))
  #save a snapshot
  filename <- paste("./videoframes/Video_Head_spin5_all_1_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  rotation = rotation + 1 
}

close3d()

#read png files and write video
av::av_encode_video(paste('videoframes/', list.files("videoframes/", '*.png'), sep = ""), 
                    framerate = 12,
                    output = 'Videos/Video2.mp4')


unlink("videoframes", recursive = T)
