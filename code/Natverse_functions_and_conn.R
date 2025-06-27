#packages, functions and CATMAID connectivity info used for the figures of the Platynereis 3d connectome paper

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.

#create directory for R-generated pictures for figure panels (ignored by git)
if (!dir.exists("pictures2")) {
  dir.create("pictures2")
  cat("Directory created:", "pictures2", "\n")
} 
if (!dir.exists("Figures")) {
  dir.create("Figures")
  cat("Directory created:", "Figures", "\n")
} 

# load packages
{
  library(catmaid)
#natverse commands
#http://natverse.org/rcatmaid/reference/index.html
#https://rdrr.io/github/natverse/nat/man/
  library(heatmaply)
  library(plotly)
  library(magick)
  options(nat.plotengine = 'rgl')
  require("graphics")
  library(cowplot)
  library(png)
  library(tidyverse)
  library(leiden, attach.required = FALSE)
  library(RColorBrewer)
  library(igraph)
  library(networkD3)
  library(visNetwork)
  library(webshot2)
  library(patchwork)
  library(rgexf)
  library(tidygraph)
}


#for revision only, then delete and decomment public server
{
  # can run this separate file using source function
  conn <- source("~/R/conn.R")
  #for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
  #for this we configure to http/1.1
  conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
}

# catmaid connection
#conn_http1 <- catmaid_login(
#  server="https://catmaid.jekelylab.ex.ac.uk/", 
#  authname="AnonymousUser",
#  config=httr::config(ssl_verifypeer=0, http_version=1)
#)

#save session info and Rstudio version info for reproducibility
writeLines(capture.output(sessionInfo()), "code/sessionInfo.txt")
writeLines(capture.output(rstudioapi::versionInfo()), "code/versionInfo.txt")


#see the available volumes
catmaid_get_volumelist(conn = NULL, pid = 11)

#read volumes
#These volumes are 3D structures in the animal's body and provide background 
#for reconstructed neurons or help create images
{
  outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                                invertFaces = T, conn = NULL, pid = 11)
  yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                             invertFaces = T, conn = NULL, pid = 11)
  acicula <-   nlapply(read.neurons.catmaid("^acicula$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  chaeta = nlapply(read.neurons.catmaid("^chaeta$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  #these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
  bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  scalebar_50um_anterior = read.neurons.catmaid("scalebar_50um_anterior", pid=11)
  scalebar_50um_ventral = read.neurons.catmaid("scalebar_50um_ventral", pid=11)
#  scalebar_50um_ventral2 = read.neurons.catmaid("scalebar_50um_ventral2", pid=11)
  
  
}

#color palettes From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/ ------------------
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933',
               '#CC6677', '#882255', '#AA4499', '#DDDDDD')

# pie charts ------------------
blues <- brewer.pal(9, 'Blues')
bluepurple <- brewer.pal(9, 'BuPu')
oranges <- brewer.pal(9, 'YlOrRd')
greens <- brewer.pal(9, 'Greens')


#define background plotting function ------------------------
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


#adjusted background plotting function for mechanosensory circuits
plot_background_mech <- function(){plot_background_ventral()
  clear3d()
  nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=0,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.05,
         col='grey80')
  plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=F, lwd=3,
         rev = FALSE, fixup = F, add=T, alpha=0.1, col="grey50")
  par3d(zoom=0.58)
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

#same without view orientation and without yolk
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

#plotting function for ventral view with yolk and acicula
plot_background_ventral <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.05,
         col="#E2E2E2") 
  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="grey70")
  par3d(zoom=0.48)
}

plot_background_ventral_no_ac <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.05,
         col="#E2E2E2") 
  par3d(zoom=0.48)
}

#background plotting with lower z plane cut
plot_background_z2  <- function(x){
  nopen3d() # opens a pannable 3d window
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
  #we define a z clipping plane for the frontal view
  par3d(zoom=0.52)
  nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
  #z-axis clip
  clipplanes3d(0, 0, -1, 105000)
  #y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
}

#plotting function for ventral view with yolk (no acicula)
plot_background_ventral2 <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.05,
         col="#E2E2E2") 
  par3d(zoom=0.38)
}

#function to retrieve skids based on two annotations
{
  skids_by_2annotations <- function(annotation1,annotation2){
    skids1 <- catmaid_skids(annotation1, pid = 11)
    skids2 <- catmaid_skids(annotation2, pid = 11)
    intersect <- intersect(skids1, skids2)
    return(intersect)
  }
  
  #function to retrieve skids based on three annotations
  skids_by_3annotations <- function(annotation1, annotation2, annotation3){
    skids1 <- catmaid_skids(annotation1, pid = 11)
    skids2 <- catmaid_skids(annotation2, pid = 11)
    skids3 <- catmaid_skids(annotation3, pid = 11)
    intersect1_2 <- intersect(skids1, skids2)
    intersect1_2_3 <- intersect(intersect1_2, skids3)
    return(intersect1_2_3)
  }
}



############################################
#read and plot neurons function
read_plot_neurons <- function(annotation1, annotation2,annotation3,
                              annotation4,annotation5,annotation6){
  neuron1 = nlapply(read.neurons.catmaid(annotation1, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  neuron2 = nlapply(read.neurons.catmaid(annotation2, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  if (hasArg(annotation3)){
    neuron3 = nlapply(read.neurons.catmaid(annotation3, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  if (hasArg(annotation4)){
    neuron4 = nlapply(read.neurons.catmaid(annotation4, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  if (hasArg(annotation5)){
    neuron5 = nlapply(read.neurons.catmaid(annotation5, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  if (hasArg(annotation6)){
    neuron6 = nlapply(read.neurons.catmaid(annotation6, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  plot_background()
  plot3d(neuron1, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="#E69F00") 
  plot3d(neuron2, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="#0072B2")
  if (hasArg(annotation3)){
    plot3d(neuron3, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#CC79A7")
  }
  if (hasArg(annotation4)){
    plot3d(neuron4, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#009E73")
  }
  if (hasArg(annotation5)){
    plot3d(neuron5, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#D55E00") 
  }
  if (hasArg(annotation6)){
    plot3d(neuron6, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="grey50")
  }
}

#read and plot a neuron function
read_plot_neuron <- function(annotation, color){
  neuron1 = nlapply(read.neurons.catmaid(annotation, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  plot_background()
  plot3d(neuron1, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=4,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col=color)
}

#read and plot a neuron function
read_plot_neuron_ventral <- function(annotation, color){
  neuron1 = nlapply(read.neurons.catmaid(annotation, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  plot_background_ventral()
  plot3d(neuron1, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=4,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col=color)
}


#read and plot a neuron function - add to existing panel
read_plot_neuron_add <- function(annotation, color, alpha){
  neuron1 = nlapply(read.neurons.catmaid(annotation, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  plot3d(neuron1, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=4,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=alpha,
         col=color)
}

#plot neurons with a z cut below all head CR neurons

read_plot_neurons_z2 <- function(annotation1, annotation2,annotation3,
                                 annotation4,annotation5,annotation6){
  neuron1 = nlapply(read.neurons.catmaid(annotation1, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  neuron2 = nlapply(read.neurons.catmaid(annotation2, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  if (hasArg(annotation3)){
    neuron3 = nlapply(read.neurons.catmaid(annotation3, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  if (hasArg(annotation4)){
    neuron4 = nlapply(read.neurons.catmaid(annotation4, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  if (hasArg(annotation5)){
    neuron5 = nlapply(read.neurons.catmaid(annotation5, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  if (hasArg(annotation6)){
    neuron6 = nlapply(read.neurons.catmaid(annotation6, pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  }
  plot_background_z2()
  plot3d(neuron1, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="#E69F00") 
  plot3d(neuron2, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="#0072B2")
  if (hasArg(annotation3)){
    plot3d(neuron3, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#CC79A7")
  }
  if (hasArg(annotation4)){
    plot3d(neuron4, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#009E73")
  }
  if (hasArg(annotation5)){
    plot3d(neuron5, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="#D55E00") 
  }
  if (hasArg(annotation6)){
    plot3d(neuron6, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=2,
           rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
           col="grey50")
  }
}



#define function to retrieve skids from a neuron list based on one to three annotations
skids_by_annotation <- function(neuron_list,annotation1,annotation2,annotation3){
  skids1 <- unlist(lapply(neuron_list,function(x) x[x$annotation==annotation1,1]))
  if(missing(annotation2)){return(skids1) #if annotation2 is missing, will return skids matching the first annotation
  } else {
    skids2 <- unlist(lapply(neuron_list,function(x) x[x$annotation==annotation2,1]))
  }
  if(missing(annotation3)){return(intersect(skids1,skids2)) #if annotation3 is missing, will return skids matching annotations 1 and 2
  }   else  {
    skids3 <- unlist(lapply(neuron_list,function(x) x[x$annotation==annotation3,1]))
  }
  skids1_2 <- intersect(skids1,skids2)
  return(intersect(skids1_2,skids3)) #will return the shared skids between the three annotations  
}

