#Generate the anatomical overview Video1 of the Platynereis 3d connectome paper
#Gaspar Jekely  2021-2023

# load nat and all associated packages, incl catmaid
source("code/Natverse_functions_and_conn.R")

#create temp dir to store video frames
mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)

# read volumes -----------
{
  outline <- catmaid_get_volume(1, rval = c("mesh3d", "catmaidmesh", "raw"),
                                invertFaces = T, conn = NULL, pid = 11)
  yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                             invertFaces = T, conn = NULL, pid = 11)
}

# read cells
--------------
{
  Sensoryneuron = nlapply(read.neurons.catmaid("^connectome_Sensory_neuron$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  Motorneuron = nlapply(read.neurons.catmaid("^connectome_motorneuron$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  Interneuron = nlapply(read.neurons.catmaid("^connectome_interneuron$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  chaeta = nlapply(read.neurons.catmaid("^chaeta$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  acicula = nlapply(read.neurons.catmaid("^acicula$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  muscle = nlapply(read.neurons.catmaid("^muscle$", pid=11, conn = conn_http1,
                                        fetch.annotations = FALSE),
                   function(x) smooth_neuron(x, sigma=6000))
  endoderm = nlapply(read.neurons.catmaid("^endoderm?$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  epithelia = nlapply(read.neurons.catmaid("^epithelia_cell$", pid=11, conn = conn_http1,
                                           fetch.annotations = FALSE),
                      function(x) smooth_neuron(x, sigma=6000))
  gland = nlapply(read.neurons.catmaid("^gland cell$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  Ciliary_band_cell = nlapply(read.neurons.catmaid("^Ciliary_band_cell$", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
  glia = nlapply(read.neurons.catmaid("^glia cell$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  pnb = nlapply(read.neurons.catmaid("^pnb$", pid=11, conn = conn_http1,
                                     fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
  notopodium = nlapply(read.neurons.catmaid("^notopodium$", pid=11, conn = conn_http1,
                                            fetch.annotations = FALSE),
                       function(x) smooth_neuron(x, sigma=6000))
  neuropodium = nlapply(read.neurons.catmaid("^neuropodium$", pid=11, conn = conn_http1,
                                             fetch.annotations = FALSE),
                        function(x) smooth_neuron(x, sigma=6000))
  
   follicle = nlapply(read.neurons.catmaid("^follicle cell$", pid=11, conn = conn_http1,
                                             fetch.annotations = FALSE),
                        function(x) smooth_neuron(x, sigma=6000))
 
   pigment = nlapply(read.neurons.catmaid("^pigment cell$", pid=11, conn = conn_http1,
                                             fetch.annotations = FALSE),
                        function(x) smooth_neuron(x, sigma=6000))
   
   all_cells = nlapply(read.neurons.catmaid("^with_soma$", pid=11, conn = conn_http1,
                                          fetch.annotations = FALSE),
                     function(x) smooth_neuron(x, sigma=6000))
  #these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
  bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
}


# Export images for Video1 ---------------

# extract connectors for left only to be able to plot them by unique colours ------------
{
SN_conn <- connectors(Sensoryneuron)
str(SN_conn)
presyn_SN_conn <- SN_conn[SN_conn$prepost == 0,]
postsyn_SN_conn <- SN_conn[SN_conn$prepost == 1,]

MN_conn <- connectors(Motorneuron)
presyn_MN_conn <- subset(MN_conn, prepost == 0)
postsyn_MN_conn <- subset(MN_conn, prepost == 1)

IN_conn <- connectors(Interneuron)
presyn_IN_conn <- IN_conn[IN_conn$prepost == 0,]
postsyn_IN_conn <- IN_conn[IN_conn$prepost == 1,]
}

{
  skids_to_plot_left <- skids_by_2annotations('Sensory neuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  SN_conn_left <- connectors(skeletons_to_plot_left)
  str(SN_conn_left)
  presyn_SN_conn_left <- SN_conn_left[SN_conn_left$prepost == 0,]
  postsyn_SN_conn_left <- SN_conn_left[SN_conn_left$prepost == 1,]
  
  skids_to_plot_left <- skids_by_2annotations('interneuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  
  IN_conn_left <- connectors(skeletons_to_plot_left)
  presyn_IN_conn_left <- IN_conn_left[IN_conn_left$prepost == 0,]
  postsyn_IN_conn_left <- IN_conn_left[IN_conn_left$prepost == 1,]
  
  skids_to_plot_left <- skids_by_2annotations('motorneuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  #we can also subset with the subset shorthand function  
  MN_conn_left <- connectors(skeletons_to_plot_left)
  presyn_MN_conn_left <- subset(MN_conn_left, prepost == 0)
  postsyn_MN_conn_left <- subset(MN_conn_left, prepost == 1)
}

# function to retrieve skids based on two annotations ------------

skids_by_2annotations <- function(annotation1,annotation2){
  annotations_cells = list()
  annotation1 <- paste("^", annotation1, "$", sep="")
  annotations_cells[[1]] <- catmaid_get_annotations_for_skeletons(annotation1, pid = 11)
  #we retrieve those skeletons that are also annotated with right_side
  return(unlist(lapply(annotations_cells,function(x) x[x$annotation==annotation2,1])))
}


# define two windows, plot background -------------
{
  nopen3d() # opens a pannable 3d window
  mfrow3d(1, 2)  #defines the two scenes
  par3d(windowRect = c(20, 30, 1200, 800)) #to define the size of the rgl window
  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
  par3d(zoom=0.48)
  next3d(clear=F)
  nview3d("left", extramat=rotationMatrix(-pi/2, pi, -0.2, 0))
  #clipplanes3d(1, 0, 0.16, -75700)
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.1,
         col="#E2E2E2") 
  par3d(zoom=0.48)
}

texts3d(99000,16000, 5000, text = "yolk", col='black', cex = 2)

for(i in 10:20){
  rgl.snapshot(paste("videoframes/Video10_", i, ".png", sep = ""))
}
rgl.pop()



# plot ventral views and side views of left side cells only ------
{
  next3d(clear=F)
  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="grey35")
  skids_to_plot_left <- skids_by_2annotations('acicula','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1, col="grey35")
  
  texts3d(99000,20000, 5000, text = "aciculae", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video11_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col="grey70") 
  skids_to_plot_left <- skids_by_2annotations('chaeta','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         add=T, alpha=1, col="grey70")
  
  texts3d(99000,20000, 5000, text = "chaetae", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video12A_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(follicle, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col=hcl.colors(length(follicle), palette='Oranges')) 
  skids_to_plot_left <- skids_by_2annotations('follicle cell','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col=hcl.colors(length(skeletons_to_plot_left), palette='Oranges'))
  
  texts3d(99000,26000, 5000, text = "follicle cells", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video12B_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(endoderm, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="grey65")
  skids_to_plot_left <- skids_by_2annotations('endoderm','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="grey65")
  
  texts3d(99000,23000, 5000, text = "endoderm", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video13_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="#4477AA")
  skids_to_plot_left <- skids_by_2annotations('stomodeum','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="#4477AA")
  
  texts3d(99000,25000, 5000, text = "stomodeum", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video14_", i, ".png", sep = ""))
  }
   rgl.pop()
}


# plot synapses ---------------
{
  next3d(clear=F)
  plot3d(presyn_SN_conn$x, presyn_SN_conn$y, presyn_SN_conn$z, size=5, alpha=1, col="#E69F00", add=T)
  next3d(clear=F)
  plot3d(presyn_SN_conn_left$x, presyn_SN_conn_left$y, presyn_SN_conn_left$z, size=5, alpha=1, col="#E69F00", add=T)
  
  texts3d(99000,30000, 5000, text = "SN synapses", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video15_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(presyn_MN_conn$x, presyn_MN_conn$y, presyn_MN_conn$z, size=5, alpha=1, col="#0072B2", add=T)
  next3d(clear=F)
  plot3d(presyn_MN_conn_left$x, presyn_MN_conn_left$y, presyn_MN_conn_left$z, size=5, alpha=1, col="#0072B2", add=T)
  
  texts3d(99000,30000, 5000, text = "MN synapses", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video16_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(presyn_IN_conn$x+1, presyn_IN_conn$y, presyn_IN_conn$z, size=5, alpha=1, col="#CC79A7", add=T)
  next3d(clear=F)
  plot3d(presyn_IN_conn_left$x+1, presyn_IN_conn_left$y, presyn_IN_conn_left$z, size=5, alpha=1, col="#CC79A7", add=T)
  
  texts3d(99000,30000, 5000, text = "IN synapses", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video17_", i, ".png", sep = ""))
  }
  
  rgl.pop()
}

# continue plotting cells -------
{
  next3d(clear=F)
  plot3d(Sensoryneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         add=T, alpha=1, col="#E69F00")
  skids_to_plot_left <- skids_by_2annotations('Sensory neuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         add=T, alpha=1, col="#E69F00")
  par3d(zoom=0.48)
  
  texts3d(99000,36000, 5000, text = "sensory neurons", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video18_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(Motorneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col="#0072B2",)
  skids_to_plot_left <- skids_by_2annotations('motorneuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col="#0072B2",)
  
  texts3d(99000,32000, 5000, text = "motoneurons", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video19_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(Interneuron, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col="#CC79A7")
  skids_to_plot_left <- skids_by_2annotations('interneuron','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=1, col="#CC79A7")
  
  texts3d(99000,30000, 5000, text = "interneurons", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video20_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(gland, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="#0072B2")
  skids_to_plot_left <- skids_by_2annotations('gland cell','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="#0072B2")
  
  texts3d(99000,16000, 5000, text = "glands", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video21_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(Ciliary_band_cell, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="#E69F00") 
  skids_to_plot_left <- skids_by_2annotations('Ciliary_band_cell','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="#E69F00") 
  
  texts3d(99000,27000, 5000, text = "ciliary bands", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video22_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(muscle, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(muscle), palette='Reds'))
  skids_to_plot_left <- skids_by_2annotations('muscle','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(skeletons_to_plot_left), palette='Reds'))
  
  texts3d(99000,20000, 5000, text = "muscles", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video23_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(glia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(glia), palette='Peach'))
  skids_to_plot_left <- skids_by_2annotations('glia cell','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(skeletons_to_plot_left), palette='Peach'))
  
  texts3d(99000,16000, 5000, text = "glia", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video24_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(pigment, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="black")
  skids_to_plot_left <- skids_by_2annotations('pigment cell','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col="black")
  
  texts3d(99000,28000, 5000, text = "pigment cell", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video25_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(pnb, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(pnb), palette='Mint', rev=T))
  skids_to_plot_left <- skids_by_2annotations('pnb','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(skeletons_to_plot_left), palette='Mint', rev=T))
  
  texts3d(99000,25000, 5000, text = "developing", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video26_", i, ".png", sep = ""))
  }
  
  rgl.pop()

  
  
  next3d(clear=F)
  plot3d(notopodium, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="#009E73")
  skids_to_plot_left <- skids_by_2annotations('notopodium','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="#009E73")
  
  texts3d(99000, 24000, 5000, text = "notopodium", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video27_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(neuropodium, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="#CC79A7")
  skids_to_plot_left <- skids_by_2annotations('neuropodium','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=1, col="#CC79A7")
  
  texts3d(99000, 27000, 5000, text = "neuropodium", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video28_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  next3d(clear=F)
  plot3d(epithelia, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(epithelia), palette='Blues')) 
  skids_to_plot_left <- skids_by_2annotations('epithelia_cell','left_side')
  skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                        fetch.annotations = FALSE),
                                   function(x) smooth_neuron(x, sigma=6000))
  next3d(clear=F)
  plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         rev = FALSE, fixup = F, add=T, forceClipregion = TRUE, alpha=1,
         col=hcl.colors(length(skeletons_to_plot_left), palette='Blues')) 
  
  texts3d(99000, 20000, 5000, text = "epidermis", col='black', cex = 2)
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video29_", i, ".png", sep = ""))
  }
  
  rgl.pop()
  
  #rotation
  

}

# add all cells -----------

# first all neurons with soma
next3d(clear=F)
plot3d(Sensoryneuron, WithConnectors = F, WithNodes = T, soma=F, lwd=1,
       add=T, alpha=1, col="#E69F00")
skids_to_plot_left <- skids_by_2annotations('Sensory neuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="#E69F00")
par3d(zoom=0.48)

next3d(clear=F)
plot3d(Motorneuron, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#0072B2",)
skids_to_plot_left <- skids_by_2annotations('motorneuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#0072B2",)

next3d(clear=F)
plot3d(Interneuron, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#CC79A7")
skids_to_plot_left <- skids_by_2annotations('interneuron','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       rev = FALSE, fixup = F, add=T, alpha=1, col="#CC79A7")


#then all remaining cells

next3d(clear=F)
plot3d(all_cells, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="grey95")
skids_to_plot_left <- skids_by_2annotations('with_soma','left_side')
skeletons_to_plot_left = nlapply(read.neurons.catmaid(skids_to_plot_left, pid=11, conn = conn_http1,
                                                      fetch.annotations = FALSE),
                                 function(x) smooth_neuron(x, sigma=6000))
next3d(clear=F)
plot3d(skeletons_to_plot_left, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=1, col="grey95")

texts3d(99000, 18000, 5000, text = "all cells", col='black', cex = 2)

for(i in 10:20){
  rgl.snapshot(paste("videoframes/Video30_", i, ".png", sep = ""))
}

rgl.pop()



# export rotation by frame for video -----

next3d(clear=F)
um1 <- par3d()$userMatrix
next3d(clear=F)
um2 <- par3d()$userMatrix
next3d(clear=F)

rotation <- 100

for (l in 1:180){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um1 %*%rotationMatrix(pi*l/90, 0, 0, 1))
  next3d(clear=F)
  nview3d(userMatrix = um2 %*%rotationMatrix(pi*l/90, 0, 0, 1))
  next3d(clear=F)
  print (l)
  #save a snapshot
  filename <- paste("./videoframes/Video31_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  rotation = rotation + 1 
}


close3d()

#read png files and write video
av::av_encode_video(paste('videoframes/', list.files("videoframes/", '*.png'), sep = ""), 
                    framerate = 10,
                    output = 'Videos/Video1.mp4')

unlink("videoframes", recursive = T)


