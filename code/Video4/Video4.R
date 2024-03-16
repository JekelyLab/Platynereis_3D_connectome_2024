#Generate Video4 of the Platynereis 3d connectome paper
#Gaspar Jekely  2024

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

# read cells --------------

stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
#these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))


INsplitPBant = nlapply(read.neurons.catmaid("^celltype78$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))

INsplitPB_RFYa = nlapply(read.neurons.catmaid("^celltype77$", pid=11),
                           function(x) smooth_neuron(x, sigma=6000))



mechanosensory_girdle_SN = nlapply(
  read.neurons.catmaid(
    skids_by_2annotations("mechanosensory_girdle", "Sensory neuron"), 
    pid=11
    ),
  function(x) smooth_neuron(x, sigma=6000)
  )

mechanosensory_girdle_IN = nlapply(
  read.neurons.catmaid(
    skids_by_2annotations("mechanosensory_girdle", "interneuron"), 
    pid=11
  ),
  function(x) smooth_neuron(x, sigma=6000)
)

mechanosensory_girdle_MN = nlapply(
  read.neurons.catmaid(
    skids_by_2annotations("mechanosensory_girdle", "motorneuron"), 
    pid=11
  ),
  function(x) smooth_neuron(x, sigma=6000)
)


# plot background -------------

plot_background_ventral()
par3d(windowRect = c(0, 0, 800, 800))
par3d(zoom=0.62)
for(i in 10:20){
  rgl.snapshot(paste("videoframes/Video10_", i, ".png", sep = ""))
}

plot3d(stomodeum, soma=T, lwd=1,
       add=T, alpha=0.2, col="#cccccc")
texts3d(115000,42000, 12000, text = "stomodeum", col='black', cex = 2)
  
for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video11B_", i, ".png", sep = ""))
  }
rgl.pop()
  
plot3d(INsplitPB_RFYa, soma=T, lwd=c(4,5),
       add=T, alpha=1, col=bluepurple[c(7,8)])
texts3d(105000,42000,10000, text = "INsplitPB-RF/Ya", col='black', cex = 2)

for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video11C_", i, ".png", sep = ""))
  }
rgl.pop()
  
plot3d(INsplitPBant, soma=T, lwd=c(4,6),
       add=T, alpha=1, col=bluepurple[c(6,9)])
texts3d(105000,42000, 12000, text = "INsplitPBant", col='black', cex = 2)
  
for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video11D_", i, ".png", sep = ""))
  }
rgl.pop()

plot3d(mechanosensory_girdle_IN, soma=T, lwd=1,
         add=T, alpha=1,
         col="#CC79A7") 
texts3d(125000,42000, 12000, text = "girdle IN", col='black', cex = 2)
  
for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video12A_", i, ".png", sep = ""))
  }
rgl.pop()

plot3d(mechanosensory_girdle_MN, soma=T, lwd=2,
         add=T, alpha=1, col="#0072B2")
texts3d(125000,42000, 12000, text = "girdle MN", col='black', cex = 2)
  
for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video13_", i, ".png", sep = ""))
  }
  
rgl.pop()
  

plot3d(mechanosensory_girdle_SN, soma=T, lwd=1,
       add=T, alpha=1,
       col="#E69F00")
texts3d(125000,42000, 12000, text = "girdle SN", col='black', cex = 2)

for(i in 10:20){
  rgl.snapshot(paste("videoframes/Video14_", i, ".png", sep = ""))
}
rgl.pop()

# export rotation by frame for video -----

um1 <- par3d()$userMatrix
rotation <- 100

for (l in 1:180){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um1 %*%rotationMatrix(pi*l/90, 0, 0, 1))
  print (l)
  #save a snapshot
  filename <- paste("./videoframes/Video31_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  rotation = rotation + 1 
}

# zoom in to head ------------
for (i in 1:36){
  nview3d(userMatrix = um1 %*%rotationMatrix(pi*180/90, 0, 0, 1))
  #progressive z clipping plane
  clipplanes3d(0, 0, -1, (230004-4091*i))

  filename <- paste("./videoframes/Video32_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  par3d(zoom=0.62-i/360)
  print(i)
}  

#rotation towards anterior
for(i in c(0:21)){print (i)
  nview3d(userMatrix = um1 %*%rotationMatrix(pi*180/90, 0, 0, 1)
                               %*%rotationMatrix(i/(21/(pi/2)), 1, 0, 0))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  
  filename <- paste("./videoframes/Video33_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}

for(i in 10:20){
  rgl.snapshot(paste("videoframes/Video34_", i, ".png", sep = ""))
}

close3d()

#read png files and write video
av::av_encode_video(paste('videoframes/', list.files("videoframes/", '*.png'), sep = ""), 
                    framerate = 10,
                    output = 'Videos/Video4.mp4')

unlink("videoframes", recursive = T)


