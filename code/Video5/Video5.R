#code to generate Video5 of the Platynereis 3d connectome paper
# showing a close-up of mechanosensory cells in the 2nd segment
#Gaspar Jekely 2022-2023

#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# create temp dir for videoframes
mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)

# read mech neuron types and muscle outlines ----------------
{
acicula_sg2l = nlapply(read.neurons.catmaid("^acicula_sg2l$", pid=11), 
                       function(x) smooth_neuron(x, sigma=6000))
chaeta_sg2l = nlapply(read.neurons.catmaid("^chaeta_sg2l$", pid=11), 
                      function(x) smooth_neuron(x, sigma=6000))
muscle_outline = read.neurons.catmaid("^muscle_outline$", pid=11)

muscle_outline = skids_by_3annotations("^muscle_outline$", "segment_2", "left_side")
muscle_outline = read.neurons.catmaid(muscle_outline, pid=11)

paratroch = skids_by_3annotations("^ciliated cell$", "segment_2", "left_side")
paratroch = read.neurons.catmaid(paratroch, pid=11)

spinGland_outline = read.neurons.catmaid("^OUTLINE spinning gland$", pid=11)

girdle = nlapply(read.neurons.catmaid("^mechanosensory_girdle$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))


PB2l = skids_by_3annotations("^celltype80$", "segment_2", "left_side")
PB2l = nlapply(read.neurons.catmaid(PB2l, pid=11),
        function(x) smooth_neuron(x, sigma=6000))
hCirrus_interparaPU2l = skids_by_3annotations("^celltype94$", "segment_2", "left_side")
hCirrus_interparaPU2l = nlapply(read.neurons.catmaid(hCirrus_interparaPU2l, pid=11),
               function(x) smooth_neuron(x, sigma=6000))
neuro_notoPU2l = skids_by_3annotations("^celltype95$", "segment_2", "left_side")
neuro_notoPU2l = nlapply(read.neurons.catmaid(neuro_notoPU2l, pid=11),
               function(x) smooth_neuron(x, sigma=6000))
paraPU2l = skids_by_3annotations("^celltype91$", "segment_2", "left_side")
paraPU2l = nlapply(read.neurons.catmaid(paraPU2l, pid=11),
               function(x) smooth_neuron(x, sigma=6000))

spinPU2l = skids_by_3annotations("^celltype89$", "segment_2", "left_side")
spinPU2l = nlapply(read.neurons.catmaid(spinPU2l, pid=11),
               function(x) smooth_neuron(x, sigma=6000))
SNFVa2l = skids_by_3annotations("^celltype143$", "segment_2", "left_side")
SNFVa2l = nlapply(read.neurons.catmaid(SNFVa2l, pid=11),
               function(x) smooth_neuron(x, sigma=6000))
interparaPM2l = skids_by_3annotations("^celltype81$", "segment_2", "left_side")
interparaPM2l = nlapply(read.neurons.catmaid(interparaPM2l, pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
CR2l = skids_by_3annotations("^celltype101$", "segment_2", "left_side")
CR2l = nlapply(read.neurons.catmaid(CR2l, pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
chaeMech2l = skids_by_3annotations("^celltype71$", "segment_2", "left_side")
chaeMech2l = nlapply(read.neurons.catmaid(chaeMech2l, pid=11),
               function(x) smooth_neuron(x, sigma=6000))
}
  
#define greyscale color list
greys <- grep('^grey', colours(), value=TRUE)

nopen3d() # opens a pannable 3d window
par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
par3d(zoom=0.7)
nview3d("ventral", extramat=rotationMatrix(pi/2, 0, 0, 1))


#z-axis clip
clipplanes3d(0, 0, -1, 155000)
#z-axis clip from top
clipplanes3d(0, 0, 1, -40000)
#x-axis clip
clipplanes3d(1, 0, 0.16, -80700)
#y-axis clip
clipplanes3d(0, 1, 0, -32700)
#x-axis clip
clipplanes3d(-1, 0, 0, 126000)

#y-axis clip
clipplanes3d(0, -1, 0, 149700)

# plot landmarks
{
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       add=T, alpha=0.07,
       col="#E2E2E2")
  plot3d(paratroch, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=0.8,
         col="grey")
  
  texts3d(125000, 50000, 125000, text = "paratroch", col='black', cex = 2)
  
  texts3d(125000, 55000, 150000, text = "dorsal", col='black', cex = 2)
  texts3d(125000, 125000, 150000, text = "ventral", col='black', cex = 2)
  texts3d(125000, 60000, 55000, text = "yolk", col='black', cex = 2)
  
  
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_1_1_", i, ".png", sep = ""))
  }
  rgl.pop()
  rgl.pop()
  rgl.pop()
  rgl.pop()
  
  plot3d(acicula_sg2l, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
       add=T, alpha=1,
       col="grey20")

  texts3d(125000, 60000, 55000, text = "aciculae", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_1_2_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(chaeta_sg2l, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       add=T, alpha=1,
       col="grey70")
  
  texts3d(125000, 60000, 55000, text = "chaetae", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_1_3_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(spinGland_outline, soma=T, lwd=25,
       add=T, alpha=0.2,
       col=blues[7])
  
  texts3d(125000, 60000, 55000, text = "spinGland", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_1_4_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  
  plot3d(muscle_outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
       add=T, alpha=1,
       col=sample(greys[50:90],length(muscle_outline), replace=TRUE) )
  
  texts3d(125000, 60000, 55000, text = "muscle", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_1_5_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=0.6, col="grey80")

  texts3d(125000, 70000, 55000, text = "nerve cord", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_1_6_", i, ".png", sep = ""))
  }
  rgl.pop()
}

#plot neurons
{
  plot3d(PB2l, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=0.8, col=oranges[8])
  
  texts3d(125000, 50000, 55000, text = "PB", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_1_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(spinPU2l, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=0.7, col=oranges[4])
  
  texts3d(125000, 53000, 55000, text = "spinPU", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_2_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(hCirrus_interparaPU2l, WithConnectors = F, WithNodes = F, soma=T, lwd=4,
         add=T, alpha=0.8, col=oranges[3])
  
  texts3d(125000, 60000, 55000, text = "interparaPU", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_3_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(neuro_notoPU2l, WithConnectors = F, WithNodes = F, soma=T, lwd=5,
         add=T, alpha=0.7, col=oranges[9])
  
  texts3d(125000, 60000, 55000, text = "neuro_notoPU", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_4_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(paraPU2l, WithConnectors = F, WithNodes = F, soma=T, lwd=5,
         add=T, alpha=0.8, col=oranges[2])
  
  texts3d(125000, 50000, 55000, text = "paraPU", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_5_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(SNFVa2l, WithConnectors = F, WithNodes = F, soma=T, lwd=3,
         add=T, alpha=1, col=bluepurple[6])
  
  texts3d(125000, 50000, 55000, text = "SNFVa", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_6_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(interparaPM2l, WithConnectors = F, WithNodes = F, soma=T, lwd=6,
         add=T, alpha=0.5, col=bluepurple[8])
  
  texts3d(125000, 60000, 55000, text = "interparaPM", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_7_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(chaeMech2l, WithConnectors = F, WithNodes = F, soma=T, lwd=6,
         add=T, alpha=0.8, col=oranges[5])
  
  texts3d(125000, 55000, 55000, text = "chaeMech", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_8_", i, ".png", sep = ""))
  }
  rgl.pop()
  
  plot3d(CR2l, WithConnectors = F, WithNodes = F, soma=T, lwd=6,
         add=T, alpha=0.8, col='red')
  
  texts3d(125000, 60000, 55000, text = "interparapodCR1", col='black', cex = 2)
  for(i in 10:20){
    rgl.snapshot(paste("videoframes/Video5_2_9_", i, ".png", sep = ""))
  }
  rgl.pop()
  
}




um <- par3d()$userMatrix
rotation <- 100

for (l in 1:180){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um %*%rotationMatrix(pi*l/90, 0, 0, -1))
  filename <- paste("videoframes/Video5_3_1_", rotation, ".png", sep = "")
  rgl.snapshot(filename)
  print (l)
  rotation = rotation + 1
}
close3d()


# read png files and write video -----------------
av::av_encode_video(paste('videoframes/', list.files("videoframes/", '*.png'), sep = ""), 
                    framerate = 10,
                    output = 'Videos/Video5.mp4')

# delete temp files
unlink("videoframes", recursive = T)


