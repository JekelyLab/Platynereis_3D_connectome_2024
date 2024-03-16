#code to generate the mushroom body video for the Platynereis 3d connectome paper
#Gaspar Jekely 2022-2023

#clear memory, load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

#create temp dir to store video frames
mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)

# read neuron groups ------------------------------------------------------

{
  INMBtype7 = nlapply(read.neurons.catmaid("^celltype122$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype8 = nlapply(read.neurons.catmaid("^INMBtype8$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype10 = nlapply(read.neurons.catmaid("^INMBtype10$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  INpara = nlapply(read.neurons.catmaid("^celltype196$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INZ = nlapply(read.neurons.catmaid("^celltype127$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))

  INMBtype6 = nlapply(read.neurons.catmaid("^celltype120$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INlasso = nlapply(read.neurons.catmaid("^celltype21$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INW = nlapply(read.neurons.catmaid("^celltype117$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  INfoot = nlapply(read.neurons.catmaid("^celltype116$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  
  INMBdev = nlapply(read.neurons.catmaid("^INMBdev$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype9 = nlapply(read.neurons.catmaid("^celltype197$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  INMBtype4 = nlapply(read.neurons.catmaid("^INMBtype4$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype5 = nlapply(read.neurons.catmaid("^celltype121$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INW = nlapply(read.neurons.catmaid("^celltype117$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  Ser_h1 = nlapply(read.neurons.catmaid("^celltype8$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  
  INMBtype1 = nlapply(read.neurons.catmaid("^INMBtype1$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype2 = nlapply(read.neurons.catmaid("^celltype198$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype3 = nlapply(read.neurons.catmaid("^INMBtype3$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INdecusshook = nlapply(read.neurons.catmaid("^celltype153$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  
  
  SNhook = nlapply(read.neurons.catmaid("^celltype18$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SNMB2 = nlapply(read.neurons.catmaid("^celltype32$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  
  MS1 = nlapply(read.neurons.catmaid("^celltype35$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  INsn = nlapply(read.neurons.catmaid("^celltype57$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  INpreMN = nlapply(read.neurons.catmaid("^celltype23$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  SNlasso = nlapply(read.neurons.catmaid("^celltype20$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  
  SNlasso = nlapply(read.neurons.catmaid("^celltype20$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  SNhorn = nlapply(read.neurons.catmaid("^celltype17$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SNstiff = nlapply(read.neurons.catmaid("^celltype167$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INtorii = nlapply(read.neurons.catmaid("^celltype191$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SNPDF_pyg = nlapply(read.neurons.catmaid("^celltype115$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INleucoPU = nlapply(read.neurons.catmaid("^celltype150$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  MNladder = nlapply(read.neurons.catmaid("^celltype151$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  INsplitPUh = nlapply(read.neurons.catmaid("^celltype83$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  INdecussM = nlapply(read.neurons.catmaid("^celltype199$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  SNbronto = nlapply(read.neurons.catmaid("^celltype168$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  SN_NS19 = nlapply(read.neurons.catmaid("^celltype133$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SN_NS1 = nlapply(read.neurons.catmaid("^celltype110$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SN_NS18 = nlapply(read.neurons.catmaid("^celltype138$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SN_NS17 = nlapply(read.neurons.catmaid("^celltype131$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SN_NS27 = nlapply(read.neurons.catmaid("^celltype137$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  SNmus = nlapply(read.neurons.catmaid("^celltype25$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SNtrpa = nlapply(read.neurons.catmaid("^celltype190$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  
  
  SNMBdev = nlapply(read.neurons.catmaid("^SNMBdev$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  SNtorii = nlapply(read.neurons.catmaid("^celltype187$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INhorn = nlapply(read.neurons.catmaid("^celltype183$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INMBdescFMRF = nlapply(read.neurons.catmaid("^celltype185$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  IMBPDF = nlapply(read.neurons.catmaid("^celltype184$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INrope = nlapply(read.neurons.catmaid("^celltype58$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INMBdesc3 = nlapply(read.neurons.catmaid("^celltype194$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBdesc2 = nlapply(read.neurons.catmaid("^celltype195$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  INMBPDF = nlapply(read.neurons.catmaid("^celltype184$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INbigloop = nlapply(read.neurons.catmaid("^celltype140$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INUturnMB = nlapply(read.neurons.catmaid("^celltype192$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  MNant = nlapply(read.neurons.catmaid("^celltype19$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  SNgolden = nlapply(read.neurons.catmaid("^celltype16$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  SN_DLSO1.2 = nlapply(read.neurons.catmaid("^celltype30$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  SN_DLSO1.1NP = nlapply(read.neurons.catmaid("^celltype172$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  cPRC = nlapply(read.neurons.catmaid("^celltype5$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  palp = nlapply(read.neurons.catmaid("^palp sensory neuron$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  antenna = nlapply(read.neurons.catmaid("^antennae_cell$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  eye = nlapply(read.neurons.catmaid("^celltype1$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  
  MBintrIN   = nlapply(read.neurons.catmaid("^MBintrIN$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  
  
  MBSN = nlapply(read.neurons.catmaid("^MBSN$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  MBprojMouth = nlapply(read.neurons.catmaid("^MBprojMouth$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  MBprojIN = nlapply(read.neurons.catmaid("^MBprojIN$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  AOSN = nlapply(read.neurons.catmaid("^apical sensory organ$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))  
  girdle = nlapply(read.neurons.catmaid("^mechanosensory_girdle$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  MB = nlapply(read.neurons.catmaid("^mushroom body$", pid=11),
             function(x) smooth_neuron(x, sigma=6000)) 
  eye_pigment = nlapply(read.neurons.catmaid("^celltype_non_neuronal11$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
eyespot_pigment = nlapply(read.neurons.catmaid("^celltype_non_neuronal10$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
}


# background function for the MB video ------------------------------------

plot_background_MB <- function(x){
  nopen3d() # opens a pannable 3d window
  #we define a z clipping plane for the frontal view
  par3d(zoom=0.52)
  nview3d("frontal")
  #z-axis clip
  clipplanes3d(0, 0, -1, 65000)
  #y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
  par3d(windowRect = c(0, 0, 1200, 676)) #resize for frontal view
  plot3d(bounding_dots, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         add=T, alpha=1,
         col="white") 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.1,
         col="#E2E2E2") 
}


# define color palettes ---------------------------------------------------

display.brewer.all(colorblindFriendly = TRUE)
blues <- brewer.pal(9, 'Blues')
pie(rep(1,9),col=blues)
bluepurple <- brewer.pal(9, 'BuPu')
pie(rep(1,9),col=bluepurple)
oranges <- brewer.pal(9, 'YlOrRd')
pie(rep(1,9),col=oranges)


# overview of larva and zoom in to MB -------------------------------------

plot_background_MB()
{

for (i in 1:44){
  clear3d()
  nview3d("frontal", extramat=(rotationMatrix(pi/2, -1, 0, 0)
                               %*%rotationMatrix(pi/44, 0, 0, 1)))
  par3d(zoom=1.52-i/44)
  #progressive z clipping plane
  clipplanes3d(0, 0, -1, (240004-4091*i))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)

  plot3d(acicula, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
      add=T, alpha=i/3,
      col="grey40")

  plot3d(eye, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
       add=T, alpha=i/30, col=Okabe_Ito[8])

  plot3d(eyespot_pigment, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
       add=T, alpha=i/30, col=Okabe_Ito[8])

  plot3d(chaeta, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=i/10,
       col="grey70") 
  plot3d(MB, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
       add=T, alpha=i/10,
       col=oranges[4]) 
  plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
       add=T, alpha=i/20,
      col=blues[5]) 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
       add=T, alpha=0.05,
       col='grey70')
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
       add=T, alpha=0.05-0.0011*i,
       col='grey80')
  filename <- paste("videoframes/Video_MB_spin1A_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
 }  
}


#rotation towards anterior
for(i in c(0:21)){print (i)
  nview3d("frontal", extramat=(rotationMatrix(pi/2, -1, 0, 0)
                               %*%rotationMatrix(i/(21/(pi/2)), 1, 0, 0)
                               %*%rotationMatrix(pi/44, 0, 0, 1)))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  
  filename <- paste("videoframes/Video_MB_spin1B_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}


# add labels

for(i in c(22:30)){print (i)
  texts3d(105000, 32000, 5000, text = "eye", col='black', cex = 2.5)
  texts3d(135000, 45000, 5000, text = "mushroom body", col='black', cex = 2.5)
  filename <- paste("videoframes/Video_MB_spin1B_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}
rgl.pop()
rgl.pop()



#fade out MB - anterior view
for(i in c(31:41)){print (i)
  nview3d("frontal", extramat=(rotationMatrix(pi/2, -1, 0, 0)
                               %*%rotationMatrix(21/(21/(pi/2)), 1, 0, 0)
                               %*%rotationMatrix(pi/44, 0, 0, 1)))
  clear3d()
  
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  plot3d(eye, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         add=T, alpha=(41-i)/10, col=Okabe_Ito[8])
 
  plot3d(eyespot_pigment, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
         add=T, alpha=(41-i)/10, col=Okabe_Ito[8])
  plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         add=T, 
         alpha=(41-i)/10,
         col=blues[5]) 
  plot3d(MB, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, 
         alpha=(41-i)/10,
         col=oranges[4]) 
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col='grey70')
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  filename <- paste("videoframes/Video_MB_spin1B_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}

extramat_to_multiply <- rotationMatrix(21/(21/(pi/2)), -1, 0, 0)


#add INMBtype1 and 2
for(i in c(0:4)){print (i)
  clear3d()
#  plot3d(MB, WithConnectors = F, WithNodes = F, soma=T, lwd=0,
#         add=T, alpha=(0.1), col="grey80")
  plot3d(INMBtype1, WithConnectors = T, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[9])
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         add=T, 
         alpha=0.01,
         col=blues[5]) 
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype1", col='black', cex = 2)
  texts3d(99000,32000, 5000, text = "INMBtype2", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1C_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
  rgl.pop()

}

#add INMBtype3
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[5])
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype3", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1D_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBtype4
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[5])
  plot3d(INMBtype4, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[2])
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype4", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1E_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBtype5
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[5])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[2])
  plot3d(INMBtype5, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[4])
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype5", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1F_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBtype6
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[5])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[2])
  plot3d(INMBtype5, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype6, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=Okabe_Ito[2])
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype6", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1G_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBtype7
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype5, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype6, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[2])
  plot3d(INMBtype7, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[7])
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype7", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1H_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBtype8
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype5, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype6, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[2])
  plot3d(INMBtype7, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[7])
  plot3d(INMBtype8, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[8])
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype8", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1I_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}


#add INMBtype9
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype5, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype6, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[2])
  plot3d(INMBtype7, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[7])
  plot3d(INMBtype8, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[8])
  plot3d(INMBtype9, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col="blue")
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype9", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1J_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}



#add INMBtype10
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype5, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype6, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[2])
  plot3d(INMBtype7, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(1), col=blues[7])
  plot3d(INMBtype8, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(1), col=blues[8])
  plot3d(INMBtype9, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(1), col="blue")
  plot3d(INMBtype10, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5)-0.3, col="black")
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBtype10", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1K_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBdev
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[5])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[2])
  plot3d(INMBtype5, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype6, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[2])
  plot3d(INMBtype7, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(1), col=blues[7])
  
  plot3d(INMBtype8, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1, col=blues[8])
  plot3d(INMBtype9, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1, col="blue")
  plot3d(INMBtype10, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1-0.3, col="black")
  
  plot3d(INMBdev, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(0.2+i/5), col=blues[3])
  
  plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.001,
         col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMB developing", col='black', cex = 2)
  
  filename <- paste("videoframes/Video_MB_spin1L_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

# add MB cells ----------------------

plot_INint <- function(){
  plot3d(INMBtype1, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[5])
  plot3d(INMBtype2, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[9])
  plot3d(INMBtype3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[5])
  plot3d(INMBtype4, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[2])
  plot3d(INMBtype5, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=blues[4])
  plot3d(INMBtype6, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(1), col=Okabe_Ito[2])
  plot3d(INMBtype7, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=(1), col=blues[7])
  
  plot3d(INMBtype8, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1, col=blues[8])
  plot3d(INMBtype9, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1, col="blue")
  plot3d(INMBtype10, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1-0.3, col="black")
  
  plot3d(INMBdev, WithConnectors = F, WithNodes = F, soma=T, lwd=2,
         add=T, alpha=1, col=blues[3])
  
}

#MBSN types SNMB2, SNtorii, SNhook, SNlasso, 
#SN_NS1, SN_NS17, SN_NS18, SN_NS19, SN_NS27, 
#SNmus, SNhorn, SNtrpa
#SNMBdev

#add SNMB2
for(i in c(0:4)){print (i)
  clear3d()

  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(SNMB2, soma=T, WithConnectors = T, lwd=4, add=T, 
         alpha=(0.2+i/5), col=oranges[2])

  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNMB2", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1M_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SNtorii
for(i in c(0:4)){print (i)
  clear3d()
  
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(SNtorii, soma=T, WithConnectors = T, lwd=2, add=T, 
         alpha=(0.2+i/5), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNtorii", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1N_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SNhook
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(SNhook, soma=T, WithConnectors = T, lwd=5, add=T, 
         alpha=(0.2+i/5), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNhook", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1O_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SNlasso
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(SNlasso, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col=oranges[8])
  plot3d(SNhook, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNlasso", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1P_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SNmus
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(SNmus, soma=T, WithConnectors = T, lwd=6, add=T, 
         alpha=(0.2+i/5), col=oranges[4])
  plot3d(SNlasso, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[8])
  plot3d(SNhook, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNmus", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1Q_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SNhorn
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(SNhorn, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col=oranges[7])
  plot3d(SNmus, soma=T, lwd=6, add=T, 
         alpha=(1), col=oranges[4])
  plot3d(SNlasso, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[8])
  plot3d(SNhook, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNhorn", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1R_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SNtrpa
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(SNtrpa, soma=T, WithConnectors = T, lwd=2, add=T, 
         alpha=(0.2+i/5), col=oranges[6])
  plot3d(SNhorn, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[7])
  plot3d(SNmus, soma=T, lwd=6, add=T, 
         alpha=(1), col=oranges[4])
  plot3d(SNlasso, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[8])
  plot3d(SNhook, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNtrpa", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1S_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SN_NS
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  
  plot3d(SN_NS1, soma=T, WithConnectors = T, lwd=5, add=T, 
         alpha=(0.2+i/5), col=oranges[2])
  plot3d(SN_NS17, soma=T, WithConnectors = T, lwd=4, add=T, 
         alpha=(0.2+i/5), col=oranges[3])
  plot3d(SN_NS18, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col=oranges[4])
  plot3d(SN_NS19, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col=oranges[5])
  plot3d(SN_NS27, soma=T, WithConnectors = T, lwd=2, add=T, 
         alpha=(0.2+i/5), col=oranges[6])
  
  plot3d(SNtrpa, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[6])
  plot3d(SNhorn, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[7])
  plot3d(SNmus, soma=T, lwd=6, add=T, 
         alpha=(1), col=oranges[4])
  plot3d(SNlasso, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[8])
  plot3d(SNhook, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SN neurosecretory1, 17-19, 27", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1T_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add SNdev
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  
  plot3d(SNMBdev, soma=T, WithConnectors = T, lwd=2, add=T, 
         alpha=(0.2+i/5), col=oranges[5])
  
  plot3d(SN_NS1, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(SN_NS17, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[3])
  plot3d(SN_NS18, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[4])
  plot3d(SN_NS19, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SN_NS27, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[6])
  
  plot3d(SNtrpa, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[6])
  plot3d(SNhorn, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[7])
  plot3d(SNmus, soma=T, lwd=6, add=T, 
         alpha=(1), col=oranges[4])
  plot3d(SNlasso, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[8])
  plot3d(SNhook, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "SNMB developing", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin1U_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

plot_SNMB <- function(){
  plot3d(SNMBdev, soma=T, lwd=2, add=T, 
         alpha=1, col=oranges[5])
  plot3d(SN_NS1, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[2])
  plot3d(SN_NS17, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[3])
  plot3d(SN_NS18, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[4])
  plot3d(SN_NS19, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SN_NS27, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[6])
  plot3d(SNtrpa, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[6])
  plot3d(SNhorn, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[7])
  plot3d(SNmus, soma=T, lwd=6, add=T, 
         alpha=(1), col=oranges[4])
  plot3d(SNlasso, soma=T, lwd=3, add=T, 
         alpha=(1), col=oranges[8])
  plot3d(SNhook, soma=T, lwd=5, add=T, 
         alpha=(1), col=oranges[5])
  plot3d(SNtorii, soma=T, lwd=2, add=T, 
         alpha=(1), col=oranges[9])
  plot3d(SNMB2, soma=T, lwd=4, add=T, 
         alpha=(1), col=oranges[2])
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
}

# add INMBproj interneurons ------------------

#INrope, INtorii, INbigloop, INUturnMB,
#INhorn, INMBPDF
#INMBdesc2, INMBdesc3, INMBdescFMRF

#plot_background_MB()

#add INrope
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INrope, soma=T, WithConnectors = T, lwd=5, add=T, 
         alpha=(0.2+i/5), col=bluepurple[3])
  
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INrope", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2A_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INtorii
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INtorii, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INtorii", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2B_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INbigloop
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INbigloop, soma=T, WithConnectors = T, lwd=2, add=T, 
         alpha=(0.2+i/5), col=bluepurple[9])
  
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INbigloop", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2C_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INUturnMB
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INUturnMB, soma=T, WithConnectors = T, lwd=4, add=T, 
         alpha=(0.2+i/5), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INUturnMB", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2D_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INhorn
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INhorn, soma=T, WithConnectors = T, lwd=2, add=T, 
         alpha=(0.2+i/5), col=bluepurple[8])
  plot3d(INUturnMB, soma=T, lwd=4, add=T, 
         alpha=(1), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INhorn", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2E_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBPDF
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INMBPDF, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col=bluepurple[6])
  plot3d(INhorn, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[8])
  plot3d(INUturnMB, soma=T, lwd=4, add=T, 
         alpha=(1), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBPDF", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2F_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBdesc2
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  
  plot3d(INMBdesc2, soma=T, WithConnectors = T, lwd=1, add=T, 
         alpha=(0.2+i/5), col=bluepurple[7])
  plot3d(INMBPDF, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[6])
  plot3d(INhorn, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[8])
  plot3d(INUturnMB, soma=T, lwd=4, add=T, 
         alpha=(1), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBdesc2", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2G_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBdesc3
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INMBdesc3, soma=T, WithConnectors = T, lwd=5, add=T, 
         alpha=(0.2+i/5), col="black")
  
  plot3d(INMBdesc2, soma=T, lwd=1, add=T, 
         alpha=(1), col=bluepurple[7])
  plot3d(INMBPDF, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[6])
  plot3d(INhorn, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[8])
  plot3d(INUturnMB, soma=T, lwd=4, add=T, 
         alpha=(1), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBdesc3", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2H_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add INMBdescFMRF
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(INMBdescFMRF, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col="blue")
  plot3d(INMBdesc3, soma=T, lwd=5, add=T, 
         alpha=(1), col="black")
  
  plot3d(INMBdesc2, soma=T, lwd=1, add=T, 
         alpha=(1), col=bluepurple[7])
  plot3d(INMBPDF, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[6])
  plot3d(INhorn, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[8])
  plot3d(INUturnMB, soma=T, lwd=4, add=T, 
         alpha=(1), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "INMBdescFMRF", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2I_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

#add MBprojMouth
for(i in c(0:4)){print (i)
  clear3d()
  plot3d(yolk, soma=F, lwd=0, add=T, alpha=0.05, col="#E2E2E2") 
  plot3d(MBprojMouth, soma=T, WithConnectors = T, lwd=3, add=T, 
         alpha=(0.2+i/5), col=Okabe_Ito[7])
  
  plot3d(INMBdesc3, soma=T, lwd=5, add=T, 
         alpha=(1), col="black")
  
  plot3d(INMBdesc2, soma=T, lwd=1, add=T, 
         alpha=(1), col=bluepurple[7])
  plot3d(INMBPDF, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[6])
  plot3d(INhorn, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[8])
  plot3d(INUturnMB, soma=T, lwd=4, add=T, 
         alpha=(1), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.001, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  texts3d(99000,26000, 5000, text = "MBprojMouth", col='black', cex = 2)
  filename <- paste("videoframes/Video_MB_spin2J_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  rgl.pop()
}

plot_INproj <- function(){
  plot3d(MBprojMouth, soma=T, lwd=3, add=T, 
         alpha=(1), col=Okabe_Ito[7])
  plot3d(INMBdesc3, soma=T, lwd=5, add=T, 
         alpha=(1), col="black")
  
  plot3d(INMBdesc2, soma=T, lwd=1, add=T, 
         alpha=(1), col=bluepurple[7])
  plot3d(INMBPDF, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[6])
  plot3d(INhorn, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[8])
  plot3d(INUturnMB, soma=T, lwd=4, add=T, 
         alpha=(1), col=bluepurple[4])
  plot3d(INbigloop, soma=T, lwd=2, add=T, 
         alpha=(1), col=bluepurple[9])
  plot3d(INtorii, soma=T, lwd=3, add=T, 
         alpha=(1), col=bluepurple[5])
  plot3d(INrope, soma=T, lwd=5, add=T, 
         alpha=(1), col=bluepurple[3])
}

plot_INint()
plot_SNMB()
plot_INproj()

# rotation towards ventral ---------------


for(i in c(0:21)){print (i)
  nview3d('frontal',extramat=(rotationMatrix(i/(21/(pi/2)), -1, 0, 0)))
  filename <- paste("videoframes/Video_MB_spin3_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}
extramat_to_multiply <- rotationMatrix(21/(21/(pi/2)), -1, 0, 0)

# rotation to right --------------
for(i in c(0:21)){print (i)
  nview3d('frontal',extramat=(extramat_to_multiply%*%rotationMatrix(i/(21/(pi/2)), 0, 0, -1)))
  filename <- paste("videoframes/Video_MB_spin4_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}

extramat_to_multiply2=(extramat_to_multiply%*%rotationMatrix(i/(21/(pi/2)), 0, 0, -1))



# rotation to anterior ---------------
for(i in c(0:28)){print (i)
  nview3d("frontal",extramat=(extramat_to_multiply2%*%rotationMatrix(i/(21/(pi/2)), 1, 1, 1)))
  filename <- paste("videoframes/Video_MB_spin5_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}

#nview3d("frontal")

# zoom in on AO and show MBSN synapses ----------------
for (i in 0:21){
  clear3d()
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=F, lwd=0,
         add=T, alpha=0.05,
         col="#E2E2E2") 
  
  plot_INint()
  plot_SNMB()
  plot_INproj()
  
  plot3d(AOSN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
         add=T, alpha=(i/21), col=Okabe_Ito[2])
  plot3d(outline, soma=F, lwd=0, add=T, alpha=0.02, col='grey80')
  #z clipping plane same as in cycle 44
  clipplanes3d(0, 0, -1, (240004-4091*44))
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 120000)
  par3d(zoom=0.52-i/70)
  filename <- paste("videoframes/Video_MB_spin6_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}



# zoom out from AO and show MBmouth with stomodeum -----------
for (i in 0:30){
  #rotation towards ventral
  nview3d(extramat=(rotationMatrix(i/25, -1, 0, 0)%*%rotationMatrix(0.2, 1, 0.1, 0.5)))
  par3d(zoom=0.22+i/100)
  # plot3d(MBintrIN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
  #       add=T, alpha=(0.6-i/36), col=Okabe_Ito[5])
  #plot3d(MBprojMouth, WithConnectors = T, WithNodes = F, soma=T, lwd=2,
  #       add=T, alpha=(0.4+i/20), col=Okabe_Ito[7])
  #plot3d(MBSN, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
  #       add=T, alpha=(0.6-i/36), col=Okabe_Ito[1])
  plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=F, lwd=3,
         add=T, alpha=0.2, col="grey60")
  filename <- paste("videoframes/Video_MB_spin7_", formatC(i, digits = 1, flag = "0"), ".png", sep = "")
  rgl.snapshot(filename)
  print(i)
}



#read png files and write video
av::av_encode_video(paste('videoframes/', list.files("videoframes/", '*.png'), sep = ""), 
                    framerate = 10,
                    output = 'Videos/Video3.mp4')

