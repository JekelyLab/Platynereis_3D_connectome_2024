#code to generate Figure14 figure supplement 1 of the Platynereis 3d connectome paper
#Gaspar Jekely Feb-Dec 2022


#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

#read neuron groups
{
  doCRunp = nlapply(read.neurons.catmaid("^celltype102$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  h_antCR = nlapply(read.neurons.catmaid("^celltype100$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  cirrusCR = nlapply(read.neurons.catmaid("^celltype101$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  pygPBunp = nlapply(read.neurons.catmaid("^celltype54$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  PB = nlapply(read.neurons.catmaid("^celltype80$", pid=11),
               function(x) smooth_neuron(x, sigma=6000))
  antPU = nlapply(read.neurons.catmaid("^celltype88$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  hCirrus_interparaPU = nlapply(read.neurons.catmaid("^celltype94$", pid=11),
                                 function(x) smooth_neuron(x, sigma=6000))
  
  pygCirrusPU = nlapply(read.neurons.catmaid("^celltype98$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  hPU = nlapply(read.neurons.catmaid("^celltype96$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  neuro_notoPU = nlapply(read.neurons.catmaid("^celltype95$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  
  
  dsoPU = nlapply(read.neurons.catmaid("^celltype97$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  
  paraPU = nlapply(read.neurons.catmaid("^celltype91$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  
  spinPU = nlapply(read.neurons.catmaid("^celltype89$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  hPU2l_asymPDF = nlapply(read.neurons.catmaid("^celltype99$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  pygCirrusPUSM = nlapply(read.neurons.catmaid("^celltype93$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  
  SNblunt = nlapply(read.neurons.catmaid("^celltype148$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SNpygM = nlapply(read.neurons.catmaid("^celltype170$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SNantler = nlapply(read.neurons.catmaid("^celltype26$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  SNPDF_pyg = nlapply(read.neurons.catmaid("^celltype115$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  
  SNFVa = nlapply(read.neurons.catmaid("^celltype143$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  SNbronto = nlapply(read.neurons.catmaid("^celltype168$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  interparaPM = nlapply(read.neurons.catmaid("^celltype81$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  chaeMech = nlapply(read.neurons.catmaid("^celltype71$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  palp = nlapply(read.neurons.catmaid("^palp sensory neuron$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  antenna = nlapply(read.neurons.catmaid("^antennae_cell$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  
  INsplitPUh = nlapply(read.neurons.catmaid("^celltype83$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  INMC3 = nlapply(read.neurons.catmaid("^celltype75$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  INrope = nlapply(read.neurons.catmaid("^celltype58$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  
  
  INsplitCR = nlapply(read.neurons.catmaid("^celltype73$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  INsplitCRATO = nlapply(read.neurons.catmaid("^celltype74$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  INsplitBronto = nlapply(read.neurons.catmaid("^celltype149$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  INsplitVent = nlapply(read.neurons.catmaid("^celltype157$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  INasc_pyg = nlapply(read.neurons.catmaid("^celltype147$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  INsplitPBant = nlapply(read.neurons.catmaid("^celltype78$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  INsplitPB = nlapply(read.neurons.catmaid("^celltype79$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INsplitPB_RFYa = nlapply(read.neurons.catmaid("^celltype77$", pid=11),
                           function(x) smooth_neuron(x, sigma=6000))
  
  INCM = nlapply(read.neurons.catmaid("^celltype60$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  
  MNant = nlapply(read.neurons.catmaid("^celltype19$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  Loop = nlapply(read.neurons.catmaid("^celltype59$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  MC3cover = nlapply(read.neurons.catmaid("^celltype87$", pid=11),
                     function(x) smooth_neuron(x, sigma=6000))
  
  
}

#load other cell clusters
{
  
  stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))  

  MUSring_pyg = nlapply(read.neurons.catmaid("^celltype_non_neuronal78$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  gut = nlapply(read.neurons.catmaid("^gut$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  girdle = nlapply(read.neurons.catmaid("^mechanosensory_girdle$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  
}

#adjusted background plotting function
plot_background_mech <- function(){
  plot_background_ventral()
  clear3d()
  nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
  plot3d(yolk, WithConnectors = F, WithNodes = F, soma=TRUE, lwd=0,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.05,
         col='grey80')
  plot3d(stomodeum, WithConnectors = F, WithNodes = F, soma=F, lwd=3,
         rev = FALSE, fixup = F, add=T, alpha=0.1, col="grey50")
  plot3d(girdle, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
         rev = FALSE, fixup = F, add=T, alpha=0.15, col="grey50")
}

#plot Mechanosensory IN types
mech_IN <- list(
  INsplitPUh, INMC3, INrope, INsplitCR, 
  INsplitCRATO, INsplitBronto, INsplitVent, INasc_pyg, 
  INsplitPBant, INsplitPB, INsplitPB_RFYa, INCM
)

mech_IN_names <- c(
  'INsplitPUh', 'INMC3', 'INrope', 'INsplitCR', 
  'INsplitCRATO', 'INsplitBronto', 'INsplitVent', 'INasc_pyg', 
  'INsplitPBant', 'INsplitPB', 'INsplitPB_RFYa', 'INCM'
)
  
for (i in 1:length(mech_IN)){
    plot_background_mech()
    plot3d(
      mech_IN[i][[1]], WithConnectors = F, WithNodes = F, soma=T, lwd=5,
      add=T, alpha=1, 
      col=sample(bluepurple[5:9], length(mech_IN[i][[1]]), replace = TRUE)
      )
    par3d(zoom=0.55)
    filename <- paste("pictures/", mech_IN_names[i], ".png", sep= "")
    rgl.snapshot(filename)
    close3d()
}

#plot all IN
plot_background_mech()
for (i in 1:length(mech_IN)){
  plot3d(
    mech_IN[i][[1]], WithConnectors = F, WithNodes = F, soma=T, lwd=5,
    add=T, alpha=1, 
    col=sample(bluepurple[5:9], length(mech_IN[i][[1]]), replace = TRUE)
  )
}
par3d(zoom=0.55)
rgl.snapshot("pictures/mech_inter_all.png")
close3d()


#plot Mechanosensory SN types
mech_sensory <- list(
  h_antCR, doCRunp, cirrusCR, pygPBunp, 
  PB, antPU, hCirrus_interparaPU, hPU,
  pygCirrusPU, neuro_notoPU, dsoPU, paraPU, 
  spinPU, hPU2l_asymPDF, pygCirrusPUSM, 
  SNblunt, SNpygM, SNantler, SNPDF_pyg, 
  SNFVa, SNbronto, interparaPM, chaeMech
)

mech_sensory_names <- c(
  'h_antCR', 'doCRunp', 'cirrusCR', 'pygPBunp', 
  'PB', 'antPU', 'hCirrus_interparaPU', 'hPU', 
  'pygCirrusPU', 'neuro_notoPU', 'dsoPU', 'paraPU', 
  'spinPU', 'hPU2l_asymPDF', 'pygCirrusPUSM', 'SNblunt', 
  'SNpygM', 'SNantler', 'SNPDF_pyg', 'SNFVa', 
  'SNbronto', 'interparaPM', 'chaeMech'
)

for (i in 1:length(mech_sensory)){
  plot_background_mech()
  plot3d(mech_sensory[i][[1]], WithConnectors = F, WithNodes = F, soma=T, lwd=4,
         rev = FALSE, fixup = F, add=T, alpha=1, 
         col=sample(oranges[4:9], length(mech_sensory[i][[1]]), replace=TRUE))
  par3d(zoom=0.55)
  if (i==1){  plot3d(scalebar_50um_ventral, WithConnectors = F, WithNodes = F, soma=F, lwd=5,
                      rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
                      col="black")
    }
  filename <- paste("pictures/", mech_sensory_names[i], ".png", sep= "")
  rgl.snapshot(filename)
  close3d()
}

#plot all SN
plot_background_mech()
for (i in 1:length(mech_sensory)){
  plot3d(mech_sensory[i][[1]], WithConnectors = F, WithNodes = F, soma=T, lwd=4,
         rev = FALSE, fixup = F, add=T, alpha=1, 
         col=sample(oranges[4:9], length(mech_sensory[i][[1]]), replace=TRUE))
}

par3d(zoom=0.55)
rgl.snapshot("pictures/mech_sensory_all.png")
close3d()


#plot Mechanosensory MN types
mech_MN <- list(Loop, MC3cover, MNant)
mech_MN_names <- c('Loop', 'MC3cover', 'MNant')
  
for (i in 1:length(mech_MN)){
    plot_background_mech()
    plot3d(mech_MN[i][[1]], WithConnectors = F, WithNodes = F, soma=T, lwd=4,
           rev = FALSE, fixup = F, add=T, alpha=1, 
           col=sample(blues[5:9], length(mech_MN[i][[1]]), replace=TRUE))
    par3d(zoom=0.55)
    filename <- paste("pictures/", mech_MN_names[i], ".png", sep= "")
    rgl.snapshot(filename)
    close3d()
}



# assemble figure ---------------------------------------------------------




{
  imgAll <- readPNG("pictures/mech_sensory_all.png")
  imgAll_IN <- readPNG("pictures/mech_inter_all.png")
  imgA <- readPNG("pictures/doCRunp.png")
  imgB <- readPNG("pictures/h_antCR.png")
  imgC <- readPNG("pictures/cirrusCR.png")
  imgD <- readPNG("pictures/pygPBunp.png")
  imgE <- readPNG("pictures/PB.png")
  imgF <- readPNG("pictures/antPU.png")
  imgG <- readPNG("pictures/hCirrus_interparaPU.png")
  imgH <- readPNG("pictures/hPU.png")
  imgI <- readPNG("pictures/pygCirrusPU.png")
  imgJ <- readPNG("pictures/neuro_notoPU.png")

  imgA2 <- readPNG("pictures/dsoPU.png")
  imgB2 <- readPNG("pictures/paraPU.png")
  imgC2 <- readPNG("pictures/spinPU.png")
  imgD2 <- readPNG("pictures/hPU2l_asymPDF.png")
  imgE2 <- readPNG("pictures/pygCirrusPUSM.png")
  imgF2 <- readPNG("pictures/SNblunt.png")
  imgG2 <- readPNG("pictures/SNpygM.png")
  imgH2 <- readPNG("pictures/SNantler.png")
  imgI2 <- readPNG("pictures/SNPDF_pyg.png")
  imgJ2 <- readPNG("pictures/SNFVa.png")
  
  imgA3 <- readPNG("pictures/SNbronto.png")
  imgB3 <- readPNG("pictures/interparaPM.png")
  imgC3 <- readPNG("pictures/chaeMech.png")

  
#read IN pictures    
  imgK <- readPNG("pictures/INsplitPUh.png")
  imgL <- readPNG("pictures/INMC3.png")
  imgM <- readPNG("pictures/INrope.png")
  imgN <- readPNG("pictures/INsplitCR.png")
  imgO <- readPNG("pictures/INsplitCRATO.png")
  imgP <- readPNG("pictures/INsplitBronto.png")
  imgQ <- readPNG("pictures/INsplitVent.png")
  imgQ2 <- readPNG("pictures/INasc_pyg.png")
  imgQ3 <- readPNG("pictures/INsplitPBant.png")
  imgR <- readPNG("pictures/INsplitPB.png")
  imgS <- readPNG("pictures/INsplitPB_RFYa.png")
  imgT <- readPNG("pictures/INCM.png")
  

  imgU <- readPNG("pictures/Loop.png")
  imgV <- readPNG("pictures/MC3cover.png")
  imgW <- readPNG("pictures/MNant.png")

  #convert png to image panel
  panelAll_SN <- ggdraw() + draw_image(imgAll, scale = 1) + 
    draw_label("all sensory", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelAll_IN <- ggdraw() + draw_image(imgAll_IN, scale = 1) + 
    draw_label("all inter", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelA <- ggdraw() + draw_image(imgA, scale = 1) + 
    draw_label("doCRunp", x = 0.5, y = 0.98,
               color = "black", size = 11)

  
  panelB <- ggdraw() + draw_image(imgB, scale = 1) + 
    draw_label("antCR/hCR", x = 0.5, y = 0.98,
               color = "black", size = 11) +
    draw_label(expression(paste("50 ", mu, "m")), x = 0.76, y = 0.072,
               color = "black", size = 10)
  
  panelC <- ggdraw() + draw_image(imgC, scale = 1) + 
    draw_label("cirrusCR", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelD <- ggdraw() + draw_image(imgD, scale = 1) + 
    draw_label("pygPBunp", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelE <- ggdraw() + draw_image(imgE, scale = 1) + 
    draw_label("PB", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelF <- ggdraw() + draw_image(imgF, scale = 1) + 
    draw_label("antPU", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelG <- ggdraw() + draw_image(imgG, scale = 1) + 
    draw_label("hCirrusPU/interparaPU", x = 0.55, y = 0.98,
               color = "black", size = 11)
  panelH <- ggdraw() + draw_image(imgH, scale = 1) + 
    draw_label("hPU", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelI <- ggdraw() + draw_image(imgI, scale = 1) + 
    draw_label("pygCirrusPU", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelJ <- ggdraw() + draw_image(imgJ, scale = 1) + 
    draw_label("neuro_notoPU", x = 0.5, y = 0.98,
               color = "black", size = 11)
  
  panelA2 <- ggdraw() + draw_image(imgA2, scale = 1) + 
    draw_label("dsoPU", x = 0.5, y = 0.98,
               color = "black", size = 11) 
  panelB2 <- ggdraw() + draw_image(imgB2, scale = 1) + 
    draw_label("paraPU", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelC2 <- ggdraw() + draw_image(imgC2, scale = 1) + 
    draw_label("spinPU", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelD2 <- ggdraw() + draw_image(imgD2, scale = 1) + 
    draw_label("hPU2l-asymPDF", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelE2 <- ggdraw() + draw_image(imgE2, scale = 1) + 
    draw_label("pygCirrusPUSM", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelF2 <- ggdraw() + draw_image(imgF2, scale = 1) + 
    draw_label("SNblunt", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelG2 <- ggdraw() + draw_image(imgG2, scale = 1) + 
    draw_label("SNpygM", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelH2 <- ggdraw() + draw_image(imgH2, scale = 1) + 
    draw_label("SNantler", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelI2 <- ggdraw() + draw_image(imgI2, scale = 1) + 
    draw_label("SNPDF-pyg", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelJ2 <- ggdraw() + draw_image(imgJ2, scale = 1) + 
    draw_label("SNFV", x = 0.5, y = 0.98,
               color = "black", size = 11)
  
  panelA3 <- ggdraw() + draw_image(imgA3, scale = 1) + 
    draw_label("SNbronto", x = 0.5, y = 0.98,
               color = "black", size = 11) 
  panelB3 <- ggdraw() + draw_image(imgB3, scale = 1) + 
    draw_label("interparaPM", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelC3 <- ggdraw() + draw_image(imgC3, scale = 1) + 
    draw_label("chaeMech", x = 0.5, y = 0.98,
               color = "black", size = 11)
  
  

  panelK <- ggdraw() + draw_image(imgK, scale = 1) + 
    draw_label("INsplitPUh", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelL <- ggdraw() + draw_image(imgL, scale = 1) + 
    draw_label("INMC3", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelM <- ggdraw() + draw_image(imgM, scale = 1) + 
    draw_label("INrope", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelN <- ggdraw() + draw_image(imgN, scale = 1) + 
    draw_label("INsplitCR", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelO <- ggdraw() + draw_image(imgO, scale = 1) + 
    draw_label("INsplitCRATO", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelP <- ggdraw() + draw_image(imgP, scale = 1) + 
    draw_label("INsplitBronto", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelQ <- ggdraw() + draw_image(imgQ, scale = 1) + 
    draw_label("INsplitVent", x = 0.5, y = 0.98,
               color = "black", size = 11)
  
  panelQ2 <- ggdraw() + draw_image(imgQ2, scale = 1) + 
    draw_label("INasc_pyg", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelQ3 <- ggdraw() + draw_image(imgQ3, scale = 1) + 
    draw_label("INsplitPBant", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelR <- ggdraw() + draw_image(imgR, scale = 1) + 
    draw_label("INsplitPB", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelS <- ggdraw() + draw_image(imgS, scale = 1) + 
    draw_label("INsplitPB-RF/Ya", x = 0.6, y = 0.98,
               color = "black", size = 11)
  panelT <- ggdraw() + draw_image(imgT, scale = 1) + 
    draw_label("INCM", x = 0.5, y = 0.98,
               color = "black", size = 11)
  
  panelU <- ggdraw() + draw_image(imgU, scale = 1) + 
    draw_label("Loop", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelV <- ggdraw() + draw_image(imgV, scale = 1) + 
    draw_label("MC3cover", x = 0.5, y = 0.98,
               color = "black", size = 11)
  panelW <- ggdraw() + draw_image(imgW, scale = 1) + 
    draw_label("MNant", x = 0.5, y = 0.98,
               color = "black", size = 11)

}


{
#define layout with textual representation
layout <- "
aAbBcCdD
########
eEfFgGhH
########
iIjJkKlL
########
mMnNoOpP
########
qQrRsStT
"
  
Fig_mech_circuits <- panelAll_SN + panelAll_IN + panelA + panelB + panelC + panelD + panelE + panelF + 
  panelG + panelH + panelI + panelJ + panelA2 + panelB2+ panelC2 + panelD2 + 
  panelE2 + panelF2 + panelG2 + panelH2 + panelI2 + panelJ2 + panelA3 + panelB3 + 
  panelC3 + panelK + panelL + panelM + panelN + panelO + panelP + panelQ + 
  panelQ2 + panelQ3 + panelR + panelS + panelT + panelU + panelV + panelW + 
  plot_layout(design = layout, heights = c(1,0.05,1,0.05,1,0.05,1,0.05,1)) +
  plot_annotation(tag_levels = 'i') & 
  theme(plot.tag = element_text(size = 12, face='plain'))
  
ggsave("Figures/Figure14_fig_suppl1.pdf", limitsize = FALSE, 
         units = c("px"), Fig_mech_circuits, width = 4800, height = 4400)
  
ggsave("Figures/Figure14_fig_suppl1.png", limitsize = FALSE, 
         units = c("px"), Fig_mech_circuits, width = 4800, height = 4400, bg='white')
  
}


