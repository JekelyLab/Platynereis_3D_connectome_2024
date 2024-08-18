#R/natverse code to generate Figure MB anatomy overview for the Platynereis 3d connectome paper
#Gaspar Jekely 2022

# load natverse and other packages, custom natverse functions and  --------
source("code/Natverse_functions_and_conn.R")


# read MB neuron types ----------------------------------------------------

{
  INMBtype1 = nlapply(read.neurons.catmaid("^INMBtype1$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype2 = nlapply(read.neurons.catmaid("^celltype198$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype3 = nlapply(read.neurons.catmaid("^INMBtype3$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype4 = nlapply(read.neurons.catmaid("^INMBtype4$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype5 = nlapply(read.neurons.catmaid("^celltype121$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype6 = nlapply(read.neurons.catmaid("^celltype120$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype7 = nlapply(read.neurons.catmaid("^celltype122$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype8 = nlapply(read.neurons.catmaid("^INMBtype8$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype9 = nlapply(read.neurons.catmaid("^celltype197$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBtype10 = nlapply(read.neurons.catmaid("^INMBtype10$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  INMBPDF = nlapply(read.neurons.catmaid("^celltype184$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SNhook = nlapply(read.neurons.catmaid("^celltype18$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  #plot3d(SNhorn, soma=T, lwd=3,
  #       add=T, alpha=0.8, col=oranges[2])
  SNhorn = nlapply(read.neurons.catmaid("^celltype17$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SNlasso = nlapply(read.neurons.catmaid("^celltype20$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SNtorii = nlapply(read.neurons.catmaid("^celltype187$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SN_NS1 = nlapply(read.neurons.catmaid("^celltype110$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  SN_NS17 = nlapply(read.neurons.catmaid("^celltype131$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SN_NS18 = nlapply(read.neurons.catmaid("^celltype138$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SN_NS19 = nlapply(read.neurons.catmaid("^celltype133$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  SN_NS27 = nlapply(read.neurons.catmaid("^celltype137$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  INhorn = nlapply(read.neurons.catmaid("^celltype183$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INtorii = nlapply(read.neurons.catmaid("^celltype191$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INMBdescFMRF = nlapply(read.neurons.catmaid("^celltype185$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  INMBPDF = nlapply(read.neurons.catmaid("^celltype184$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  INrope = nlapply(read.neurons.catmaid("^celltype58$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INMBdesc2 = nlapply(read.neurons.catmaid("^celltype195$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INMBdesc3 = nlapply(read.neurons.catmaid("^celltype194$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  INbigloop = nlapply(read.neurons.catmaid("^celltype140$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INUturnMB = nlapply(read.neurons.catmaid("^celltype192$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  
  MBprojMouth = nlapply(read.neurons.catmaid("^MBprojMouth$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
  
  stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))  
  MB = nlapply(read.neurons.catmaid("^mushroom body$", pid=11),
               function(x) smooth_neuron(x, sigma=6000)) 
  SNMBdev = nlapply(read.neurons.catmaid("^SNMBdev$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
  INMBdev = nlapply(read.neurons.catmaid("^INMBdev$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000)) 
  MBintrIN   = read.neurons.catmaid("^MBintrIN$", pid=11)

  
  INMBtype7 = nlapply(read.neurons.catmaid("^celltype122$", pid=11),
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
  
  
  INMBtype9 = nlapply(read.neurons.catmaid("^celltype197$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INproT2 = nlapply(read.neurons.catmaid("^celltype125$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  INMBtype5 = nlapply(read.neurons.catmaid("^celltype121$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INW = nlapply(read.neurons.catmaid("^celltype117$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  Ser_h1 = nlapply(read.neurons.catmaid("^celltype8$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  
  
  INMBtype2 = nlapply(read.neurons.catmaid("^celltype198$", pid=11),
                      function(x) smooth_neuron(x, sigma=6000))
  INdecusshook = nlapply(read.neurons.catmaid("^celltype153$", pid=11),
                         function(x) smooth_neuron(x, sigma=6000))
  
  
  MS1 = nlapply(read.neurons.catmaid("^celltype35$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
  INsn = nlapply(read.neurons.catmaid("^celltype57$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  INpreMN = nlapply(read.neurons.catmaid("^celltype23$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  MNring_INproT2_target = nlapply(read.neurons.catmaid("^MNring_INproT2_target$", pid=11),
                                  function(x) smooth_neuron(x, sigma=6000))
  MNring_INproT2_target_MUS_targets = nlapply(read.neurons.catmaid("^MNring_INproT2_target_MUS_targets$", pid=11),
                                              function(x) smooth_neuron(x, sigma=6000))
  
  INhook = nlapply(read.neurons.catmaid("^celltype119$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INbiax = nlapply(read.neurons.catmaid("^celltype76$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  
  
  INDLSO = nlapply(read.neurons.catmaid("^celltype186$", pid=11),
                   function(x) smooth_neuron(x, sigma=6000))
  INhook = nlapply(read.neurons.catmaid("^celltype119$", pid=11),
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
  
  MNant = nlapply(read.neurons.catmaid("^celltype19$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  cPRC = nlapply(read.neurons.catmaid("^celltype5$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  palp = nlapply(read.neurons.catmaid("^palp sensory neuron$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  antenna = nlapply(read.neurons.catmaid("^antennae_cell$", pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  
  MBON = nlapply(read.neurons.catmaid("^MBON$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  SNMBinputs = nlapply(read.neurons.catmaid("^MBSNinputs$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  SNMB = nlapply(read.neurons.catmaid("^MBSN$", pid=11),
                 function(x) smooth_neuron(x, sigma=6000))
  MBcentralINoutput = nlapply(read.neurons.catmaid("^MBcentralINoutput$", pid=11),
                              function(x) smooth_neuron(x, sigma=6000))
  MBprojINoutput = nlapply(read.neurons.catmaid("^MBprojINoutput$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))
  MNant = nlapply(read.neurons.catmaid("^celltype19$", pid=11),
                  function(x) smooth_neuron(x, sigma=6000))
  
}

# load other cell clusters ------------------------------------------------

{
eye = nlapply(read.neurons.catmaid("^celltype1$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
AOSN = nlapply(read.neurons.catmaid("^apical sensory organ$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))
stomodeum = nlapply(read.neurons.catmaid("^stomodeum$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))  
EC = nlapply(read.neurons.catmaid("^celltype_non_neuronal90$", pid=11),
                function(x) smooth_neuron(x, sigma=6000))

#read in two blocks because curl often breaks for larger sets
connectome_left = nlapply(read.neurons.catmaid(skids_by_2annotations("connectome", 
                "left_side"), pid=11, conn = conn_http1,
                fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
connectome_right = nlapply(read.neurons.catmaid(skids_by_2annotations("connectome", 
                "right_side"), pid=11, conn = conn_http1,
                fetch.annotations = FALSE),
                function(x) smooth_neuron(x, sigma=6000))
}


# plot MBintrIN -----------------------------------------------------------

{
plot_background()
plot3d(INMBtype1, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[5])
plot3d(INMBtype2, soma=T, lwd=2,
       add=T, alpha=1, col=blues[9])
plot3d(INMBtype3, soma=T, lwd=3,
       add=T, alpha=1, col=blues[5])
plot3d(INMBtype4, soma=T, lwd=2,
       add=T, alpha=1, col=blues[8])
plot3d(INMBtype5, soma=T, lwd=3,
       add=T, alpha=1, col=blues[4])
plot3d(INMBtype6, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[2])
plot3d(INMBtype7, soma=T, lwd=3,
       add=T, alpha=1, col=blues[7])
plot3d(INMBdev, soma=T, lwd=2,
       add=T, alpha=1, col=blues[3])

plot3d(eye, soma=T, lwd=1,
       add=T, alpha=0.2, col=Okabe_Ito[8])
texts3d(36000,42000, 24000, text = "eye", col='black', cex = 2.5)
texts3d(96000,22000, 25000, text = "yolk", col='grey30', cex = 2.5)

plot3d(
  scalebar_50um_anterior, lwd = 4,
  add = T, alpha = 1, col = "black"
)
texts3d(104000, 116000, 300, text = "50 μm", col = "black", cex = 2)


filename <- paste("pictures/Figure_MB_overviewA.png")
rgl.snapshot(filename)
close3d()
}

# plot MBSN ---------------------------------------------------------------

{
plot_background()
plot3d(SNtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNlasso, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SNhook, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[5])
plot3d(SNhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS19, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[7])
plot3d(SN_NS1, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[8])
plot3d(SN_NS17, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SN_NS18, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS27, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNMBdev, soma=T, lwd=3,
       add=T, alpha=0.4, col=oranges[4])
plot3d(eye, soma=T, lwd=1,
       add=T, alpha=0.2, col=Okabe_Ito[8])
texts3d(36000,42000, 24000, text = "eye", col='black', cex = 2.5)
texts3d(96000,22000, 25000, text = "yolk", col='grey30', cex = 2.5)

filename <- paste("pictures/Figure_MB_overviewB.png")
rgl.snapshot(filename)

clear3d()
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
plot3d(yolk, soma=TRUE, lwd=0,
       add=T, alpha=0.05,
       col='grey70')

#z-axis clip
clipplanes3d(0, 0, -1, 65000)
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)))


plot3d(SNtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNlasso, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SNhook, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[5])
plot3d(SNhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS19, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[7])
plot3d(SN_NS1, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[8])
plot3d(SN_NS17, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SN_NS18, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS27, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNMBdev, soma=T, lwd=3,
       add=T, alpha=0.4, col=oranges[4])

plot3d(stomodeum, soma=F, lwd=2,
       add=T, alpha=0.2, col="grey50")
plot3d(AOSN, soma=T, lwd=1,
       add=T, alpha=0.3, col=Okabe_Ito[2])
texts3d(78000,132000, 24000, text = "stomodeum", col='black', cex = 2.5)
texts3d(72000,62000, 2000, text = "apical sensory organ", col='black', cex = 2.5)

par3d(zoom=0.62)
filename <- paste("pictures/Figure_MB_overviewBventr.png")
rgl.snapshot(filename)
close3d()
}


# plot MBON ---------------------------------------------------------------

{
plot_background()
plot3d(INMBdescFMRF, soma=T, lwd=3,
         add=T, alpha=0.8, col=blues[8])
plot3d(INMBdesc2, soma=T, lwd=2,
       add=T, alpha=0.8, col=blues[7])
plot3d(INMBdesc3, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[7])
plot3d(INbigloop, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[5])
plot3d(INMBPDF, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[9])
plot3d(INtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=blues[9])

plot3d(INhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[7])
plot3d(INUturnMB, soma=T, lwd=3,
       add=T, alpha=0.8, col=Okabe_Ito[5])
plot3d(INrope, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[2])
plot3d(MBprojMouth, soma=T, lwd=2,
       add=T, alpha=0.8, col=bluepurple[9])
plot3d(eye, soma=T, lwd=1,
       add=T, alpha=0.2, col=Okabe_Ito[8])
texts3d(36000,42000, 24000, text = "eye", col='black', cex = 2.5)
texts3d(96000,22000, 25000, text = "yolk", col='grey30', cex = 2.5)

filename <- paste("pictures/Figure_MB_overviewCant.png")
rgl.snapshot(filename)
close3d()
}  


# plot all MB neurons -----------------------------------------------------

{
plot_background()
plot3d(INMBtype1, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[5])
plot3d(INMBtype2, soma=T, lwd=2,
       add=T, alpha=1, col=blues[9])
plot3d(INMBtype3, soma=T, lwd=3,
       add=T, alpha=1, col=blues[5])
plot3d(INMBtype4, soma=T, lwd=2,
       add=T, alpha=1, col=blues[8])
plot3d(INMBtype5, soma=T, lwd=3,
       add=T, alpha=1, col=blues[4])
plot3d(INMBtype6, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[2])
plot3d(INMBtype7, soma=T, lwd=3,
       add=T, alpha=1, col=blues[7])
plot3d(INMBdev, soma=T, lwd=2,
       add=T, alpha=1, col=blues[3])

plot3d(eye, soma=T, lwd=1,
       add=T, alpha=0.2, col=Okabe_Ito[8])
texts3d(36000,42000, 24000, text = "eye", col='black', cex = 2.5)


plot3d(SNtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNlasso, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SNhook, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[5])
plot3d(SNhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS19, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[7])
plot3d(SN_NS1, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[8])
plot3d(SN_NS17, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SN_NS18, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS27, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNMBdev, soma=T, lwd=3,
       add=T, alpha=0.4, col=oranges[4])



plot3d(INMBdescFMRF, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[8])
plot3d(INMBdesc2, soma=T, lwd=2,
       add=T, alpha=0.8, col=blues[7])
plot3d(INMBdesc3, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[7])
plot3d(INbigloop, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[5])
plot3d(INMBPDF, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[9])
plot3d(INtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=blues[9])

plot3d(INhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[7])
plot3d(INUturnMB, soma=T, lwd=3,
       add=T, alpha=0.8, col=Okabe_Ito[5])
plot3d(INrope, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[2])
plot3d(MBprojMouth, soma=T, lwd=2,
       add=T, alpha=0.8, col=bluepurple[9])
texts3d(96000,22000, 25000, text = "yolk", col='grey30', cex = 2.5)

filename <- paste("pictures/Figure_MB_overviewDant.png")
rgl.snapshot(filename)

clear3d()
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
plot3d(yolk, soma=TRUE, lwd=0,
       add=T, alpha=0.05,
       col='grey70')
#z-axis clip
clipplanes3d(0, 0, -1, 95000)
par3d(zoom=0.62)

plot3d(INMBtype1, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[5])
plot3d(INMBtype2, soma=T, lwd=2,
       add=T, alpha=1, col=blues[9])
plot3d(INMBtype3, soma=T, lwd=3,
       add=T, alpha=1, col=blues[5])
plot3d(INMBtype4, soma=T, lwd=2,
       add=T, alpha=1, col=blues[8])
plot3d(INMBtype5, soma=T, lwd=3,
       add=T, alpha=1, col=blues[4])
plot3d(INMBtype6, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[2])
plot3d(INMBtype7, soma=T, lwd=3,
       add=T, alpha=1, col=blues[7])
plot3d(INMBdev, soma=T, lwd=2,
       add=T, alpha=1, col=blues[3])



plot3d(SNtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNlasso, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SNhook, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[5])
plot3d(SNhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS19, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[7])
plot3d(SN_NS1, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[8])
plot3d(SN_NS17, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SN_NS18, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS27, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNMBdev, soma=T, lwd=3,
       add=T, alpha=0.4, col=oranges[4])


plot3d(INMBdescFMRF, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[8])
plot3d(INMBdesc2, soma=T, lwd=2,
       add=T, alpha=0.8, col=blues[7])
plot3d(INMBdesc3, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[7])
plot3d(INbigloop, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[5])
plot3d(INMBPDF, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[9])
plot3d(INtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=blues[9])

plot3d(INhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[7])
plot3d(INUturnMB, soma=T, lwd=3,
       add=T, alpha=0.8, col=Okabe_Ito[5])
plot3d(INrope, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[2])
plot3d(MBprojMouth, soma=T, lwd=2,
       add=T, alpha=0.8, col=bluepurple[9])
plot3d(stomodeum, soma=F, lwd=2,
       add=T, alpha=0.2, col="grey50")
texts3d(78000,132000, 24000, text = "stomodeum", col='black', cex = 2.5)

filename <- paste("pictures/Figure_MB_overviewDventr.png")
rgl.snapshot(filename)

close3d()
plot_background()
clear3d()
nview3d('ventral', extramat=(rotationMatrix(1.43, 0, 0, -1)%*%rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
plot3d(yolk, soma=TRUE, lwd=0,
       add=T, alpha=0.05,
       col='grey70')

plot3d(INMBtype1, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[5])
plot3d(INMBtype2, soma=T, lwd=2,
       add=T, alpha=1, col=blues[9])
plot3d(INMBtype3, soma=T, lwd=3,
       add=T, alpha=1, col=blues[5])
plot3d(INMBtype4, soma=T, lwd=2,
       add=T, alpha=1, col=blues[8])
plot3d(INMBtype5, soma=T, lwd=3,
       add=T, alpha=1, col=blues[4])
plot3d(INMBtype6, soma=T, lwd=2,
       add=T, alpha=1, col=Okabe_Ito[2])
plot3d(INMBtype7, soma=T, lwd=3,
       add=T, alpha=1, col=blues[7])
plot3d(INMBdev, soma=T, lwd=2,
       add=T, alpha=1, col=blues[3])



plot3d(SNtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNlasso, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SNhook, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[5])
plot3d(SNhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS19, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[7])
plot3d(SN_NS1, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[8])
plot3d(SN_NS17, soma=T, lwd=2,
       add=T, alpha=0.8, col=oranges[4])
plot3d(SN_NS18, soma=T, lwd=4,
       add=T, alpha=0.8, col=oranges[2])
plot3d(SN_NS27, soma=T, lwd=3,
       add=T, alpha=0.8, col=oranges[3])
plot3d(SNMBdev, soma=T, lwd=3,
       add=T, alpha=0.4, col=oranges[4])


plot3d(INMBdescFMRF, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[8])
plot3d(INMBdesc2, soma=T, lwd=2,
       add=T, alpha=0.8, col=blues[7])
plot3d(INMBdesc3, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[7])
plot3d(INbigloop, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[5])
plot3d(INMBPDF, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[9])
plot3d(INtorii, soma=T, lwd=4,
       add=T, alpha=0.8, col=blues[9])

plot3d(INhorn, soma=T, lwd=3,
       add=T, alpha=0.8, col=blues[7])
plot3d(INUturnMB, soma=T, lwd=3,
       add=T, alpha=0.8, col=Okabe_Ito[5])
plot3d(INrope, soma=T, lwd=4,
       add=T, alpha=0.8, col=Okabe_Ito[2])
plot3d(MBprojMouth, soma=T, lwd=2,
       add=T, alpha=0.8, col=bluepurple[9])
plot3d(stomodeum, soma=F, lwd=3,
       add=T, alpha=0.3, col="grey80")
texts3d(58000,125000, 29000, text = "stomodeum", col='black', cex = 2.2)



#z-axis clip
clipplanes3d(0, 0, -1, 95000)
par3d(zoom=0.62)
#y-axis clip
clipplanes3d(0, 1, 0, -35000)


filename <- paste("pictures/Figure_MB_overviewDlat.png")
rgl.snapshot(filename)

close3d()
}

# plot all MB neurons full body, ventral ----------------------------------

{
plot_background()
clear3d()
nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
#z-axis clip
par3d(zoom=0.72)
  
plot3d(INMBtype1, soma=T, lwd=2,
         add=T, alpha=1, col=Okabe_Ito[5])
  plot3d(INMBtype2, soma=T, lwd=2,
         add=T, alpha=1, col=blues[9])
  plot3d(INMBtype3, soma=T, lwd=3,
         add=T, alpha=1, col=blues[5])
  plot3d(INMBtype4, soma=T, lwd=2,
         add=T, alpha=1, col=blues[8])
  plot3d(INMBtype5, soma=T, lwd=3,
         add=T, alpha=1, col=blues[4])
  plot3d(INMBtype6, soma=T, lwd=2,
         add=T, alpha=1, col=Okabe_Ito[2])
  plot3d(INMBtype7, soma=T, lwd=3,
         add=T, alpha=1, col=blues[7])
  plot3d(INMBdev, soma=T, lwd=2,
         add=T, alpha=1, col=blues[3])
  

  plot3d(SNtorii, soma=T, lwd=4,
         add=T, alpha=0.8, col=oranges[3])
  plot3d(SNlasso, soma=T, lwd=3,
         add=T, alpha=0.8, col=oranges[4])
  plot3d(SNhook, soma=T, lwd=3,
         add=T, alpha=0.8, col=oranges[5])
  plot3d(SNhorn, soma=T, lwd=3,
         add=T, alpha=0.8, col=oranges[2])
  plot3d(SN_NS19, soma=T, lwd=3,
         add=T, alpha=0.8, col=oranges[7])
  plot3d(SN_NS1, soma=T, lwd=2,
         add=T, alpha=0.8, col=oranges[8])
  plot3d(SN_NS17, soma=T, lwd=2,
         add=T, alpha=0.8, col=oranges[4])
  plot3d(SN_NS18, soma=T, lwd=4,
         add=T, alpha=0.8, col=oranges[2])
  plot3d(SN_NS27, soma=T, lwd=3,
         add=T, alpha=0.8, col=oranges[3])
  plot3d(SNMBdev, soma=T, lwd=3,
         add=T, alpha=0.4, col=oranges[4])
  
  plot3d(INMBdescFMRF, soma=T, lwd=3,
         add=T, alpha=0.8, col=blues[8])
  plot3d(INMBdesc2, soma=T, lwd=2,
         add=T, alpha=0.8, col=blues[7])
  plot3d(INMBdesc3, soma=T, lwd=4,
         add=T, alpha=0.8, col=Okabe_Ito[7])
  plot3d(INbigloop, soma=T, lwd=3,
         add=T, alpha=0.8, col=blues[5])
  plot3d(INMBPDF, soma=T, lwd=3,
         add=T, alpha=0.8, col=blues[9])
  plot3d(INtorii, soma=T, lwd=4,
         add=T, alpha=0.8, col=blues[9])
  
  plot3d(INhorn, soma=T, lwd=3,
         add=T, alpha=0.8, col=blues[7])
  plot3d(INUturnMB, soma=T, lwd=3,
         add=T, alpha=0.8, col=Okabe_Ito[5])
  plot3d(INrope, soma=T, lwd=4,
         add=T, alpha=0.8, col=Okabe_Ito[2])
  plot3d(MBprojMouth, soma=T, lwd=2,
         add=T, alpha=0.8, col=bluepurple[9])
  plot3d(stomodeum, soma=F, lwd=2,
         add=T, alpha=0.2, col="grey50")
  texts3d(78000,132000, 24000, text = "stomodeum", col='black', cex = 2.2)


  plot3d(yolk, soma=TRUE, lwd=0,
         add=T, alpha=0.05,
         col='grey70')
  plot3d(acicula, soma=T, lwd=3,
         add=T, alpha=1,
         col="grey50")
  texts3d(34000,162000, 84000, text = "aciculae", col='black', cex = 2.5)
  
  plot3d(chaeta, soma=F, lwd=1,
         add=T, alpha=1,
         col="grey70") 
  texts3d(144000,152000, 54000, text = "chaetae", col='grey40', cex = 2.5)
  texts3d(10000,72000, 12000, text = "MB", col='black', cex = 2.5)
  plot3d(connectome_left, soma=F, lwd=1,
         add=T, alpha=0.3,
         col="#56B4E9") 
  plot3d(connectome_right, soma=F, lwd=1,
         add=T, alpha=0.3,
         col="#56B4E9")
  plot3d(EC, soma=T, lwd=1,
         add=T, alpha=0.1,
         col="grey") 
  filename <- paste("pictures/Figure_MB_overview_full.png")
  rgl.snapshot(filename)
  close3d()
}

# MB circuits ----------------------------------------------

# plot overview images for summary graph

#plot and save MBintrIN, MBON
{
  plot_background()
  plot3d(MBintrIN, WithConnectors = F, soma=T, lwd=2,
         add=T, alpha=0.6, col="#CC79A7")
  plot3d(MBON, soma=TRUE, lwd=3,
         add=T, alpha=0.8,
         col="#56B4E9")
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  rgl.snapshot("pictures/Figure_MBMBintrIN_MBON.png")
  close3d()
}

#plot and save SNMBinputs, SNMB
{
  plot_background()
  plot3d(SNMBinputs, soma=TRUE, lwd=2,
         add=T, alpha=0.6,
         col="#D55E00")
  plot3d(SNMB, soma=TRUE, lwd=3,
         add=T, alpha=0.8,
         col="#E69F00")
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  rgl.snapshot("pictures/Figure_SNMBinputs_SNMB.png")
  close3d()
}

#plot and save MBcentralINoutput, MBONoutput
{
  plot_background()
  plot3d(MBcentralINoutput, soma=TRUE, lwd=2,
         add=T, alpha=0.6,
         col="#0072B2")
  plot3d(MBprojINoutput, soma=TRUE, lwd=3,
         add=T, alpha=1,
         col="#0072B2")
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  rgl.snapshot("pictures/Figure_MBcentralINoutput_MBONoutput.png")
  close3d()
}

#plot and save MNant
{
  plot_background()
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  plot3d(MNant, soma=TRUE, lwd=2,
         add=T, alpha=1,
         col=c("grey20","grey40"))
  rgl.snapshot("pictures/Figure_MB_MNant.png")
  close3d()
}


# get connectivity between MB cell clusters and their input and output --------

{
cell_groups <- list(MBintrIN, MBON, SNMBinputs, SNMB, MBcentralINoutput, 
                    MBprojINoutput, MNant)
N_cell_groups <- length(cell_groups)
N_cell_groups

cell_group_attr <- data.frame (
  cell_group_names  = c("MBintrIN", "MBON", "SN inputs", "SNMB", 
                        "centralIN", "projIN", "MNant"),
  type = c("MBintrIN", "MBON", "SN", "MBSN", 
           "IN", "IN", "MN"),
  level = c("2", "2", "1", "1", 
            "3", "3", "4")
)
dim(cell_group_attr)

}

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
{
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=N_cell_groups)
rownames(synapse_matrix) <- cell_group_attr$cell_group_names
colnames(synapse_matrix) <- cell_group_attr$cell_group_names
synapse_matrix
}

write.csv(as.data.frame(synapse_matrix), "source_data/Figure9_source_data1.txt")

# graph conversion --------------------------------------------------------

{
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


#use visNetwork to plot the network

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
    visHierarchicalLayout(levelSeparation=550, 
                          direction='LR',
                          sortMethod='hubsize',
                          shakeTowards='roots') %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0.3),
             scaling=list(min=2, max=12),
             color = list(inherit=TRUE, opacity=0.5),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 1.2, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(background=Conn_graph.visn$nodes$color, border='black'),
             opacity=0.9,
             shape='dot', 
             font=list(color='black', size=52),
             scaling = list(label=list(enabled=TRUE, min=44, max=52)),
             level= Conn_graph.visn$nodes$level) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, algorithm='hierarchical',labelOnly=FALSE)) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    visGroups(groupname = "MBSN", shape = "square", 
              opacity=1, color="#E69F00") %>%
    visGroups(groupname = "MBintrIN", shape = "square", 
              opacity=1, color="#CC79A7") %>%
    visGroups(groupname = "MBON", shape = "square", 
              opacity=1, color="#56B4E9") %>%
    visGroups(groupname = "IN", shape = "dot", 
              opacity=1, color="#0072B2") %>%
    visGroups(groupname = "SN", shape = "dot", 
              opacity=1, color="#D55E00") %>%
    visGroups(groupname = "MN", shape = "dot", 
              opacity=1, color="#cccccc")  %>%
    addFontAwesome()
  
  
  visNet
}

#save as html
saveNetwork(visNet, "pictures/visNetwork_MBgrouped_circuit.html")
webshot2::webshot(url="pictures/visNetwork_MBgrouped_circuit.html",
                  file="pictures/visNetwork_MBgrouped_circuit.png",
                  vwidth=1000, vheight=300, #define the size of the browser window
                  cliprect = c(40,110,950, 280), zoom=10)

}

# get connectivity between MB cell types and their partners and ma --------

{
#these are the MBintr cell types and their partners
cell_groups <-  list(INMBtype7, INpara, INZ, INMBtype6, 
                     INlasso, INW, INfoot, 
                     INMBtype5, 
                     INMBtype2,
                     SNlasso, SNtorii, 
                     SNhook, SNhorn,INMBPDF, INrope,  
                     INMBdesc3, INMBdesc2, INtorii, INhorn,
                     INhook, INbigloop, INUturnMB,  
                     SNstiff, SNPDF_pyg, SNbronto,  
                     palp, antenna, 
                     INleucoPU, INsplitPUh, INdecussM, INdecusshook,  
                     INDLSO,  MNladder, MNant, Ser_h1  
)
N_cell_groups <- length(cell_groups)
N_cell_groups

#add attributes to nodes
{
  cell_group_attr <- data.frame (
    cell_group_names  = c("INMBtype7", "INpara", "INZ", "INMBtype6", 
                          "INlasso", "INW", "INfoot",
                          "INMBtype5", 
                          "INMBtype2",
                          " SNlasso", "SNtorii", 
                          "SNhook","SNhorn","INMBPDF","INrope",  
                          "INMBdesc3","INMBdesc2","INtorii","INhorn",
                          "INhook","INbigloop","INUturnMB", 
                          "SNstiff","SNPDFpyg","SNbronto", 
                          "palp","antenna",
                          "INleucoPU","INsplitPUh","INdecussM","INdecusshook",  
                          "INDLSO","MNladder","MNant","Ser-h1"),
    type = c("MBintrIN", "IN", "IN", "MBintrIN", 
             "IN", "IN", "IN",
             "MBintrIN",
             "MBintrIN",
             "MBSN","MBSN", 
             "MBSN","MBSN","MBON","MBON",  
             "MBON","MBON","MBON","MBON",
             "IN","MBON","MBON", 
             "SN","SN","SN", 
             "SN","SN", 
             "IN","IN","IN","IN",  
             "IN","MN","MN", "MN"),
    level = c("3", "5", "5", "3", 
              "5", "6", "5", 
              "3", 
              "3",
              "2","2", 
              "2","2","4","4",  
              "4","4","4","4",
              "5","4","4", 
              "1","9","1", 
              "1","1", 
              "7","7","6","6",  
              "5","8","8","8")
  )
  dim(cell_group_attr)
}


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
max(synapse_matrix)

write.csv(as.data.frame(synapse_matrix), "source_data/Figure9_source_data2.txt")

#plot with ggplot
{
  as.data.frame((synapse_matrix)) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
    ggplot() +
    geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction), 
               stroke = 0)  + 
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
  ggsave("pictures/MBproj_syn_matrix.pdf", 
         width = nrow(synapse_matrix)/1.3, 
         height = ncol(synapse_matrix)/1.6, limitsize = TRUE, 
         units = c("cm"))
  
  # Saving R ggplot with R ggsave Function
  ggsave("pictures/MBproj_syn_matrix.png", 
         width = 1700, 
         height = 1300, limitsize = TRUE, 
         units = c("px"))
  
  write.csv2(synapse_matrix, file = "supplements/Figure_MB_circuits_synapse_matrix.csv")
}


# graph conversion --------------------------------------------------------

#edge weight filtering on the matrix to remove weak edges
synapse_matrix[synapse_matrix < 3] <- 0

#with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
Conn_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE,
  diag = TRUE
)

#calculate node weighted degree
degree=degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)

# use visNetwork to plot the network --------------------------------------

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
    visIgraphLayout(layout = "layout_nicely", physics = TRUE, 
                    randomSeed = 42, type="square") %>%
    visHierarchicalLayout(levelSeparation=250, 
                          nodeSpacing=10,
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
             font=list(color='black', size=44),
             scaling = list(label=list(enabled=TRUE, min=48, max=56)),
             level= Conn_graph.visn$nodes$level) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, algorithm='hierarchical',labelOnly=FALSE)) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    visGroups(groupname = "MBSN", shape = "square", 
              opacity=1, color="#E69F00") %>%
    visGroups(groupname = "MBintrIN", shape = "square", 
              opacity=1, color="#CC79A7") %>%
    visGroups(groupname = "MBON", shape = "square", 
              opacity=1, color="#56B4E9") %>%
    visGroups(groupname = "IN", shape = "dot", 
              opacity=1, color="#0072B2") %>%
    visGroups(groupname = "SN", shape = "dot", 
              opacity=1, color="#D55E00") %>%
    visGroups(groupname = "MN", shape = "dot", 
              opacity=1, color="#cccccc")  %>%
    addFontAwesome() %>%
    visLegend(
      #addNodes = list(
      # list(label = "66 syn", shape = "icon", 
      #     icon = list(code = "f2d1", size = 30, color = "#D55E00"))), 
      useGroups = TRUE,  width=0.1,ncol = 1,
      position='right', stepY=70)
  
  
  visNet
}

#save as html
saveNetwork(visNet, "pictures/visNetwork_MBproj_circuit.html", selfcontained = TRUE)
webshot2::webshot(url="pictures/visNetwork_MBproj_circuit.html",
                  file="pictures/visNetwork_MBproj_circuit.png",
                  vwidth=1200, vheight=600, #define the size of the browser window
                  cliprect = c(40,70,950, 440), zoom=5, delay = 2)

}

# MB inputs outputs part of figure ----------------------------------------

{
  plot_background()
  plot3d(INpara, soma=TRUE, lwd=3,
         add=T, alpha=1,
         col=bluepurple[5])
  plot3d(INZ, soma=TRUE, lwd=4,
         add=T, alpha=0.8,
         col=blues[8])
  plot3d(INlasso, soma=TRUE, lwd=2,
         add=T, alpha=0.7,
         col=Okabe_Ito[6])
  plot3d(INW, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[7])
  plot3d(INfoot, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(INproT2, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=blues[4])
  plot3d(INdecusshook, soma=TRUE, lwd=3,
         add=T, alpha=0.8,
         col=Okabe_Ito[6])
  plot3d(MB, soma=T, lwd=0,
         add=T, alpha=0.15, col="grey50")  
  plot3d(Ser_h1, soma=TRUE, lwd=3,
         add=T, alpha=0.6,
         col=bluepurple[6])
  texts3d(54000,104000, 22000, text = "INpara", col='black', cex = 2.5)
  texts3d(52000,46000, 27000, text = "INZ", col='black', cex = 2.5) 
  texts3d(69000,36000, 22000, text = "INlasso", col='black', cex = 2.5)
  texts3d(73000,76000, 8000, text = "INW", col='black', cex = 2.5)
  texts3d(56000,63000, 7000, text = "INfoot", col='black', cex = 2.5)
  texts3d(66000,90000, 12000, text = "INproT2", col='black', cex = 2.5)
  texts3d(38000,39000, 10000, text = "Ser-h1", col='black', cex = 2.5)
  texts3d(66000,53000, 6000, text = "INdecusshook", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_outputs.png")
  close3d()
}

#plot and save INMBtype7
{
  plot_background()
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  plot3d(INMBtype7, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(INpara, soma=TRUE, lwd=3,
         add=T, alpha=1,
         col=Okabe_Ito[6])
  plot3d(INZ, soma=TRUE, lwd=4,
         add=T, alpha=0.8,
         col=Okabe_Ito[7])
  texts3d(34000,64000, 12000, text = "INMBtype7", col='black', cex = 2.5)
  texts3d(54000,104000, 27000, text = "INpara", col='black', cex = 2.5)
  texts3d(52000,46000, 27000, text = "INZ", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_type7.png")
  close3d()
}

#plot and save INMBtype6
{
  plot_background()
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  
  plot3d(INMBtype6, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(INlasso, soma=TRUE, lwd=2,
         add=T, alpha=0.7,
         col=Okabe_Ito[6])
  plot3d(INW, soma=TRUE, lwd=3,
         add=T, alpha=0.8,
         col=Okabe_Ito[7])
  plot3d(INfoot, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[5])
  texts3d(30000,67000, 12000, text = "INMBtype6", col='black', cex = 2.5)
  texts3d(50000,44000, 22000, text = "INlasso", col='black', cex = 2.5)
  texts3d(73700,76000, 8000, text = "INW", col='black', cex = 2.5)
  texts3d(56000,64000, 7000, text = "INfoot", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_type6.png")
  close3d()
}

#plot and save INMBtype9
{
  plot_background()
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  plot3d(INMBtype9, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(INproT2, soma=TRUE, lwd=3,
         add=T, alpha=0.8,
         col=Okabe_Ito[6])
  texts3d(38000,69000, 12000, text = "INMBtype9", col='black', cex = 2.5)
  texts3d(56000,92000, 11000, text = "INproT2", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_type9.png")
  close3d()
}

#plot and save INMBtype5
{
  plot_background()
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  plot3d(INMBtype5, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(INW, soma=TRUE, lwd=4,
         add=T, alpha=0.7,
         col=Okabe_Ito[6])
  plot3d(Ser_h1, soma=TRUE, lwd=3,
         add=T, alpha=0.7,
         col=Okabe_Ito[7])
  texts3d(32000,70000, 5000, text = "INMBtype5", col='black', cex = 2.5)
  texts3d(66000,85000, 6000, text = "INW", col='black', cex = 2.5)
  texts3d(40000,39000, 10000, text = "Ser-h1", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_type5.png")
  close3d()
}

#plot and save INMBtype2
{
  plot_background()
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  plot3d(INMBtype2, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(INdecusshook, soma=TRUE, lwd=3,
         add=T, alpha=0.8,
         col=Okabe_Ito[6])
  texts3d(28000,54000, 12000, text = "INMBtype2", col='black', cex = 2.5)
  texts3d(66000,53000, 6000, text = "INdecusshook", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_type2.png")
  close3d()
}

#plot INdecusshook output systems and save
{
  plot_background()
  plot3d(INdecusshook, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[6])
  plot3d(INMBtype2, soma=TRUE, lwd=2,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(SNhook, soma=TRUE, lwd=3,
         add=T, alpha=1,
         col=Okabe_Ito[1])
  plot3d(INfoot, soma=TRUE, lwd=2,
         add=T, alpha=0.6,
         col=Okabe_Ito[3])
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  plot3d(Ser_h1, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[7])
  
  texts3d(32000,54000, 12000, text = "INMBtype2", col='black', cex = 2.5)
  texts3d(73000,62000, 6000, text = "INdecusshook", col='black', cex = 2.5)
  texts3d(34000,73000, 8000, text = "SNhook", col='black', cex = 2.5)
  texts3d(57000,53000, 8000, text = "INfoot", col='black', cex = 2.5)
  texts3d(40000,40000, 10000, text = "Ser-h1", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_INdecusshook.png")
  close3d()
}

#plot INproT2 output systems and save
{
  plot_background()
  clear3d()
  plot3d(yolk, soma=TRUE, lwd=0,
         add=T, alpha=0.05,
         col='grey70')
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  plot3d(INproT2, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[6])
  plot3d(INMBtype9, soma=TRUE, lwd=2,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(MS1, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[1])
  plot3d(INsn, soma=TRUE, lwd=2,
         add=T, alpha=0.6,
         col=Okabe_Ito[3])
  plot3d(INpreMN, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[7])
  plot3d(MNring_INproT2_target, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[6])
  plot3d(MNring_INproT2_target_MUS_targets, soma=TRUE, lwd=2,
         add=T, alpha=0.3,
         col=Okabe_Ito[8])
  
  texts3d(40000,68000, 10000, text = "INMBtype9", col='black', cex = 2.5)
  texts3d(78000,98000, 12000, text = "INproT2", col='black', cex = 2.5)
  texts3d(74000,57000, 4000, text = "MS1", col='black', cex = 2.5)
  texts3d(91000,71000, 7000, text = "INsn", col='black', cex = 2.5)
  texts3d(64000,42000, 28000, text = "INpreMN", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_INproT2_ant.png")
  
  nview3d("ventral")
  par3d(zoom=0.55)
  #z-axis clip
  clipplanes3d(0, 0, -1, 95000)
  texts3d(42000,85000, 72000, text = "MNring MUS t.", col='black', cex = 2.5)
  texts3d(64000,84000, 79000, text = "MNring", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_INproT2_ventr.png")
  close3d()
}

#plot INW output systems and save
{
  plot_background()
  plot3d(INW, soma=TRUE, lwd=4,
         add=T, alpha=1,
         col=Okabe_Ito[6])
  plot3d(INMBtype5, soma=TRUE, lwd=2,
         add=T, alpha=1,
         col=Okabe_Ito[2])
  plot3d(INMBtype6, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[3])
  plot3d(SNlasso, soma=TRUE, lwd=2,
         add=T, alpha=0.6,
         col=Okabe_Ito[1])
  plot3d(INhook, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[7])
  plot3d(INbiax, soma=TRUE, lwd=2,
         add=T, alpha=0.8,
         col=Okabe_Ito[5])
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.15, col="grey50")
  texts3d(40000,82000, 10000, text = "INMBtype5", col='black', cex = 2.5)
  texts3d(40000,68000, 10000, text = "INMBtype6", col='black', cex = 2.5)
  texts3d(75000,74000, 12000, text = "INW", col='black', cex = 2.5)
  texts3d(34000,54000, 4000, text = "SNlasso", col='black', cex = 2.5)
  texts3d(58000,52000, 7000, text = "INhook", col='black', cex = 2.5)
  texts3d(54000,37000, 28000, text = "INbiax", col='black', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_circuits_INW.png")
  close3d()
}


#plot all SN inputs to MB
{
  plot_background()
  clear3d()
  plot3d(yolk, soma=TRUE, lwd=0,
         add=T, alpha=0.05,
         col='grey70')
  #z-axis clip
  clipplanes3d(0, 0, -1, 125000)
  plot3d(SNstiff, soma=T, lwd=3,
         add=T, alpha=1, col=oranges[6])
  plot3d(SNbronto, soma=T, lwd=3,
         add=T, alpha=1, col=oranges[5])
  plot3d(antenna, soma=T, lwd=3,
         add=T, alpha=0.2, col=oranges[8])
  plot3d(palp, soma=T, lwd=3,
         add=T, alpha=0.2, col=oranges[5])
  par3d(zoom=0.6)
  texts3d(100000,66000, 100, text = "antenna", col='black', cex = 2.5)
  texts3d(42000,102000, 16000, text = "palp", col='black', cex = 2.5)
  texts3d(74000,128000, 21000, text = "SNstiff", col='black', cex = 2.5)
  texts3d(74000,82000, 8000, text = "SNbronto", col='black', cex = 2.5)
  texts3d(58000,42000, 88000, text = "yolk", col='grey30', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_all_SN_inputs_ant.png")
  
  clear3d()
  nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
  plot3d(yolk, soma=TRUE, lwd=0,
         add=T, alpha=0.05,
         col='grey70')
  plot3d(SNPDF_pyg, soma=T, lwd=3,
         add=T, alpha=1, col=oranges[4])
  plot3d(SNstiff, soma=T, lwd=3,
         add=T, alpha=1, col=oranges[6])
  plot3d(SNbronto, soma=T, lwd=3,
         add=T, alpha=1, col=oranges[5])
  plot3d(antenna, soma=T, lwd=3,
         add=T, alpha=0.2, col=oranges[8])
  plot3d(palp, soma=T, lwd=3,
         add=T, alpha=0.2, col=oranges[5])
  par3d(zoom=0.76)
  texts3d(129000,76000, 7000, text = "antenna", col='black', cex = 2.5)
  texts3d(42000,122000, 16000, text = "palp", col='black', cex = 2.5)
  texts3d(74000,122000, 41000, text = "SNstiff", col='black', cex = 2.5)
  texts3d(68000,142000, 156000, text = "SNPDF-pyg", col='black', cex = 2.5)
  texts3d(74000,120000, 8000, text = "SNbronto", col='black', cex = 2.5)
  texts3d(58000,142000, 88000, text = "yolk", col='grey30', cex = 2.5)
  rgl.snapshot("pictures/Figure_MB_all_SN_inputs_ventr.png")
  
  close3d()
}

#plot palp and partners
{
  plot_background()
  plot3d(palp, soma=T, lwd=3,
         add=T, alpha=0.5, col=oranges[5])
  plot3d(
    scalebar_50um_anterior, lwd = 4,
    add = T, alpha = 1, col = "black"
  )
  texts3d(104000, 116000, 300, text = "50 μm", col = "black", cex = 2)
  
  plot3d(INtorii, soma=T, lwd=3,
         add=T, alpha=0.8, col=blues[9])
  plot3d(INhorn, soma=T, lwd=5,
         add=T, alpha=0.8, col=blues[5])
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.1, col="grey50")
  texts3d(42000,102000, 16000, text = "palp", col='black', cex = 2.5)
  texts3d(42000,59000, 21000, text = "INhorn", col='black', cex = 2.5)
  texts3d(32000,71000, 8000, text = "INtorii", col='black', cex = 2.5)
  
  rgl.snapshot("pictures/Figure_MB_palp_and_partners.png")
  close3d()
}

#plot antenna, SNbronto and partners
{
  plot_background()
  plot3d(antenna, soma=T, lwd=3,
         add=T, alpha=0.2, col=oranges[8])
  plot3d(INtorii, soma=T, lwd=3,
         add=T, alpha=0.8, col=blues[9])
  plot3d(INrope, soma=T, lwd=5,
         add=T, alpha=0.8, col=blues[4])
  plot3d(SNbronto, soma=T, lwd=3,
         add=T, alpha=1, col=oranges[5])
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.1, col="grey50")
  
  texts3d(74000,82000, 8000, text = "SNbronto", col='black', cex = 2.5)
  texts3d(104000,72000, 100, text = "antenna", col='black', cex = 2.5)
  texts3d(35000,60000, 6000, text = "INrope", col='black', cex = 2.5)
  texts3d(32000,80000, 8000, text = "INtorii", col='black', cex = 2.5)
  
  rgl.snapshot("pictures/Figure_MB_antenna_SNbronto_partners.png")
  close3d()
}

#plot SNstiff SNhorn and partners
{
  plot_background()
  plot3d(SNstiff, soma=T, lwd=3,
         add=T, alpha=0.8, col=oranges[8])
  plot3d(SNhorn, soma=T, lwd=2,
         add=T, alpha=1, col=oranges[9])
  par3d(zoom=0.54)
  plot3d(INhorn, soma=T, lwd=3,
         add=T, alpha=0.6, col=blues[9])
  plot3d(INMBPDF, soma=T, lwd=5,
         add=T, alpha=0.7, col=blues[4])
  plot3d(MB, soma=T, lwd=1,
         add=T, alpha=0.1, col="grey50")
  texts3d(77000,120000, 21000, text = "SNstiff", col='black', cex = 2.5)
  texts3d(20000,57000, 8000, text = "SNhorn", col='black', cex = 2.5)
  texts3d(44000,57500, 6000, text = "INhorn", col='black', cex = 2.5)
  texts3d(32000,74000, 8000, text = "INMBPDF", col='black', cex = 2.5)
  
  rgl.snapshot("pictures/Figure_MB_SNstiff_SNhorn_and_partners.png")
  close3d()
}

# assemble figure -------------------------------------------------------------

#read images
{
img_MBover <- readPNG("pictures/Figure_MB_overview_full.png")

imgC <- readPNG("pictures/Figure_MB_overviewA.png")
imgD <- readPNG("pictures/Figure_MB_overviewB.png")
imgE <- readPNG("pictures/Figure_MB_overviewBventr.png")
imgF <- readPNG("pictures/Figure_MB_overviewCant.png")
imgG <- readPNG("pictures/Figure_MB_overviewCventr.png")

imgH <- readPNG("pictures/Figure_MB_overviewDant.png")
imgI <- readPNG("pictures/Figure_MB_overviewDventr.png")
imgJ <- readPNG("pictures/Figure_MB_overviewDlat.png")
}

#convert png to image panel
{
panel_MBover <- ggdraw() + draw_image(img_MBover, scale = 1)

panelMBintrIN <- ggdraw() + draw_image(imgC, scale = 1) + 
  draw_label("MBintrIN, anterior", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1) 

panelMBSN <- ggdraw() + draw_image(imgD, scale = 1) + 
  draw_label("MBSN, anterior", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
panelMBON <- ggdraw() + draw_image(imgF, scale = 1) + 
  draw_label("MBON, anterior", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)

panelAllMB <- ggdraw() + draw_image(imgH, scale = 1) + 
  draw_label("all MB neurons, anterior", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
panelAllMBventr <- ggdraw() + draw_image(imgI, scale = 1) + 
  draw_label("all MB neurons, ventral", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
panelAllMBright <- ggdraw() + draw_image(imgJ, scale = 1) + 
  draw_label("all MB neurons, right", x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)

}

#second half of figure
#read images
{
img1b <- readPNG('pictures/visNetwork_MBproj_circuit.png')
img2b <- readPNG('pictures/visNetwork_MBproj_circuit_sel1.png')
img3b <- readPNG('pictures/visNetwork_MBproj_circuit_sel2.png')
img4b <- readPNG('pictures/visNetwork_MBproj_circuit_sel3.png')
img5b <- readPNG('pictures/visNetwork_MBproj_circuit_sel4.png')

img_group <- readPNG('pictures/visNetwork_MBgrouped_circuit.png')

img_small1 <- readPNG("pictures/Figure_SNMBinputs_SNMB.png")
img_small2 <- readPNG("pictures/Figure_MBMBintrIN_MBON.png")
img_small3 <- readPNG("pictures/Figure_MBcentralINoutput_MBONoutput.png")
img_small4 <- readPNG("pictures/Figure_MB_MNant.png")
}


#convert png to image panel
{
panelMBcelltypes <- ggdraw() + draw_image(img1b, scale = 1) + 
  draw_label("MB cell types and partners connectivity", x = 0.22, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)  + 
  draw_label("# of synapses", x = 0.89, y = 0.9, size = 8, hjust = 1) +
  draw_label("6", x = 0.82, y = 0.85, size = 8, hjust = 1) + 
  draw_label("29", x = 0.82, y = 0.8, size = 8, hjust = 1) +
  draw_line(x = c(0.83, 0.88), y = c(0.85, 0.85), size = 0.3, color = 'grey') +
  draw_line(x = c(0.83, 0.88), y = c(0.8, 0.8), size = 1.2, color = 'grey')

panelSNinputs <- ggdraw() + draw_image(img2b, scale = 1) + 
  draw_label("partners of SN inputs", x = 0.15, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)

panelSNMBpartners<- ggdraw() + draw_image(img3b, scale = 1) + 
  draw_label("partners of SNMB", x = 0.15, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)

panelINMB_intr_partners <- ggdraw() + draw_image(img4b, scale = 1) + 
  draw_label("partners of INMBintr", x = 0.15, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)

panelINMBproj_partners <- ggdraw() + draw_image(img5b, scale = 1) + 
  draw_label("partners of INMBproj", x = 0.15, y = 0.98, fontfamily = "sans", fontface = "plain",
             color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)

graph_grouped  <- ggdraw() + draw_image(img_group, scale = 1) + 
  draw_label("MB neuron-category connectivity", x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
             size = 11,lineheight = 0.9)  + 
  draw_label("# of synapses", x = 0.99, y = 0.9, size = 8, hjust = 1) +
  draw_label("4", x = 0.92, y = 0.8, size = 8, hjust = 1) + 
  draw_label("69", x = 0.92, y = 0.72, size = 8, hjust = 1) +
  draw_line(x = c(0.93, 0.98), y = c(0.8, 0.8), size = 0.3, color = 'grey') +
  draw_line(x = c(0.93, 0.98), y = c(0.72, 0.72), size = 1.2, color = 'grey')

panel_small1  <- ggdraw() + draw_image(img_small1, scale = 1)
panel_small2  <- ggdraw() + draw_image(img_small2, scale = 1)
panel_small3  <- ggdraw() + draw_image(img_small3, scale = 1)
panel_small4  <- ggdraw() + draw_image(img_small4, scale = 1)
}

#merge panels
{
panel_grouped_plots <- plot_grid(panel_small1,panel_small2,panel_small3,panel_small4,
                                 ncol=4,
                                 align="left")
}

#bottom part of fig - MB inputs outputs
#read images
{
  imgSN1 <- readPNG("pictures/Figure_MB_all_SN_inputs_ant.png")
  imgSN2 <- readPNG("pictures/Figure_MB_all_SN_inputs_ventr.png")
  imgSN3 <- readPNG("pictures/Figure_MB_palp_and_partners.png")
  imgSN4 <- readPNG("pictures/Figure_MB_antenna_SNbronto_partners.png")
  imgSN5 <- readPNG("pictures/Figure_MB_SNstiff_SNhorn_and_partners.png")
  
  img_MBout <- readPNG("pictures/Figure_MB_outputs.png")
  
  img_type2 <- readPNG("pictures/Figure_MB_circuits_type2.png")
  img_type5 <- readPNG("pictures/Figure_MB_circuits_type5.png")
  img_type6 <- readPNG("pictures/Figure_MB_circuits_type6.png")
  img_type7 <- readPNG("pictures/Figure_MB_circuits_type7.png")
  img_type9 <- readPNG("pictures/Figure_MB_circuits_type9.png")
  
  
  img_T2a <- readPNG("pictures/Figure_MB_circuits_INproT2_ant.png")
  img_T2v <- readPNG("pictures/Figure_MB_circuits_INproT2_ventr.png")
  img_INW <- readPNG("pictures/Figure_MB_circuits_INW.png")
  img_INde <- readPNG("pictures/Figure_MB_circuits_INdecusshook.png")
}

#convert png to image panel
{
  panelSNA <- cowplot::ggdraw() + draw_image(imgSN1, scale = 1) + 
    draw_label("MB SN inputs", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  panelSNB <- ggdraw() + draw_image(imgSN2, scale = 1) + 
    draw_label("MB SN inputs, ventral view", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panelSNC <- ggdraw() + draw_image(imgSN3, scale = 1) + 
    draw_label("palp partners", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panelSND <- ggdraw() + draw_image(imgSN4, scale = 1) + 
    draw_label("antenna and SNbronto partners", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panelSNE <- ggdraw() + draw_image(imgSN5, scale = 1) + 
    draw_label("SNhorn, SNstiff partners", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  
  panel_type2 <- ggdraw() + draw_image(img_type2, scale = 1) + 
    draw_label("INMBtype2", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panel_type5 <- ggdraw() + draw_image(img_type5, scale = 1) + 
    draw_label("INMBtype5", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panel_type6 <- ggdraw() + draw_image(img_type6, scale = 1) + 
    draw_label("INMBtype6", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panel_type7 <- ggdraw() + draw_image(img_type7, scale = 1) + 
    draw_label("INMBtype7", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panel_type9 <- ggdraw() + draw_image(img_type9, scale = 1) + 
    draw_label("INMBtype9", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  panel_T2a <- ggdraw() + draw_image(img_T2a, scale = 1) + 
    draw_label("INproT2 partners, anterior", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panel_T2v <- ggdraw() + draw_image(img_T2v, scale = 1) + 
    draw_label("INproT2 partners, ventral", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  panel_INW <- ggdraw() + draw_image(img_INW, scale = 1) + 
    draw_label("INW partners", x = 0.4, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  panel_INde <- ggdraw() + draw_image(img_INde, scale = 1) + 
    draw_label("INdecussHook partners", x = 0.5, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
  panel_MBout <- ggdraw() + draw_image(img_MBout, scale = 1) + 
    draw_label("MB IN outputs", x = 0.3, y = 0.97, fontfamily = "sans", fontface = "plain",
               color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1)
  
}

#assemble fig
{

layout1 <- "
AAAA
BBBB"

panel_grouped_graph <-  graph_grouped + panel_grouped_plots +
  plot_layout(design=layout1, heights = c(1,1))

panel_grouped_graph <- plot_grid(panel_grouped_graph)

layout3 <- 'AB
CD'

Figure9_fig_suppl1 <- panelSNinputs + panelSNMBpartners+ panelINMB_intr_partners + panelINMBproj_partners +
  plot_layout(design=layout3, heights = c(1,1,1,1), widths = c(1,1,1,1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=12, face='plain'))

ggsave("Figures/Figure9_fig_suppl1.pdf", limitsize = FALSE, 
       units = c("px"), Figure9_fig_suppl1, width = 4000, height = 2400)  

ggsave("Figures/Figure9_fig_suppl1.png", limitsize = FALSE, 
       units = c("px"), Figure9_fig_suppl1, width = 4000, height = 2400, bg='white')  

#define layout with textual representation
layout <- "
AAAABBCCDD
AAAA######
AAAAEEFFGG
"

Fig_MB_overview <-  panel_MBover + panelMBSN +  panelMBintrIN + 
   panelMBON + panelAllMB + panelAllMBventr + panelAllMBright +
  plot_layout(design = layout, heights = c(1, 0.05, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))


ggsave("Figures/Figure8.pdf", limitsize = FALSE, 
       units = c("px"), Fig_MB_overview, width = 4000, height = 1800)

ggsave("Figures/Figure8.png", limitsize = FALSE, 
       units = c("px"), Fig_MB_overview, width = 4000, height = 1900, bg='white')


layout <- "
JJJJKKKKKK
##########
LLMMNNOOPP
##########
QQRRSSTTUU
##########
VVWWXXYYZZ
"

Fig_MB_overview2 <-  
  panel_grouped_graph + panelMBcelltypes  +
  panelSNA + panelSNB +panelSNC + panelSND +panelSNE +
  panel_MBout + panel_type2 + panel_type5 + panel_type6 + panel_type7 + panel_type9 +
  panel_T2a + panel_T2v + panel_INW + panel_INde +
  plot_layout(design = layout, heights = c(1.30, 0.05, 1, 0.05, 1, 0.05,1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 12, face='plain'))


ggsave("Figures/Figure9.pdf", limitsize = FALSE, 
       units = c("px"), Fig_MB_overview2, width = 4000, height = 3500)

ggsave("Figures/Figure9.png", limitsize = FALSE, 
       units = c("px"), Fig_MB_overview2, width = 4000, height = 3500, bg='white')

}

