# Code to generate Figure6 and 7 of the Platynereis 3d connectome paper
# Gaspar Jekely
# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# read sensory cell clusters

Dorsal_SO <- nlapply(
  read.neurons.catmaid("^Dorsal_sensory_organ$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
Adult_eye <- nlapply(
  read.neurons.catmaid("^celltype1$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
nuchal <- nlapply(
  read.neurons.catmaid("^nuchal organ$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
DLSO <- nlapply(
  read.neurons.catmaid("Dorsolateral sense organs", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
palp <- nlapply(
  read.neurons.catmaid("^palp sensory neuron$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
antenna <- nlapply(
  read.neurons.catmaid("^antennae_cell$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MBSN <- nlapply(
  read.neurons.catmaid("^MBSN$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
eyespot <- nlapply(
  read.neurons.catmaid("^eyespot_PRC$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
ANS_SN <- nlapply(
  read.neurons.catmaid("^apical sensory organ$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
HSG3 <- nlapply(
  read.neurons.catmaid("head sensory group 3", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
cMS_cell <- nlapply(
  read.neurons.catmaid("central_MS_cell", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
headPU <- nlapply(
  read.neurons.catmaid("headPU", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
antler <- nlapply(
  read.neurons.catmaid("celltype26", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
vMIP <- nlapply(
  read.neurons.catmaid("celltype24", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNbicil <- nlapply(
  read.neurons.catmaid("celltype50", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNcirri <- nlapply(
  read.neurons.catmaid("anterior pair of tentacular cirri", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)


# read head interneuron clusters

visualIN <- nlapply(
  read.neurons.catmaid("visual_interneuron", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
ANS_IN <- nlapply(
  read.neurons.catmaid("ns_plexusIN", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INdesc_decus <- nlapply(
  read.neurons.catmaid("INdesc-decuss-head", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INcentrHead <- nlapply(
  read.neurons.catmaid("central head IN cluster", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INMBproj <- nlapply(
  read.neurons.catmaid("MBprojIN", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INMBintr <- nlapply(
  read.neurons.catmaid("MBintrIN", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MBprojMouth <- nlapply(
  read.neurons.catmaid("MBprojMouth", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)


# read head MN clusters

vMN <- nlapply(
  read.neurons.catmaid("vMN", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
prototroch  <- nlapply(
  read.neurons.catmaid("celltype_non_neuronal3", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
ciliomotor <- nlapply(
  read.neurons.catmaid("ciliomotor_head", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MNgland <- nlapply(
  read.neurons.catmaid("celltype166", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MNheadV <- nlapply(
  read.neurons.catmaid("celltype178", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)



# read MB cells and minicircuit cell groups

SNhorn <- nlapply(
  read.neurons.catmaid("^celltype17$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNhook <- nlapply(
  read.neurons.catmaid("^celltype18$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNlasso <- nlapply(
  read.neurons.catmaid("^celltype20$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNtorii <- nlapply(
  read.neurons.catmaid("^celltype187$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INtorii <- nlapply(
  read.neurons.catmaid("^celltype191$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INlasso <- nlapply(
  read.neurons.catmaid("^celltype21$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INfoot <- nlapply(
  read.neurons.catmaid("^celltype116$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INdecussfoot <- nlapply(
  read.neurons.catmaid("^celltype152$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INMBPDF <- nlapply(
  read.neurons.catmaid("^celltype184$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INhorn <- nlapply(
  read.neurons.catmaid("^celltype183$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INdecusshook <- nlapply(
  read.neurons.catmaid("^celltype153$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INMBtype5 <- nlapply(
  read.neurons.catmaid("^celltype121$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INRGWa <- nlapply(
  read.neurons.catmaid("^celltype6$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW <- nlapply(
  read.neurons.catmaid("^celltype117$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INZ <- nlapply(
  read.neurons.catmaid("^celltype127$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INproT2 <- nlapply(
  read.neurons.catmaid("^celltype125$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INsn <- nlapply(
  read.neurons.catmaid("^celltype57$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INUturn <- nlapply(
  read.neurons.catmaid("^celltype144$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INUturnMB <- nlapply(
  read.neurons.catmaid("^celltype190$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INbigloop <- nlapply(
  read.neurons.catmaid("^celltype140$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INMBintr <- nlapply(
  read.neurons.catmaid("^MBintrIN$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INMBtype7 <- nlapply(
  read.neurons.catmaid("^celltype122$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MNant <- nlapply(
  read.neurons.catmaid("^celltype19$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
vMN1_2 <- nlapply(
  read.neurons.catmaid("^vMN1-2$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MBmouth <- nlapply(
  read.neurons.catmaid("^celltype189$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MNant_ciliary_target <- nlapply(
  read.neurons.catmaid("^MNant_ciliary_target$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)


# read vMN and presyn cells
SNantler <- nlapply(
  read.neurons.catmaid("^celltype26$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MS <- nlapply(
  read.neurons.catmaid("^MS_pre_vMN$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
eyespotPRCR3 <- nlapply(
  read.neurons.catmaid("^celltype33$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INpreMN <- nlapply(
  read.neurons.catmaid("^celltype23$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
vMN <- nlapply(
  read.neurons.catmaid("^vMN$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MUSlong <- nlapply(
  read.neurons.catmaid("^vMN1-2_mus$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
vMN_ciliated_targets <- nlapply(
  read.neurons.catmaid("^vMN_ciliated_targets$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
acicula <- nlapply(
  read.neurons.catmaid("^acicula$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)


# read eye-eyespot integration celltypes
PRC <- nlapply(
  read.neurons.catmaid("^celltype1$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
eyespotPRCR1 <- nlapply(
  read.neurons.catmaid("^celltype34$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
IN1 <- nlapply(
  read.neurons.catmaid("^celltype2$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INR <- nlapply(
  read.neurons.catmaid("^celltype171$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNDLSO1.1NP <- nlapply(
  read.neurons.catmaid("^celltype172$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INton <- nlapply(
  read.neurons.catmaid("^celltype3$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

# read descending neuron types

SNantler <- nlapply(
  read.neurons.catmaid("celltype26", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNbronto <- nlapply(
  read.neurons.catmaid("celltype168", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
SNblunt <- nlapply(
  read.neurons.catmaid("celltype148", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
hCRdesc <- nlapply(
  read.neurons.catmaid("hCRdesc", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
hPUdesc <- nlapply(
  read.neurons.catmaid("hPUdesc", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INdesc_decus <- nlapply(
  read.neurons.catmaid("INdesc-decuss-head", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INMBproj <- nlapply(
  read.neurons.catmaid("MBprojIN", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INdesc_decuss_head_girdle <- nlapply(
  read.neurons.catmaid("INdesc-decuss-head-girdle", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)


# plot sensory neurons in head --------------------------------------------

plot_background()


# plot sensory cells without soma  --------------------------------------------
plot_SN_No_soma <- function(x) {
  plot3d(Dorsal_SO,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[2]
  )
  plot3d(MBSN,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[1]
  )
  plot3d(Adult_eye,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[3]
  )
  plot3d(nuchal,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[4]
  )
  plot3d(DLSO,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[7]
  )
  plot3d(palp,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[6]
  )
  plot3d(antenna,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[5]
  )
  plot3d(eyespot,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 2,
    rev = FALSE, fixup = F, add = T, alpha = 0.5, col = "black"
  )
  plot3d(ANS_SN,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Tol_muted[3]
  )
}
  # plot3d(SNcirri, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
  #       rev = FALSE, fixup = F, add=T, alpha=1, col=Tol_muted[4])
  # plot3d(HSG3, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
  #       rev = FALSE, fixup = F, add=T, alpha=1, col=Tol_muted[4], depth_mask = FALSE)
  # plot3d(cMS_cell, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
  #       rev = FALSE, fixup = F, add=T, alpha=1, col='red', depth_mask = FALSE)
  # plot3d(headPU, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
  #       rev = FALSE, fixup = F, add=T, alpha=1, col=Tol_muted[6], depth_mask = FALSE)
  # plot3d(antler, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
  #       rev = FALSE, fixup = F, add=T, alpha=1, col=Tol_muted[7], depth_mask = FALSE)
  # plot3d(vMIP, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
  #       rev = FALSE, fixup = F, add=T, alpha=1, col=Tol_muted[8], depth_mask = FALSE)
  # plot3d(SNbicil, WithConnectors = F, WithNodes = F, soma=F, lwd=1,
  #       rev = FALSE, fixup = F, add=T, alpha=1, col='#44BB99', depth_mask = FALSE)

  plot_SN_No_soma()

  # add a text label
  add_text <- function(x) {
    texts3d(69000, 16000, 30000, text = "DSO", col = "black", cex = 2)
    texts3d(35000, 49000, 21000, text = "eye", col = "black", cex = 2)
    texts3d(23000, 48000, 42000, text = "nuchal", col = "black", cex = 2)
    texts3d(46000, 20000, 30000, text = "DLSO", col = "black", cex = 2)
    texts3d(42000, 117000, 7000, text = "palp", col = "black", cex = 2)
    texts3d(87500, 83000, 4000, text = "antenna", col = "black", cex = 2)
    texts3d(20500, 80500, 15000, text = "SNMB", col = "black", cex = 2)
    texts3d(72500, 45500, 5000, text = "AO_SN", col = "black", cex = 2)
    texts3d(24000, 62000, 10000, text = "eyespot", col = "black", cex = 2)
    # texts3d(27300,114500, 47000, text = "cirrus", col='black', cex = 2)
  }
  add_text()
  # plot sensory cells with soma
  plot_SN_soma <- function(x, alpha) {
    alpha <- alpha
    plot3d(Dorsal_SO,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[2]
    )
    plot3d(MBSN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(Adult_eye,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[3]
    )
    plot3d(nuchal,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[4]
    )
    plot3d(DLSO,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[7]
    )
    plot3d(palp,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[6]
    )
    plot3d(antenna,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[5]
    )
    plot3d(eyespot,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 0.6, col = "black"
    )
    plot3d(ANS_SN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Tol_muted[3]
    )
    #  plot3d(SNcirri, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
    #         rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[4])
    #  plot3d(HSG3, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
    #         rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[4], depth_mask = FALSE)
    #  plot3d(cMS_cell, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
    #         rev = FALSE, fixup = F, add=T, alpha=alpha, col='red', depth_mask = FALSE)
    #  plot3d(headPU, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
    #         rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[6], depth_mask = FALSE)
    #  plot3d(antler, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
    #         rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[7], depth_mask = FALSE)
    #  plot3d(vMIP, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
    #         rev = FALSE, fixup = F, add=T, alpha=alpha, col=Tol_muted[8], depth_mask = FALSE)
    #  plot3d(SNbicil, WithConnectors = F, WithNodes = F, soma=T, lwd=1,
    #         rev = FALSE, fixup = F, add=T, alpha=alpha, col='#44BB99', depth_mask = FALSE)
  }
  plot_SN_soma(x, alpha = 0.4)

  rgl.snapshot("pictures/Figure_head_sensory_anterior.png")
  close3d()


  # plot SN neuropils without soma in two separate pictures

  plot_background()
  plot3d(Dorsal_SO,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[2]
  )
  texts3d(69000, 17000, 26000, text = "DSO", col = "black", cex = 4)
  plot3d(Adult_eye,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[3]
  )
  texts3d(27000, 49000, 27000, text = "eye", col = "black", cex = 4)
  plot3d(eyespot,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 2,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = "black"
  )
  texts3d(25000, 62000, 20000, text = "eyespot", col = "black", cex = 3.5)
  plot3d(palp,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[6]
  )
  texts3d(42000, 117000, 7000, text = "palp", col = "black", cex = 4)
  plot3d(DLSO,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[7]
  )
  texts3d(36000, 30000, 30000, text = "DLSO", col = "black", cex = 4)
  rgl.snapshot("pictures/Figure_head_SN_neuropil1.png")
  close3d()

  plot_background()
  plot3d(MBSN,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[1]
  )
  texts3d(25000, 74500, 15000, text = "SNMB", col = "black", cex = 4)
  texts3d(89500, 75000, 4000, text = "antenna", col = "black", cex = 4)
  texts3d(23000, 48000, 42000, text = "nuchal", col = "black", cex = 4)
  texts3d(72500, 45500, 5000, text = "AO_SN", col = "black", cex = 4)
  plot3d(antenna,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 0.6, col = Okabe_Ito[5]
  )
  plot3d(nuchal,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 3,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[4]
  )
  plot3d(ANS_SN,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Tol_muted[3]
  )
  rgl.snapshot("pictures/Figure_head_SN_neuropil2.png")
  close3d()


  # unique head SNs --------------------------------------------

  plot_background()

  plot_SN_unique_soma <- function(x, alpha) {
    alpha <- alpha
    plot3d(HSG3,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Tol_muted[4], depth_mask = FALSE
    )
    plot3d(cMS_cell,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = "red", depth_mask = FALSE
    )
    plot3d(headPU,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Tol_muted[6], depth_mask = FALSE
    )
    plot3d(antler,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Tol_muted[7], depth_mask = FALSE
    )
    plot3d(vMIP,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Tol_muted[8], depth_mask = FALSE
    )
    plot3d(SNbicil,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = "#44BB99", depth_mask = FALSE
    )
  }
  plot_SN_unique_soma(x, alpha = 0.6)


  texts3d(47200, 70000, 9000, text = "HSG3", col = "black", cex = 2)
  texts3d(71300, 114500, 7000, text = "SNbicil", col = "black", cex = 2)
  texts3d(84500, 81000, 4000, text = "SNantler", col = "black", cex = 2)
  texts3d(55000, 43300, 2000, text = "headPU", col = "black", cex = 2)
  texts3d(74000, 58500, 11000, text = "MS1,2", col = "black", cex = 2)
  texts3d(79500, 94000, 4000, text = "SNMIP-vc", col = "black", cex = 2)


  rgl.snapshot("pictures/Figure_head_SN_unique.png")
  close3d()


  # plot interneurons in head  --------------------------------------------

  plot_background()
  # plot IN without soma
  plot_IN_no_soma <- function(x) {
    alpha <- 1
    plot3d(visualIN,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[2]
    )
    plot3d(ANS_IN,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[8]
    )
    plot3d(INdesc_decus,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[4]
    )
    plot3d(INMBintr,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[1]
    )
    plot3d(INMBproj,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[5]
    )
    plot3d(INcentrHead,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[3]
    )
  }
  plot_IN_no_soma()

  add_text_IN <- function(x) {
    texts3d(59000, 32000, 15000, text = "INcentr", col = "black", cex = 2)
    texts3d(35500, 44500, 42000, text = "INeye", col = "black", cex = 2)
    texts3d(87500, 78000, 4000, text = "INdesc-decuss", col = "black", cex = 2)
    texts3d(22500, 83500, 20000, text = "INMBintr", col = "black", cex = 2)
    texts3d(25500, 64500, 25000, text = "INMBproj", col = "black", cex = 2)
    texts3d(74500, 36000, 14000, text = "AO_IN", col = "black", cex = 2)
  }
  add_text_IN()

  # plot cells with soma
  plot_IN_soma <- function(x, alpha) {
    alpha <- alpha
    plot3d(visualIN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha + 0.2, col = Okabe_Ito[2]
    )
    plot3d(ANS_IN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha + 0.2, col = Okabe_Ito[8]
    )
    plot3d(INdesc_decus,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[4]
    )
    plot3d(INMBintr,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(INMBproj,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[5]
    )
    plot3d(INcentrHead,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[3]
    )
  }
  plot_IN_soma(x, alpha = 0.4)


  rgl.snapshot("pictures/Figure_head_IN_anterior.png")
  close3d()

  # plot IN neuropils without soma in two separate pictures

  plot_background()
  plot3d(visualIN,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[2]
  )
  plot3d(ANS_IN,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[8]
  )
  plot3d(INMBintr,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[1]
  )
  plot3d(INMBproj,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[5]
  )
  texts3d(35500, 44000, 32000, text = "INeye", col = "black", cex = 4)
  texts3d(35500, 73500, 20000, text = "INMB", col = "black", cex = 4)
  texts3d(74500, 36000, 20000, text = "AO_IN", col = "black", cex = 4)

  rgl.snapshot("pictures/Figure_head_IN_neuropil1.png")
  close3d()

  plot_background()
  plot3d(INcentrHead,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[3]
  )
  plot3d(INdesc_decus,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[4]
  )
  texts3d(59000, 32000, 20000, text = "INcentr", col = "black", cex = 4)
  texts3d(87500, 78000, 4000, text = "INdesc-decuss", col = "black", cex = 4)

  rgl.snapshot("pictures/Figure_head_IN_neuropil2.png")
  close3d()


  # plot MNs in head --------------------------------------------

  plot_background()

  # plot cells without soma
  plot_MN_no_soma <- function(x) {
    plot3d(vMN,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[1]
    )
    plot3d(ciliomotor,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[5]
    )
    plot3d(MNheadV,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[7]
    )
    plot3d(MNgland,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[8]
    )
  }
  plot_MN_no_soma()

  add_text_MN <- function(x) {
    texts3d(52500, 53800, 30000, text = "MNant", col = "black", cex = 2)
    texts3d(38500, 44700, 33000, text = "Ser-h1", col = "black", cex = 2)
    texts3d(62000, 95300, 7000, text = "vMN", col = "black", cex = 2)
    texts3d(75000, 85000, 9000, text = "MNheadV", col = "black", cex = 2)
    texts3d(38000, 120500, 30000, text = "cMN", col = "black", cex = 2)
    texts3d(85500, 38000, 20000, text = "MNgland", col = "black", cex = 2)
    texts3d(78500, 51400, 20000, text = "MC", col = "black", cex = 2)
    texts3d(68000, 100300, 7000, text = "MNakro", col = "black", cex = 2)
  }


  # plot cells with soma
  plot_MN_soma <- function(x, alpha) {
    alpha <- alpha
    plot3d(vMN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(ciliomotor,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[5]
    )
    plot3d(MNheadV,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[7]
    )
    plot3d(MNgland,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[8]
    )
  }
  plot_MN_soma(x, alpha = 0.6)
  add_text_MN()
  rgl.snapshot("pictures/Figure_head_MN_anterior.png")
  
  nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
  
  close3d()

  # plot MB neurons in head --------------------------------------------

  plot_background()

  # plot cells without soma
  plot_MB_no_soma <- function(x) {
    plot3d(INMBintr,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[1]
    )
    plot3d(INMBproj,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[2]
    )
    plot3d(MBmouth,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[3]
    )
    plot3d(MBSN,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[7]
    )
  }
  plot_MB_no_soma()
  # plot cells with soma
  plot_MB_soma <- function(x, alpha) {
    alpha <- alpha
    plot3d(INMBintr,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(INMBproj,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha + 0.3, col = Okabe_Ito[2]
    )
    plot3d(MBmouth,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[3]
    )
    plot3d(MBSN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[7]
    )
  }
  plot_MB_soma(x, alpha = 0.5)

  add_text_MB <- function(x) {
    texts3d(52500, 72000, 17000, text = "INMBintr", col = "black", cex = 2)
    texts3d(49500, 64500, 15000, text = "INMBproj", col = "black", cex = 2)
    texts3d(32500, 89000, 18000, text = "SNMB", col = "black", cex = 2)
    texts3d(68500, 102000, 18000, text = "MBmouth", col = "black", cex = 2)
  }
  add_text_MB()

  rgl.snapshot("pictures/Figure_head_MB.png")
  close3d()

  # plot postural control system of vMN1-2  --------------------------------------------
  # plot ventral view

  plot_postural <- function(x, alpha) {
    plot3d(SNantler,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(MS,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[4]
    )
    plot3d(eyespotPRCR3,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[3]
    )
    plot3d(INpreMN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[7]
    )
    plot3d(vMN,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[5]
    )
    plot3d(MUSlong,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[6]
    )
    plot3d(vMN_ciliated_targets,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 0.5, col = "grey50"
    )
    plot3d(vMIP,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 0.5, col = Okabe_Ito[2]
    )
  }
  plot_background_ventral()
  plot_postural(x, alpha = 0.7)

  # plot labels at the position of the soma
  {
    texts3d(52075, 79230, 4000, text = "SNantler", col = "black", cex = 2)
    texts3d(76261, 122581, 2000, text = "MS", col = "black", cex = 2)
    texts3d(84261, 152581, 10000, text = "SNMIP-vc", col = "black", cex = 2)
    texts3d(21368, 71912, 22000, text = "eyespotPRCR3", col = "black", cex = 2)
    texts3d(13368, 55912, 40000, text = "ciliated cells", col = "black", cex = 2)
    texts3d(73747, 88903, 36160, text = "INpreMN", col = "black", cex = 2)
    texts3d(103988, 87558.45, 18400, text = "vMN", col = "black", cex = 2)
    texts3d(99988.34, 155558.45, 92400, text = "MUSlong", col = "black", cex = 2)
  }
  rgl.snapshot("pictures/Figure_postural_control_ventr.png")
  close3d()

  # plot anterior view
  plot_background_lower_clip <- function() {
    nopen3d() # opens a pannable 3d window
    # plot3d(outline, WithConnectors = F, WithNodes = F, soma=F, lwd=2,
    #      rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=0.07,
    #     col="#E2E2E2")
    plot3d(bounding_dots,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
      rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 1,
      col = "white"
    )
    #  par3d(windowRect = c(20, 30, 600, 800)) #to define the size of the rgl window
    #  nview3d("ventral", extramat=rotationMatrix(0, 1, 0, 0))
    plot3d(yolk,
      WithConnectors = F, WithNodes = F, soma = F, lwd = 2,
      rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
      col = "#E2E2E2"
    )
    # we define a z clipping plane for the frontal view
    par3d(zoom = 0.52)
    nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
    # z-axis clip
    clipplanes3d(0, 0, -1, 75000)
    # y-axis clip
    clipplanes3d(1, 0, 0.16, 7000)
    # x-axis clip
    clipplanes3d(0, -1, 0.16, 130000)
    par3d(windowRect = c(0, 0, 800, 800)) # resize for frontal view
  }
  plot_background_lower_clip()
  plot_postural(x, alpha = 0.8)
  # plot labels at the position of the soma
  {
    texts3d(57075, 72230, 2000, text = "SNantler", col = "black", cex = 2)
    texts3d(74261, 56581, 2000, text = "MS", col = "black", cex = 2)
    texts3d(77261, 103581, 2000, text = "SNMIP-vc", col = "black", cex = 2)
    texts3d(30368, 75368, 33368, text = "eyespotPRCR3", col = "black", cex = 2)
    texts3d(35368, 123368, 12368, text = "ciliated cells", col = "black", cex = 2)
    texts3d(63747, 36581, 2000, text = "INpreMN", col = "black", cex = 2)
    texts3d(65988, 98558.45, 18400, text = "vMN", col = "black", cex = 2)
    texts3d(115988, 98558.45, 18400, text = "MUSlong", col = "black", cex = 2)
  }

  plot3d(
    scalebar_50um_anterior, lwd = 4,
    add = T, alpha = 1, col = "black"
  )
  texts3d(104000, 116000, 300, text = "50 μm", col = "black", cex = 2)
  
  rgl.snapshot("pictures/Figure_postural_control_ant.png")
  close3d()

  # plot MNant system with direct head sensory inputs  --------------------------------------------

  plot_MNant_circ <- function(x, alpha) {
    plot3d(MNant,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(SNhook,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(SNhorn,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[3]
    )
    plot3d(SNantler,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[7]
    )
    plot3d(SNtorii,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[5]
    )
    plot3d(MNant_ciliary_target,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = 0.4, col = "grey50"
    )
  }

  plot_background_lower_clip()

  plot_MNant_circ(x, alpha = 1)
  # plot labels at the position of the soma
  {
    texts3d(66075, 72230, 1000, text = "SNantler", col = "black", cex = 2)
    texts3d(54261, 41581, 2000, text = "MNant", col = "black", cex = 2)
    texts3d(40261, 71581, 2000, text = "SNhook", col = "black", cex = 2)
    texts3d(15368, 66368, 33368, text = "SNhorn", col = "black", cex = 2)
    texts3d(27368, 83368, 23368, text = "SNtorii", col = "black", cex = 2)
    texts3d(32368, 125368, 12368, text = "ciliated cells", col = "black", cex = 2)
  }
  rgl.snapshot("pictures/Figure_MNant.png")
  close3d()

  # plot eye-eyespot integration --------------------------------------------

  plot_eyes <- function(x, alpha) {
    plot3d(PRC,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(eyespotPRCR1,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[8]
    )
    plot3d(IN1,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[3]
    )
    plot3d(INR,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = 1, col = Okabe_Ito[6]
    )
    plot3d(SNDLSO1.1NP,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[5]
    )
    plot3d(INton,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[7]
    )
  }
  plot_background()
  plot_eyes(x, alpha = 0.8)

  # plot labels at the position of the soma
  {
    texts3d(24075, 46230, 4000, text = "PRC", col = "black", cex = 2)
    texts3d(33261, 34581, 2000, text = "INR", col = "black", cex = 2)
    texts3d(30261, 74581, 2000, text = "eyespotPRCR1", col = "black", cex = 2)
    texts3d(15368, 66368, 33368, text = "SNhorn", col = "black", cex = 2)
    texts3d(32768, 67368, 33368, text = "IN1", col = "black", cex = 2)
    texts3d(48768, 60068, 23368, text = "INton", col = "black", cex = 2)
    texts3d(60768, 34368, 23368, text = "SNDLSO1.1NP", col = "black", cex = 2)
  }
  rgl.snapshot("pictures/Figure_eyes.png")
  close3d()

  # plot head descending pathways --------------------------------------------

  plot_head_desc <- function(x, alpha) {
    plot3d(SNantler,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[1]
    )
    plot3d(SNblunt,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[8]
    )
    plot3d(SNbronto,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[3]
    )
    plot3d(hCRdesc,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[6]
    )
    plot3d(hPUdesc,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[5]
    )
    plot3d(INdesc_decus,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[7]
    )
    plot3d(INMBproj,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[2]
    )
    plot3d(INdesc_decuss_head_girdle,
      WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
      rev = FALSE, fixup = F, add = T, alpha = alpha, col = Okabe_Ito[4]
    )
  }
  plot_background()
  plot_head_desc(x, alpha = 0.6)

  # plot labels at the position of the soma
  texts3d(61261, 17581, 2000, text = "hCR", col = "black", cex = 2)
  texts3d(65368, 77368, 1368, text = "SNantler", col = "black", cex = 2)
  texts3d(57768, 87368, 1368, text = "SNbronto", col = "black", cex = 2)
  texts3d(90768, 67368, 1368, text = "SNblunt", col = "black", cex = 2)
  texts3d(60768, 66368, 1368, text = "hPU", col = "black", cex = 2)
  texts3d(63768, 54368, 3368, text = "INdesc-decuss", col = "black", cex = 2)
  texts3d(33768, 56368, 3368, text = "INMBproj", col = "black", cex = 2)
  texts3d(50261, 98581, 2000, text = "INdesc_decuss_head_girdle", col = "black", cex = 2)

  rgl.snapshot("pictures/Figure_head_desc.png")
  close3d()

  
  # plot all cell groups in head --------------------------------------------
  
  plot_background()
  plot3d(
    scalebar_50um_anterior, lwd = 4,
    add = T, alpha = 1, col = "black"
  )
  texts3d(104000, 116000, 300, text = "50 μm", col = "black", cex = 2)
  
  plot_SN_soma(alpha=0.5)
  plot_SN_unique_soma(alpha=0.5)
  plot_IN_soma(alpha=0.5)
  plot_MB_soma(alpha=0.5)
  plot_MN_soma(alpha=0.5)
  plot_MNant_circ(alpha=0.5)
  plot_eyes(alpha=0.5)
  plot_head_desc(alpha=0.5)
  plot_postural(alpha=0.5)
  
  rgl.snapshot("pictures/Figure_head_all_cells.png")
  close3d()
  
  
  
  # plot synaptic connectivity among cell groups in the head --------------------------------------------


  # these are the SN IN and MN cell groups in the figure
  cell_groups <- list(
    Dorsal_SO, Adult_eye, nuchal, DLSO, MBSN, eyespot, ANS_SN,
    visualIN, ANS_IN, INdesc_decus,
    INcentrHead, INMBproj, INMBintr, vMN, ciliomotor
  )
  N_cell_groups <- length(cell_groups)

  # iterate through cell group neuron lists and get connectivity for all agains all
  # define empty synapse list with the right dimensions
  synapse_list <- vector("list", N_cell_groups * N_cell_groups)
  for (i in 1:N_cell_groups) {
    for (j in 1:N_cell_groups) {
      # get connectors between two cell groups
      presyn_skids <- attr(cell_groups[i][[1]], "df")$skid
      postsyn_skids <- attr(cell_groups[j][[1]], "df")$skid
      connectivity <- catmaid_get_connectors_between(
        pre = presyn_skids,
        post = postsyn_skids, pid = 11
      )
      # check the number of synapses from group1 -> group2
      N_synapses <- dim(connectivity)[1]
      print((i * N_cell_groups - N_cell_groups) + j)
      # add value to synapse list
      synapse_list[[(i * N_cell_groups - N_cell_groups) + j]] <- N_synapses
    }
  }

  # change "NULL" to 0
  synapse_list[synapse_list == "NULL"] <- 0
  # convert synapse list into a matrix of appropriate dimensions
  synapse_matrix <- matrix(unlist(synapse_list), byrow = TRUE, nrow = N_cell_groups)

  cell_group_names <- list(
    "Dorsal_SO", "Adult_eye", "nuchal", "DLSO", "SNMB", "eyespot", "AO_SN",
    "visualIN", "AO_IN", "INdesc_decus",
    "INcentrHead", " INMBproj", "INMBintr", "vMN", "ciliomotor"
  )
  rownames(synapse_matrix) <- as.character(cell_group_names)
  colnames(synapse_matrix) <- as.character(cell_group_names)
  synapse_matrix

  # check matrix with simple heatmap plot
  heatmap((synapse_matrix), # show matrix
    Rowv = NA, Colv = NA,
    cexRow = 0.7, cexCol = 0.7, revC = T,
    scale = "none",
    col = hcl.colors(300, "Oslo", alpha = 1, rev = FALSE, fixup = TRUE),
    symm = F, margins = c(3, 2)
  )

  # plot with ggplot
  as.data.frame((synapse_matrix)) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses") %>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE)) %>%
    ggplot() +
    geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction),
      stroke = 0
    ) +
    scale_x_discrete(limits = (c(
      "Dorsal_SO", "Adult_eye", "nuchal", "DLSO", "palp", "antenna", "SNMB", "eyespot", "AO_SN",
      "visualIN", "AO_IN", "INdesc_decus",
      "INcentrHead", " INMBproj", "INMBintr", "vMN", "ciliomotor", "MNgland", "MNheadV"
    ))) +
    #  coord_flip() +
    scale_y_discrete(limits = rev(c(
      "Dorsal_SO", "Adult_eye", "nuchal", "DLSO", "palp", "antenna", "SNMB", "eyespot", "AO_SN",
      "visualIN", "AO_IN", "INdesc_decus",
      "INcentrHead", " INMBproj", "INMBintr", "vMN", "ciliomotor", "MNgland", "MNheadV"
    ))) +
    #  coord_flip() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 10),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    labs(x = "postsynaptic cell groups", y = "presynaptic cell groups", title = " ") +
    scale_size_area(max_size = 5) +
    guides(color = "legend") +
    scale_colour_gradient2(
      low = "#0072B2",
      mid = "#D55E00",
      high = "#D55E00",
      midpoint = 0.5,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    ) +
    #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
    theme(panel.background = element_rect(fill = "grey98", color = "white"))

  # Saving R ggplot with R ggsave Function
  ggsave("pictures/Head_cellgroups_syn_matrix.pdf",
    width = nrow(synapse_matrix) / 1.3,
    height = ncol(synapse_matrix) / 1.6, limitsize = TRUE,
    units = c("cm")
  )

  # Saving R ggplot with R ggsave Function
  ggsave("pictures/Head_cellgroups_syn_matrix.png",
    width = 1700,
    height = 1300, limitsize = TRUE,
    units = c("px")
  )

  # we do edge weight filtering on the matrix to remove weak edges
  synapse_matrix[synapse_matrix < 10] <- 0

  # with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
  Conn_graph <- graph_from_adjacency_matrix(
    synapse_matrix,
    mode = c("directed"),
    weighted = TRUE,
    diag = TRUE
  )

  # Convert to object suitable for networkD3
  Connectome_graph_d3 <- igraph_to_networkD3(Conn_graph)
  # weights of the edges
  E(Conn_graph)$weight

  # calculate node weighted degree
  degree <- degree(Conn_graph, v = V(Conn_graph), mode = c("all"), loops = TRUE, normalized = FALSE)

  # Convert to object suitable for networkD3
  Connectome_graph_d3 <- igraph_to_networkD3(Conn_graph)

  # assign weights to the nodes
  Connectome_graph_d3$nodes$weight <- degree

  Connectome_graph_d3$links$group <- Connectome_graph_d3$links$source

  # Define 'group' based on SN/IN type:
  Connectome_graph_d3$nodes$group <- as.factor(c(
    "SN", "SN", "SN", "SN", "SN", "SN", "SN",
    "INmid", "INmid", "INmid", "INmid", "INMBintr", "INMBintr", "MN", "MN"
  ))


  # Give a color for each group:
  my_color <- 'd3.scaleOrdinal() .domain(["SN", "INmid","INMBintr", "MN"]) .range(["#E69F00", "#0072B2", "#99DDFF", "#000000"])'

  
  # Plot sankeyNetwork
  SN_head <- networkD3::sankeyNetwork(
    Links = Connectome_graph_d3$links, Nodes = Connectome_graph_d3$nodes, Source = "source",
    Target = "target", NodeID = "name", Value = "value",
    LinkGroup = NULL, units = "", NodeGroup = "group",
    colourScale = my_color, fontSize = 56,
    fontFamily = "Arial", nodeWidth = 30, nodePadding = 150, margin = NULL,
    height = 700, width = 2000, iterations = 1000, sinksRight = TRUE
  )
  SN_head
  saveNetwork(SN_head, "pictures/Sankey_Head_regions.html")
  
  #do snapshot in browser after fine tuning of node positions
  
  #webshot::webshot(
 #   url = "pictures/Sankey_Head_regions.html",
 #   file = "pictures/Sankey_Head_regions.png",
 #  vwidth = 2000, vheight = 750, # define the size of the browser window
 #   cliprect = c(100, 0, 1800, 750), zoom = 5, delay = 2
#)
  
  write.csv(as.data.frame(synapse_matrix), "source_data/Figure6_source_data1.txt")
  read_csv("source_data/Figure6_source_data1.txt")
  
  # make label

  label_conn <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0)
  dim(label_conn) <- c(4, 4)
  rownames(label_conn) <- c("SN", "INmid", "INMB", "MN")
  colnames(label_conn) <- c("SN", "INmid", "INMB", "MN")
  label_conn
  # with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
  label_graph <- graph_from_adjacency_matrix(label_conn, mode = c("directed"), weighted = TRUE, diag = TRUE)

  # Convert to object suitable for networkD3
  label_graph_d3 <- igraph_to_networkD3(label_graph)

  # Define 'group' based on SN/IN type:
  label_graph_d3$nodes$group <- as.factor(c("SN", "INmid", "INMBintr", "MN"))
  # Define 'group' based on SN/IN type:
  label_graph_d3$links$group <- as.factor(c("link"))

  # Give a color for each group:
  my_color_label <- 'd3.scaleOrdinal() .domain(["SN", "INmid","INMB", "MN","link"]) .range(["#E69F00", "#0072B2", "#99DDFF", "#000000","#FFFFFF"])'

  # Plot sankeyNetwork
  label <- networkD3::sankeyNetwork(
    Links = label_graph_d3$links, Nodes = label_graph_d3$nodes, Source = "source",
    Target = "target", NodeID = "name", Value = "value",
    LinkGroup = "group", units = "", NodeGroup = "group",
    colourScale = my_color_label, fontSize = 38,
    fontFamily = "Arial", nodeWidth = 30, nodePadding = 88, margin = NULL,
    height = 100, width = 700, iterations = 10, sinksRight = F
  )
  label

  saveNetwork(label, "pictures/Sankey_label.html")
  webshot2::webshot(
    url = "pictures/Sankey_label.html",
    file = "pictures/Sankey_label.png",
    vwidth = 2000, vheight = 110, # define the size of the browser window
    cliprect = c(50, 20, 610, 90), zoom = 5
  )

 

  # plot synaptic connectivity among postural control --------------------------------------------

  # these are the differentiated MB cell types
  postural_celltypes <- list(SNantler, eyespotPRCR3, INpreMN, INsn, MS, vMIP, vMN_ciliated_targets, MUSlong, vMN)
  N_postural_celltypes <- length(postural_celltypes)
  postural_celltypes_names <- list("SNantler", "eyespotPRCR3", "INpreMN", "INsn", "MS", "vMIP", "vMN_ciliated_targets", "MUSlong", "vMN")


  # iterate through cell group neuron lists and get connectivity for all against all
  # define empty synapse list with the right dimensions
  postural_synapse_list <- vector("list", N_postural_celltypes * N_postural_celltypes)

  for (i in 1:N_postural_celltypes) {
    for (j in 1:N_postural_celltypes) {
      MNant_ciliary_target # get connectors between two cell groups
      presyn_skids <- attr(postural_celltypes[i][[1]], "df")$skid
      postsyn_skids <- attr(postural_celltypes[j][[1]], "df")$skid
      connectivity <- catmaid_get_connectors_between(
        pre = presyn_skids,
        post = postsyn_skids, pid = 11
      )
      # check the number of synapses from group1 -> group2
      N_synapses <- dim(connectivity)[1]
      print((i * N_postural_celltypes - N_postural_celltypes) + j)
      # add value to synapse list
      postural_synapse_list[[(i * N_postural_celltypes - N_postural_celltypes) + j]] <- N_synapses
    }
  }


  postural_synapse_list

  # change "NULL" to 0
  postural_synapse_list[postural_synapse_list == "NULL"] <- 0
  # convert synapse list into a matrix of appropriate dimensions
  postural_synapse_matrix <- matrix(unlist(postural_synapse_list), byrow = TRUE, nrow = N_postural_celltypes)
  postural_synapse_matrix

  rownames(postural_synapse_matrix) <- as.character(postural_celltypes_names)
  colnames(postural_synapse_matrix) <- as.character(postural_celltypes_names)
  postural_synapse_matrix

  # we do edge weight filtering on the matrix to remove weak edges
  postural_synapse_matrix[postural_synapse_matrix < 2] <- 0
  postural_synapse_matrix

  # with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
  Conn_graph_postural <- graph_from_adjacency_matrix(
    postural_synapse_matrix,
    mode = c("directed"),
    weighted = TRUE, diag = TRUE
  )

  # Convert to object suitable for networkD3
  Conn_graph_postural_d3 <- igraph_to_networkD3(Conn_graph_postural)


  # Define 'group' based on SN/IN type:
  Conn_graph_postural_d3$nodes$group <- as.factor(c(
    "SN", "SN", "IN", "IN", "SN",
    "SN", "MN", "effector", "effector"
  ))
  Conn_graph_postural_d3$nodes$group
  # Give a color for each group:
  my_color_post <- 'd3.scaleOrdinal() .domain(["SN", "IN","effector", "MN"]) .range(["#E69F00", "#0072B2", "#99DDFF", "#000000"])'

  # Plot sankeyNetwork
  SN_postural <- networkD3::sankeyNetwork(
    Links = Conn_graph_postural_d3$links, Nodes = Conn_graph_postural_d3$nodes, Source = "source",
    Target = "target", NodeID = "name", Value = "value",
    LinkGroup = NULL, units = "", NodeGroup = "group",
    colourScale = my_color_post, fontSize = 56,
    fontFamily = "Arial", nodeWidth = 30, nodePadding = 88, margin = NULL,
    height = 500, width = 1800, iterations = 1000, sinksRight = F
  )

  SN_postural
  saveNetwork(SN_postural, "pictures/Sankey_vMN_circuit.html")
  webshot2::webshot(
    url = "pictures/Sankey_vMN_circuit.html",
    file = "pictures/Sankey_vMN_circuit.png",
    vwidth = 3000, vheight = 600, # define the size of the browser window
    cliprect = c(0, 0, 1800, 520), zoom = 1, delay = 2
  )

  write.csv(as.data.frame(postural_synapse_matrix), "source_data/Figure7_source_data1.txt")
  read_csv("source_data/Figure7_source_data1.txt")
  
  # plot synaptic connectivity among MNant circuit --------------------------------------------


  # these are the differentiated MNant cell types
  MNant_circ_celltypes <- list(MNant_ciliary_target, SNantler, SNhook, SNtorii, SNhorn, MNant)
  N_MNant_circ_celltypes <- length(MNant_circ_celltypes)
  MNant_circ_celltypes_names <- list("ciliary band cells", "SNantler", "SNhook", "SNtorii", "SNhorn", "MNant")


  # iterate through cell group neuron lists and get connectivity for all against all
  # define empty synapse list with the right dimensions
  MNant_circ_synapse_list <- vector("list", N_MNant_circ_celltypes * N_MNant_circ_celltypes)

  for (i in 1:N_MNant_circ_celltypes) {
    for (j in 1:N_MNant_circ_celltypes) {
      MNant_ciliary_target # get connectors between two cell groups
      presyn_skids <- attr(MNant_circ_celltypes[i][[1]], "df")$skid
      postsyn_skids <- attr(MNant_circ_celltypes[j][[1]], "df")$skid
      connectivity <- catmaid_get_connectors_between(
        pre = presyn_skids,
        post = postsyn_skids, pid = 11
      )
      # check the number of synapses from group1 -> group2
      N_synapses <- dim(connectivity)[1]
      print((i * N_MNant_circ_celltypes - N_MNant_circ_celltypes) + j)
      # add value to synapse list
      MNant_circ_synapse_list[[(i * N_MNant_circ_celltypes - N_MNant_circ_celltypes) + j]] <- N_synapses
    }
  }

  # change "NULL" to 0
  MNant_circ_synapse_list[MNant_circ_synapse_list == "NULL"] <- 0
  # convert synapse list into a matrix of appropriate dimensions
  MNant_circ_synapse_matrix <- matrix(unlist(MNant_circ_synapse_list), byrow = TRUE, nrow = N_MNant_circ_celltypes)
  MNant_circ_synapse_matrix

  rownames(MNant_circ_synapse_matrix) <- as.character(MNant_circ_celltypes_names)
  colnames(MNant_circ_synapse_matrix) <- as.character(MNant_circ_celltypes_names)
  MNant_circ_synapse_matrix

  # we do edge weight filtering on the matrix to remove weak edges
  MNant_circ_synapse_matrix[MNant_circ_synapse_matrix < 1] <- 0
  MNant_circ_synapse_matrix

  # with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
  Conn_graph_MNant_circ <- graph_from_adjacency_matrix(
    MNant_circ_synapse_matrix,
    mode = c("directed"),
    weighted = TRUE, diag = TRUE
  )

  # calculate weighted degrees https://igraph.org/r/doc/strength.html
  weights <- degree(Conn_graph_MNant_circ,
    mode = c("all"),
    loops = TRUE
  )

  # Convert to object suitable for networkD3
  Conn_graph_MNant_circ_d3 <- igraph_to_networkD3(Conn_graph_MNant_circ)

  # assign weights to the nodes
  Conn_graph_MNant_circ_d3$nodes$weight <- weights

  # Define 'group' based on SN/IN type:
  Conn_graph_MNant_circ_d3$nodes$group <- as.factor(c(
    "effector", "SN", "SN", "SN",
    "SN", "MN"
  ))
  Conn_graph_MNant_circ_d3$nodes$group

  # Give a color for each group:
  my_color_MB <- 'd3.scaleOrdinal() .domain(["SN", "effector", "MN"]) .range(["#E69F00", "#0072B2", "#99DDFF", "#000000"])'

  # Plot sankeyNetwork
  SN_ant <- networkD3::sankeyNetwork(
    Links = Conn_graph_MNant_circ_d3$links, Nodes = Conn_graph_MNant_circ_d3$nodes, Source = "source",
    Target = "target", NodeID = "name", Value = "value",
    LinkGroup = NULL, units = "", NodeGroup = "group",
    colourScale = my_color_MB, fontSize = 56,
    fontFamily = "Arial", nodeWidth = 30, nodePadding = 88, margin = NULL,
    height = 230, width = 1800, iterations = 1000, sinksRight = F
  )

  SN_ant
  saveNetwork(SN_ant, "pictures/Sankey_MNant_circuit.html")
  webshot2::webshot(
    url = "pictures/Sankey_MNant_circuit.html",
    file = "pictures/Sankey_MNant_circuit.png",
    vwidth = 3000, vheight = 350, # define the size of the browser window
    cliprect = c(180, 0, 1460, 260), zoom = 5, delay = 2
  )

  write.csv(as.data.frame(MNant_circ_synapse_matrix), "source_data/Figure7_source_data2.txt")
  read_csv("source_data/Figure7_source_data2.txt")

  # plot synaptic connectivity among eye circuits --------------------------------------------
  
  
  # these are the differentiated MNant cell types
  eye_circ_celltypes <- list(PRC, IN1, eyespotPRCR1, INR, SNDLSO1.1NP)
  N_eye_circ_celltypes <- length(eye_circ_celltypes)
  eye_circ_celltypes_names <- list("PRC", "IN1", "eyespot-PRCR1", "INR", 
                                   "SNDLSO1.1NP")
  
  # iterate through cell group neuron lists and get connectivity for all against all
  # define empty synapse list with the right dimensions
  eye_circ_synapse_list <- c()
  
  for (i in 1:N_eye_circ_celltypes) {
    for (j in 1:N_eye_circ_celltypes) {
       # get connectors between two cell groups
      presyn_skids <- attr(eye_circ_celltypes[i][[1]], "df")$skid
      postsyn_skids <- attr(eye_circ_celltypes[j][[1]], "df")$skid
      connectivity <- catmaid_get_connectors_between(
        pre = presyn_skids,
        post = postsyn_skids, pid = 11
      )
      # check the number of synapses from group1 -> group2
      N_synapses <- dim(connectivity)[1]
      if(length(connectivity) == 0) {N_synapses = 0}
      print((i * N_eye_circ_celltypes - N_eye_circ_celltypes) + j)
      # add value to synapse list
      eye_circ_synapse_list <- c(eye_circ_synapse_list, N_synapses)
    }
  }
  eye_circ_synapse_list
  
  # convert synapse list into a matrix of appropriate dimensions
  eye_circ_synapse_matrix <- matrix(unlist(eye_circ_synapse_list), byrow = TRUE, nrow = N_eye_circ_celltypes)
  eye_circ_synapse_matrix
  
  rownames(eye_circ_synapse_matrix) <- as.character(eye_circ_celltypes_names)
  colnames(eye_circ_synapse_matrix) <- as.character(eye_circ_celltypes_names)
  eye_circ_synapse_matrix
  
  # we do edge weight filtering on the matrix to remove weak edges
  eye_circ_synapse_matrix[eye_circ_synapse_matrix < 1] <- 0
  eye_circ_synapse_matrix
  
  # with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
  Conn_graph_eye_circ <- graph_from_adjacency_matrix(
    eye_circ_synapse_matrix,
    mode = c("directed"),
    weighted = TRUE, diag = TRUE
  )
  
  # calculate weighted degrees https://igraph.org/r/doc/strength.html
  weights <- degree(Conn_graph_eye_circ,
                    mode = c("all"),
                    loops = TRUE
  )
  
  # Convert to object suitable for networkD3
  Conn_graph_eye_circ_d3 <- igraph_to_networkD3(Conn_graph_eye_circ)
  
  # assign weights to the nodes
  Conn_graph_eye_circ_d3$nodes$weight <- weights
  
  # Define 'group' based on SN/IN type:
  Conn_graph_eye_circ_d3$nodes$group <- as.factor(c(
    "SN", "IN", "SN", 
    "IN", "SN"  )
    )
  Conn_graph_eye_circ_d3$nodes$group
  
  # Give a color for each group:
  my_color_eye <- 'd3.scaleOrdinal() .domain(["SN", "IN", "MN"]) .range(["#E69F00", "#0072B2", "#99DDFF", "#000000"])'
  
  # Plot sankeyNetwork
  eye_sankey <- networkD3::sankeyNetwork(
    Links = Conn_graph_eye_circ_d3$links, Nodes = Conn_graph_eye_circ_d3$nodes, Source = "source",
    Target = "target", NodeID = "name", Value = "value",
    LinkGroup = NULL, units = "", NodeGroup = "group",
    colourScale = my_color_eye, fontSize = 56,
    fontFamily = "Arial", nodeWidth = 30, nodePadding = 88, margin = NULL,
    height = 360, width = 1500, iterations = 1000, sinksRight = F
  )
  
  eye_sankey
  saveNetwork(eye_sankey, "pictures/Sankey_eye_circuit.html")
  webshot2::webshot(
    url = "pictures/Sankey_eye_circuit.html",
    file = "pictures/Sankey_eye_circuit.png",
    vwidth = 1500, vheight = 360, # define the size of the browser window
    cliprect = c(180, 0, 1120, 360), zoom = 5, delay = 2
  )
  
  write.csv(as.data.frame(eye_circ_synapse_matrix), "source_data/Figure7_source_data3.txt")
  read_csv("source_data/Figure7_source_data3.txt")
  
  
# retrieve head celltype connectivity matrix and plot number of postsyn - presyn partners --------------------------------------------

# function to retrieve skids based on two annotations

  skids_by_2annotations <- function(annotation1, annotation2) {
    annotations_cells <- list()
    annotations_cells[[1]] <- catmaid_get_annotations_for_skeletons(annotation1, pid = 11, conn = conn_http1)
    # we retrieve those skeletons that are also annotated with right_side
    return(unlist(lapply(annotations_cells, function(x) x[x$annotation == annotation2, 1])))
  }
  # function to retrieve skids based on three annotations
  skids_by_3annotations <- function(annotation1, annotation2, annotation3) {
    annotations_cells <- list()
    annotations_cells[[1]] <- catmaid_get_annotations_for_skeletons(annotation1, pid = 11, conn = conn_http1)
    # we retrieve those skeletons that are also annotated with right_side
    skids1 <- unlist(lapply(annotations_cells, function(x) x[x$annotation == annotation2, 1]))
    skids2 <- unlist(lapply(annotations_cells, function(x) x[x$annotation == annotation3, 1]))
    return(intersect(skids1, skids2))
  }


  # retrieve skids

  # define empty SN skids list
  episphere_SN_skids <- list()
  counts <- 0
  # iterate through the celltypes, check if there are skids in the episphere and retrive those head skids in a list
  for (i in c(1:202)) {
    annotation <- paste("annotation:^celltype", i, "$", sep = "")
    # we retrieve those skids that are annotated with celltype_i and 'episphere'
    skids <- skids_by_3annotations(annotation, "episphere", "Sensory neuron")
    # if no skids are returned, next
    if (identical(skids, integer(0))) {
      next
    }
    counts <- counts + 1
    episphere_SN_skids[[counts]] <- skids
    print(i)
    print(skids)
  }

  # define empty IN skids list
  episphere_IN_skids <- list()
  counts <- 0
  # iterate through the celltypes, check if there are skids in the episphere and retrive those head skids in a list
  for (i in c(1:202)) {
    annotation <- paste("annotation:^celltype", i, "$", sep = "")
    # we retrieve those skids that are annotated with celltype_i and 'episphere'
    skids <- skids_by_3annotations(annotation, "episphere", "interneuron")
    # if no skids are returned, next
    if (identical(skids, integer(0))) {
      next
    }
    counts <- counts + 1
    episphere_IN_skids[[counts]] <- skids
    print(i)
    print(skids)
  }

  # define empty MN skids list
  episphere_MN_skids <- list()
  counts <- 0
  # iterate through the celltypes, check if there are skids in the episphere and retrive those head skids in a list
  for (i in c(1:202)) {
    annotation <- paste("annotation:^celltype", i, "$", sep = "")
    # we retrieve those skids that are annotated with celltype_i and 'episphere'
    skids <- skids_by_3annotations(annotation, "episphere", "motorneuron")
    # if no skids are returned, next
    if (identical(skids, integer(0))) {
      next
    }
    counts <- counts + 1
    episphere_MN_skids[[counts]] <- skids
    print(i)
    print(skids)
  }

  episphere_skids <- c(episphere_SN_skids, episphere_IN_skids, episphere_MN_skids)


  # retrieve synaptic connectivity and generate matrix

  # define empty synapse list with the right dimensions
  head_synapse_list <- vector("list", length(episphere_skids) * length(episphere_skids))
  length(episphere_skids)

  # retrieve all against all number of synapses
  for (i in 1:length(episphere_skids)) {
    for (j in 1:length(episphere_skids)) {
      # get connectors between two cell lists
      presyn_skids <- unlist(episphere_skids[i])
      postsyn_skids <- unlist(episphere_skids[j])
      connectivity <- catmaid_get_connectors_between(
        pre = presyn_skids,
        post = postsyn_skids, pid = 11
      )
      # check the number of synapses from group1 -> group2
      N_synapses <- dim(connectivity)[1]
      counter <- ((i * length(episphere_skids) - length(episphere_skids)) + j)
      print(counter)
      # add value to synapse list
      print(N_synapses)
      # change NULL to 0
      if (is.null(N_synapses)) {
        N_synapses <- 0
      }
      head_synapse_list[[counter]] <- N_synapses
    }
  }

  # convert synapse list into a matrix of appropriate dimensions
  head_synapse_matrix <- matrix(unlist(head_synapse_list), byrow = TRUE, nrow = length(episphere_skids))

  # get the neuron name of the first skid in the skids list
  celltype_names <- list()
  for (i in 1:length(episphere_skids)) {
    name <- catmaid_get_neuronnames(episphere_skids[[i]][1], pid = 11)
    name <- sub("_.*$", "", name)
    celltype_names[i] <- name
  }

  # add row and column names to matrix
  rownames(head_synapse_matrix) <- unlist(celltype_names)
  colnames(head_synapse_matrix) <- unlist(celltype_names)
  # check duplicated names
  celltype_names[duplicated(celltype_names)]


  # plot and save matrix

  # plot with heatmaply
  library(heatmaply)
  heatmaply(sqrt(head_synapse_matrix),
    column_text_angle = 90, row_text_angle = 0,
    fontsize_row = 12,
    fontsize_col = 12,
    Rowv = TRUE, Colv = TRUE,
    hclust_method = "ward.D",
    dist_method = "manhattan",
    point_size_mat = sqrt(head_synapse_matrix),
    node_type = "scatter",
    show_dendrogram = c(F, F),
    grid_gap = 0,
    hide_colorbar = F,
    plot_method = c("ggplot"),
    revC = TRUE,
    # col=brewer.pal(9, 'Spectral'),
    col = c("white", "#0072B2", "#0072B2", "#D55E00", "#D55E00", "#D55E00"),
    # col=hcl.colors(20, "Oslo", alpha = 1, rev = F, fixup = F),
    xlab = "postsynaptic cell types",
    ylab = "presynaptic cell types",
    main = "connectivity matrix of head celltypes (sqrt of summed synapses)",
  )

  # save matrix
  write.csv2(head_synapse_matrix, file = "data/head_celltypes_syn_matrix.csv")

  # Saving R ggplot with R ggsave Function
  ggsave("pictures/head_celltypes_syn_matrix.pdf",
    width = 7000,
    height = 6600, limitsize = TRUE,
    units = "px",
    scale = 1
  )
  # Saving R ggplot with R ggsave Function
  ggsave("pictures/head_celltypes_syn_matrix.png",
    width = 8000,
    height = 6400, limitsize = TRUE,
    units = "px"
  )

  # plot with ggplot
  as.data.frame(head_synapse_matrix) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses") %>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE)) %>%
    ggplot() +
    geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction),
      stroke = 0,
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 3),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    ) +
    # order data as in the celltype_names list
    scale_x_discrete(limits = as.character(celltype_names)) +
    scale_y_discrete(limits = as.character(rev(celltype_names))) +
    labs(x = "postsynaptic cell types", y = "presynaptic cell types", title = " ") +
    scale_size_area(max_size = 1) +
    guides(color = "legend") +
    scale_colour_gradient2(
      low = "#0072B2",
      mid = "#D55E00",
      high = "#D55E00",
      midpoint = 0.5,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    ) +
    #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
    theme(panel.background = element_rect(fill = "grey98", color = "white"))

  # Saving R ggplot with R ggsave Function
  ggsave("pictures/head_celltypes_syn_matrix.pdf",
    width = 2000,
    height = 1600, limitsize = TRUE,
    units = c("px")
  )
  # Saving R ggplot with R ggsave Function
  ggsave("pictures/head_celltypes_syn_matrix.png",
    width = 1700,
    height = 1400, limitsize = TRUE,
    units = c("px")
  )

  
# load cell type connectivity ---------------

syn_tb <- readRDS("source_data/Figure4_source_data1.rds")
syn_tb

celltypes_table <- read_csv("supplements/Supplementary_Table1.csv")

head_only_celltypes <- celltypes_table %>% 
  filter(`soma position` == "head") %>%
  select(`name of cell type`) %>%
  pull()
  
# binarize and only consider >2 synapses
syn_tb_bin <- syn_tb %>%
  activate(edges) %>%
  mutate(syn_filtered = ifelse(synapses <= 3, 0, 1))

syn_tb_bin_summed <- syn_tb_bin %>%
  activate(edges) %>%
  select(from, syn_filtered) %>%
  group_by(from) %>%
  mutate(postsyn_partners = sum(syn_filtered)) %>%
  as_tibble() %>%
  ungroup() %>%
  group_by(to) %>%
  mutate(presyn_partners = sum(syn_filtered)) %>%
  ungroup()

#add names of from and to nodes
syn_tb_bin_summed_named <- syn_tb_bin_summed %>%
  mutate(presyn_name = syn_tb %>%
           activate(nodes) %>%
           select(name) %>%
           as_tibble() %>%
           slice(from) %>%
           pull()
  ) %>%
  mutate(postsyn_name = syn_tb %>%
           activate(nodes) %>%
           select(name) %>%
           as_tibble() %>%
           slice(to) %>%
           pull()
  ) %>%
  mutate(presyn_type = syn_tb %>%
           activate(nodes) %>%
           select(group) %>%
           as_tibble() %>%
           slice(from) %>%
           pull()
  ) %>%
  mutate(postsyn_type = syn_tb %>%
           activate(nodes) %>%
           select(group) %>%
           as_tibble() %>%
           slice(to) %>%
           pull()
  ) 

syn_tb_bin_summed_named

syn_tb_bin_summed_named_pre_head <- syn_tb_bin_summed_named %>%
  filter(presyn_name %in% head_only_celltypes)

syn_tb_bin_summed_named_post_head <- syn_tb_bin_summed_named %>%
  filter(postsyn_name %in% head_only_celltypes)


#plot theme

plot_theme <-  theme(
  text = element_text(size = 10),
  axis.text.y = element_text(angle = 90, hjust = 1, size = 10),
  axis.title = element_text(size = 12),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.length.y = unit(1, "mm"),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.key.size = unit(3, "mm"),
  legend.position = "top",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 8),
  axis.line = element_blank()
)

#filter cells to plot
presyn_ordered_names <- syn_tb_bin_summed_named_pre_head %>%
  select(presyn_name, presyn_type, postsyn_partners) %>%
  filter(postsyn_partners > 3) %>%
  arrange(desc(postsyn_partners)) %>%
  pull(presyn_name) %>%
  unique()

postsyn_ordered_names <- syn_tb_bin_summed_named_post_head %>%
  select(postsyn_name, postsyn_type, presyn_partners) %>%
  filter(presyn_partners > 3) %>%
  arrange(desc(presyn_partners)) %>%
  pull(postsyn_name) %>%
  unique()

# plot by number of postsyn cell types ----------
plot_pre_head <- syn_tb_bin_summed_named_pre_head %>%
  select(presyn_name, presyn_type, postsyn_partners) %>%
  unique() %>%
  ggplot(aes(presyn_name, postsyn_partners, fill = presyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "interneuron" = "#0072B2",
    "motoneuron" = "#CC79A7"
  )) +
  guides(fill = "none") +
  geom_col() +
  scale_x_discrete(limits = presyn_ordered_names) +
  labs(
    y = "# postsyn cell types", 
    x = "head neuronal cell types"
    ) +
  plot_theme +
  geom_text(aes(label = (presyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(
    limits = c(0, 18
                                )
  )
plot_pre_head

# plot by number of presyn cell types -------

plot_post_head <- syn_tb_bin_summed_named_post_head %>%
  select(postsyn_name, postsyn_type, presyn_partners) %>%
  unique() %>%
  ggplot(aes(postsyn_name, presyn_partners, fill = postsyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "interneuron" = "#0072B2",
    "motoneuron" = "#CC79A7",
    "effector" = "#DDDDDD"
  ),
  labels = c(
    "sensory neuron" = "SN",
    "interneuron" = "IN",
    "motoneuron" = "MN",
    "effector" = "effector"
  )
  ) +
  geom_col() +
  scale_x_discrete(limits = postsyn_ordered_names) +
  labs(
    y = "# presyn neuron types", 
    x = "head cell types",
    fill = ""
  ) +
  plot_theme +
  geom_text(aes(label = (postsyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(
    limits = c(0, 18
  )
  )

plot_post_head


# plot histogram of number of cells per head cell type -----

hist_cells_per_head_type <- celltypes_table %>% 
  filter(`soma position` == "head") %>%
  filter(`Sensory/inter/motor neuron` == "Sensory neuron" |
           `Sensory/inter/motor neuron` == "Interneuron" |
           `Sensory/inter/motor neuron` == "Motorneuron") %>%
  ggplot(aes(`number of cells`)) +
  geom_histogram(fill = "grey30", color = '#EEEEEE') +
  labs(y = "# head neuronal celltypes", x = "# cells per celltype") +
  scale_x_continuous(breaks = c(2, 5, 10, 20, 28), labels = c("2", "5", "10", "20", "28")) +
  theme(
      text = element_text(size = 10),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      axis.text.y = element_text(angle = 90, hjust = 1, size = 10),
      axis.title = element_text(size = 12),
      axis.ticks.length.x = unit(1, "mm"),
      axis.ticks.length.y = unit(1, "mm"),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.1)
    )
hist_cells_per_head_type


# Assemble Figure 6 ---------------------------------------------

# Figure panel A

# read png convert to image panel
panelA <- ggdraw() + draw_image(readPNG("pictures/Figure_head_all_cells.png"), scale = 1) +
    draw_label("all neurons",
      x = 0.22, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )  +
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
    

# create panel A
Fig1A <- plot_grid(panelA,
    ncol = 1,
    align = "h",
    # A negative rel_height shrinks space between elements
    rel_widths = c(1),
    rel_heights = c(1)
  )
rm(panelA)


# Figure panel B

img_SN <- readPNG("pictures/Figure_head_sensory_anterior.png")
img2 <- readPNG("pictures/Figure_head_SN_neuropil1.png")
img3 <- readPNG("pictures/Figure_head_SN_neuropil2.png")

panel_SN <- ggdraw() + draw_image(img_SN, scale = 1) +
    draw_label("Sensory",
               x = 0.22, y = 0.95, fontfamily = "sans", fontface = "plain",
               color = "black", size = 12, alpha = 1
    ) 
panelB1 <- ggdraw() + draw_image(img2, scale = 1)
panelB2 <- ggdraw() + draw_image(img3, scale = 1)

Fig1B <- plot_grid(panelB1, panelB2,
                     ncol = 1,
                     align = "h",
                     # A negative rel_height shrinks space between elements
                     rel_widths = c(1),
                     rel_heights = c(1)
  )
  
# plot in a multi panel figure
Fig1B <- plot_grid(panelB1, panelB2,
    ncol = 1,
    align = "h",
    # A negative rel_height shrinks space between elements
    rel_widths = c(1),
    rel_heights = c(1)
  )
rm(img2, img3, panelB1, panelB2)


# Figure panel C
# read png
img4 <- readPNG("pictures/Figure_head_SN_unique.png")
# convert png to image panel
panelC <- ggdraw() + draw_image(img4, scale = 1) +
    draw_label("other sensory types",
      x = 0.35, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )


# create panel C
Fig1C <- plot_grid(panelC,
    ncol = 1,
    align = "h",
    # A negative rel_height shrinks space between elements
    rel_widths = c(1),
    rel_heights = c(1)
  )
rm(img4, panelC)


# Figure panel D
# read png
img5 <- readPNG("pictures/Figure_head_IN_anterior.png")
# convert png to image panel
panelD <- ggdraw() + draw_image(img5, scale = 1) +
    draw_label("Inter",
      x = 0.22, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )


# create panel C
Fig1D <- plot_grid(panelD,
    ncol = 1,
    align = "h",
    # A negative rel_height shrinks space between elements
    rel_widths = c(1),
    rel_heights = c(1)
  )
rm(img5, panelD)


# Figure panel E
img6 <- readPNG("pictures/Figure_head_IN_neuropil1.png")
img7 <- readPNG("pictures/Figure_head_IN_neuropil2.png")

panelE1 <- ggdraw() + draw_image(img6, scale = 1)
panelE2 <- ggdraw() + draw_image(img7, scale = 1)

# plot in a multi panel figure
Fig1E <- plot_grid(panelE1, panelE2,
    ncol = 1,
    align = "h",
    # A negative rel_height shrinks space between elements
    rel_widths = c(1),
    rel_heights = c(1)
  )
rm(img6, img7, panelE1, panelE2)


# combine panels  A-E
Fig1A_E <- plot_grid(Fig1A, panel_SN, Fig1B, Fig1D, Fig1E, 
    ncol = 5,
    rel_widths = c(2, 2, 1, 2, 1),
    labels = c("A", "B", "C", "D", "E", "F"),
    label_size = 14, label_y = 1, label_x = 0,
    label_fontfamily = "sans", label_fontface = "plain"
  ) +
    theme(plot.margin = unit(c(1, 1, 1, 1), units = "pt"), )

# Figure panel F-J
# read png
img10 <- readPNG("pictures/Figure_head_MN_anterior.png")
# convert pdf to image panel
panelF <- ggdraw() + draw_image(img10, scale = 1) +
    draw_label("Motor",
      x = 0.22, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )
img11 <- readPNG("pictures/Sankey_Head_regions.png")
# convert pdf to image panel
panelG <- ggdraw() + draw_image(img11, scale = 1) +
    draw_label("connectivity between head regions",
      x = 0.3, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )

# create panel
Fig1F_H <- plot_grid(Fig1C, panelF, panelG,
    ncol = 3,
    align = "h",
    labels = c("F", "G", "H"),
    label_size = 14, label_y = 1, label_x = 0,
    label_fontfamily = "sans", label_fontface = "plain",
    # A negative rel_height shrinks space between elements
    rel_widths = c(2, 2, 4),
    rel_heights = c(1)
  )
rm(img10, img11, panelF, panelG)

Fig6 <- plot_grid(Fig1A_E, Fig1F_H,
    ncol = 1,
    rel_heights = c(1, 1)
  ) +
    theme(plot.margin = unit(c(1, 1, 1, 1), units = "pt"))

  
ggsave("Figures/Figure6.png",
    limitsize = FALSE,
    units = c("px"), Fig6,
    width = 4800, height = 2400, bg = "white"
  )
  
ggsave("Figures/Figure6.pdf",
         limitsize = FALSE,
         units = c("px"), Fig6, width = 4800, height = 2400
  )
  

# Assemble Figure 7 ---------------------------------------------

# read png
img14 <- readPNG("pictures/Postural_control_schematic.png")
img15 <- readPNG("pictures/Figure_postural_control_ventr.png")
img16 <- readPNG("pictures/Figure_postural_control_ant.png")
img17 <- readPNG("pictures/Sankey_vMN_circuit.png")
img18 <- readPNG("pictures/Sankey_MNant_circuit.png")
img19 <- readPNG("pictures/Figure_MNant.png")

# convert png to imagae panel
panelA <- ggdraw() + draw_image(img14, scale = 1) +
    draw_label("postural control",
      x = 0.5, y = 0.99, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )
panelB <- ggdraw() + draw_image(img15, scale = 1) +
    draw_label("ventral view",
      x = 0.35, y = 0.99, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    ) +
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

panelC <- ggdraw() + draw_image(img16, scale = 1) +
    draw_label("anterior view",
      x = 0.3, y = 0.99, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    ) +
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
    draw_label("d", x = 0.1, y = 0.93, size = 8) +
    draw_label("v", x = 0.1, y = 0.79, size = 8) 

panelD <- ggdraw() + draw_image(img17, scale = 1) +
    draw_label("postural control circuit",
      x = 0.3, y = 0.99, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )
panelE <- ggdraw() + draw_image(img18, scale = 1) +
    draw_label("MNant circuit",
      x = 0.25, y = 0.93, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )
panelF <- ggdraw() + draw_image(img19, scale = 1) +
    draw_label("MNant ciliomotor circuit",
      x = 0.5, y = 0.96, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )


# read png
img20 <- readPNG("pictures/Figure_eyes.png")
img21 <- readPNG("pictures/Sankey_eyes.png")

# convert png to image panel
panelG <- ggdraw() + draw_image(img20, scale = 1) +
    draw_label("eyes and eyespot",
      x = 0.4, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )
# convert png to imagae panel
panelH <- ggdraw() + draw_image(img21, scale = 1) +
    draw_label("eyes and eyespot",
      x = 0.37, y = 0.94, fontfamily = "sans", fontface = "plain",
      color = "black", size = 12, alpha = 1
    )

#define layout

layout = "
AAABBBCCCCDDDDDFFFF
AAABBBCCCCEEEEEFFFF
GGGGHHHH###########
GGGGHHHHJJJJJKKKKKK
GGGGIIIIJJJJJKKKKKK
"

Figure7 <- panelA + panelB + panelC +
  panelD + panelE + panelF + panelG + panelH +
  hist_cells_per_head_type +
  plot_post_head + plot_pre_head +
  plot_layout(
  design = layout, 
  heights = c(1, 1, 0.1, 0.2, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure7.png",
    limitsize = FALSE,
    units = c("px"), Figure7,
    width = 4800, height = 2500, bg = "white"
  )

ggsave("Figures/Figure7.pdf",
       limitsize = FALSE,
       units = c("px"), Figure7, width = 4800, height = 2500
)

