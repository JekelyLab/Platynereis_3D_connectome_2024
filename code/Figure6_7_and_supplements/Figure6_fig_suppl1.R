# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# load and plot neuron groups ---------------------------------------------

# load stomodeum and MB cells as reference
stomodeum <- nlapply(
  read.neurons.catmaid("stomodeum", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
mushroom_body <- nlapply(
  read.neurons.catmaid("mushroom body", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)


# plot neurons and export fig panels
{
  read_plot_neurons("^celltype1$", "^celltype2$", "celltype171")
  plot3d(scalebar_50um_anterior,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 5,
    rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 1,
    col = "black"
  )
  # use expression to define mu with its latex code
  texts3d(25000, 49000, 27000, text = "eye", col = "black", cex = 2)
  texts3d(32000, 71000, 52000, text = "IN1", col = "black", cex = 2)
  texts3d(34000, 38000, 30000, text = "INR", col = "black", cex = 2)
  rgl.snapshot("pictures/Figure_head_Suppl_eye.png")
  close3d()

  read_plot_neurons("^celltype168$", "^celltype149$", "celltype87", "celltype58")
  texts3d(60000, 79000, 7000, text = "SNbronto", col = "black", cex = 2)
  texts3d(45000, 120000, 3000, text = "INsplitBronto", col = "black", cex = 2)
  texts3d(76000, 120000, 3000, text = "MC3cover", col = "black", cex = 2)
  texts3d(35000, 60000, 3000, text = "INrope", col = "black", cex = 2)
  plot3d(mushroom_body,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
    rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
    col = "grey20"
  )
  rgl.snapshot("pictures/Figure_head_Suppl_SNbronto.png")
  close3d()

  read_plot_neurons("^celltype50$", "^celltype11$")
  texts3d(72000, 104000, 4000, text = "SNbicil", col = "black", cex = 2)
  texts3d(41000, 117000, 3000, text = "cMNATO", col = "black", cex = 2)
  rgl.snapshot("pictures/Figure_head_Suppl_SNbicil.png")
  close3d()


  read_plot_neurons("^celltype13$", "^celltype55$", "^celltype56$")
  texts3d(32000, 45000, 4000, text = "SNnuch", col = "black", cex = 2)
  texts3d(57000, 68000, 3000, text = "INarc1", col = "black", cex = 2)
  texts3d(61000, 48000, 5000, text = "INarc2", col = "black", cex = 2)
  rgl.snapshot("pictures/Figure_head_Suppl_SNnuch.png")
  close3d()


  read_plot_neurons("^celltype130$", "^celltype6$")
  texts3d(59000, 34000, 4000, text = "SN-NS5", col = "black", cex = 2)
  texts3d(50000, 53000, 52000, text = "INRGW", col = "black", cex = 2)
  rgl.snapshot("pictures/Figure_head_Suppl_SN-NS5.png")
  close3d()


  read_plot_neurons("^celltype5$", "^celltype6$", "^celltype7$")
  texts3d(54000, 34000, 4000, text = "cPRC", col = "black", cex = 2)
  texts3d(50000, 53000, 5000, text = "INRGW", col = "black", cex = 2)
  texts3d(70000, 35000, 5000, text = "INNOS", col = "black", cex = 2)
  rgl.snapshot("pictures/Figure_head_Suppl_cPRC.png")
  close3d()


  read_plot_neurons("^celltype34$", "^celltype171$")
  texts3d(34000, 64000, 4000, text = "eyespot-PRCR1", col = "black", cex = 2)
  texts3d(34000, 38000, 30000, text = "INR", col = "black", cex = 2)
  rgl.snapshot("pictures/Figure_head_Suppl_eyespot-PRCR1.png")
  close3d()



  read_plot_neurons_z2(
    "^CRhead$", "^celltype60$", "^celltype55$",
    "^celltype73$", "^celltype58$"
  )
  texts3d(52000, 25000, 4000, text = "CRhead", col = "black", cex = 2)
  texts3d(68000, 121000, 3000, text = "INCM", col = "black", cex = 2)
  texts3d(57000, 58000, 3000, text = "INarc1", col = "black", cex = 2)
  texts3d(42000, 120000, 5000, text = "INsplitCR", col = "black", cex = 2)
  texts3d(28000, 68000, 3000, text = "INrope", col = "black", cex = 2)
  plot3d(mushroom_body,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
    rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
    col = "grey20"
  )
  rgl.snapshot("pictures/Figure_head_Suppl_CRhead.png")
  close3d()

  read_plot_neurons("^celltype52$", "^celltype6$", "^celltype144$")
  texts3d(54000, 31000, 4000, text = "SN47Ach", col = "black", cex = 2)
  texts3d(50000, 53000, 5000, text = "INRGW", col = "black", cex = 2)
  texts3d(41000, 84000, 5000, text = "INUturn", col = "black", cex = 2)
  rgl.snapshot("pictures/Figure_head_Suppl_SN47Ach.png")
  close3d()
}


# MB plotting -------------------------------------------------------------


# load and plot different MB cell types and their postsyn partners

read_plot_neurons("^celltype187$", "^celltype191$", "^celltype183$")
texts3d(30000, 70000, 4000, text = "SNtorii", col = "black", cex = 2)
texts3d(38000, 79000, 3000, text = "INtorii", col = "black", cex = 2)
texts3d(42000, 57000, 5000, text = "INhorn", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_SNtorii.png")
close3d()

read_plot_neurons("^celltype17$", "^celltype183$", "^celltype19$")
texts3d(22000, 57000, 4000, text = "SNhorn", col = "black", cex = 2)
texts3d(41000, 57000, 3000, text = "INhorn", col = "black", cex = 2)
texts3d(52000, 50000, 52000, text = "MNant", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_SNhorn.png")
close3d()

read_plot_neurons("^celltype18$", "^celltype19$", "celltype153")
texts3d(36000, 72000, 7000, text = "SNhook", col = "black", cex = 2)
texts3d(52000, 50000, 52000, text = "MNant", col = "black", cex = 2)
texts3d(65000, 53000, 3000, text = "INdecusshook", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_SNhook.png")
close3d()

read_plot_neurons("^celltype20$", "^celltype160$", "^celltype6$", "^celltype192$", "^INW_mature$")
texts3d(32000, 55000, 4000, text = "SNlasso", col = "black", cex = 2)
texts3d(60000, 34000, 3000, text = "INlasso_postSN", col = "black", cex = 2)
texts3d(50000, 53000, 52000, text = "INRGW", col = "black", cex = 2)
texts3d(42000, 73000, 5000, text = "INUturnMB", col = "black", cex = 2)
texts3d(64000, 77000, 2000, text = "INW", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_SNlasso.png")
close3d()


# load and plot different MB cell types and cell groups
read_plot_neurons("^celltype133$", "^celltype110$", "^celltype138$", "^celltype131$", "^celltype137$")
texts3d(24000, 71000, 14000, text = "SN-NS19", col = "black", cex = 2)
texts3d(50000, 84000, 12000, text = "SN-NS1", col = "black", cex = 2)
texts3d(32000, 60000, 11000, text = "SN-NS18", col = "black", cex = 2)
texts3d(48500, 70000, 18000, text = "SN-NS17", col = "black", cex = 2)
texts3d(48000, 96000, 18000, text = "SN-NS27", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_MBSN-NS.png")
close3d()

# load and plot different MB cell types and cell groups
read_plot_neurons("^INMBtype1$", "^INMBtype2$", "^INMBtype3$", "^INMBtype4$")
texts3d(26000, 67000, 27000, text = "INMBtype1", col = "black", cex = 2)
texts3d(52000, 58000, 27000, text = "INMBtype2", col = "black", cex = 2)
texts3d(52000, 65000, 17000, text = "INMBtype3", col = "black", cex = 2)
texts3d(42000, 85000, 23000, text = "INMBtype4", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_MBtype1_4.png")
close3d()



# load and plot different MB cell types and cell groups
read_plot_neurons("^INMBtype5$", "^INMBtype6$", "^INMBtype7$")
texts3d(23000, 54000, 12000, text = "INMBtype5", col = "black", cex = 2)
texts3d(23000, 71000, 17000, text = "INMBtype6", col = "black", cex = 2)
texts3d(52000, 57000, 22000, text = "INMBtype7", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_MBtype5_7.png")
close3d()

# load and plot different MB cell types and cell groups
read_plot_neurons("^celltype25$", "^celltype190$", "^celltype16$")
texts3d(23000, 59000, 22000, text = "SNmus", col = "black", cex = 2)
texts3d(38000, 86000, 27000, text = "SNtrpa", col = "black", cex = 2)
texts3d(49000, 62000, 22000, text = "SNgolden", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_SNMBtypes.png")
close3d()

# load and plot different MB cell types and cell groups
read_plot_neurons("^celltype140$", "^celltype184$", "^celltype58$")
texts3d(23000, 75000, 22000, text = "INbigloop", col = "black", cex = 2)
texts3d(44000, 63000, 22000, text = "INMBPDF", col = "black", cex = 2)
texts3d(47000, 64500, 3000, text = "INrope", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_INproj1.png")
close3d()


# load and plot different MB cell types and cell groups
read_plot_neurons("^celltype194$", "^celltype185$", "^celltype195$")
texts3d(25000, 71000, 12000, text = "INMBdesc2", col = "black", cex = 2)
texts3d(102000, 51000, 12000, text = "INMBdesc3", col = "black", cex = 2)
texts3d(50000, 58000, 22000, text = "INMBdescFMRF", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
rgl.snapshot("pictures/Figure_head_Suppl_INproj2.png")
close3d()


# load and plot different MB cell types and cell groups
read_plot_neurons("^SNMBdev$", "^INMBdev$", "^celltype189$")
texts3d(22000, 49000, 22000, text = "SNMBdev", col = "black", cex = 2)
texts3d(50000, 59000, 22000, text = "INMBdev", col = "black", cex = 2)
texts3d(45000, 52000, 12000, text = "MBmouth", col = "black", cex = 2)
plot3d(mushroom_body,
  WithConnectors = F, WithNodes = F, soma = T, lwd = 0.1,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.1,
  col = "grey20"
)
plot3d(stomodeum,
  WithConnectors = F, WithNodes = F, soma = F, lwd = 4,
  rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.04,
  col = "grey10"
)
rgl.snapshot("pictures/Figure_head_Suppl_MB_SNdev_INMBdev.png")
close3d()



# Make multi-panel figure -------------------------------------------------

{
  panelA <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_eye.png")) +
    draw_label("eye PRC",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    # use expression to define mu with its latex code
    draw_label(expression(paste("50 ", mu, "m")),
      x = 0.73, y = 0.08, fontfamily = "sans", fontface = "plain",
      color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1
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
  panelB <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_eyespot-PRCR1.png")) +
    draw_label("eyespot-PRCR1",
      x = 0.35, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelC <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_cPRC.png")) +
    draw_label("ciliary PRC",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelD <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNnuch.png")) +
    draw_label("nuchal organ",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelE <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SN47Ach.png")) +
    draw_label("SN47Ach",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelF <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SN-NS5.png")) +
    draw_label("SN-NS5",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelG <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNbicil.png")) +
    draw_label("SNbicil",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelH <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNbronto.png")) +
    draw_label("SNbronto",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelI <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_CRhead.png")) +
    draw_label("head CR",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelJ <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNhorn.png")) +
    draw_label("SNhorn",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelK <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNlasso.png")) +
    draw_label("SNlasso",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelL <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNhook.png")) +
    draw_label("SNhook",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelM <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNtorii.png")) +
    draw_label("SNtorii",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )


  panelN <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_MBSN-NS.png")) +
    draw_label("SNMB-NS",
      x = 0.25, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelO <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_SNMBtypes.png")) +
    draw_label("other SNMB types",
      x = 0.35, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelP <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_MBtype1_4.png")) +
    draw_label("MBIN morph. types 1-4",
      x = 0.45, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelQ <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_MBtype5_7.png")) +
    draw_label("MBIN morph. types 5-7",
      x = 0.45, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelR <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_INproj1.png")) +
    draw_label("MB projection",
      x = 0.35, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelS <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_INproj2.png")) +
    draw_label("MB projection",
      x = 0.4, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
  panelT <- ggdraw() + draw_image(readPNG("pictures/Figure_head_Suppl_MB_SNdev_INMBdev.png")) +
    draw_label("MB developing",
      x = 0.40, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11, angle = 0, lineheight = 0.9, alpha = 1
    )
}

{
  Fig <- plot_grid(panelA, panelB, panelC, panelD, panelE, panelF, panelG, panelH, panelI, panelJ, panelK, panelL, panelM,
    panelN, panelO, panelP, panelQ, panelR, panelS, panelT,
    ncol = 5,
    rel_heights = c(1, 1, 1),
    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T"),
    label_size = 12, label_y = 1, label_x = -0.04,
    label_fontfamily = "sans", label_fontface = "plain"
  ) +
    theme(plot.margin = unit(c(1, 1, 1, 3), units = "pt"))

  ggsave("Figures/Figure6_fig_suppl1.pdf",
    limitsize = FALSE,
    units = c("px"), Fig, width = 4000, height = 3400
  )
  ggsave("Figures/Figure6_fig_suppl1.png",
    limitsize = FALSE,
    units = c("px"), Fig, width = 4000, height = 3400, bg = "white"
  )
}
