# R code to generate Figure 3 fig suppl 4of the 3d Platynereis connectome paper
# Uses Natverse/catmaid and accesses the data on CATMAID
# Gaspar Jekely 2022

# load packages, functions and anatomical references
source("code/Natverse_functions_and_conn.R")

# read and plot example cell types -------------------------------------
{
  # motoneurons
  read_plot_neuron("celltype164", "#0072B2")
  texts3d(30000, 8000, 6000, text = "MNakro", cex = 4)
  plot3d(scalebar_50um_anterior, lwd = 3, color = "black")
  rgl.snapshot("pictures/MNakro.png")
  close3d()

  read_plot_neuron("celltype166", "#0072B2")
  texts3d(30000, 8000, 6000, text = "MNgland", cex = 4)
  rgl.snapshot("pictures/MNgland.png")
  close3d()

  read_plot_neuron("celltype178", "#0072B2")
  texts3d(30000, 8000, 6000, text = "MNheadV", cex = 4)
  rgl.snapshot("pictures/MNheadV.png")
  close3d()

  read_plot_neuron_ventral("celltype151", "#0072B2")
  texts3d(30000, 8000, 20000, text = "MNladder", cex = 3)
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
  rgl.snapshot("pictures/MNladder.png")
  close3d()



  # interneurons
  read_plot_neuron("celltype199", "#CC79A7")
  texts3d(33000, 8000, 6000, text = "INdecussM", cex = 4)
  rgl.snapshot("pictures/INdecussM.png")
  close3d()

  read_plot_neuron("celltype179", "#CC79A7")
  texts3d(35000, 8000, 6000, text = "INsqNSasym", cex = 4)
  rgl.snapshot("pictures/INsqNSasym.png")
  close3d()

  read_plot_neuron("celltype196", "#CC79A7")
  texts3d(30000, 8000, 6000, text = "INpara", cex = 4)
  rgl.snapshot("pictures/INpara.png")
  close3d()

  read_plot_neuron("celltype191", "#CC79A7")
  texts3d(30000, 8000, 6000, text = "INtorii", cex = 4)
  rgl.snapshot("pictures/INtori.png")
  close3d()

  read_plot_neuron_ventral("celltype150", "#CC79A7")
  texts3d(33000, 8000, 20000, text = "INleucoPU", cex = 3)
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
  rgl.snapshot("pictures/INleucoPU.png")
  close3d()

  read_plot_neuron("celltype75", "#CC79A7")
  texts3d(30000, 8000, 6000, text = "INMC3", cex = 4)
  rgl.snapshot("pictures/INMC.png")
  close3d()

  read_plot_neuron("celltype22", "#CC79A7")
  texts3d(35000, 8000, 6000, text = "INdecussPre", cex = 4)
  rgl.snapshot("pictures/INdecussPre.png")
  close3d()

  read_plot_neuron("celltype144", "#CC79A7")
  texts3d(30000, 8000, 6000, text = "INUturn", cex = 4)
  rgl.snapshot("pictures/INUturn.png")
  close3d()


  read_plot_neuron("celltype186", "#CC79A7")
  texts3d(30000, 8000, 6000, text = "INDLSO", cex = 4)
  rgl.snapshot("pictures/INDLSO.png")
  close3d()

  read_plot_neuron_ventral("celltype145", "#CC79A7")
  texts3d(33000, 8000, 20000, text = "INATOpyg", cex = 3)
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
  rgl.snapshot("pictures/INATOpyg.png")
  close3d()

  # Sensory neurons
  read_plot_neuron("celltype16", "#E69F00")
  texts3d(30000, 8000, 6000, text = "SNgolden", cex = 4)
  rgl.snapshot("pictures/SNgold.png")
  close3d()

  read_plot_neuron("celltype187", "#E69F00")
  texts3d(30000, 8000, 6000, text = "SNtorii", cex = 4)
  rgl.snapshot("pictures/SNtori.png")
  close3d()

  read_plot_neuron("celltype26", "#E69F00")
  texts3d(35000, 8000, 6000, text = "SNantlerPDF", cex = 4)
  rgl.snapshot("pictures/SNantlerPDF.png")
  close3d()

  read_plot_neuron("celltype48", "#E69F00")
  texts3d(30000, 8000, 6000, text = "SN-WLD", cex = 4)
  rgl.snapshot("pictures/SN-WLD.png")
  close3d()

  read_plot_neuron_ventral("celltype170", "#E69F00")
  texts3d(30000, 8000, 20000, text = "SNpygM", cex = 3)
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
  rgl.snapshot("pictures/SNpyg.png")
  close3d()
  
  read_plot_neuron_ventral("celltype67", "#CC79A7")
  texts3d(30000, 8000, 20000, text = "MNbow", cex = 3)
  plot3d(scalebar_50um_ventral, lwd = 3, color = "black")
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  par3d(zoom = 0.53)
  rgl.snapshot("pictures/MNbow.png")
  close3d()
}


# assemble figure ------------------------------

MNakro <- ggdraw() + draw_image(readPNG("pictures/MNakro.png"), scale = 1) +
  draw_label(expression(paste("50 ", mu, " m")), x = 0.68, y = 0.07, size = 10) +
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

MNgland <- ggdraw() + draw_image(readPNG("pictures/MNgland.png"), scale = 1)
MNheadV <- ggdraw() + draw_image(readPNG("pictures/MNheadV.png"), scale = 1)
MNladder <- ggdraw() + draw_image(readPNG("pictures/MNladder.png"), scale = 1)

INsqNSasym <- ggdraw() + draw_image(readPNG("pictures/INsqNSasym.png"), scale = 1)
INdecussM <- ggdraw() + draw_image(readPNG("pictures/INdecussM.png"), scale = 1)
INpara <- ggdraw() + draw_image(readPNG("pictures/INpara.png"), scale = 1)
INtorii <- ggdraw() + draw_image(readPNG("pictures/INtori.png"), scale = 1)
INleucoPU <- ggdraw() + draw_image(readPNG("pictures/INleucoPU.png"), scale = 1)
INMC3 <- ggdraw() + draw_image(readPNG("pictures/INMC.png"), scale = 1)
INdecussPre <- ggdraw() + draw_image(readPNG("pictures/INdecussPre.png"), scale = 1)
INUturn <- ggdraw() + draw_image(readPNG("pictures/INUturn.png"), scale = 1)
INDLSO <- ggdraw() + draw_image(readPNG("pictures/INDLSO.png"), scale = 1)
INATOpyg <- ggdraw() + draw_image(readPNG("pictures/INATOpyg.png"), scale = 1)

MNbow <- ggdraw() + draw_image(readPNG("pictures/MNbow.png"), scale = 1) +
  draw_label(expression(paste("50 ", mu, " m")), x = 0.83, y = 0.1, size = 10)  +
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


SNgolden <- ggdraw() + draw_image(readPNG("pictures/SNgold.png"), scale = 1)
SNtorii <- ggdraw() + draw_image(readPNG("pictures/SNtori.png"), scale = 1)
SNantlerPDF <- ggdraw() + draw_image(readPNG("pictures/SNantlerPDF.png"), scale = 1)
SNWLD <- ggdraw() + draw_image(readPNG("pictures/SN-WLD.png"), scale = 1)
SNpygM <- ggdraw() + draw_image(readPNG("pictures/SNpyg.png"), scale = 1)

neurons_ant <- plot_grid(MNakro, MNgland, INdecussPre, INUturn, INsqNSasym, 
  INpara, SNWLD, SNgolden, MNheadV, INMC3, 
  SNantlerPDF, INDLSO, INdecussM, INtorii, SNtorii, 
  ncol = 5, rel_heights = c(1, 1, 1)
)

neurons_vent <- plot_grid(
  MNbow, INATOpyg, SNpygM, INleucoPU, MNladder,
                         ncol = 5, rel_heights = c(1, -0.2, 1, -0.2, 1, -0.2, 1.5)
)

# layout and save ---------------------------------------------------------

# define layout with textual representation for pathcwork assembly of figure
layout <- "
aaaa
AAAA
"

Figure3_fig_suppl4 <- neurons_ant + neurons_vent +
  plot_layout(design = layout, heights = c(
    3, 1.66
  )) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(
    size = 12, face = "plain"
  ))


ggsave("Figures/Figure3_fig_suppl4.png",
  limitsize = FALSE,
  units = c("px"), Figure3_fig_suppl4, 
  width = 4000, height = 3200, bg = "white"
)

ggsave("Figures/Figure3_fig_suppl4.pdf",
  limitsize = FALSE,
  units = c("px"), Figure3_fig_suppl4, width = 4000, height = 1800
)
