# R code to generate Fig 1 fig suppl anatomy of the 3d Platynereis connectome paper
# Gaspar Jekely 2024

# load natverse and other packages, some custom natverse functions and catmaid connectivity info -----------
source("code/Natverse_functions_and_conn.R")

# load skeletons ----------------

Ciliary_band_cell <- nlapply(
  read.neurons.catmaid("^Ciliary_band_cell$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
stomodeum <- nlapply(
  read.neurons.catmaid("^stomodeum$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

pyg_SN_skids <- skids_by_2annotations("pygidium", "Sensory neuron")
pyg_SN <- nlapply(
  read.neurons.catmaid(pyg_SN_skids, pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

# plot overview ----------------

plot_background_ventral2()
par3d(windowRect = c(0, 0, 800, 800))
par3d(zoom=0.65)

plot3d(
  acicula, soma = TRUE, lwd = 2, add = TRUE, 
  alpha = 1, col = "black"
)
plot3d(
  chaeta, soma = TRUE, lwd = 3, add = TRUE, 
  alpha = 0.7, col = "#cccccc"
)
plot3d(
  stomodeum, soma = TRUE, lwd = 1, add = TRUE, 
  alpha = 0.5, col = sample(bluepurple[1:9], 676, replace = TRUE)
)
plot3d(
  pyg_SN, soma = TRUE, lwd = 4, add = TRUE, 
  alpha = 0.2, col = "#0072B2"
)
plot3d(Ciliary_band_cell,
       soma = TRUE, lwd = 1,
       add = TRUE, alpha = 0.6,
       col = "#E69F00"
)
plot3d(
  scalebar_50um_ventral, lwd = 2, add = TRUE, 
  alpha = 1, col = "black"
)
plot3d(
  outline, lwd = 2, add = TRUE, 
  alpha = 0.02, col = "black"
)

# make snapshot
rgl.snapshot("pictures/Naomi_overview.png")
close3d()

# assemble figure ---------------

imgSEM <- readPNG("images_notR/Platynereis_SEM_ventral_d9628_291um.png")
img_larva <- readPNG("images_notR/Platynereis_3d_LM_360um.png")
imgNaomi <- readPNG("pictures/Naomi_overview.png")

panelSEM <- ggdraw() + draw_image(imgSEM) +
  draw_label("head", x = 0.5, y = 0.87, size = 12, color = "white") +
  draw_label("peristomium", x = 0.52, y = 0.74, size = 12, color = "white") +
  draw_label("sg0", x = 0.53, y = 0.7, size = 12, color = "white") +
  draw_label("sg1", x = 0.53, y = 0.63, size = 12, color = "white") +
  draw_label("sg2", x = 0.54, y = 0.54, size = 12, color = "white") +
  draw_label("sg3", x = 0.55, y = 0.38, size = 12, color = "white") +
  draw_label("pygidium", x = 0.56, y = 0.24, size = 12, color = "white") +
  draw_label("chaetae", x = 0.1, y = 0.63, size = 12, color = "white") +
  draw_label("ciliary \nband", x = 0.77, y = 0.8, size = 12, color = "white") +
  draw_label("Platynereis", x = 0.3, y = 0.99, size = 10, fontface = "italic", color = "white") + 
  draw_label("larva", x = 0.55, y = 0.99, size = 10, fontface = "plain", color = "white") +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
             x = 0.18, y = 0.08, size = 10, color = "white") +
  draw_line(
    x = c(0.1, 0.271),
    y = c(0.06, 0.06), size = 0.3, color = "white") 

panelLarva <- ggdraw() + draw_image(img_larva) +
  draw_label("head", x = 0.48, y = 0.87, size = 12, color = "white") +
  draw_label("peristomium", x = 0.5, y = 0.71, size = 12, color = "white") +
  draw_label("sg1", x = 0.51, y = 0.63, size = 12, color = "white") +
  draw_label("sg2", x = 0.52, y = 0.54, size = 12, color = "white") +
  draw_label("sg3", x = 0.52, y = 0.36, size = 12, color = "white") +
  draw_label("pygidium", x = 0.52, y = 0.22, size = 12, color = "white") +
  draw_label("chaetae", x = 0.1, y = 0.7, size = 12, color = "white") +
  draw_label("ciliary \nband", x = 0.71, y = 0.83, size = 12, color = "white") +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
             x = 0.165, y = 0.08, size = 10, color = "white") +
  draw_line(
    x = c(0.1, 0.238),
    y = c(0.06, 0.06), size = 0.3, color = "white") 

panelNaomi <- ggdraw() + draw_image(imgNaomi) +
  draw_label("head", x = 0.52, y = 0.93, size = 12) +
  draw_label("stomodeum", x = 0.78, y = 0.91, size = 12) +
  draw_line(
    x = c(0.54, 0.75),
    y = c(0.77, 0.89), size = 0.3) +
  draw_label("sg1", x = 0.56, y = 0.58, size = 12) +
  draw_label("sg2", x = 0.55, y = 0.46, size = 12) +
  draw_label("sg3", x = 0.5, y = 0.32, size = 12) +
  draw_label("pygidium", x = 0.4, y = 0.07, size = 12) +
  draw_label("chaetae", x = 0.17, y = 0.43, size = 12) +
  draw_label("ciliary \nband", x = 0.89, y = 0.76, size = 12) +
  draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
             x = 0.8, y = 0.1, size = 10) +
  draw_label("aciculae", x = 0.83, y = 0.26, size = 12) +
  draw_line(
    x = c(0.58, 0.74, 0.76),
    y = c(0.15, 0.25, 0.41), size = 0.3) +
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
  
layout <-
  "A#B#C"

Fig1_fig_suppl1 <- panelSEM + panelLarva + panelNaomi +
  plot_layout(design = layout, widths = c(1, 0.02, 1, 0.02, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure1_fig_suppl1.png",
       limitsize = FALSE,
       units = c("px"), Fig1_fig_suppl1, width = 4000, height = 1420, bg = "white"
)

ggsave("Figures/Figure1_fig_suppl1.pdf",
       limitsize = FALSE,
       units = c("px"), Fig1_fig_suppl1, width = 4000, height = 1420, bg = "white"
)
