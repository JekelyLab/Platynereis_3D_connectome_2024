# R/natverse code to generate Figure 16 for the Platynereis 3d connectome paper
# Gaspar Jekely March 2022-2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# load cell clusters ------------------------------------------------------

{
  CR <- nlapply(
    read.neurons.catmaid("^CRneurons$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  PU <- nlapply(
    read.neurons.catmaid("^PUneurons$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  PB <- nlapply(
    read.neurons.catmaid("^Biciliated_penetrating_cell$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INsplit <- nlapply(
    read.neurons.catmaid("^INsplit$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MUSlongV <- nlapply(
    read.neurons.catmaid("^celltype_non_neuronal76$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  MUStrans <- nlapply(
    read.neurons.catmaid("^celltype_non_neuronal74$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  stomodeum <- nlapply(
    read.neurons.catmaid("^stomodeum$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  girdle <- nlapply(
    read.neurons.catmaid("^mechanosensory_girdle$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
}

segmental_colors <- brewer.pal(6, "Paired")
pie(rep(1, 6), col = segmental_colors, segmental_colors)

# adjusted background plotting function
plot_background_mech <- function() {
  plot_background_ventral()
  clear3d()
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1)))
  plot3d(yolk,
    WithConnectors = F, WithNodes = F, soma = TRUE, lwd = 0,
    add = T, forceClipregion = F, alpha = 0.05,
    col = "grey80"
  )
  plot3d(stomodeum,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 3,
    add = T, alpha = 0.1, col = "grey50"
  )
  par3d(zoom = 0.58)
}

body_regions <- c("episphere", "segment_0", "segment_1", "segment_2", "segment_3", "pygidium")
body_region_names <- c("head", "sg0", "sg1", "sg2", "sg3", "pyg")
# matrix to count cell numbers per segment
N_cells_seg <- matrix(rep(0, 36), nrow = 6, ncol = 6)
rownames(N_cells_seg) <- c("head", "sg0", "sg1", "sg2", "sg3", "pyg")
colnames(N_cells_seg) <- c("CR", "PU", "PB", "INsplit", "MUSlongV", "MUStrans")

# plot CR neurons by segment
{
  plot_background_mech()

  for (i in 1:6) {
    skids <- skids_by_2annotations("CRneurons", as.character(body_regions[i]))
    print(skids)
    N_cells_seg[i, 1] <- length(skids)
    neurons <- nlapply(
      read.neurons.catmaid(skids, pid = 11),
      function(x) smooth_neuron(x, sigma = 6000)
    )
    lapply(neurons, function(x) {
      plot3d(x,
        WithConnectors = F, WithNodes = F, soma = T, lwd = 3,
        add = T, alpha = 1,
        col = segmental_colors[i]
      )
    })
  }

  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(acicula,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.3, col = "grey50"
  )
  plot3d(scalebar_50um_ventral,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 3,
    add = T, forceClipregion = F, alpha = 1,
    col = "black"
  )

  rgl.snapshot("pictures/Figure_mec_seg_CR.png")
  close3d()
}

# plot PU neurons by segment
{
  plot_background_mech()

  for (i in 1:6) {
    skids <- skids_by_2annotations("PUneurons", as.character(body_regions[i]))
    print(skids)
    N_cells_seg[i, 2] <- length(skids)
    neurons <- nlapply(
      read.neurons.catmaid(skids, pid = 11),
      function(x) smooth_neuron(x, sigma = 6000)
    )
    lapply(neurons, function(x) {
      plot3d(x,
        WithConnectors = F, WithNodes = F, soma = T, lwd = 3,
        add = T, alpha = 1,
        col = segmental_colors[i]
      )
    })
  }
  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(acicula,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.3, col = "grey50"
  )

  rgl.snapshot("pictures/Figure_mec_seg_PU.png")
  close3d()
}

# plot PB neurons by segment
{
  plot_background_mech()

  for (i in c(1, 3, 4, 5, 6)) {
    skids <- skids_by_2annotations("Biciliated_penetrating_cell", as.character(body_regions[i]))
    print(skids)
    N_cells_seg[i, 3] <- length(skids)
    neurons <- nlapply(
      read.neurons.catmaid(skids, pid = 11),
      function(x) smooth_neuron(x, sigma = 6000)
    )
    lapply(neurons, function(x) {
      plot3d(x,
        WithConnectors = F, WithNodes = F, soma = T, lwd = 3,
        add = T, alpha = 1,
        col = segmental_colors[i]
      )
    })
  }
  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(acicula,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.3, col = "grey50"
  )

  rgl.snapshot("pictures/Figure_mec_seg_PB.png")
  close3d()
}

# plot INsplit by segment
{
  plot_background_mech()

  for (i in 1:5) {
    skids <- skids_by_2annotations("INsplit", as.character(body_regions[i]))
    print(skids)
    N_cells_seg[i, 4] <- length(skids)
    neurons <- nlapply(
      read.neurons.catmaid(skids, pid = 11),
      function(x) smooth_neuron(x, sigma = 6000)
    )
    lapply(neurons, function(x) {
      plot3d(x,
        WithConnectors = F, WithNodes = F, soma = T, lwd = 3,
        add = T, alpha = 1,
        col = segmental_colors[i]
      )
    })
  }
  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(acicula,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.3, col = "grey50"
  )

  rgl.snapshot("pictures/Figure_mec_seg_INsplit.png")
  close3d()
}

# plot MUSlongV by segment
{
  plot_background_mech()

  for (i in 1:6) {
    skids <- skids_by_2annotations("celltype_non_neuronal76", as.character(body_regions[i]))
    print(skids)
    N_cells_seg[i, 5] <- length(skids)
    neurons <- nlapply(
      read.neurons.catmaid(skids, pid = 11),
      function(x) smooth_neuron(x, sigma = 6000)
    )
    lapply(neurons, function(x) {
      plot3d(x,
        WithConnectors = F, WithNodes = F, soma = T, lwd = 3,
        add = T, alpha = 1,
        col = segmental_colors[i]
      )
    })
  }
  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(acicula,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.3, col = "grey50"
  )


  rgl.snapshot("pictures/Figure_mec_seg_MUSlongV.png")
  close3d()
}

# plot MUStrans by segment
{
  plot_background_mech()

  for (i in 2:6) {
    skids <- skids_by_2annotations("celltype_non_neuronal74", as.character(body_regions[i]))
    print(skids)
    N_cells_seg[i, 6] <- length(skids)
    neurons <- nlapply(
      read.neurons.catmaid(skids, pid = 11),
      function(x) smooth_neuron(x, sigma = 6000)
    )
    lapply(neurons, function(x) {
      plot3d(x,
        WithConnectors = F, WithNodes = F, soma = T, lwd = 3,
        add = T, alpha = 1,
        col = segmental_colors[i]
      )
    })
  }

  plot3d(girdle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 0.15, col = "grey50"
  )
  plot3d(acicula,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.3, col = "grey50"
  )

  rgl.snapshot("pictures/Figure_mec_seg_MUStrans.png")
  close3d()
}

# plot all cell groups
{
  plot_background_mech()

  plot3d(CR,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.5, col = segmental_colors[1]
  )
  plot3d(PU,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.5, col = segmental_colors[2]
  )
  plot3d(PB,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.5, col = segmental_colors[3]
  )
  plot3d(INsplit,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.5, col = segmental_colors[4]
  )
  plot3d(MUSlongV,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 3,
    add = T, alpha = 0.5, col = segmental_colors[6]
  )
  plot3d(MUStrans,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 3,
    add = T, alpha = 0.8, col = segmental_colors[5]
  )


  rgl.snapshot("pictures/Figure_mec_seg_all.png")
  nview3d("ventral", extramat = (rotationMatrix(0.35, 1, 0, 0) %*% rotationMatrix(0.05, 0, 0, 1) %*% rotationMatrix(-pi / 2, 0, 0, 1)))
  rgl.snapshot("pictures/Figure_mec_seg_all_lat.png")

  close3d()
}

# plot number of cells per segment ---------------------------------------

{
  as.data.frame((N_cells_seg)) %>%
    rownames_to_column(var = "cell_group") %>%
    pivot_longer(-cell_group, names_to = "body_region", values_to = "number_of_cells") %>%
    group_by(body_region) %>%
    ggplot(aes(x = body_region, y = cell_group)) +
    geom_raster(aes(fill = sqrt(number_of_cells))) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      panel.grid = element_blank()
    ) +
    labs(x = "cell class", y = "segment", title = " ") +
    scale_x_discrete(limits = c("CR", "PU", "PB", "INsplit", "MUSlongV", "MUStrans")) +
    scale_y_discrete(limits = rev(body_region_names)) +
    scale_fill_gradientn(colours = c("white", "#0072B2")) +
    geom_text(aes(label = number_of_cells)) +
    guides(size = "none", fill = "none")

  # Saving R ggplot with R ggsave Function
  ggsave("pictures/mech_girdle_seg_homol_number_of_cells.png",
    width = 1100,
    height = 1300, limitsize = TRUE,
    units = c("px")
  )
}



# assmeble figure -------------------------------------------------------------

{
  imgCR <- readPNG("pictures/Figure_mec_seg_CR.png")
  imgPU <- readPNG("pictures/Figure_mec_seg_PU.png")
  imgPB <- readPNG("pictures/Figure_mec_seg_PB.png")
  imgINsplit <- readPNG("pictures/Figure_mec_seg_INsplit.png")
  imgMUSlongV <- readPNG("pictures/Figure_mec_seg_MUSlongV.png")
  imgMUStrans <- readPNG("pictures/Figure_mec_seg_MUStrans.png")
  imgAll <- readPNG("pictures/Figure_mec_seg_all.png")
  imgAll_lat <- readPNG("pictures/Figure_mec_seg_all_lat.png")
  imgTable <- readPNG("pictures/mech_girdle_seg_homol_number_of_cells.png")

  panelCR <- ggdraw() + draw_image(imgCR, scale = 1) +
    draw_label("CR",
      x = 0.2, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    ) +
    draw_label(expression(paste("50 ", mu, " m")),
      x = 0.78, y = 0.07, fontfamily = "sans", fontface = "plain",
      color = "black", size = 10
    )

  panelPU <- ggdraw() + draw_image(imgPU, scale = 1) +
    draw_label("PU",
      x = 0.2, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    )
  panelPB <- ggdraw() + draw_image(imgPB, scale = 1) +
    draw_label("PB",
      x = 0.2, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    )
  panelNsplit <- ggdraw() + draw_image(imgINsplit, scale = 1) +
    draw_label("INsplit all types",
      x = 0.42, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    )
  panelMUSlongV <- ggdraw() + draw_image(imgMUSlongV, scale = 1) +
    draw_label("MUSlongV",
      x = 0.3, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    )
  panelMUStrans <- ggdraw() + draw_image(imgMUStrans, scale = 1) +
    draw_label("MUStrans",
      x = 0.3, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    )

  panelAll <- ggdraw() + draw_image(imgAll, scale = 1) +
    draw_label("all",
      x = 0.4, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    )

  panelTable <- ggdraw() + draw_image(imgTable) +
    draw_label("cells per segment",
      x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 11
    )


layout <- "
ABCD
####
EFGH
"

Figure_mech_seg_homol <- panelCR + panelPU + panelPB +
    panelNsplit + panelMUSlongV + panelMUStrans +
    panelAll + panelTable +
    plot_layout(design = layout, heights = c(1, 0.05, 1, 0.05, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12, face = "plain"))

  ggsave(
    "Figures/Figure13.png",
    limitsize = FALSE,
    units = c("px"), Figure_mech_seg_homol,
    width = 2400, height = 1640, bg = "white"
  )

  ggsave(
    "Figures/Figure13.pdf",
    limitsize = FALSE,
    units = c("px"), Figure_mech_seg_homol,
    width = 2400, height = 1640
  )
}
