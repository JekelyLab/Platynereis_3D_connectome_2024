# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/libraries_functions_and_CATMAID_conn.R")

# functions ------------------
load_neuron <- function(annotation) {
  nlapply(
    read.neurons.catmaid(
      annotation,
      pid = 11
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
}

# define background plotting function
plot_background <- function(x) {
  nopen3d() # opens a pannable 3d window
  plot3d(bounding_dots, lwd = 1, add = T, alpha = 1, col = "white")
  # we define a z clipping plane for the frontal view
  par3d(zoom = 0.42)
  nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
  # z-axis clip
  clipplanes3d(0, 0, -1, 165000)
  # y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  # x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
  par3d(windowRect = c(0, 0, 800, 800)) # resize for frontal view
}

# read background skeletons -----------------

yolk <- catmaid_get_volume(4,
  rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 11
)
bounding_dots <- nlapply(
  read.neurons.catmaid("^bounding_dots$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
scalebar_50um_anterior <- read.neurons.catmaid("^scalebar_50um_anterior$", pid = 11)
scalebar_50um_ventral = read.neurons.catmaid("^scalebar_50um_ventral$", pid=11)

acicula <- nlapply(
  read.neurons.catmaid("^acicula$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

# read data ------------------------
celltypes_table <- read_csv("supplements/Supplementary_Table1.csv")
celltypes_RD_in <- read.csv2("data/all_celltypes_radial_density_in.csv", sep = ",")
celltypes_RD_out <- read.csv2("data/all_celltypes_radial_density_out.csv", sep = ",")


celltypes_RD_in_tb <- celltypes_RD_in %>%
  rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Radial.density.of.input.synapses", "neuron", .x, fixed = TRUE)) %>%
  pivot_longer(matches("0"),
    names_to = c("nm"),
    values_to = "synapses"
  ) %>%
  mutate(synapse_type = "post") %>%
  mutate(nm = as.integer(nm)) %>%
  mutate(celltype_name = gsub("^(.*?)_.*$", "\\1", neuron)) %>%
  group_by(nm, celltype_name) %>%
  mutate(mean = mean(synapses)) %>%
  select(nm, synapse_type, synapses, celltype_name, mean)

celltypes_RD_in_tb

celltypes_RD_out_tb <- celltypes_RD_out %>%
  rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Radial.density.of.output.synapses", "neuron", .x, fixed = TRUE)) %>%
  pivot_longer(matches("0"),
    names_to = c("nm"),
    values_to = "synapses"
  ) %>%
  mutate(synapse_type = "pre") %>%
  mutate(nm = as.integer(nm)) %>%
  mutate(celltype_name = gsub("^(.*?)_.*$", "\\1", neuron)) %>%
  group_by(nm, celltype_name) %>%
  mutate(mean = mean(synapses)) %>%
  select(nm, synapse_type, synapses, celltype_name, mean)

celltypes_in_out_tb <- full_join(celltypes_RD_in_tb, celltypes_RD_out_tb)

# save source data --------------
write.table(celltypes_in_out_tb, "source_data/Figure3_fig_suppl7_source_data1.txt", sep = "\t")
celltypes_in_out_tb <- read.table("source_data/Figure3_fig_suppl7_source_data1.txt", sep = "\t")

# plot for all cell types ---------------

for (i in 47:47) {
  Catmaid_annot <- paste("celltype", i, sep = "")
  neuron_name <- celltypes_table %>%
    filter(`CATMAID annotation` == Catmaid_annot) %>%
    pull(`name of cell type`) %>%
    unique()

  n_pre <- celltypes_in_out_tb %>%
    filter(celltype_name == neuron_name) %>%
    group_by(synapse_type) %>%
    mutate(total_syn = sum(synapses)) %>%
    filter(synapse_type == "pre") %>%
    pull(total_syn) %>%
    unique()
  n_post <- celltypes_in_out_tb %>%
    filter(celltype_name == neuron_name) %>%
    group_by(synapse_type) %>%
    mutate(total_syn = sum(synapses)) %>%
    filter(synapse_type == "post") %>%
    pull(total_syn) %>%
    unique()

  if (n_pre < 12 | n_post < 12) {
    print("too few synapes")
    next()
  }

  # define scale of plot
  row_number_end_of_scale <- celltypes_in_out_tb %>%
    filter(celltype_name == neuron_name) %>%
    group_by(celltype_name) %>%
    mutate(last_synapse = max(which(mean != 0))) %>%
    pull(last_synapse) %>%
    unique()

  end_of_scale <- celltypes_in_out_tb %>%
    filter(celltype_name == neuron_name) %>%
    group_by(celltype_name) %>%
    slice(row_number_end_of_scale) %>%
    pull(nm)

  # scale to 50, 100, 150 um
  end_of_scale <- ifelse(end_of_scale < 50000, 50000,
    ifelse(end_of_scale < 100000, 100000,
      ifelse(end_of_scale < 150000, 150000, end_of_scale)
    )
  )

  celltypes_in_out_plot <- celltypes_in_out_tb %>%
    filter(celltype_name == neuron_name) %>%
    ggplot() +
    geom_smooth(aes(
      x = nm / 1000, y = synapses, group = synapse_type,
      color = synapse_type
    ), show.legend = FALSE, span = 0.1, level = 0.95) +
    geom_point(aes(
      x = nm / 1000, y = mean, shape = synapse_type,
      color = synapse_type
    ), show.legend = FALSE, size = 1) +
    theme(panel.background = element_rect(fill = "grey95", color = "grey")) +
    labs(x = "distance from soma (µm)", y = "radial synapse density") +
    scale_color_manual(values = c(
      Okabe_Ito[6],
      Okabe_Ito[5]
    )) +
    scale_x_reverse(limits = c(end_of_scale / 1000, 0)) +
    scale_shape_manual(values = c(1, 2)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    coord_flip()

  neuron_skeleton <- load_neuron(Catmaid_annot)
  synapses <- connectors(neuron_skeleton)
  presyn <- synapses[synapses$prepost == 0, ]
  postsyn <- synapses[synapses$prepost == 1, ]

  segment <- celltypes_table %>%
    filter(`CATMAID annotation` == Catmaid_annot) %>%
    pull(`soma position`)

  if (segment == "head") {
    plot_background()
    plot3d(neuron_skeleton,
      soma = TRUE, lwd = 4,
      add = T, alpha = 0.7, col = Okabe_Ito[1]
    )
    plot3d(presyn$x, presyn$y, presyn$z,
      size = 6,
      add = T, alpha = 1, col = Okabe_Ito[5]
    )
    plot3d(postsyn$x, postsyn$y, postsyn$z,
      size = 8,
      add = T, alpha = 1, col = Okabe_Ito[6]
    )
    plot3d(yolk, add = T, alpha = 0.05, col = "#E2E2E2")
    texts3d(32000, 10000, 2000, text = neuron_name, col = "black", cex = 3)
    if(neuron_name == "PRC"){
      plot3d(
        scalebar_50um_anterior, lwd = 4,
        add = T, alpha = 1, col = "black"
      )
      texts3d(104000, 116000, 300, text = "50 μm", col = "black", cex = 2)
    }
  } else {
    plot_background_ventral()
    plot3d(neuron_skeleton,
      soma = TRUE, lwd = 4,
      add = T, alpha = 0.7, col = Okabe_Ito[1]
    )
    plot3d(presyn$x, presyn$y, presyn$z,
      size = 6,
      add = T, alpha = 1, col = Okabe_Ito[5]
    )
    plot3d(postsyn$x, postsyn$y, postsyn$z,
      size = 8,
      add = T, alpha = 1, col = Okabe_Ito[6]
    )
    texts3d(32000, 20000, 34000, text = neuron_name, col = "black", cex = 3 / 8 * 6)
    }
  

  rgl.snapshot("pictures/celltypes_in_out_plot.png")
  close3d()

  # make plot of neuron and synapse plot

  panel_neurons <- ggdraw() + draw_image(readPNG("pictures/celltypes_in_out_plot.png"))

  panel_neuron_plot <- plot_grid(panel_neurons, celltypes_in_out_plot,
    rel_widths = c(1, 0.7)
  )

  file_name <- paste("pictures/Fig1_suppl2_panel_", neuron_name, ".png", sep = "")
  ggsave(file_name,
    units = c("px"), panel_neuron_plot,
    width = 1800, height = 1200, bg = "white"
  )
}


# assemble figure ------------------

panel_MC <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MC.png"
))
panel_MNant <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MNant.png"
))
panel_SNlasso <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_SNlasso.png"
))
panel_MNspinning <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MNspinning.png"
))
panel_INarc1 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INarc1.png"
))
panel_INarc2 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INarc2.png"
))
panel_INrope <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INrope.png"
))
panel_Loop <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_Loop.png"
))
panel_INCM <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INCM.png"
))
panel_MNspider_a <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MNspider-ant.png"
))
panel_MNspider_p <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MNspider-post.png"
))
panel_MNbiramous <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MNbiramous.png"
))
panel_chaeMech <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_chaeMech.png"
))
panel_INsplitCR <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INsplitCR.png"
))
panel_INsplitPB <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INsplitPB.png"
))
panel_cioMNcover <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_cioMNcover.png"
))
panel_MN1 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MN1.png"
))
panel_MN2 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MN2.png"
))
panel_MN3 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MN3.png"
))
panel_MC3cover <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MC3cover.png"
))
panel_INW <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INW.png"
))
panel_INproT2 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INproT2.png"
))
panel_INZ <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INZ.png"
))
panel_Ser_tr1 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_Ser-tr1.png"
))
panel_SNFVa <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_SNFVa.png"
))
panel_INUturn <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INUturn.png"
))
panel_INasc_pyg <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INasc-pyg.png"
))
panel_INleucoPU <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INleucoPU.png"
))
panel_INdecussfoot <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INdecussfoot.png"
))
panel_INsplitVent <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INsplitVent.png"
))
panel_MNgland_head <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MNgland-head.png"
))
panel_SNpygM <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_SNpygM.png"
))
panel_INsplitBronto <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INsplitBronto.png"
))
panel_MNladder <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_MNladder.png"
))
panel_INdecusshook <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INdecusshook.png"
))
panel_INlasso_postSN <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INlasso-postSN.png"
))
panel_PRC <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_PRC.png"
))
panel_IN1 <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_IN1.png"
))
panel_INton <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INton.png"
))
panel_INcomm_Upcross <- ggdraw() + draw_image(readPNG(
  "pictures/Fig1_suppl2_panel_INcomm-Upcross.png"
))

layout_suppl2 <- "
ABCD
EFGH
IJKL
MNOP
QRST
"


Figure3_fig_suppl7 <-
  panel_PRC + panel_IN1 + panel_INton + panel_MC +
    panel_MNant + panel_SNlasso + panel_INarc1 + panel_INarc2 +
    panel_INrope + panel_INW + panel_INproT2 + panel_INlasso_postSN +
    panel_INZ + panel_INUturn + panel_INdecussfoot + panel_INdecusshook +
    panel_MNgland_head + panel_MN1 + panel_MN2 + panel_MN3 +
    plot_layout(design = layout_suppl2) +
    plot_annotation(tag_levels = ("A")) &
    theme(plot.tag = element_text(size = 14, face = "plain"))

ggsave("Figures/Figure3_fig_suppl7.png",
  limitsize = FALSE,
  units = c("px"), Figure3_fig_suppl7,
  width = 7200, height = 6000, bg = "white"
)


ggsave("Figures/Figure3_fig_suppl7.pdf",
  limitsize = FALSE,
  units = c("px"), Figure3_fig_suppl7,
  width = 7200, height = 6000, bg = "black"
)


layout_suppl3 <- "
ABCD
EFGH
IJKL
MNOP
QRST
"

Figure3_fig_suppl8 <-
  panel_MNspinning + panel_Loop + panel_MNspider_a + panel_MNspider_p +
    panel_MC3cover + panel_Ser_tr1 + panel_MNbiramous + panel_chaeMech +
    panel_INsplitCR + panel_INsplitPB + panel_INCM + panel_INsplitVent +
    panel_INsplitBronto + panel_cioMNcover + panel_SNFVa + panel_INasc_pyg +
    panel_INleucoPU + panel_SNpygM + panel_MNladder + panel_INcomm_Upcross +
    plot_layout(design = layout_suppl3) +
    plot_annotation(tag_levels = ("A")) &
    theme(plot.tag = element_text(size = 14, face = "plain"))

ggsave("Figures/Figure3_fig_suppl8.png",
  limitsize = FALSE,
  units = c("px"), Figure3_fig_suppl8,
  width = 7200, height = 6000, bg = "white"
)

ggsave("Figures/Figure3_fig_suppl8.pdf",
  limitsize = FALSE,
  units = c("px"), Figure3_fig_suppl8,
  width = 7200, height = 6000, bg = "black"
)
