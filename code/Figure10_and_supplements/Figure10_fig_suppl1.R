# R code to generate Figure 10 fig suppl of the 3d Platynereis connectome paper
# Gaspar Jekely 2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")


#import data of radial density of input and iutput synapses ------------------------
#generated in CATMAID with the morphology plot widget with a radius of 1000 nm and basis spline interpolation, with the root node (soma) as centre
desc_RD_in <- read.csv2("data/Radial_density_of_input_synapses_descending_head.csv", sep = ",")
desc_RD_out <- read.csv2("data/Radial_density_of_output_synapses_descending_head.csv", sep = ",")
decuss_RD_in <- read.csv2("data/Radial_density_of_input_synapses_decussating.csv", sep = ",")
decuss_RD_out <- read.csv2("data/Radial_density_of_output_synapses_decussating.csv", sep = ",")

# tidy radial density data -----------------------------------------------------------

desc_RD_in_tb <- desc_RD_in %>%
  rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Radial.density.of.input.synapses", "neuron", .x, fixed = TRUE)) %>%
  pivot_longer(matches("0"),
               names_to = c("nm"),
               values_to = "synapses"
  ) %>%
  mutate(neuron_type = "descending") %>%
  mutate(synapse_type = "post") %>%
  mutate(nm = as.integer(nm)) %>%
  group_by(nm) %>%
  mutate(mean = mean(synapses))
desc_RD_out_tb <- desc_RD_out %>%
  rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Radial.density.of.output.synapses", "neuron", .x, fixed = TRUE)) %>%
  pivot_longer(matches("0"),
               names_to = c("nm"),
               values_to = "synapses"
  ) %>%
  mutate(neuron_type = "descending") %>%
  mutate(synapse_type = "pre") %>%
  mutate(nm = as.integer(nm)) %>%
  group_by(nm) %>%
  mutate(mean = mean(synapses))
decuss_RD_in_tb <- decuss_RD_in %>%
  rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Radial.density.of.input.synapses", "neuron", .x, fixed = TRUE)) %>%
  pivot_longer(matches("0"),
               names_to = c("nm"),
               values_to = "synapses"
  ) %>%
  mutate(neuron_type = "decussating") %>%
  mutate(synapse_type = "post") %>%
  mutate(nm = as.integer(nm)) %>%
  group_by(nm) %>%
  mutate(mean = mean(synapses))
decuss_RD_out_tb <- decuss_RD_out %>%
  rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Radial.density.of.output.synapses", "neuron", .x, fixed = TRUE)) %>%
  pivot_longer(matches("0"),
               names_to = c("nm"),
               values_to = "synapses"
  ) %>%
  mutate(neuron_type = "decussating") %>%
  mutate(synapse_type = "pre") %>%
  mutate(nm = as.integer(nm)) %>%
  group_by(nm) %>%
  mutate(mean = mean(synapses))


# combine tibbles
desc_tb <- full_join(desc_RD_in_tb, desc_RD_out_tb)
decuss_tb <- full_join(decuss_RD_in_tb, decuss_RD_out_tb)
desc_decuss_in_out_tb <- full_join(desc_tb, decuss_tb)

write.table(desc_decuss_in_out_tb, "source_data/Figure10_fig_suppl1_source_data1.txt", sep = "\t")
desc_decuss_in_out_tb <- read.table("source_data/Figure10_fig_suppl1_source_data1.txt", sep = "\t")
# plot syn distirbution ----------------------------

# plot for descending
desc_in_out_plot <- desc_decuss_in_out_tb %>%
  filter(neuron_type == "descending") %>%
  ggplot() +
  geom_smooth(aes(
    x = nm / 1000, y = synapses, group = synapse_type,
    color = synapse_type
  ), show.legend = FALSE) +
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
  scale_x_reverse() +
  scale_shape_manual(values = c(1, 2)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.text = element_text(size = 10)
  ) +
  coord_flip()
desc_in_out_plot

# plot for decussating
decuss_in_out_plot <- desc_decuss_in_out_tb %>%
  filter(neuron_type == "decussating") %>%
  ggplot() +
  geom_smooth(aes(
    x = nm / 1000, y = synapses, group = synapse_type,
    color = synapse_type
  ), show.legend = FALSE) +
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
  scale_x_reverse() +
  scale_shape_manual(values = c(1, 2)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.text = element_text(size = 10)
  ) +
  coord_flip()
decuss_in_out_plot

# table of cell type names -----------------

decuss_celltype_skids <- skids_by_3annotations("decussating", "episphere", "celltype")
decuss_names <- as_tibble(catmaid_get_neuronnames(decuss_celltype_skids, pid = 11)) %>%
  mutate(name = sub("_.*$", "", value)) %>%
  select(name) %>% arrange(name) %>%
  pull() %>% unique() %>% 
  paste0(collapse = ", ")


desc_celltype_skids <- skids_by_3annotations("descending", "episphere", "celltype")
desc_names <- as_tibble(catmaid_get_neuronnames(desc_celltype_skids, pid = 11)) %>%
  mutate(name = sub("_.*$", "", value)) %>%
  select(name) %>% arrange(name) %>%
  pull() %>% unique() %>% 
  paste0(collapse = ", ")

# table of desc and decuss names
table_desc_decuss_names <- plot_ly(
  type = "table",
  columnwidth = c(1), 
  header = list(
    values = c("decussating")
  ),
  cells = list(
    values = list(c(decuss_names, "descending", desc_names)),
    fill = list(
      color = "white"  # Specify background colors
    )
  )
)
table_desc_decuss_names


saveNetwork(table_desc_decuss_names, "pictures/table_desc_decuss_names.html")
webshot::webshot(
  url = "pictures/table_desc_decuss_names.html",
  file = "pictures/table_desc_decuss_names.png",
  vwidth = 250, vheight = 490, # define the size of the browser window
  zoom = 10, cliprect = c(10, 45, 205, 400)
)

#read skids and neurons ----------------------------
decuss_skids <- skids_by_2annotations("decussating", "episphere")
decuss <- nlapply(
  read.neurons.catmaid(decuss_skids, pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

desc_skids <- skids_by_2annotations("descending", "episphere")
desc <- nlapply(
  read.neurons.catmaid(desc_skids, pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

# get and parse connectors ------------------
decuss_conn <- connectors(decuss)
presyn_decuss_conn <- decuss_conn[decuss_conn$prepost == 0, ]
postsyn_decuss_conn <- decuss_conn[decuss_conn$prepost == 1, ]
desc_conn <- connectors(desc)
presyn_desc_conn <- desc_conn[desc_conn$prepost == 0, ]
postsyn_desc_conn <- desc_conn[desc_conn$prepost == 1, ]

# plot desc and decuss ---------------
plot_background_ventral_no_ac()
plot3d(
  decuss,
  lwd = 2, alpha = 0.6,
  col = Okabe_Ito[1], add = T
)
plot3d(
  desc,
  lwd = 2, alpha = 0.5,
  col = Okabe_Ito[5], add = T
)
plot3d(scalebar_50um_ventral, lwd = 3, color = "black")
texts3d(129000, 116000, 182500, text = "50 μm", col = "black", cex = 2)
par3d(zoom = 0.53)

rgl.snapshot("pictures/desc_decuss_neurons_ventral.png")
nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
rgl.pop()
rgl.pop()
# z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/desc_decuss_neurons_frontal.png")
close3d()


# plot decussating synapses ---------------------
plot_background_ventral_no_ac()
plot3d(
  presyn_decuss_conn$x,
  presyn_decuss_conn$y,
  presyn_decuss_conn$z,
  size = 6, alpha = 0.7,
  col = Okabe_Ito[5], add = T
)
plot3d(
  postsyn_decuss_conn$x,
  postsyn_decuss_conn$y,
  postsyn_decuss_conn$z,
  size = 6, alpha = 0.7,
  col = Okabe_Ito[6], add = T
)
par3d(zoom = 0.53)

rgl.snapshot("pictures/decuss_post_pre_ventral.png")
nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
# z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/decuss_post_pre_frontal.png")
close3d()


# plot descending synapses ---------------------
plot_background_ventral_no_ac()
plot3d(
  presyn_desc_conn$x,
  presyn_desc_conn$y,
  presyn_desc_conn$z,
  size = 6, alpha = 0.7,
  col = Okabe_Ito[5], add = T
)
plot3d(
  postsyn_desc_conn$x,
  postsyn_desc_conn$y,
  postsyn_desc_conn$z,
  size = 6, alpha = 0.7,
  col = Okabe_Ito[6], add = T
)
par3d(zoom = 0.53)

rgl.snapshot("pictures/desc_post_pre_ventral.png")
nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
# z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/desc_post_pre_frontal.png")
close3d()


# plot postsyn of decuss and desc ---------------
plot_background_ventral_no_ac()
# plot only the presyn connectors
plot3d(
  postsyn_decuss_conn$x,
  postsyn_decuss_conn$y,
  postsyn_decuss_conn$z,
  size = 7, alpha = 0.8,
  col = Okabe_Ito[1], add = T
)
plot3d(
  postsyn_desc_conn$x,
  postsyn_desc_conn$y,
  postsyn_desc_conn$z,
  size = 7, alpha = 0.8,
  col = Okabe_Ito[5], add = T
)
par3d(zoom = 0.53)

rgl.snapshot("pictures/desc_decuss_postsyn_ventral.png")
nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
# z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/desc_decuss_postsyn_frontal.png")
close3d()



# plot presyn of decuss and desc ---------------
plot_background_ventral_no_ac()
# plot only the presyn connectors
plot3d(
  presyn_decuss_conn$x,
  presyn_decuss_conn$y,
  presyn_decuss_conn$z,
  size = 7, alpha = 0.8,
  col = Okabe_Ito[1], add = T
)
plot3d(
  presyn_desc_conn$x,
  presyn_desc_conn$y,
  presyn_desc_conn$z,
  size = 7, alpha = 0.8,
  col = Okabe_Ito[5], add = T
)
par3d(zoom = 0.53)

rgl.snapshot("pictures/desc_decuss_presyn_ventral.png")
nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
# z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/desc_decuss_presyn_frontal.png")
close3d()


# retrieve partners --------------------

trunk_connectome_skids <- skids_by_2annotations("torso", "connectome")
decuss_skids
desc_skids

decuss_trunk_connections <- catmaid_get_connectors_between(
  pre_skids = decuss_skids, 
  post_skids = trunk_connectome_skids, pid = 11
)

desc_trunk_connections <- catmaid_get_connectors_between(
  pre_skids = desc_skids, 
  post_skids = trunk_connectome_skids, pid = 11
)

desc_trunk_connections_skids <- desc_trunk_connections %>%
  as_tibble() %>%
  pull(post_skid) %>% unique()
decuss_trunk_connections_skids <- decuss_trunk_connections %>%
  as_tibble() %>%
  pull(post_skid) %>% unique()

decuss_trunk_partners <- nlapply(
  read.neurons.catmaid(decuss_trunk_connections_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
)
desc_trunk_partners <- nlapply(
  read.neurons.catmaid(desc_trunk_connections_skids, pid=11),
  function(x) smooth_neuron(x, sigma=6000)
)

# plot trunk partners ----------------

plot_background_ventral_no_ac()
plot3d(
  decuss_trunk_partners, soma = FALSE,
  lwd = 1, alpha = 0.6,
  col = Okabe_Ito[1], add = T
)
plot3d(
  desc_trunk_partners, soma = FALSE,
  lwd = 2, alpha = 0.6,
  col = Okabe_Ito[5], add = T
)
par3d(zoom = 0.53)
segments3d(
  x = as.vector(c(500, 145000)),
  y = as.vector(c(110000, 110000)),
  z = as.vector(c(78000, 78000))
)
rgl.snapshot("pictures/desc_decuss_trunk_partners.png")
close3d()

# VNC cross sections -------------------
# VNC cross-section plotting function
plot_background_VNC <- function(z_top, z_bot) {
  plot_background_ventral()
  clipplanes3d(0, 0, 1, z_top)
  clipplanes3d(0, 0, -1, z_bot)
  nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
  par3d(windowRect = c(0, 0, 800, 400)) # resize for frontal view
  par3d(zoom = 0.15)
  # x-axis clip
  clipplanes3d(0, 1, 0.16, -90000)
  # y-axis clip
  clipplanes3d(1, 0, 0.16, -13000)
  plot3d(outline, add = T, alpha = 0.14, col = "#E2E2E2")
}

plot_background_VNC(-73000, 78000)
plot3d(
  decuss, soma = FALSE,
  lwd = 5, alpha = 0.7,
  col = Okabe_Ito[1], add = T
)
plot3d(
  desc, soma = FALSE,
  lwd = 5, alpha = 0.6,
  col = Okabe_Ito[5], add = T
)
rgl.snapshot("pictures/neurites_desc_decuss_VNC.png")
close3d()

plot_background_VNC(-73000, 78000)
plot3d(
  decuss_trunk_partners, soma = FALSE,
  lwd = 4, alpha = 0.7,
  col = Okabe_Ito[1], add = T
)
plot3d(
  desc_trunk_partners, soma = FALSE,
  lwd = 4, alpha = 0.6,
  col = Okabe_Ito[5], add = T
)
rgl.snapshot("pictures/neurites_desc_decuss_VNC_targets.png")
close3d()

# assemble figure ------------------


panel_decuss_desc_f <- ggdraw() + draw_image(
  readPNG("pictures/desc_decuss_neurons_frontal.png")
) +
  draw_label("decussating",
    x = 0.1, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[1], fontface = "bold"
  ) +
  draw_label("descending",
    x = 0.1, y = 0.94, hjust = 0,
    size = 10, color = Okabe_Ito[5], fontface = "bold"
  ) +
  draw_label("neurites",
    x = 0.7, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[8], fontface = "plain"
  )

panel_decuss_desc_post_v <- ggdraw() + draw_image(
  readPNG("pictures/desc_decuss_postsyn_ventral.png")
) +
  draw_label("decussating",
    x = 0.1, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[1], fontface = "bold"
  ) +
  draw_label("descending",
    x = 0.1, y = 0.94, hjust = 0,
    size = 10, color = Okabe_Ito[5], fontface = "bold"
  ) +
  draw_label("postsynapse",
    x = 0.6, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[8], fontface = "plain"
  )

panel_decuss_desc_pre_v <- ggdraw() + draw_image(
  readPNG("pictures/desc_decuss_presyn_ventral.png")
) +
  draw_label("decussating",
    x = 0.1, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[1], fontface = "bold"
  ) +
  draw_label("descending",
    x = 0.1, y = 0.94, hjust = 0,
    size = 10, color = Okabe_Ito[5], fontface = "bold"
  ) +
  draw_label("presynapse",
    x = 0.6, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[8], fontface = "plain"
  )


panel_decuss_post_pre_v <- ggdraw() + draw_image(
  readPNG("pictures/decuss_post_pre_ventral.png")
)

panel_decuss_post_pre_f <- ggdraw() + draw_image(
  readPNG("pictures/decuss_post_pre_frontal.png")
) +
  draw_label("decussating",
    x = 0.1, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[8], fontface = "plain"
  ) +
  draw_label("presynapse",
    x = 0.6, y = 0.94, hjust = 0,
    size = 10, color = Okabe_Ito[5], fontface = "bold"
  ) +
  draw_label("postsynapse",
    x = 0.6, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[6], fontface = "bold"
  )

panel_desc_post_pre_v <- ggdraw() + draw_image(
  readPNG("pictures/desc_post_pre_ventral.png")
)

panel_desc_post_pre_f <- ggdraw() + draw_image(
  readPNG("pictures/desc_post_pre_frontal.png")
) +
  draw_label("descending",
    x = 0.1, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[8], fontface = "plain"
  ) +
  draw_label("presynapse",
    x = 0.6, y = 0.94, hjust = 0,
    size = 10, color = Okabe_Ito[5], fontface = "bold"
  ) +
  draw_label("postsynapse",
    x = 0.6, y = 0.99, hjust = 0,
    size = 10, color = Okabe_Ito[6], fontface = "bold"
  )

panel_decuss_desc_v <- ggdraw() + draw_image(
  readPNG("pictures/desc_decuss_neurons_ventral.png")
)  +
  draw_label("decussating",
             x = 0.1, y = 0.99, hjust = 0,
             size = 10, color = Okabe_Ito[1], fontface = "bold"
  ) +
  draw_label("descending",
             x = 0.1, y = 0.94, hjust = 0,
             size = 10, color = Okabe_Ito[5], fontface = "bold"
  ) +
  draw_label("neurites",
             x = 0.7, y = 0.99, hjust = 0,
             size = 10, color = Okabe_Ito[8], fontface = "plain"
  ) 


panel_desc_decuss_partners <- ggdraw() + draw_image(readPNG(
  "pictures/desc_decuss_trunk_partners.png"
  )
  ) +
  draw_label("decussating",
             x = 0.1, y = 0.99, hjust = 0,
             size = 10, color = Okabe_Ito[1], fontface = "bold"
  ) +
  draw_label("descending",
             x = 0.1, y = 0.94, hjust = 0,
             size = 10, color = Okabe_Ito[5], fontface = "bold"
  ) +
  draw_label("targets",
             x = 0.7, y = 0.99, hjust = 0,
             size = 10, color = Okabe_Ito[8], fontface = "plain"
  )

panel_desc_decuss_cross <- ggdraw() + draw_image(readPNG(
  "pictures/neurites_desc_decuss_VNC.png"
)
) +
  draw_label("decussating",
             x = 0.1, y = 0.99, hjust = 0,
             size = 10, color = Okabe_Ito[1], fontface = "bold"
  ) +
  draw_label("descending",
             x = 0.1, y = 0.89, hjust = 0,
             size = 10, color = Okabe_Ito[5], fontface = "bold"
  )

panel_desc_decuss_partners_cross <- ggdraw() + draw_image(readPNG(
  "pictures/neurites_desc_decuss_VNC_targets.png"
)
) +
  draw_label("decussating targets",
             x = 0.1, y = 0.99, hjust = 0,
             size = 10, color = Okabe_Ito[1], fontface = "bold"
  ) +
  draw_label("descending targets",
             x = 0.1, y = 0.89, hjust = 0,
             size = 10, color = Okabe_Ito[5], fontface = "bold"
  )


panel_table <- ggdraw() + draw_image(readPNG("pictures/table_desc_decuss_names.png"))

panel_in_out <- plot_grid(
  panel_decuss_post_pre_f, panel_decuss_post_pre_v, decuss_in_out_plot, 
  panel_desc_post_pre_f, panel_desc_post_pre_v, desc_in_out_plot, 
  rel_widths = c(1, 1.3, 0.9, 1, 1.3, 0.9), ncol = 6, rel_heights = c(1, 1, 1, 1, 1, 1),
  labels = c("H", "I", "J", "K", "L", "M"),
  label_size = 12, label_y = 1.03, label_x = 0,
  label_fontfamily = "sans", label_fontface = "plain"
  )

layout <- "
ABCDEF
ABCDEG
HHHHHH
HHHHHH
"

Figure10_fig_suppl1 <- panel_decuss_desc_v + panel_table +
  panel_decuss_desc_post_v +
  panel_decuss_desc_pre_v + 
  panel_desc_decuss_partners + panel_desc_decuss_cross + panel_desc_decuss_partners_cross +
  panel_in_out +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "G"))) &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure10_fig_suppl1.png",
  limitsize = FALSE,
  units = c("px"), Figure10_fig_suppl1,
  width = 4800, height = 2200, bg = "white"
)


ggsave("Figures/Figure10_fig_suppl1.pdf",
  limitsize = FALSE,
  units = c("px"), Figure10_fig_suppl1,
  width = 4800, height = 2200, bg = "black"
)


