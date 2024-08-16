# R code to generate the anatomical overview images of the Platynereis connectome in Figure 1 of the connectome paper
# Gaspar Jekely 2021-2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info -----------
source("code/Natverse_functions_and_conn.R")

# read graph
conn.tb <- readRDS("supplements/connectome_graph_tibble.rds")

# select subsets
connectome_neuron <- conn.tb %>%
  filter(group == "Sensory neuron" | group == "interneuron" | group == "motorneuron") %>%
  select(skids) %>%
  pull()

writeLines(connectome_neuron, "data/connectome_neuron_skids.csv")

connectome <- conn.tb %>%
  select(skids) %>%
  pull()

writeLines(connectome, "data/connectome_skids.csv")

connectome_no_frag <- conn.tb %>%
  filter(group != "fragmentum") %>%
  select(skids) %>%
  pull()

connectome_frag <- conn.tb %>%
  filter(group == "fragmentum") %>%
  select(skids) %>%
  pull()
writeLines(connectome_frag, "data/connectome_frag.csv")

Sensory_neuron <- conn.tb %>%
  filter(group == "Sensory neuron") %>%
  select(skids) %>%
  pull()

writeLines(Sensory_neuron, "data/connectome_SN_skids.csv")

motorneuron <- conn.tb %>%
  filter(group == "motorneuron") %>%
  select(skids) %>%
  pull()

writeLines(motorneuron, "data/connectome_MN_skids.csv")

interneuron <- conn.tb %>%
  filter(group == "interneuron") %>%
  select(skids) %>%
  pull()

writeLines(interneuron, "data/connectome_IN_skids.csv")

# read cells (I don't recommend running this in a batch because it is
# memory intensive and can run into curl fetch errors)

neurons <- nlapply(
  read.neurons.catmaid(connectome_neuron,
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

connectome <- nlapply(
  read.neurons.catmaid(connectome_no_frag,
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

Sensoryneuron <- nlapply(
  read.neurons.catmaid(Sensory_neuron, pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

Motorneuron <- nlapply(
  read.neurons.catmaid(motorneuron, pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

Interneuron <- nlapply(
  read.neurons.catmaid(interneuron, pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

muscle <- nlapply(
  read.neurons.catmaid("^muscle$",
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

endoderm <- nlapply(
  read.neurons.catmaid("^endoderm$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

# read in two blocks because curl often breaks for larger sets
epithelia_left <- nlapply(
  read.neurons.catmaid(
    skids_by_2annotations(
      "epithelia_cell",
      "left_side"
    ),
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

epithelia_right <- nlapply(
  read.neurons.catmaid(
    skids_by_2annotations(
      "epithelia_cell",
      "right_side"
    ),
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

skids_by_2annotations("epithelia_cell", "left_side")

gland <- nlapply(
  read.neurons.catmaid("^gland cell$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

Ciliary_band_cell <- nlapply(
  read.neurons.catmaid("^Ciliary_band_cell$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

connectome_pigment_skids <- skids_by_2annotations("pigment cell", "connectome")

connectome_pigment <- nlapply(
  read.neurons.catmaid(connectome_pigment_skids,
                       pid = 11, conn = conn_http1,
                       fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

glia <- nlapply(
  read.neurons.catmaid("^glia cell$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

stomodeum <- nlapply(
  read.neurons.catmaid("^stomodeum$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

# read in two blocks because curl often breaks for larger sets
pnb_left <- nlapply(
  read.neurons.catmaid(
    skids_by_2annotations(
      "pnb",
      "left_side"
    ),
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

pnb_right <- nlapply(
  read.neurons.catmaid(
    skids_by_2annotations(
      "pnb",
      "right_side"
    ),
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

notopodium <- nlapply(
  read.neurons.catmaid("^notopodium$",
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)

neuropodium <- nlapply(
  read.neurons.catmaid("^neuropodium$",
    pid = 11, conn = conn_http1,
    fetch.annotations = FALSE
  ),
  function(x) smooth_neuron(x, sigma = 6000)
)


# check if there are any cells with two or more tagged somas
sum <- summary(epithelia_left)
attributes(sum)
sum[sum$nsoma != 1, ]
as.numeric(rownames(sum[sum$nsoma == 2, ]))

# 3D plotting  ------------------------------------------------------------

# extract connectors to be able to plot them by unique colours
{
  SN_conn <- connectors(Sensoryneuron)
  str(SN_conn)
  presyn_SN_conn <- SN_conn[SN_conn$prepost == 0, ]
  postsyn_SN_conn <- SN_conn[SN_conn$prepost == 1, ]

  IN_conn <- connectors(Interneuron)
  presyn_IN_conn <- IN_conn[IN_conn$prepost == 0, ]
  postsyn_IN_conn <- IN_conn[IN_conn$prepost == 1, ]

  # we can also subset with the subset shorthand function
  MN_conn <- connectors(Motorneuron)
  presyn_MN_conn <- subset(MN_conn, prepost == 0)
  postsyn_MN_conn <- subset(MN_conn, prepost == 1)
}

# plot only the synapses coloured by SN MN IN
# cb friendly colour codes interneuron = "#CC79A7", motoneuron = "#0072B2",  `sensory neuron` = "#E69F00"
{
  plot_background_ventral()
  
  # plot only the presyn connectors
  plot3d(presyn_SN_conn$x, presyn_SN_conn$y, presyn_SN_conn$z, size = 5, alpha = 0.5, col = "#E69F00", add = T)
  plot3d(presyn_MN_conn$x, presyn_MN_conn$y, presyn_MN_conn$z, size = 5, alpha = 0.5, col = "#0072B2", add = T)
  plot3d(presyn_IN_conn$x + 1, presyn_IN_conn$y, presyn_IN_conn$z, size = 5, alpha = 0.5, col = "#CC79A7", add = T)
  plot3d(scalebar_50um_ventral2, lwd = 3, color = "black")
  texts3d(97000, 109000, 182500, text = "50 μm", col = "black", cex = 2)

  }
rgl.snapshot("pictures/connectome_SN_IN_MN_synapses_ventral.png")
# we define a z clipping plane for the frontal view
{
  par3d(windowRect = c(20, 30, 800, 800)) # resize for frontal view
  nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
  clipplanes3d(0, 0, -1, 60000)
  par3d(zoom = 0.63)
}
rgl.snapshot("pictures/connectome_SN_IN_MN_synapses_frontal.png")
close3d()


# plot presyn and postsyn sites for SN, IN, MN in separate colours
{
  plot_background_ventral()
  par3d(windowRect = c(20, 30, 500, 800))
  par3d(zoom = 0.34)
  plot3d(presyn_SN_conn$x, presyn_SN_conn$y, presyn_SN_conn$z,
    size = 4, alpha = 0.3, col = "#0072B2", add = T
  )
  plot3d(postsyn_SN_conn$x, postsyn_SN_conn$y, postsyn_SN_conn$z,
    size = 5, alpha = 0.5, col = "#D55E00", add = T
  )

  # make snapshot
  rgl.snapshot("pictures/SN_pre_postsynaptic_sites.png")
  close3d()

  plot_background_ventral()
  par3d(windowRect = c(20, 30, 500, 800))
  par3d(zoom = 0.34)
  plot3d(presyn_IN_conn$x, presyn_IN_conn$y, presyn_IN_conn$z,
    size = 4, alpha = 0.3, col = "#0072B2", add = T
  )
  plot3d(postsyn_IN_conn$x, postsyn_IN_conn$y, postsyn_IN_conn$z,
    size = 5, alpha = 0.5, col = "#D55E00", add = T
  )

  # make snapshot
  rgl.snapshot("pictures/IN_pre_postsynaptic_sites.png")
  close3d()

  plot_background_ventral()
  par3d(windowRect = c(20, 30, 500, 800))
  par3d(zoom = 0.34)
  plot3d(presyn_MN_conn$x, presyn_MN_conn$y, presyn_MN_conn$z,
    size = 4, alpha = 0.3, col = "#0072B2", add = T
  )
  plot3d(postsyn_MN_conn$x, postsyn_MN_conn$y, postsyn_MN_conn$z,
    size = 5, alpha = 0.5, col = "#D55E00", add = T
  )

  # make snapshot
  rgl.snapshot("pictures/MN_pre_postsynaptic_sites.png")
  close3d()
}

# VNC cross sections ------------------------------------------------------

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

{
  plot_background_VNC(-73000, 78000)
  plot3d(Sensoryneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 5,
    add = T, alpha = 0.8, col = "#E69F00"
  )
  plot3d(Motorneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 4,
    rev = FALSE, fixup = F, add = T, alpha = 0.9, col = "#0072B2"
  )
  plot3d(connectome,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 7,
    rev = FALSE, fixup = F, add = T, alpha = 0.25, col = "#cccccc"
  )
  rgl.snapshot("pictures/neurites_VNC_SN_MN1.png")
  close3d()

  plot_background_VNC(-108300, 111000)
  clipplanes3d(0, 1, 0.16, -105000)
  plot3d(Sensoryneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 5,
    add = T, alpha = 0.9, col = "#E69F00"
  )
  plot3d(Motorneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 5,
    rev = FALSE, fixup = F, add = T, alpha = 0.9, col = "#0072B2"
  )
  plot3d(connectome,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 7,
    rev = FALSE, fixup = F, add = T, alpha = 0.25, col = "#cccccc"
  )
  rgl.snapshot("pictures/neurites_VNC_SN_MN2.png")
  close3d()

  plot_background_VNC(-108300, 111000)
  clipplanes3d(0, 1, 0.16, -105000)
  plot3d(Interneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 4,
    add = T, alpha = 0.9, col = "#CC79A7"
  )
  plot3d(Motorneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 4,
    rev = FALSE, fixup = F, add = T, alpha = 0.9, col = "#0072B2"
  )
  plot3d(connectome,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 7,
    rev = FALSE, fixup = F, add = T, alpha = 0.25, col = "#cccccc"
  )
  rgl.snapshot("pictures/neurites_VNC_IN_MN2.png")
  close3d()

  plot_background_VNC(-133000, 138000)
  clipplanes3d(0, 1, 0.16, -107000)
  clipplanes3d(-1, 0, 0.16, 136000)
  plot3d(Sensoryneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 5,
    add = T, alpha = 0.9, col = "#E69F00"
  )
  plot3d(Motorneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 5,
    rev = FALSE, fixup = F, add = T, alpha = 0.9, col = "#0072B2"
  )
  plot3d(connectome,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 7,
    rev = FALSE, fixup = F, add = T, alpha = 0.25, col = "#cccccc"
  )
  rgl.snapshot("pictures/neurites_VNC_SN_MN3.png")
  close3d()
}

# plot lateral view
{
  plot_background_ventral()
  # plot only the presyn connectors
  plot3d(presyn_SN_conn$x, presyn_SN_conn$y, presyn_SN_conn$z, size = 5, alpha = 0.5, col = "#E69F00", add = T)
  plot3d(presyn_MN_conn$x, presyn_MN_conn$y, presyn_MN_conn$z, size = 5, alpha = 0.5, col = "#0072B2", add = T)
  plot3d(presyn_IN_conn$x + 1, presyn_IN_conn$y, presyn_IN_conn$z, size = 5, alpha = 0.5, col = "#CC79A7", add = T)
  nview3d("left", extramat = rotationMatrix(-pi / 2, pi, -0.2, 0))
  clipplanes3d(1, 0, 0.16, -75700)
  plot3d(outline,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 2,
    rev = FALSE, fixup = F, add = T, forceClipregion = F, alpha = 0.07,
    col = "#E2E2E2"
  )
  par3d(zoom = 0.52)
}

rgl.snapshot("pictures/connectome_SN_IN_MN_synapses_left.png")

close3d()


# plot SN IN MN cells without a soma
{
  plot_background_ventral()
  plot3d(Sensoryneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 1, col = "#E69F00"
  )
  plot3d(Motorneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = "#0072B2",
  )
  plot3d(Interneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = "#CC79A7"
  )
}
rgl.snapshot("pictures/connectome_SN_IN_MN_cells_nosoma_ventral.png")
close3d()

# plot effectors only
{
  plot_background_ventral()
  plot3d(gland,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = "#0072B2"
  )
  plot3d(Ciliary_band_cell,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = "#E69F00"
  )
  plot3d(connectome_pigment,
         WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
         rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
         col = "black"
  )
  plot3d(muscle,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 2,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = hcl.colors(853, palette = "Reds")
  )
}
rgl.snapshot("pictures/connectome_effectors_frontal.png")

# add all other cells
{
  plot3d(Sensoryneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    add = T, alpha = 1, col = "#E69F00"
  )
  plot3d(Motorneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = "#0072B2",
  )
  plot3d(Interneuron,
    WithConnectors = F, WithNodes = F, soma = F, lwd = 1,
    rev = FALSE, fixup = F, add = T, alpha = 1, col = "#CC79A7"
  )
  plot3d(glia,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = hcl.colors(78, palette = "Peach")
  )
  plot3d(pnb_left,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = hcl.colors(length(pnb_left), palette = "Mint", rev = T)
  )
  plot3d(pnb_right,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = hcl.colors(length(pnb_right), palette = "Mint", rev = T)
  )
  plot3d(epithelia_left,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = hcl.colors(length(epithelia_left), palette = "Blues")
  )
  plot3d(epithelia_right,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 1,
    rev = FALSE, fixup = F, add = T, forceClipregion = TRUE, alpha = 1,
    col = hcl.colors(length(epithelia_right), palette = "Blues")
  )
}
rgl.snapshot("pictures/connectome_body_all_cells_ventral.png")
close3d()

# plot SN IN and MN only on left side with soma ---------------------------

# plot SN IN MN on left side only
# cb friendly colour codes interneuron = "#CC79A7", motoneuron = "#0072B2",  `sensory neuron` = "#E69F00"
{
  plot_background_ventral()
  skids_to_plot_left <- skids_by_2annotations("motorneuron", "left_side")
  skeletons_to_plot_left <- nlapply(
    read.neurons.catmaid(skids_to_plot_left,
      pid = 11, conn = conn_http1,
      fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  plot3d(skeletons_to_plot_left,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    rev = FALSE, fixup = F, add = T, alpha = 0.6, col = "#0072B2"
  )
  rgl.snapshot("pictures/MN_left.png")

  skids_to_plot_left <- skids_by_2annotations("Sensory neuron", "left_side")
  skeletons_to_plot_left <- nlapply(
    read.neurons.catmaid(skids_to_plot_left,
      pid = 11, conn = conn_http1,
      fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  plot3d(skeletons_to_plot_left,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    add = T, alpha = 0.6, col = "#E69F00"
  )
  par3d(zoom = 0.47)
  rgl.snapshot("pictures/SN_MN_left.png")

  segments3d(
    x = as.vector(c(500, 145000)),
    y = as.vector(c(110000, 110000)),
    z = as.vector(c(78000, 78000))
  )
  segments3d(
    x = as.vector(c(500, 145000)),
    y = as.vector(c(110000, 110000)),
    z = as.vector(c(111000, 111000))
  )
  segments3d(
    x = as.vector(c(500, 145000)),
    y = as.vector(c(110000, 110000)),
    z = as.vector(c(138000, 138000))
  )
  rgl.snapshot("pictures/SN_MN_left_segments.png")
  close3d()
}

# plot SN and IN left side only
{
  IN_left <- nlapply(
    read.neurons.catmaid(
      skids_by_3annotations(
        "connectome_neuron", "interneuron",
        "left_side"
      ),
      pid = 11, conn = conn,
      fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  plot_background_ventral()
  plot3d(skeletons_to_plot_left,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    rev = FALSE, fixup = F, add = T, alpha = 0.6, col = "#E69F00"
  )
  plot3d(IN_left,
    WithConnectors = F, WithNodes = F, soma = T, lwd = 2,
    rev = FALSE, fixup = F, add = T, alpha = 0.6, col = "#CC79A7"
  )
  par3d(zoom = 0.47)
  rgl.snapshot("pictures/SN_IN_left.png")
  close3d()
}

# SN IN MN statistics plots -----------------------------------------------

# plot the ratio of input to output synapses (I-O)/(I+O)
# for sensory, inter and motoneurons, as a function of cable length
# the data were downloaded from the Catmaid measurements widget


# synapse distribution plots --------------------------------

{
  # data exported from CATMAID Morhology plot with a resolution of 1000 nm
  # and root node as start point
  # annoations: connectome_Sensory_neuron, connectome_interneuron and connectome_motorneuron
  SN_in <- read.csv2("data/SN_Radial_density_of_input_synapses.csv", sep = ",")
  SN_out <- read.csv2("data/SN_Radial_density_of_output_synapses.csv", sep = ",")
  IN_in <- read.csv2("data/IN_Radial_density_of_input_synapses.csv", sep = ",")
  IN_out <- read.csv2("data/IN_Radial_density_of_output_synapses.csv", sep = ",")
  MN_in <- read.csv2("data/MN_Radial_density_of_input_synapses.csv", sep = ",")
  MN_out <- read.csv2("data/MN_Radial_density_of_output_synapses.csv", sep = ",")
}

# tidy the data -----------------------------------------------------------

{
  MN_in_tb <- MN_in %>%
    rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    rename_with(~ gsub("Radial.density.of.input.synapses", "neuron", .x, fixed = TRUE)) %>%
    pivot_longer(matches("0"),
      names_to = c("nm"),
      values_to = "synapses"
    ) %>%
    mutate(neuron_type = "MN") %>%
    mutate(synapse_type = "post") %>%
    mutate(nm = as.integer(nm)) %>%
    group_by(nm) %>%
    mutate(mean = mean(synapses))
  MN_in_tb

  MN_out_tb <- MN_out %>%
    rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    rename_with(~ gsub("Radial.density.of.output.synapses", "neuron", .x, fixed = TRUE)) %>%
    pivot_longer(matches("0"),
      names_to = c("nm"),
      values_to = "synapses"
    ) %>%
    mutate(neuron_type = "MN") %>%
    mutate(synapse_type = "pre") %>%
    mutate(nm = as.integer(nm)) %>%
    group_by(nm) %>%
    mutate(mean = mean(synapses))

  IN_in_tb <- IN_in %>%
    rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    rename_with(~ gsub("Radial.density.of.input.synapses", "neuron", .x, fixed = TRUE)) %>%
    pivot_longer(matches("0"),
      names_to = c("nm"),
      values_to = "synapses"
    ) %>%
    mutate(neuron_type = "IN") %>%
    mutate(synapse_type = "post") %>%
    mutate(nm = as.integer(nm)) %>%
    group_by(nm) %>%
    mutate(mean = mean(synapses))

  IN_out_tb <- IN_out %>%
    rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    rename_with(~ gsub("Radial.density.of.output.synapses", "neuron", .x, fixed = TRUE)) %>%
    pivot_longer(matches("0"),
      names_to = c("nm"),
      values_to = "synapses"
    ) %>%
    mutate(neuron_type = "IN") %>%
    mutate(synapse_type = "pre") %>%
    mutate(nm = as.integer(nm)) %>%
    group_by(nm) %>%
    mutate(mean = mean(synapses))

  SN_in_tb <- SN_in %>%
    rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    rename_with(~ gsub("Radial.density.of.input.synapses", "neuron", .x, fixed = TRUE)) %>%
    pivot_longer(matches("0"),
      names_to = c("nm"),
      values_to = "synapses"
    ) %>%
    mutate(neuron_type = "SN") %>%
    mutate(synapse_type = "post") %>%
    mutate(nm = as.integer(nm)) %>%
    group_by(nm) %>%
    mutate(mean = mean(synapses))

  SN_out_tb <- SN_out %>%
    rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    rename_with(~ gsub("Radial.density.of.output.synapses", "neuron", .x, fixed = TRUE)) %>%
    pivot_longer(matches("0"),
      names_to = c("nm"),
      values_to = "synapses"
    ) %>%
    mutate(neuron_type = "SN") %>%
    mutate(synapse_type = "pre") %>%
    mutate(nm = as.integer(nm)) %>%
    group_by(nm) %>%
    mutate(mean = mean(synapses))

  # combine tibbles
  MN_syn <- full_join(MN_in_tb, MN_out_tb)
  SN_syn <- full_join(SN_in_tb, SN_out_tb)
  IN_syn <- full_join(IN_in_tb, IN_out_tb)
  SN_MN_syn <- full_join(MN_syn, SN_syn)
  All_syn <- full_join(SN_MN_syn, IN_syn)
}

# generate synapse distribution plots
{
  # plot for SN
  syn1 <- All_syn %>%
    filter(neuron_type == "SN") %>%
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
    labs(x = "distance from soma (µm)", y = "mean radial synapse density") +
    scale_color_manual(values = c(
      "#D55E00",
      "#0072B2"
    )) +
    scale_shape_manual(values = c(1, 2)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 10)
    )
  syn1

  d <- ggplot(diamonds, aes(carat, price))
  d + geom_point(alpha = 1 / 10)

  # plot for IN
  syn2 <- All_syn %>%
    filter(neuron_type == "IN") %>%
    ggplot(aes(x = nm / 1000, y = synapses, group = synapse_type, color = synapse_type)) +
    geom_smooth(show.legend = FALSE) +
    geom_point(aes(
      x = nm / 1000, y = mean, shape = synapse_type,
      color = synapse_type
    ), show.legend = FALSE, size = 1) +
    theme(panel.background = element_rect(fill = "grey95", color = "grey")) +
    labs(x = "distance from soma (µm)", y = "mean radial synapse density") +
    scale_color_manual(values = c(
      "#D55E00",
      "#0072B2"
    )) +
    scale_shape_manual(values = c(1, 2)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 10)
    )
  syn2

  # plot for MN
  syn3 <- All_syn %>%
    filter(neuron_type == "MN") %>%
    ggplot(aes(x = nm / 1000, y = synapses, group = synapse_type, color = synapse_type)) +
    geom_smooth(show.legend = FALSE) +
    geom_point(aes(
      x = nm / 1000, y = mean, shape = synapse_type,
      color = synapse_type
    ), show.legend = FALSE, size = 1) +
    scale_color_manual(values = c(
      "#D55E00",
      "#0072B2"
    )) +
    theme(panel.background = element_rect(fill = "grey95", color = "grey")) +
    labs(x = "distance from soma (µm)", y = "mean radial synapse density") +
    scale_shape_manual(values = c(1, 2)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.key.size = unit(4, "mm")
    ) + # Apply guides function
    guides(color = guide_legend("synapse"))
  syn3

  ggsave("plots/SN_syn_distr.png",
    width = 800, height = 800, limitsize = FALSE,
    units = c("px"), syn1
  )
  ggsave("plots/IN_syn_distr.png",
    width = 800, height = 800, limitsize = FALSE,
    units = c("px"), syn2
  )
  ggsave("plots/MN_syn_distr.png",
    width = 800, height = 800, limitsize = FALSE,
    units = c("px"), syn3
  )

  write.table(All_syn, "source_data/Figure1_source_data1.txt", sep = "\t")
  All_syn <- read.table("source_data/Figure1_source_data1.txt", sep = "\t")
}
# combine bottom multi-panels


# make table about cells and synapses -------------------------------------

{
  table <- plot_ly(
    type = "table",
    columnwidth = c(10, 7),
    columnorder = c(0, 1),
    header = list(
      values = c("Type", "Number"),
      align = c("center", "center"),
      line = list(width = 1, color = "black"),
      fill = list(color = c("#E69F00", "#0072B2")),
      font = list(family = "Arial", size = 14, color = "white")
    ),
    cells = list(
      values = rbind(
        c(
          "total cells in the body", "effector cells in the body", "neurons in connectome", "brain neurons in connectome",
          "trunk neurons in connectome", "sensory neurons in connectome",
          "interneurons in connectome", "motor neurons in connectome",
          "effectors in connectome"
        ),
        c(
          "9162", "1176", "1642", "724",
          "918", "473",
          "930", "239", "428"
        )
      ),
      align = c("center", "center"),
      line = list(color = "black", width = 0.3),
      font = list(family = "Arial", size = 12, color = c("black"))
    )
  )

  table
  saveNetwork(table, "pictures/connectome_stats_table.html")
  webshot2::webshot(
    url = "pictures/connectome_stats_table.html",
    file = "pictures/connectome_stats_table.png",
    vwidth = 400, vheight = 300, # define the size of the browser window
    cliprect = c(50, 20, 350, 220), zoom = 10
  )
}

# assemble figure -------------------------------------

# read png
{
  imgSEM <- readPNG("pictures/Platynereis_SEM_inverted_nolabel.png")
  
  img_all <- readPNG("pictures/connectome_body_all_cells_ventral.png")
  img_nosoma_ventr <- readPNG("pictures/connectome_SN_IN_MN_cells_nosoma_ventral.png")
  img_eff <- readPNG("pictures/connectome_effectors_frontal.png")

  img_table <- readPNG("pictures/connectome_stats_table.png")

  img_syn_ventr <- readPNG("pictures/connectome_SN_IN_MN_synapses_ventral.png")
  img_syn_ant <- readPNG("pictures/connectome_SN_IN_MN_synapses_frontal.png")
  img_syn_left <- readPNG("pictures/connectome_SN_IN_MN_synapses_left.png")


  img_SN_MN_left <- readPNG("pictures/SN_MN_left_segments.png")
  img_SN_IN_left <- readPNG("pictures/SN_IN_left.png")

  img_SN_MN_cross1 <- readPNG("pictures/neurites_VNC_SN_MN1.png")
  img_SN_MN_cross2 <- readPNG("pictures/neurites_VNC_SN_MN2.png")
  img_IN_MN_cross2 <- readPNG("pictures/neurites_VNC_IN_MN2.png")
  img_SN_MN_cross3 <- readPNG("pictures/neurites_VNC_SN_MN3.png")

  img_SN_syn_distr <- readPNG("plots/SN_syn_distr.png")
  img_IN_syn_distr <- readPNG("plots/IN_syn_distr.png")
  img_MN_syn_distr <- readPNG("plots/MN_syn_distr.png")

  img_cables_A <- readPNG("plots/SN_IN_MN_A.png")
  img_cables_B <- readPNG("plots/SN_IN_MN_B.png")
  img_cables_C <- readPNG("plots/SN_IN_MN_C.png")

  img_SN_prepost <- readPNG("pictures/SN_pre_postsynaptic_sites.png")
  img_IN_prepost <- readPNG("pictures/IN_pre_postsynaptic_sites.png")
  img_MN_prepost <- readPNG("pictures/MN_pre_postsynaptic_sites.png")
}

# convert png to image panel
{
  panelSEM <- ggdraw() + draw_image(imgSEM) +
    draw_label("head", x = 0.5, y = 0.85, size = 9) +
    draw_label("peristomium", x = 0.52, y = 0.69, size = 9) +
    draw_label("sg0", x = 0.54, y = 0.63, size = 9) +
    draw_label("sg1", x = 0.55, y = 0.54, size = 9) +
    draw_label("sg2", x = 0.56, y = 0.38, size = 9) +
    draw_label("sg3", x = 0.57, y = 0.22, size = 9) +
    draw_label("pygidium", x = 0.58, y = 0.07, size = 9) +
    draw_label("chaetae", x = 0.1, y = 0.6, size = 9) +
    draw_label("cilia", x = 0.85, y = 0.72, size = 9) +
    draw_label("Platynereis", x = 0.3, y = 0.99, size = 10, fontface = "italic") + 
    draw_label("larva", x = 0.55, y = 0.99, size = 10, fontface = "plain") +
    draw_label(paste("50 ", "\u00B5", "m", sep = ""), 
               x = 0.2, y = 0.07, size = 10)

  panel_all <- ggdraw() + draw_image(img_all) +
    draw_label("all cells",
      x = 0.2, y = 0.98, fontfamily = "sans", fontface = "plain",
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
    draw_label("a", x = 0.1, y = 0.93, size = 8) +
    draw_label("p", x = 0.1, y = 0.79, size = 8) 
  

  panel_nosoma_ventr <- ggdraw() + draw_image(img_nosoma_ventr) +
    draw_label("SN",
      x = 0.9, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("IN",
      x = 0.9, y = 0.9, fontfamily = "sans", fontface = "plain",
      color = "#CC79A7", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("MN",
      x = 0.9, y = 0.85, fontfamily = "sans", fontface = "plain",
      color = "#0072B2", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("neurons (no soma)",
      x = 0.4, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("ventral",
      x = 0.9, y = 0.1, hjust = 1, fontfamily = "sans", fontface = "plain",
      color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1
    )

  panel_eff <- ggdraw() + draw_image(img_eff) +
    draw_label("effectors",
      x = 0.3, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("gland", x = 0.85, y = 0.3, color = "#0072B2", size = 10) +
    draw_label("ciliary band", x = 0.85, y = 0.25, color = "#E69F00", size = 10) +
    draw_label("pigment cell", x = 0.83, y = 0.2, color = "black", size = 10) +
    draw_label("muscle", x = 0.85, y = 0.15, color = "red", size = 10)
    

  panel_syn_ventr <- ggdraw() + draw_image(img_syn_ventr) +
    draw_label("all presynaptic sites",
      x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("ventral",
      x = 0.3, y = 0.1, hjust = 1, fontfamily = "sans", fontface = "plain",
      color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("SN",
      x = 0.9, y = 0.95, fontfamily = "sans", fontface = "plain",
      color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("IN",
      x = 0.9, y = 0.9, fontfamily = "sans", fontface = "plain",
      color = "#CC79A7", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("MN",
      x = 0.9, y = 0.85, fontfamily = "sans", fontface = "plain",
      color = "#0072B2", size = 10, angle = 0, lineheight = 0.9, alpha = 1
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
  

  panel_table <- ggdraw() + draw_image(img_table)


  panel_syn_ant <- ggdraw() + draw_image(img_syn_ant) +
    draw_label("frontal",
      x = 0.9, y = 0.1, hjust = 1, fontfamily = "sans", fontface = "plain",
      color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1
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
    draw_label("v", x = 0.1, y = 0.79, size = 8) +
    draw_line(
      x = c(0.6, 0.7),
      y = c(0.78, 0.88), size = 0.3) +
    draw_label("ciliary band nerve", x = 0.66, y = 0.9, size = 9 )
  

  panel_syn_left <- ggdraw() + draw_image(img_syn_left) +
    draw_label("left",
      x = 0.9, y = 0.1, hjust = 1, fontfamily = "sans", fontface = "plain",
      color = "black", size = 9, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    geom_segment(aes(x = 0.1,
                     y = 0.9,
                     xend = 0.23,
                     yend = 0.9),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
    geom_segment(aes(x = 0.23,
                     y = 0.9,
                     xend = 0.1,
                     yend = 0.9),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
    draw_label("v", x = 0.07, y = 0.9, size = 8) +
    draw_label("d", x = 0.26, y = 0.9, size = 8) 
  

  panel_SN_MN_left <- ggdraw() + draw_image(img_SN_MN_left) + draw_label("SN",
    x = 0.93, y = 0.95, fontfamily = "sans", fontface = "plain",
    color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1
  ) +
    draw_label("MN",
      x = 0.93, y = 0.9, fontfamily = "sans", fontface = "plain",
      color = "#0072B2", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("neurons (with soma), left side",
      x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("1", x = 0.1, y = 0.6, size = 8) +
    draw_label("2", x = 0.1, y = 0.46, size = 8) +
    draw_label("3", x = 0.1, y = 0.34, size = 8)

  panel_SN_IN_left <- ggdraw() + draw_image(img_SN_IN_left) + draw_label("SN",
    x = 0.93, y = 0.95, fontfamily = "sans", fontface = "plain",
    color = "#E69F00", size = 10, angle = 0, lineheight = 0.9, alpha = 1
  ) +
    draw_label("IN",
      x = 0.93, y = 0.9, fontfamily = "sans", fontface = "plain",
      color = "#CC79A7", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    ) +
    draw_label("neurons (with soma), left side",
      x = 0.5, y = 0.98, fontfamily = "sans", fontface = "plain",
      color = "black", size = 10, angle = 0, lineheight = 0.9, alpha = 1
    )

  panel_SN_MN_cross1 <- ggdraw() + draw_image(img_SN_MN_cross1)
  panel_SN_MN_cross2 <- ggdraw() + draw_image(img_SN_MN_cross2)
  panel_IN_MN_cross2 <- ggdraw() + draw_image(img_IN_MN_cross2)
  panel_SN_MN_cross3 <- ggdraw() + draw_image(img_SN_MN_cross3)
  panel_SN_MN_cross <- plot_grid(panel_SN_MN_cross1, NULL, panel_SN_MN_cross2, NULL,
    panel_IN_MN_cross2, NULL, panel_SN_MN_cross3,
    ncol = 1,
    rel_heights = c(1, -0.2, 1, -0.2, 1, -0.2, 1)
  ) +
    draw_label("SN", x = 0.93, y = 0.95, color = "#E69F00", size = 10) +
    draw_label("MN", x = 0.93, y = 0.9, color = "#0072B2", size = 10) +
    draw_label("MN", x = 0.93, y = 0.45, color = "#0072B2", size = 10) +
    draw_label("IN", x = 0.93, y = 0.4, color = "#CC79A7", size = 10) +
    draw_label("1", x = 0.05, y = 0.95, size = 8) +
    draw_label("2", x = 0.05, y = 0.70, size = 8) +
    draw_label("2", x = 0.05, y = 0.45, size = 8) +
    draw_label("3", x = 0.05, y = 0.23, size = 8)

  panel_SN_prepost <- ggdraw() + draw_image(img_SN_prepost) +
    draw_label("SN presyn",
      x = 0.4, y = 0.99,
      size = 10, color = "#0072B2"
    ) +
    draw_label("postsyn",
      x = 0.8, y = 0.99,
      size = 10, color = "#D55E00"
    ) +
    draw_label("aciculae", x = 0.89, y = 0.28, size = 8) +
    draw_line(
      x = c(0.59, 0.77, 0.79),
      y = c(0.15, 0.28, 0.43), size = 0.3) +
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
  

  panel_IN_prepost <- ggdraw() + draw_image(img_IN_prepost) +
    draw_label("IN presyn",
      x = 0.4, y = 0.99,
      size = 10, color = "#0072B2"
    ) +
    draw_label("postsyn",
      x = 0.8, y = 0.99,
      size = 10, color = "#D55E00"
    )

  panel_MN_prepost <- ggdraw() + draw_image(img_MN_prepost) +
    draw_label("MN presyn",
      x = 0.4, y = 0.99,
      size = 10, color = "#0072B2"
    ) +
    draw_label("postsyn",
      x = 0.8, y = 0.99,
      size = 10, color = "#D55E00"
    )

  panel_SN_syn_distr <- ggdraw() + draw_image(img_SN_syn_distr) +
    draw_label("sensory neurons", x = 0.4, y = 0.98, size = 10) +
    draw_label("presynaptic", x = 0.58, y = 0.68, size = 10, color = "#0072B2") +
    draw_label("postsynaptic", x = 0.45, y = 0.38, size = 10, color = "#D55E00")
  panel_IN_syn_distr <- ggdraw() + draw_image(img_IN_syn_distr) +
    draw_label("interneurons", x = 0.4, y = 0.98, size = 10) +
    draw_label("presynaptic", x = 0.62, y = 0.55, size = 10, color = "#0072B2") +
    draw_label("postsynaptic", x = 0.52, y = 0.73, size = 10, color = "#D55E00")
  panel_MN_syn_distr <- ggdraw() + draw_image(img_MN_syn_distr) +
    draw_label("motor neurons", x = 0.4, y = 0.98, size = 10) +
    draw_label("presynaptic", x = 0.7, y = 0.68, size = 10, color = "#0072B2") +
    draw_label("postsynaptic", x = 0.45, y = 0.36, size = 10, color = "#D55E00")

  panel_img_cables_A <- ggdraw() + draw_image(img_cables_A)
  panel_img_cables_B <- ggdraw() + draw_image(img_cables_B)
  panel_img_cables_C <- ggdraw() + draw_image(img_cables_C)
}

# save fig -----------
{
# define layout with textual representation for pathchwork assembly of figure
layout <- "
AAAAAAAAABBBBBBBBBCCCCCCCCCDDDDDDDDDEEEEEEEEEEEEEE
##################################################
FFFFFFFFFGGGGGGGGHHHHHHHHHIIIIIIIIIJJJJJJJJJKKKKKK
##################################################
LLLLLLLMMMMMMMMM#NNNNNNNOOOOOOOOO#PPPPPPPQQQQQQQQQ
"

# arrange three synapse panels with cowplot plot_grid
panel_syn_ant_left <- plot_grid(panel_syn_ant, panel_syn_left, nrow = 2)

panel_synapses <- plot_grid(panel_syn_ventr, panel_syn_ant_left,
    nrow = 1, rel_widths = c(2, 1)
  )

Fig1 <- panelSEM + panel_all + panel_nosoma_ventr + panel_eff + panel_table +
    panel_syn_ventr + panel_syn_left + panel_syn_ant + panel_SN_IN_left +
    panel_SN_MN_left + panel_SN_MN_cross +
    panel_SN_prepost + panel_SN_syn_distr +
    panel_IN_prepost + panel_IN_syn_distr +
    panel_MN_prepost + panel_MN_syn_distr +
    plot_layout(design = layout, guides = "collect", heights = c(1, 0.02, 1, 0.02, 0.8)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure1.png",
    limitsize = FALSE,
    units = c("px"), Fig1, width = 4400, height = 3400, bg = "white"
  )
}


ggsave("Figures/Figure1.pdf",
  limitsize = FALSE,
  units = c("px"), Fig1, width = 4400, height = 3400
)

