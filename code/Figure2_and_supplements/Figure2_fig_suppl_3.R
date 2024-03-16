# R code to generate Figure 2 fig suppl 3 of connectome modules in the 3d Platynereis connectome paper
# Gaspar Jekely 2024

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# plot graph with coordinates from gephi ----------------------------------

# read graph
connect.tb <- readRDS("supplements/connectome_graph_tibble.rds")

# plot neurons by colours matching network pic -----------------------------

#background plot
plot_neuron_modules <- function(neurons, color, alpha){
  plot3d(neurons, soma=FALSE, lwd=2,
         add=T, alpha=alpha,
         col=color)
  par3d(zoom = 0.53)
}

# select modules --------

connect.tb %>%
  select(names, partition) %>%
  as_tibble() %>%
  filter(names == "MN1_l")

# eye, 11
# MB, 6
# AO, 2
# postural, 7
plot_background()
# plot modules  in the colour that was used in the network plot
for (i in c(2, 6, 11, 7)) {
  skids <- connect.tb %>%
    filter(partition == i) %>%
    pull(skids)
  
  color <- connect.tb %>%
    filter(partition == i) %>%
    pull(color) %>%
    unique()
  
  neurons = nlapply(read.neurons.catmaid(skids, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  plot_neuron_modules(
    neurons, 
    color, 
    runif(length(skids), min = 0.4, max = 0.7)
  )
}
par3d(zoom = 0.25)

rgl.snapshot("pictures/conn_module_eye_MB_AO_postural.png")
close3d()

# plot neuropils --------------

#Mechano r, 1
#Mechano l, 4
#"muscle 1", 8
#"muscle 2", 3
#"muscle 3", 13
#"muscle 4", 9
#"muscle 5", 12

#background plot
plot_neuron_modules <- function(neurons, color, alpha){
  plot3d(neurons, soma=FALSE, lwd=3,
         add=T, alpha=alpha,
         col=color)
  par3d(zoom = 0.53)
}

plot_background_ventral()

# clipping plane --------------

clipplanes3d(0, 0, 1, -73000)
clipplanes3d(0, 0, -1, 78000)
nview3d("frontal", extramat = rotationMatrix(0.2, 1, 0.1, 0.5))
par3d(windowRect = c(0, 0, 800, 800)) # resize for frontal view
par3d(zoom = 0.22)
# x-axis clip
clipplanes3d(0, 1, 0.16, -90000)
# y-axis clip
clipplanes3d(1, 0, 0.16, -13000)

# plot modules  in the colour that was used in the network plot
for (i in c(1, 4, 6, 7)) {
  skids <- connect.tb %>%
    filter(partition == i) %>%
    pull(skids)
  
  color <- connect.tb %>%
    filter(partition == i) %>%
    pull(color) %>%
    unique()
  
  neurons = nlapply(read.neurons.catmaid(skids, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  plot_neuron_modules(
    neurons, 
    color, 
    runif(length(skids), min = 0.8, max = 1)
  )
}
par3d(zoom = 0.22)

rgl.snapshot("pictures/conn_module_Mech_MB_post.png")
close3d()


# assemble figure ------------------


panel_1 <- ggdraw() + draw_image(readPNG("pictures/conn_module_eye_MB_AO_postural.png")) +
  draw_label("brain neuropils", x = 0.4, y = 0.99, size = 10) +
  draw_label("visual", x = 0.47, y = 0.44, size = 11) +
  draw_label("MB and central brain", x = 0.3, y = 0.35, size = 11) +
  draw_label("anterior NS", x = 0.5, y = 0.6, size = 11) +
  draw_label("postural control", x = 0.3, y = 0.2, size = 11)

panel_2 <- ggdraw() + draw_image(readPNG("pictures/conn_module_Mech_MB_post.png")) +
  draw_label("VNC neuropils", x = 0.4, y = 0.99, size = 10) +
  draw_label("Mech. (r)", x = 0.85, y = 0.35, size = 11) +  
  draw_label("Mech. (l)", x = 0.2, y = 0.35, size = 11) +
  draw_label("MB and central brain", x = 0.54, y = 0.52, size = 11) +
  draw_label("postural control", x = 0.5, y = 0.63, size = 11)

# define layout
layout <- "
A#B
"

Figure2_fig_suppl_neuropil <- panel_1 + panel_2  +
  plot_layout(design = layout,  widths = c(1, 0.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure2_fig_suppl_3.png",
       limitsize = FALSE,
       units = c("px"), Figure2_fig_suppl_neuropil, width = 1800, height = 960, bg = "white"
)

ggsave("Figures/Figure2_fig_suppl_3.pdf",
       limitsize = FALSE,
       units = c("px"), Figure2_fig_suppl_neuropil, width = 1800, height = 960
)

