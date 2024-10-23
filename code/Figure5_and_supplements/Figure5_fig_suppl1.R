#R/natverse code to generate Figure 5 fig suppl 1 of the Platynereis 3d connectome paper
#Gaspar Jekely 2023

#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/libraries_functions_and_CATMAID_conn.R")

load_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(
    annotation, pid=11), 
    function(x) smooth_neuron(x, sigma=6000))
}

# load neurons ---------------
dcv <- load_neuron("^dense cored vesicles$")
dcv2 <- load_neuron("^nr dense core vesicles$")

yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                           invertFaces = T, conn = NULL, pid = 11)
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
scalebar_50um_ventral = read.neurons.catmaid("^scalebar_50um_ventral$", pid=11)

# dcv neuron plotting ----------

plot_background_ventral_no_ac()
plot3d(
  dcv, soma = TRUE, lwd = c(2,4), 
  color = bluepurple[2:9], alpha = c(3:7/10)
  )
plot3d(
  dcv2, soma = TRUE, lwd = c(2,4), 
  color = blues[2:9], alpha = c(3:7/10)
)
plot3d(scalebar_50um_ventral, lwd = 3, color = "black")
par3d(zoom=0.53)

rgl.snapshot("pictures/dcv_neurons_ventral.png")
nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
rgl.pop()
#z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/dcv_neurons_frontal.png")
close3d()

# dcv across neuron classes -----------------

dcv <- load_neuron("^dense cored vesicles$")
dcv2 <- load_neuron("^nr dense core vesicles$")

skids_dcv1 <- skids_by_2annotations("^dense cored vesicles$", "with_soma")
skids_dcv2 <- skids_by_2annotations("^nr dense core vesicles$", "with_soma")

all_dcv_skids <- union(skids_dcv1, skids_dcv2)

SN_skids <- skids_by_2annotations("^Sensory neuron$", "with_soma")
IN_skids <- skids_by_2annotations("^interneuron$", "with_soma")
MN_skids <- skids_by_2annotations("^motorneuron$", "with_soma")

num_SN_dcv <- length(intersect(all_dcv_skids, SN_skids))
num_IN_dcv <- length(intersect(all_dcv_skids, IN_skids))
num_MN_dcv <- length(intersect(all_dcv_skids, MN_skids))

cell_class_dvc <- tibble(class = c("SN", "IN", "MN"),
                         num_cell_with_dvc = c(num_SN_dcv, num_IN_dcv, num_MN_dcv)
                         )
dcv_num_plot <- cell_class_dvc |>
  ggplot(aes(class, num_cell_with_dvc, fill = class)) +
  geom_col() +
  scale_fill_manual(values = c("#E69F00", "#CC79A7", "#0072B2"),
                    limits = c("SN", "IN", "MN"))  +
  scale_x_discrete(limits = c("SN", "IN", "MN")) +
  labs(x = "", y = "number of cells with dcv") +
  guides(fill = "none")

# geom_dotplot()NULL# geom_dotplot()# crop dcv images ---------------------------------------------------------

crop_catmaid("nice DCV", 800, 800, 1, 0, "pictures/", 11, 21)

#list all files starting with crop
system("ls ./pictures/crop_*")

# assemble figure ----------------

dcv1 <- ggdraw() + draw_image(magick::image_read("pictures/crop_6l53xv.tiff"))
dcv2 <- ggdraw() + draw_image(magick::image_read("pictures/crop_w7ypfl.tiff"))
dcv3 <- ggdraw() + draw_image(magick::image_read("pictures/crop_vijo5k.tiff"))
dcv4 <- ggdraw() + draw_image(magick::image_read("pictures/crop_rd2i2q.tiff"))
dcv5 <- ggdraw() + draw_image(magick::image_read("pictures/crop_qde6ql.tiff"))
dcv6 <- ggdraw() + draw_image(magick::image_read("pictures/crop_pd3d1g.tiff"))
dcv7 <- ggdraw() + draw_image(magick::image_read("pictures/crop_p0ruba.tiff"))
dcv8 <- ggdraw() + draw_image(magick::image_read("pictures/crop_oszyea.tiff"))
dcv9 <- ggdraw() + draw_image(magick::image_read("pictures/crop_fvxd95.tiff"))
dcv10 <- ggdraw() + draw_image(magick::image_read("pictures/crop_blvkgf.tiff"))
dcv11 <- ggdraw() + draw_image(magick::image_read("pictures/crop_apdckt.tiff"))
dcv12 <- ggdraw() + draw_image(magick::image_read("pictures/crop_7omp7s.tiff"))

panel_dcvv <- ggdraw() + draw_image(readPNG("pictures/dcv_neurons_ventral.png")) + 
  draw_label(
    "dcv neurons, ventral view", x = 0.5, y = 0.99, 
    color = "black", size = 11
  ) +
  draw_label(
    expression(paste("50 ", mu, " m")), 
    x = 0.73, y = 0.08, 
    color = "black", size = 10
    )  +
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

panel_dcva <-ggdraw() + draw_image(readPNG("pictures/dcv_neurons_frontal.png")) + 
  draw_label(
    "dcv neurons, anterior view", x = 0.5, y = 0.99, 
    color = "black", size = 11
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

layout <- "
AbBC#D#E
AbBG#H#I
AbBK#L#M
AbBN#O#P
"

Figure5_fig_suppl1 <- panel_dcvv + panel_dcva + 
  dcv_num_plot + dcv1 + dcv2 + dcv3 + dcv4 + dcv5 + dcv6 + dcv7 + dcv8 + dcv9 + dcv10 + dcv11 + dcv12 +
  plot_layout(design = layout, widths = c(1.2, 1, 0.25, 0.25, 0.01, 0.25, 0.01, 0.25)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(size=12, face='plain'))

ggsave("Figures/Figure5_fig_suppl1.png", limitsize = FALSE, 
       units = c("px"), Figure5_fig_suppl1, 
       width = 2600, height = 1000, bg='white'
       )  

ggsave("Figures/Figure5_fig_suppl1.pdf", limitsize = FALSE, 
       units = c("px"), Figure5_fig_suppl1, 
       width = 2600, height = 1000
       )  




