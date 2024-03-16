# R code to generate Figure 3 fig suppl 5 of the 3d Platynereis connectome paper
# Gaspar Jekely 2022-2023

# load packages, functions and anatomical references
source("code/libraries_functions_and_CATMAID_conn.R")

# Sholl analysis --------

annotation_neuronal_celltypelist <- list()
# read all neuronal cell types and all annotations
for (i in c(1:202)) {
  annotation <- paste("annotation:^celltype", i, "$", sep = "")
  # retrieve all annotations for the same neurons and create the annotations data frames
  annotation_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

list_position <- 0
cycle <- 0
sholl_left_mat <- matrix(data = NA, nrow = length(annotation_neuronal_celltypelist), ncol = 180)
sholl_right_mat <- matrix(data = NA, nrow = length(annotation_neuronal_celltypelist), ncol = 180)

for (df1 in annotation_neuronal_celltypelist) { # iterate through the cell type list  (all presyn cells per body region)
  left_skids <- df1[df1$annotation == "left_side", 1]
  right_skids <- df1[df1$annotation == "right_side", 1]
  cycle <- cycle + 1

  # skip asymmetric or midline cells (e.g. MC, pygPBunp)
  if (length(left_skids) == 0 | length(right_skids) == 0) {
    sholl_left_mat[cycle, ] <- rep(0, 180)
    sholl_right_mat[cycle, ] <- rep(0, 180)
    next
  }

  left_neurons <- nlapply(read.neurons.catmaid(left_skids, pid = 11, fetch.annotations = T), function(x) smooth_neuron(x, sigma = 6000))
  right_neurons <- nlapply(read.neurons.catmaid(right_skids, pid = 11, fetch.annotations = T), function(x) smooth_neuron(x, sigma = 6000))

  print(left_neurons)

  # do sholl analysis on left and right cells separately
  sholl_left <- sholl_analysis(
    left_neurons,
    start = soma(left_neurons),
    starting.radius = 1000,
    ending.radius = 180000,
    radius.step = 1000
  )

  print(right_neurons)

  sholl_right <- sholl_analysis(
    right_neurons,
    start = soma(right_neurons),
    starting.radius = 1000,
    ending.radius = 180000,
    radius.step = 1000
  )

  sholl_intersections_left <- lapply(sholl_left, function(x) (x$intersections))
  sholl_means_left <- rowMeans(as.data.frame(sholl_intersections_left[1:length(sholl_intersections_left)]))

  sholl_intersections_right <- lapply(sholl_right, function(x) (x$intersections))
  sholl_means_right <- rowMeans(as.data.frame(sholl_intersections_right[1:length(sholl_intersections_right)]))

  print(cycle)
  sholl_left_mat[cycle, ] <- sholl_means_left
  sholl_right_mat[cycle, ] <- sholl_means_right
}

# smoothing and averaging of Sholl plots by cell type ----------------------

# define smoothing function
lowess_fun <- function(x, f = 0) {
  return(lowess(x, f = 0.04, iter = 100))
}

lowess_sholl_left <- apply(sholl_left_mat, 1, lowess_fun, f = 0.04)
lowess_sholl_right <- apply(sholl_right_mat, 1, lowess_fun, f = 0.04)

# test plot
plot(lowess_sholl_left[[2]])

# create smoothed matrix
sholl_left_smoothed_mat <- sholl_left_mat
sholl_right_smoothed_mat <- sholl_right_mat

# add smoothed values
for (i in c(1:length(lowess_sholl_left))) {
  sholl_left_smoothed_mat[i, ] <- lowess_sholl_left[[i]]$y
  sholl_right_smoothed_mat[i, ] <- lowess_sholl_right[[i]]$y
}

# plot the smoothed sholl diagrams
matplot(t(sholl_left_smoothed_mat), type = "l")
sholl_left_smoothed_mat


# sholl - assign names ----------------------
# assign rownames
# get the neuron name of the first skid in the skids list
neuronal_celltype_names <- list()
first_skids <- lapply(annotation_neuronal_celltypelist, function(x) x$skid[1])
skids <- list()
for (i in 1:202) {
  name <- catmaid_get_neuronnames(unlist(first_skids)[i], pid = 11)
  name <- sub("_.*$", "", name)
  neuronal_celltype_names[i] <- name
  skids[i] <- (unlist(first_skids)[i])
}


rownames(sholl_left_smoothed_mat) <- neuronal_celltype_names
rownames(sholl_right_smoothed_mat) <- neuronal_celltype_names

# remove all-zero rows
# sholl_left_smoothed_mat.no0 <- sholl_left_smoothed_mat[rowSums(sholl_left_smoothed_mat) != 0, ]
# sholl_right_smoothed_mat.no0 <- sholl_right_smoothed_mat[rowSums(sholl_right_smoothed_mat) != 0, ]

# sholl data - tidy and merge left-right ----------------------------------

# left side

# convert to tibble
sholl_left_tb <- sholl_left_smoothed_mat %>%
  as_tibble(rownames = "neuron") %>%
  pivot_longer(
    cols = contains("V"),
    names_to = "radius",
    values_to = "sholl"
  ) %>%
  mutate(radius = gsub("V", "", radius)) %>%
  mutate_at(vars(radius), as.integer) %>%
  mutate(side = ("left"))
sholl_left_tb
# same for right side
# convert to tibble
sholl_right_tb <- sholl_right_smoothed_mat %>%
  as_tibble(rownames = "neuron") %>%
  pivot_longer(
    cols = contains("V"),
    names_to = "radius",
    values_to = "sholl"
  ) %>%
  mutate(radius = gsub("V", "", radius)) %>%
  mutate_at(vars(radius), as.integer) %>%
  mutate(side = ("right"))
sholl_right_tb

# join left and right tibble
sholl_tb <- full_join(sholl_left_tb, sholl_right_tb)

# write the data as supplement
write_csv2(sholl_tb, "source_data/Figure3_fig_suppl5_source_data1.csv")
sholl_tb <- read_csv2("source_data/Figure3_fig_suppl5_source_data1.csv")


# plot sholl heatmaps -----------------


#read file generated by Figure3.R
n_cells.tb <- read.table("source_data/Figure3_source_data_1.txt")

# get neuron names in SN IN MN order
SN_IN_MN_names <- n_cells.tb %>%
  filter(SN_MN_IN %in% c("Sensory neuron", "interneuron", "motorneuron")) %>%
  filter(region %in% "head") %>%
  select(celltype) %>%
  pull()
SN_IN_MN_names
asymmetric_neurons <- c(
  "MS1", "SN-IRP2-burs", "SN-IRP2-FMRF", "SNasym",
  "SN47Ach", "pygPBunp", "VentraltrunkPUunp-Glialike",
  "doCRunp", "SN-NS4", "SN-NS6", "SN-YF5cil", "MC2-biax-1", "INATOpyg",
  "cMNPDF-vcl1", "cMNATO", "cMNdc", "INsqNSasym", "MC",
  "hPU2l-asymPDF"
)

# neuron list without the asymmetric neurons
SN_IN_MN_names_new <- SN_IN_MN_names[!grepl(paste(asymmetric_neurons, collapse = "|"), SN_IN_MN_names)]


SN_IN_MN_names
hm_left <- sholl_tb %>%
  filter(side == "left") %>%
  ggplot(aes(x = radius, y = neuron, fill = sholl)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "#0072B2", "#E69F00", "#E69F00", "#E69F00", "#D55E00", "#D55E00", "#D55E00"),
    na.value = "white"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7, vjust = -1),
    axis.text.y = element_text(size = 5, hjust = 1),
    axis.text.x = element_text(size = 4.2, angle = 90, hjust = 0.5),
    legend.position = c(.88, .85),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    panel.grid.major = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    plot.title = element_text(size = 7, vjust = -8, hjust = 0),
    plot.margin = margin(t = -12, r = 2, b = 1, l = 0)
  ) +
  scale_y_discrete(limits = SN_IN_MN_names_new) +
  labs(x = "distance from soma (µm)", title = "left") +
  coord_flip()
hm_left

hm_right <- sholl_tb %>%
  filter(side == "right") %>%
  ggplot(aes(x = radius, y = neuron, fill = sholl)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_gradientn(
    colors = c("white", "#0072B2", "#E69F00", "#E69F00", "#E69F00", "#D55E00", "#D55E00", "#D55E00"),
    na.value = "white"
  ) +
  scale_x_reverse() +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 7, vjust = -1),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 5, hjust = 1),
    plot.title = element_text(size = 7, vjust = -5, hjust = 0),
    panel.grid.major = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.1),
    plot.margin = margin(t = -12, r = 2, b = 0, l = 0)
  ) +
  scale_y_discrete(position = "right", limits = SN_IN_MN_names_new) +
  labs(x = "distance from soma (µm)", title = "right") +
  coord_flip()
hm_right

# define layout to combine left and right panel with patchwork
layout_sholl <- "A
B"

# combine panels and save
sholl <- hm_left + hm_right +
  plot_layout(design = layout_sholl, widths = c(1, 1)) &
  theme(
    plot.margin = margin(t = -12, r = 1, b = 0, l = 0),
  )

ggsave("plots/sholl_heatmap_left_right.png",
  width = 3500, height = 1500, limitsize = FALSE,
  units = c("px"), sholl
)

ggsave("Figures/Figure3_fig_suppl5.png",
       width = 3500, height = 1500, limitsize = FALSE,
       units = c("px"), sholl
)
ggsave("Figures/Figure3_fig_suppl5.pdf",
       width = 3500, height = 1500, limitsize = FALSE,
       units = c("px"), sholl
)

