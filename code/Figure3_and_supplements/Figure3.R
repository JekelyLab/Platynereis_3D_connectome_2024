# R code to generate Figure 3 of the 3d Platynereis connectome paper
# Uses Natverse/catmaid and accesses the data on CATMAID
# Gaspar Jekely 2022

# load packages, functions and anatomical references
source("code/Natverse_functions_and_conn.R")

# table for neuronal cell types -------------------------------

annotation_neuronal_celltypelist <- list()
# read all neuronal cell types and all annotations
for (i in c(1:202)) {
  annotation <- paste("annotation:^celltype", i, "$", sep = "")
  # retrieve all annotations for the same neurons and create the annotations data frames
  annotation_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

# define the seven body regions, matching the catmaid annotations
regions <- c("episphere", "peristomium", "segment_0", "segment_1", "segment_2", "segment_3", "pygidium")
types <- paste("celltype", 1:202, sep = "")

# we retrieve those skids that match our annotations
# three iterated lapply functions to retrieve skids by cell type, body region
cells_per_body_region <- lapply(regions, function(r) lapply(types, function(m) skids_by_annotation(annotation_neuronal_celltypelist, m, r)))

n_cells <- matrix(nrow = 7, ncol = 202)
# count the occurrence of each type in each body region per side
for (j in 1:7) {
  for (i in 1:202) {
    n_cells[j, i] <- length(cells_per_body_region[[j]][[i]])
  }
}

# add row and column names
row.names(n_cells) <- c("head", "peristomium", "sg0", "sg1", "sg2", "sg3", "pyg")

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

# neuronal cell type names
colnames(n_cells) <- neuronal_celltype_names

# query SN IN MN annotations
contains_SN <- c()
contains_IN <- c()
contains_MN <- c()

for (i in seq_along(annotation_neuronal_celltypelist)) {
  my_list <- annotation_neuronal_celltypelist[[i]]
  annotations <- my_list[my_list$skid == skids[i], "annotation"]
  contains_SN[i] <- "Sensory neuron" %in% annotations
  contains_IN[i] <- "interneuron" %in% annotations
  contains_MN[i] <- "motorneuron" %in% annotations
}

combined_annotations <- character(length(contains_SN))
combined_annotations[contains_SN] <- "Sensory neuron"
combined_annotations[contains_IN] <- "interneuron"
combined_annotations[contains_MN] <- "motorneuron"

# convert to tibble
n_cells.tb <- as.data.frame(n_cells) %>%
  rownames_to_column(var = "region") %>%
  pivot_longer(cols = -region, names_to = "celltype", values_to = "number") %>%
  mutate(skid = rep(skids, 7)) %>%
  mutate(SN_MN_IN = rep(combined_annotations, 7))

# Reorder 'celltype' based on 'SN_MN_IN' and 'number'
n_cells.tb <- n_cells.tb %>%
  group_by(celltype) %>%
  mutate(sum_cells = sum(number)) %>%
  ungroup() %>%
  arrange(SN_MN_IN, desc(sum_cells)) %>%
  mutate(SN_MN_IN = factor(SN_MN_IN,
    levels = c("Sensory neuron", "interneuron", "motorneuron"),
    ordered = TRUE
  )) %>%
  mutate(celltype = recode_factor(celltype,
    "Sensory neuron" = "Sensory neuron",
    "interneuron" = "interneuron",
    "motorneuron" = "motorneuron",
    .ordered = TRUE
  ))

SN_names <- n_cells.tb %>%
  filter(SN_MN_IN %in% c("Sensory neuron")) %>%
  filter(region %in% "head") %>%
  select(celltype) %>%
  pull()
IN_names <- n_cells.tb %>%
  filter(SN_MN_IN %in% c("interneuron")) %>%
  filter(region %in% "head") %>%
  select(celltype) %>%
  pull()
MN_names <- n_cells.tb %>%
  filter(SN_MN_IN %in% c("motorneuron")) %>%
  filter(region %in% "head") %>%
  select(celltype) %>%
  pull()


# theme for plot ----------------------------------------------------------

theme_celltypes <- theme(
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  plot.title = element_text(size = 12, face = "plain"),
  axis.text.x = element_text(
    angle = 90, hjust = 1, vjust = 0.5,
    size = 8, color = "black"
  ),
  axis.text.y = element_text(
    angle = 0, hjust = 1, vjust = 0.5,
    size = 8, color = "black"
  ),
  axis.title = element_blank(),
  axis.ticks = element_line(linewidth = 0.2),
  axis.ticks.length = unit(0.4, "mm"),
  axis.line = element_line(linewidth = 0.2),
  legend.text = element_text(size = 8, hjust = 0, margin = margin(r = 2)),
  legend.key.height = unit(0.1, "mm"),
  legend.box.spacing = unit(0, "mm"),
  legend.position = "top"
)

# plots ----------

SN_plot <- n_cells.tb %>%
  ggplot(aes(celltype, region, color = SN_MN_IN)) +
  geom_point(
    data = n_cells.tb %>% filter(number > 0),
    aes(alpha = log(number)),
    shape = 15,
    stroke = 0,
    size = 4
  ) +
  geom_text(
    data = n_cells.tb %>% filter(number > 0), # Filter data first
    aes(label = number), color = "black", size = 2
  ) +
  theme_half_open() +
  theme_celltypes +
  scale_color_manual(
    values = c(
      "Sensory neuron" = "#E69F00",
      "interneuron" = "#CC79A7",
      "motorneuron" = "#0072B2"
    ),
    labels = c(
      "Sensory neuron" = "SN",
      "interneuron" = "IN",
      "motorneuron" = "MN"
    )
  ) +
  scale_y_discrete(limits = rev(c("head", "sg0", "sg1", "sg2", "sg3", "pyg"))) +
  scale_x_discrete(limits = SN_names) +
  guides(alpha = "none", color = "none") +
  labs(color = NULL, title = "sensory neurons")

IN_plot <- n_cells.tb %>%
  ggplot(aes(celltype, region, color = SN_MN_IN)) +
  geom_point(
    data = n_cells.tb %>% filter(number > 0),
    aes(alpha = log(number)),
    shape = 15,
    stroke = 0,
    size = 4
  ) +
  geom_text(
    data = n_cells.tb %>% filter(number > 0), # Filter data first
    aes(label = number), color = "black", size = 2
  ) +
  theme_half_open() +
  theme_celltypes +
  scale_color_manual(
    values = c(
      "Sensory neuron" = "#E69F00",
      "interneuron" = "#CC79A7",
      "motorneuron" = "#0072B2"
    ),
    labels = c(
      "Sensory neuron" = "SN",
      "interneuron" = "IN",
      "motorneuron" = "MN"
    )
  ) +
  scale_y_discrete(limits = rev(c("head", "sg0", "sg1", "sg2", "sg3", "pyg"))) +
  scale_x_discrete(limits = IN_names) +
  guides(alpha = "none", color = "none") +
  labs(color = NULL, title = "interneurons")

SN_plot
IN_plot

length(SN_names)
length(IN_names)
length(MN_names)


MN_plot <- n_cells.tb %>%
  ggplot(aes(celltype, region, color = SN_MN_IN)) +
  geom_point(
    data = n_cells.tb %>% filter(number > 0),
    aes(alpha = log(number)),
    shape = 15,
    stroke = 0,
    size = 4 
  ) +
  geom_text(
    data = n_cells.tb %>% filter(number > 0), # Filter data first
    aes(label = number), color = "black", size = 2
  ) +
  theme_half_open() +
  theme_celltypes +
  scale_color_manual(values = c(
    "Sensory neuron" = "#E69F00",
    "interneuron" = "#CC79A7",
    "motorneuron" = "#0072B2"
  )) +
  scale_y_discrete(limits = rev(c("head", "sg0", "sg1", "sg2", "sg3", "pyg"))) +
  scale_x_discrete(limits = MN_names) +
  guides(alpha = "none", color = "none") +
  labs(title = "motor neurons")


n_cells.tb <- n_cells.tb %>%
  select(-skid)

# save tibble
write.table(n_cells.tb, file = "source_data/Figure3_source_data_1.txt", sep ="\t")

# table for non-neuronal cell types -------------------------------

annotation_non_neuronal_celltypelist <- list()
# read all non-neuronal cell types and all annotations
for (i in c(1:92)) {
  annotation <- paste("annotation:^celltype_non_neuronal", i, "$", sep = "")
  # retrieve all annotations for the same neurons and create the annotations data frames
  annotation_non_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

# define the seven body regions, matching the catmaid annotations
regions <- c("episphere", "peristomium", "segment_0", "segment_1", "segment_2", "segment_3", "pygidium")
types <- paste("celltype_non_neuronal", 1:92, sep = "")

# retrieve those skids that match our annotations
# three iterated lapply functions to retrieve skids by cell type, body region
non_neuronal_cells_per_body_region <- lapply(regions, function(r) lapply(types, function(m) skids_by_annotation(annotation_non_neuronal_celltypelist, m, r)))

n_non_neuronal_cells <- matrix(nrow = 7, ncol = 92)
# count the occurrence of each type in each body region per side
for (j in 1:7) {
  for (i in 1:92) {
    n_non_neuronal_cells[j, i] <- length(non_neuronal_cells_per_body_region[[j]][[i]])
  }
}

# add row and column names
row.names(n_non_neuronal_cells) <- c("head", "peristomium", "sg0", "sg1", "sg2", "sg3", "pyg")

# get the neuron name of the first skid in the skids list
non_neuronal_celltype_names <- list()
first_skids_non_neuronal <- lapply(annotation_non_neuronal_celltypelist, function(x) x$skid[1])
skids <- list()
for (i in 1:92) {
  name <- catmaid_get_neuronnames(unlist(first_skids_non_neuronal)[i], pid = 11)
  name <- sub("_.*$", "", name)
  non_neuronal_celltype_names[i] <- name
  skids[i] <- (unlist(first_skids_non_neuronal)[i])
}

# non_neuronal cell type names
colnames(n_non_neuronal_cells) <- non_neuronal_celltype_names

# query SN IN MN annotations
contains_MUS <- c()
contains_cil <- c()
contains_FC <- c()
contains_gl <- c()
contains_pig <- c()
contains_other <- c()

for (i in seq_along(annotation_non_neuronal_celltypelist)) {
  my_list <- annotation_non_neuronal_celltypelist[[i]]
  annotations <- my_list[my_list$skid == skids[i], "annotation"]
  contains_other[i] <- "with_soma" %in% annotations
  contains_MUS[i] <- "muscle" %in% annotations
  contains_cil[i] <- "ciliated cell" %in% annotations
  contains_FC[i] <- "follicle cell" %in% annotations
  contains_gl[i] <- "gland cell" %in% annotations
  contains_pig[i] <- "pigment cell" %in% annotations
}


combined_annotations <- character(length(contains_MUS))
combined_annotations[contains_other] <- "other" # the ones not overwritten below
combined_annotations[contains_MUS] <- "muscle"
combined_annotations[contains_cil] <- "ciliated cell"
combined_annotations[contains_FC] <- "follicle cell"
combined_annotations[contains_gl] <- "gland cell"
combined_annotations[contains_pig] <- "pigment cell"


# convert to tibble
non_n_cells.tb <- as.data.frame(n_non_neuronal_cells) %>%
  rownames_to_column(var = "region") %>%
  pivot_longer(cols = -region, names_to = "celltype", values_to = "number") %>%
  mutate(skid = rep(skids, 7)) %>%
  mutate(ANNOT = rep(combined_annotations, 7))

# Reorder 'celltype' based on 'ANNOT' and 'number'
non_n_cells.tb <- non_n_cells.tb %>%
  group_by(celltype) %>%
  mutate(sum_cells = sum(number)) %>%
  ungroup() %>%
  arrange(ANNOT, desc(sum_cells)) %>%
  mutate(ANNOT = factor(ANNOT,
    levels = c(
      "muscle", "ciliated cell",
      "gland cell", "pigment cell",
      "follicle cell", "other"
    ),
    ordered = TRUE
  )) %>%
  mutate(celltype = recode_factor(
    celltype,
    "muscle" = "muscle",
    "gland cell" = "gland cell",
    "ciliated cell" = "ciliated cell",
    "pigment cell" = "pigment cell",
    "follicle cell" = "follicle cell",
    "other" = "other",
    .ordered = TRUE
  ))

names <- non_n_cells.tb %>%
  filter(region %in% "head") %>%
  select(celltype) %>%
  pull()

# non_n_cell_plot ---------------------------------------------------------

# neuron list without the asymmetric neurons
names <- names[!grepl(paste("MUSph", collapse = "|"), names)]
names
length(names)
non_n_cell_plot <- non_n_cells.tb %>%
  ggplot(aes(celltype, region, color = ANNOT)) +
  geom_point(
    data = non_n_cells.tb %>% filter(number > 0),
    aes(alpha = log(number)),
    shape = 15,
    stroke = 0,
    size = 4
  ) +
  geom_text(
    data = non_n_cells.tb %>% filter(number > 0), # Filter data first
    aes(label = number), color = "black", size = 2
  ) +
  theme_half_open() +
  theme_celltypes +
  scale_color_manual(
    values = c(
      "muscle" = "#D55E00",
      "ciliated cell" = "#F0E442",
      "gland cell" = "#56B4E9",
      "pigment cell" = "#009E73",
      "follicle cell" = "#aaB4E9",
      "other" = "#AAAAAA"
    ),
    labels = c(
      "muscle" = "muscle",
      "ciliated cell" = "ciliated",
      "gland cell" = "gland",
      "pigment cell" = "pigment",
      "follicle cell" = "follicle",
      "other" = "other"
    )
  ) +
  scale_x_discrete(limits = names) +
  scale_y_discrete(limits = rev(c("head", "sg0", "sg1", "sg2", "sg3", "pyg"))) +
  guides(alpha = "none") +
  labs(color = NULL, title = "non-neuronal cell types")
non_n_cell_plot

non_n_cells.tb <- non_n_cells.tb %>%
  select(-skid)

# save tibble
write.table(non_n_cells.tb, file = "source_data/Figure3_source_data_2.txt", sep = "\t")


# plot number of cells per cell type as histogram --------------------------

# count the number of skids per celltype
neuronal_cells_per_celltype <-
  sapply(annotation_neuronal_celltypelist, function(x) length(unique(x$skid)))

non_neuronal_cells_per_celltype <-
  sapply(annotation_non_neuronal_celltypelist, function(x) length(unique(x$skid)))

# plot histogram for neurons
hist_neuro <- as.data.frame(neuronal_cells_per_celltype[1:length(neuronal_cells_per_celltype)]) %>%
  ggplot(aes(x = neuronal_cells_per_celltype[1:length(neuronal_cells_per_celltype)])) +
  geom_histogram(binwidth = 0.1) +
  labs(x = "cells per neuron type", y = "neuron types", title = " ") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(size = 0.2),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 12, vjust = 16),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.4, "mm"),
    axis.line = element_line(linewidth = 0.2)
  ) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))
hist_neuro


hist_non_neur <- as.data.frame(non_neuronal_cells_per_celltype[1:length(non_neuronal_cells_per_celltype)]) %>%
  ggplot(aes(x = non_neuronal_cells_per_celltype[1:length(non_neuronal_cells_per_celltype)])) +
  geom_histogram(binwidth = 0.1) +
  labs(x = "cells per non-neuronal cell type", y = "cell types", title = " ") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(size = 0.2),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 12, vjust = 16),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.4, "mm"),
    axis.line = element_line(linewidth = 0.2)
  ) +
  scale_x_log10(breaks = c(1, 2, 10, 50, 1000)) +
  scale_y_continuous(breaks = c(0, 4, 8, 12))
hist_non_neur


# assemble figure ------------------------------

# define layout 
layout <- "
AAAA
BBBB
CCDE
FFFF
"

Figure3 <- SN_plot + IN_plot +
  MN_plot + hist_neuro + hist_non_neur + 
  non_n_cell_plot +
  plot_layout(design = layout, heights = c(
    1,1,1,1
  )) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(
    size = 12, face = "plain"
  ))


ggsave("Figures/Figure3.png",
  limitsize = FALSE,
  units = c("px"), Figure3, width = 4000, height = 3080, bg = "white"
)

ggsave("Figures/Figure3.pdf",
  limitsize = FALSE,
  units = c("px"), Figure3, width = 4000, height = 3080
)
