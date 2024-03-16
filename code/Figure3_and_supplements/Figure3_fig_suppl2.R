# R code to generate Figure3 figure suppl 2 in the  Platynereis 3d connectome paper
# Gaspar Jekely Feb 2021

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/libraries_functions_and_CATMAID_conn.R")

# read all annotations for cell types -------------------------------------
annotation_neuronal_celltypelist <- list(length = 202)

# define empty cell type skids list
celltype_skids <- list()

for (i in c(1:202)) {
  annotation <- paste("annotation:^celltype", i, "$", sep = "")
  annotation_neuronal_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
  # use the sourced skids_by_2annotations function to retrieve all skids per cell type that are
  # also annotated with with_soma
  skids <- skids_by_2annotations(annotation, "with_soma")
  # if no skids are returned, next
  if (identical(skids, integer(0))) {
    print("this celltype is missing!")
  }
  celltype_skids[[i]] <- skids
}

# get neuron names and rename matrix col and row ---------------------------------

# we use the Catmaid neuron name before the "_" to parse the generic cell type name
# get the neuron name of the first skid in the skids list
neuronal_celltype_names <- list()

for (i in 1:length(celltype_skids)) {
  name <- catmaid_get_neuronnames(celltype_skids[[i]][1], pid = 11)
  name <- sub("_.*$", "", name)
  neuronal_celltype_names[i] <- name
}

# check duplicated names
neuronal_celltype_names[duplicated(neuronal_celltype_names)]

# check length
length(neuronal_celltype_names)

# make annotations table for 202 neuronal cell types ------------------------------


# create a tibble to store annotations per cell type
annot_tb <- tibble(celltype = character(), annotation = character(), present = numeric())

# define a list of anatomical annotations to search for
annot_to_search <- list(
  "ectoderm",
  "episphere",
  "segment_0", "segment_1", "segment_2", "segment_3", "torso",
  "pygidium",
  "neurosecretory_plexus",
  "mushroom body", "Apical_organ", "Dorsal_sensory_organ",
  "Dorsolateral_sense_organs", "antenna",
  "eyespot",
  "mechanosensory_girdle", "nuchal organ",
  "Sensory neuron", "sensory_motor_neuron",
  "sensory_neurosecretory_neuron",
  "interneuron", "premotor", "motorneuron", "ciliomotors",
  "rhabdomeric photoreceptor", "ciliary photoreceptor",
  "sensory_cilia",
  "Uniciliated_penetrating_cell", "Biciliated_penetrating_cell",
  "Uniciliated_nonpenetrating_cell", "Biciliated_nonpenetrating_cell",
  "Multiciliated_penetrating_cell",
  "dense cored vesicles",
  "asymmetric neuron",
  "biaxonal", "decussating", "commissural", "ipsilateral",
  "contralateral",
  "pseudounipolar", "descending", "ascending",
  "global_reach",
  "head_trunk", "siGOLD", "glutamatergic",
  "serotonergic", "cholinergic",
  "adrenergic", "dopaminergic"
)

annot_to_search[duplicated(annot_to_search)]

# iterate through cell type list and look for annotations - populate the annot_tb
for (i in 1:202) {
  for (j in 1:length(annot_to_search)) {
    match <- lapply(
      annotation_neuronal_celltypelist[[i]],
      function(ch) grep(annot_to_search[[j]], ch)
    )
    if (length(match$annotation) != 0) {
      print("hit")
      annot_tb <- annot_tb %>% add_row(
        celltype = paste("celltype", i, sep = ""),
        annotation = annot_to_search[[j]], present = 1
      )
    } else {
      print("no hit")
      annot_tb <- annot_tb %>% add_row(
        celltype = paste("celltype", i, sep = ""),
        annotation = annot_to_search[[j]], present = 0
      )
    }
  }
}

# convert to matrix for heatmap plotting
m <- matrix(as.numeric(as.matrix(pivot_wider(annot_tb,
  names_from = c(celltype),
  values_from = present
))[1:length(annot_to_search), 2:203]), nrow = length(annot_to_search))

dim(m)
# add row and column names
rownames(m) <- annot_to_search
colnames(m) <- neuronal_celltype_names


# plot heatmap of annotations and celltypes
hm_annot <- heatmaply(m,
  colors = c("white", Okabe_Ito[5]),
  column_text_angle = 90, row_text_angle = 0,
  dist_method = "manhattan",
  hclustfun = hclust,
  hclust_method = "ward.D2",
  dendrogram = "column",
  show_dendrogram = c(F, T),
  k_row = NA, k_col = 12,
  grid_size = 1.1,
  grid_gap = 0.1,
  hide_colorbar = T,
  plot_method = c("ggplot"),
  fontsize_row = 11,
  fontsize_col = 11,
  revC = F,
  grid_color = "grey90"
)

write.table(
  m, "source_data/Figure3_fig_suppl2_source_data1.txt", 
  sep = "\t")
m <- read.table("source_data/Figure3_fig_suppl2_source_data1.txt", 
                sep = "\t")

# save as html
saveNetwork(hm_annot, "pictures/Suppl_Heatmap_celltypes_annotations.html")

webshot2::webshot(
  url = "pictures/Suppl_Heatmap_celltypes_annotations.html",
  file = "Figures/Figure3_fig_suppl2.png",
  vwidth = 3000, vheight = 1250, # define the size of the browser window
  cliprect = c(10, 17, 3000, 1250), zoom = 2
)

webshot::webshot(
  url = "pictures/Suppl_Heatmap_celltypes_annotations.html",
  file = "Figures/Figure3_fig_suppl2.pdf",
  vwidth = 3000, vheight = 1250, # define the size of the browser window
  zoom = 1
)

