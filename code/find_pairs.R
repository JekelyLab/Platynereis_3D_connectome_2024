library(catmaid)
library(tidyverse)

source("~/R/conn.R")

source("code/Natverse_functions_and_conn.R")

# get all skids with annotations decussating, commissural, contralateral
skids_decussating <- catmaid_skids("annotation:^decussating$", pid=11)
skids_commissural <- catmaid_skids("annotation:^commissural$", pid=11)
skids_contralateral <- catmaid_skids("annotation:^contralateral$", pid=11)

skids_crossing <- unique(c(skids_decussating, skids_commissural, skids_contralateral))

# other variables that will be needed later
skids_celltype <-  catmaid_skids("annotation:celltype$", pid=11)
midline <- read.neuron.catmaid(2428987, pid = 11)

# check if it has synapses from/to known celltypes, throw out if it doesn't
# maybe the we could check those that don't in the future, but for now let's
# stick to more promising ones
connected_skids <- catmaid_get_connectors_between(skids_crossing,
                                                  skids_celltype,
                                                  pid=11) %>%
  select(pre_skid,
         post_skid) %>%
  unlist() %>%
  unname() %>%
  unique()

# use this new filtered subset for further searches
skids_to_search <- intersect(skids_crossing,
                             connected_skids)

# separate into SN, IN, MN subsets for easier searching

skids_SN <- catmaid_skids("annotation:^Sensory neuron$",
                          pid=11)
skids_IN <- catmaid_skids("annotation:^interneuron$",
                          pid=11)
skids_MN <- catmaid_skids("annotation:^motorneuron$",
                          pid=11)

skids_SN_to_search <- intersect(skids_to_search,
                                skids_SN)
skids_IN_to_search <- intersect(skids_to_search,
                                skids_IN)
skids_MN_to_search <- intersect(skids_to_search,
                                skids_MN)

# find point where it crosses the midline

for (k in seq_along(skids_to_search)) {
  print(skids_to_search[k])
  neuron1 <- read.neuron.catmaid(skids_to_search[k],
                                 pid=11)
  #distance <- 999999999999
  #for (i in seq_along(neuron1$d$PointNo)) {
  #  xyz1 <- c(neuron1$d$X[i],
  #            neuron1$d$Y[i],
  #            neuron1$d$Z[i])
  #  for (j in seq_along(midline$d$PointNo)) {
  #    xyz_mid <- c(midline$d$X[j],
  #                 midline$d$Y[j],
  #                 midline$d$Z[j])
  #    x <- rbind(xyz1,
  #               xyz_mid)
  #    dist <- stats::dist(x)
  #    if (dist < distance) {
  #      distance <<- dist
  #      node1 <<- neuron1$d$PointNo[i]
  #      node_mid <<- midline$d$PointNo[j]
  #      xyz <<- xyz1
  #    }
  #  }
  #}
  
  distance <- 999999999999
  for (i in seq_along(neuron1$d$PointNo)) {
    z1 <- neuron1$d$Z[i] - 150
    z2 <- neuron1$d$Z[i] + 150
    for (j in seq_along(midline$d$PointNo)) {
      if ((midline$d$Z[j] > z1) && (midline$d$Z[j] < z2)) {
        dist <- abs(neuron1$d$X[i] -  midline$d$X[j])
        if (dist < distance) {
          distance <<- dist
          node1 <<- neuron1$d$PointNo[i]
          xyz <<- c(neuron1$d$X[i],
                    neuron1$d$Y[i],
                    neuron1$d$Z[i])
        }
      }
    }
  }
  
  # find all other skids (from the same set of annotations) which cross within
  # distance x
  path <- paste("11/node/list?",
                "left=", xyz[1] - 175, "&",
                "top=", xyz[2] - 175, "&",
                "z1=", xyz[3] - 175, "&",
                "right=", xyz[1] + 175, "&",
                "bottom=", xyz[2] + 175, "&",
                "z2=", xyz[3] + 175, "&",
                "&with_relation_map=none",
                sep = "")
  print(path)
  neurons_surrounding <- catmaid_fetch(path = path)
  
  skids_bb <- vector()
  for (i in seq_along(neurons_surrounding[[1]])) {
    new_skid <- neurons_surrounding[[1]][[i]][[8]]
    skids_bb <- c(skids_bb,
                  new_skid)
  }

  skids_bb_filtered <- unique(skids_bb)
  # select only crossover skids
  #skids_bb_filtered <- intersect(skids_crossing,
  #                      skids_bb_filtered)
  # exclude original search neuron
  skids_bb_filtered <- setdiff(skids_bb_filtered,
                      neuron1)
  
  if (length(skids_bb_filtered) == 0) {
    next
  }
  
  # get skids on opposite side of our neuron
  annot_neuron1 <- catmaid_get_annotations_for_skeletons(skids_to_search[k],
                                                           pid=11)
  if (nrow(
    filter(annot_neuron1, annotation == "left_side"))
    >0 ) {
    opposite_side <- "right_side"
  } else if (nrow(
    filter(annot_neuron1, annotation == "right_side"))
    >0 ) {
    opposite_side <- "left_side"
  }
  
  skids_bb_contra <- catmaid_get_annotations_for_skeletons(skids_bb_filtered,
                                                           pid=11) %>%
    filter(annotation == opposite_side) %>%
    select(skid) %>%
    unname() %>%
    unlist() %>%
    unique()
    
  if (length(skids_bb_contra) == 0) {
    next
  }
  print("skids contra")
  print(skids_bb_contra)
  neurons_bb_contra <- read.neurons.catmaid(skids_bb_contra,
                                              pid=11)
    
  neuron1_name <- catmaid_get_neuronnames(skids_to_search[k], pid=11)
  print(neuron1_name)

  plot_background_ventral()
  plot3d(neuron1, soma=T, lwd=2, col = "red")
  plot3d(neurons_bb_contra, soma=T, lwd=2, col = "Blue")
  texts3d(100,100, 100, text = neuron1_name, cex = 4)
  
  snapshot_path <- paste("pictures/potential_pairs_of_",
                         neuron1_name,
                         ".png",
                         sep = "")
  rgl.snapshot(snapshot_path); close3d()
  # check if skid 2 has synapses to/from any of the same celltypes as skid 1
}
