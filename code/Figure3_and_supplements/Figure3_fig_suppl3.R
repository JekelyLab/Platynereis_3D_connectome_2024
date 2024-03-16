# R code to generate Fig3 fig suppl3 of the Platynereis connectome paper
# Gaspar Jekely 2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/libraries_functions_and_CATMAID_conn.R")

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


# Apply the function using two nested lapply loops
results <- sapply(annot_to_search, function(ann1) {
  sapply(annot_to_search, function(ann2) {
    length(skids_by_2annotations(ann1, ann2))
  })
})

# add row and col names
colnames(results) <- annot_to_search
rownames(results) <- annot_to_search
results

# convert to tibble
results.tb <- results %>%
  as_tibble(rownames = "annotation1") %>%
  pivot_longer(
    -annotation1,
    names_to = "annotation2",
    values_to = "number"
  )

results.tb

write.table(results.tb, "source_data/Figure3_fig_suppl3_source_data1.txt", sep = "\t")
results.tb <- read.table("source_data/Figure3_fig_suppl3_source_data1.txt", sep = "\t")


results.tb %>%
  ggplot(aes(annotation2, annotation1)) +
  geom_raster(aes(fill = sqrt(number))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 10),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  labs(x = "annotation 2", y = "annotation 1", title = " ") +
  scale_x_discrete(limits = unlist(annot_to_search)) +
  scale_y_discrete(limits = rev(unlist(annot_to_search))) +
  scale_fill_gradientn(colours = c("white", "#0072B2")) +
  geom_text(aes(label = number, size = number / (number + 0.1))) +
  scale_radius(range = c(0, 2)) +
  guides(size = "none")


# Saving R ggplot with R ggsave Function
ggsave("Figures/Figure3_fig_suppl3.png",
  width = 4500,
  height = 3400, limitsize = TRUE,
  units = c("px")
)

# Saving R ggplot with R ggsave Function
ggsave("Figures/Figure3_fig_suppl3.pdf",
  width = length(annot_to_search) / 1.3,
  height = length(annot_to_search) / 1.6, limitsize = TRUE,
  units = c("cm")
)
