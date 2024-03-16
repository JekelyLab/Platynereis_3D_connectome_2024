# code to compare the connectivity of cell types on the left and right body sides of the 3 day old Platynereis larva
# Gaspar Jekely 2023

# load packages, functions and anatomical references
source("code/libraries_functions_and_CATMAID_conn.R")

# load all cells with celltype or celltype_non_neuronal annotation
celltypes <- catmaid_get_annotations_for_skeletons("^celltype$", pid = 11)

non_neuronal_celltypes <- catmaid_get_annotations_for_skeletons(
  "^celltype_non_neuronal$",
  pid = 11
)

# bind data.frames
all_celltypes <- rbind(celltypes, non_neuronal_celltypes)

# pull skids
skids <- as_tibble(all_celltypes) %>%
  select(skid) %>%
  pull()

# get neuron names (redundant but fast)
names <- catmaid_get_neuronnames(skids, pid = 11)

# assign names to skids
all_celltypes_tb <- as_tibble(all_celltypes) %>%
  mutate(name = names) %>%
  select(skid, annotation, name)

# retrieve connectivity
unique_skids <- unique(skids)
connectivity <- catmaid_get_connectors_between(
  pre = unique_skids,
  post = unique_skids, pid = 11
)

# tidy
connectivity_tb <- as_tibble(connectivity) %>%
  select(pre_skid, post_skid, connector_id)

left_skids <- all_celltypes_tb %>%
  filter(annotation == "left_side") %>%
  pull(skid)

right_skids <- all_celltypes_tb %>%
  filter(annotation == "right_side") %>%
  pull(skid)

mid_skids <- all_celltypes_tb %>%
  filter(annotation == "middle") %>%
  pull(skid)

# add side column
all_celltypes_tb_annot <- all_celltypes_tb %>%
  mutate(side = case_when(
    skid %in% left_skids ~ "left",
    skid %in% right_skids ~ "right",
    skid %in% mid_skids ~ "middle"
  ))

# add cell class column
SN_skids <- all_celltypes_tb %>%
  filter(annotation == "Sensory neuron") %>%
  pull(skid)

IN_skids <- all_celltypes_tb %>%
  filter(annotation == "interneuron") %>%
  pull(skid)

MN_skids <- all_celltypes_tb %>%
  filter(annotation == "motorneuron") %>%
  pull(skid)

# add class and colour column
all_celltypes_tb_annot <- all_celltypes_tb_annot %>%
  mutate(class = ifelse(
    skid %in% SN_skids, "sensory neuron",
    ifelse(skid %in% IN_skids, "interneuron",
      ifelse(skid %in% MN_skids, "motor neuron", "effector")
    )
  )) %>%
  mutate(celltype_color = ifelse(
    class == "sensory neuron", "#E69F00",
    ifelse(class == "interneuron", "#CC79A7",
      ifelse(class == "motor neuron", "#0072B2", "grey50")
    )
  ))

# generate connectivity tbl
left_right_connectivity <- tibble(
  presynaptic = character(),
  postsynaptic = character(),
  synapses = numeric()
)

SN_names <- c()
IN_names <- c()
MN_names <- c()
class <- c()

# query synapse connectivity for left or right to left+right by celltype
for (i in 1:202) {
  skip_outer <- FALSE

  for (j in 1:202) {
    left_skids <- all_celltypes_tb_annot %>%
      filter(annotation == paste("celltype", i, sep = "")) %>%
      filter(side == "left") %>%
      pull(skid)

    right_skids <- all_celltypes_tb_annot %>%
      filter(annotation == paste("celltype", i, sep = "")) %>%
      filter(side == "right") %>%
      pull(skid)

    left_right_skids <- all_celltypes_tb_annot %>%
      filter(annotation == paste("celltype", j, sep = "")) %>%
      filter(side == "left" | side == "right" | side == "middle") %>%
      pull(skid)

    # get neuron name by truncating at "_" the name of the first skid in each set
    left_name <- all_celltypes_tb_annot %>%
      filter(skid == left_skids[1]) %>%
      pull(name) %>%
      unique()
    left_name <- sub("_.*$", "", left_name)

    right_name <- all_celltypes_tb_annot %>%
      filter(skid == right_skids[1]) %>%
      pull(name) %>%
      unique()
    right_name <- sub("_.*$", "", right_name)

    left_right_name <- all_celltypes_tb_annot %>%
      filter(skid == left_right_skids[1]) %>%
      pull(name) %>%
      unique()
    left_right_name <- sub("_.*$", "", left_right_name)

    if (length(left_name) == 0 | length(right_name) == 0) {
      print(paste("no left-right pairs ", "celltype", i))
      skip_outer <- TRUE
      break
    }

    # number of synapses left
    num_syn_l <- length(
      connectivity_tb %>%
        filter(pre_skid %in% left_skids) %>%
        filter(post_skid %in% left_right_skids) %>%
        pull(connector_id)
    )

    new_row_l <- data.frame(
      presynaptic = paste(
        left_name, "left",
        sep = " "
      ),
      postsynaptic = left_right_name,
      synapses = num_syn_l
    )

    left_right_connectivity <- add_row(
      left_right_connectivity, new_row_l
    )

    # number of synapses right
    num_syn_r <- length(
      connectivity_tb %>%
        filter(pre_skid %in% right_skids) %>%
        filter(post_skid %in% left_right_skids) %>%
        pull(connector_id)
    )

    new_row_r <- data.frame(
      presynaptic = paste(
        right_name, "right",
        sep = " "
      ),
      postsynaptic = left_right_name,
      synapses = num_syn_r
    )

    left_right_connectivity <- add_row(
      left_right_connectivity, new_row_r
    )
  }

  # if celltype has no left or right rep (e.g. middle cells) skip
  if (skip_outer) {
    next # Skip the rest of the current iteration of the outer loop
  }

  # gather cell type names separately for SN IN MN
  class <- all_celltypes_tb_annot %>%
    filter(skid == left_skids[1]) %>%
    pull(class) %>%
    unique()

  if (class == "sensory neuron") {
    SN_names <- c(SN_names, left_name)
  } else if (class == "interneuron") {
    IN_names <- c(IN_names, left_name)
  } else if (class == "motor neuron") {
    MN_names <- c(MN_names, left_name)
  }
}

# connectivity with non_neuronal cells
effector_names <- c()

for (i in 1:202) {
  skip_outer <- FALSE

  for (j in 1:91) {
    left_skids <- all_celltypes_tb_annot %>%
      filter(annotation == paste("celltype", i, sep = "")) %>%
      filter(side == "left") %>%
      pull(skid)

    right_skids <- all_celltypes_tb_annot %>%
      filter(annotation == paste("celltype", i, sep = "")) %>%
      filter(side == "right") %>%
      pull(skid)

    left_right_skids <- all_celltypes_tb_annot %>%
      filter(annotation == paste("celltype_non_neuronal", j, sep = "")) %>%
      filter(side == "left" | side == "right" | side == "middle") %>%
      pull(skid)

    # get neuron name by truncating at "_" the name of the first skid in each set
    left_name <- all_celltypes_tb_annot %>%
      filter(skid == left_skids[1]) %>%
      pull(name) %>%
      unique()
    left_name <- sub("_.*$", "", left_name)

    right_name <- all_celltypes_tb_annot %>%
      filter(skid == right_skids[1]) %>%
      pull(name) %>%
      unique()
    right_name <- sub("_.*$", "", right_name)

    left_right_name <- all_celltypes_tb_annot %>%
      filter(skid == left_right_skids[1]) %>%
      pull(name) %>%
      unique()
    left_right_name <- sub("_.*$", "", left_right_name)

    if (length(left_name) == 0 | length(right_name) == 0) {
      print(paste("no left-right pairs ", "celltype", i))
      skip_outer <- TRUE
      break
    }

    # gather cell type names for effectors
    effector_names[j] <- left_right_name

    # number of synapses left
    num_syn_l <- length(
      connectivity_tb %>%
        filter(pre_skid %in% left_skids) %>%
        filter(post_skid %in% left_right_skids) %>%
        pull(connector_id)
    )

    new_row_l <- data.frame(
      presynaptic = paste(
        left_name, "left",
        sep = " "
      ),
      postsynaptic = left_right_name,
      synapses = num_syn_l
    )

    left_right_connectivity <- add_row(
      left_right_connectivity, new_row_l
    )

    # number of synapses right
    num_syn_r <- length(
      connectivity_tb %>%
        filter(pre_skid %in% right_skids) %>%
        filter(post_skid %in% left_right_skids) %>%
        pull(connector_id)
    )

    new_row_r <- data.frame(
      presynaptic = paste(
        right_name, "right",
        sep = " "
      ),
      postsynaptic = left_right_name,
      synapses = num_syn_r
    )

    left_right_connectivity <- add_row(
      left_right_connectivity, new_row_r
    )
  }

  # if celltype has no left or right rep (e.g. middle cells) skip
  if (skip_outer) {
    next # Skip the rest of the current iteration of the outer loop
  }
}

# define name vectors
ordered_names_pre <- c(SN_names, IN_names, MN_names)
ordered_names_post <- c(SN_names, IN_names, MN_names, effector_names)

# save connectivity ------------------------
write.table(left_right_connectivity, "source_data/Figure4_fig_suppl2_source_data1.txt", sep = "\t")
left_right_connectivity <- read.table("source_data/Figure4_fig_suppl2_source_data1.txt", sep = "\t")

# filter postsyn cells --------------------
no_syn_cells <- left_right_connectivity %>%
  group_by(postsynaptic) %>%
  mutate(total_syn = sum(synapses, na.rm = TRUE)) %>%
  filter(total_syn == 0) %>%
  select(postsynaptic) %>%
  unique() %>%
  pull()

no_presyn_cells <- left_right_connectivity %>%
  group_by(presynaptic) %>%
  mutate(total_syn = sum(synapses, na.rm = TRUE)) %>%
  filter(total_syn == 0) %>%
  select(presynaptic) %>%
  unique() %>%
  pull()

# remove postsyn cell names that receive no synapses
ordered_names_post_red <- ordered_names_post[
  !ordered_names_post %in% no_syn_cells
]

# define colours for axis text
presyn_colour <- all_celltypes_tb_annot %>%
  select(name, celltype_color) %>%
  mutate(celltype_name = sub("_.*$", "", name)) %>%
  filter(celltype_name %in% ordered_names_pre) %>%
  group_by(celltype_name) %>%
  slice(1) %>%
  arrange(match(celltype_name, ordered_names_pre)) %>%
  ungroup() %>%
  select(celltype_color) %>%
  pull()

postsyn_colour <- all_celltypes_tb_annot %>%
  select(name, celltype_color) %>%
  mutate(celltype_name = sub("_.*$", "", name)) %>%
  filter(celltype_name %in% ordered_names_post_red) %>%
  group_by(celltype_name) %>%
  slice(1) %>%
  arrange(match(celltype_name, ordered_names_post_red)) %>%
  ungroup() %>%
  select(celltype_color) %>%
  pull()

# plot ----------------

plot <- left_right_connectivity %>%
  ggplot(aes(
    x = postsynaptic, y = presynaptic, fill = sqrt(synapses)
  )) +
  geom_raster() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = 4,
      color = unlist(postsyn_colour)
    ),
    axis.text.y = element_text(
      angle = 0, hjust = 1, vjust = 0.5, size = 4,
      color = rev(unlist(presyn_colour))
    ),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.ticks = element_line(linewidth = 0.1)
  ) +
  scale_fill_gradientn(colours = c("white", "#0072B2", "#0072B2", "#0072B2", "#D55E00", "#D55E00", "#D55E00", "#D55E00", "#D55E00")) +
  guides(size = "none") +
  labs(
    x = "postsynaptic cell groups",
    y = "presynaptic cell groups",
    title = " ", fill = "sqrt\nsynapses"
  ) +
  scale_x_discrete(
    limits =
      unlist(ordered_names_post_red)
  ) +
  scale_y_discrete(
    limits =
      rev(c(
        paste(unlist(ordered_names_pre), "left", sep = " "),
        paste(unlist(ordered_names_pre), "right", sep = " ")
      ))
  )

# Saving R ggplot with R ggsave Function
ggsave("Figures/Figure4_fig_suppl2.png", plot,
  width = 4600,
  height = 6400, limitsize = TRUE,
  units = c("px")
)
ggsave("Figures/Figure4_fig_suppl2.pdf", plot,
  width = 4600,
  height = 6400, limitsize = TRUE,
  units = c("px")
)


# correlation calculations --------------------------------

left_conn <- left_right_connectivity %>%
  as_tibble() %>%
  filter(presynaptic %in% paste(ordered_names_pre, "left")) %>%
  filter(postsynaptic %in% ordered_names_post_red) %>%
  arrange(match(presynaptic, ordered_names_pre)) %>%
  pull(synapses)

left_names <- left_right_connectivity %>%
  as_tibble() %>%
  filter(presynaptic %in% paste(ordered_names_pre, "left")) %>%
  filter(postsynaptic %in% ordered_names_post_red) %>%
  arrange(match(presynaptic, ordered_names_pre)) %>%
  pull(presynaptic)

right_conn <- left_right_connectivity %>%
  as_tibble() %>%
  filter(presynaptic %in% paste(ordered_names_pre, "right")) %>%
  filter(postsynaptic %in% ordered_names_post_red) %>%
  arrange(match(presynaptic, ordered_names_pre)) %>%
  pull(synapses)

right_names <- left_right_connectivity %>%
  as_tibble() %>%
  filter(presynaptic %in% paste(ordered_names_pre, "right")) %>%
  filter(postsynaptic %in% ordered_names_post_red) %>%
  arrange(match(presynaptic, ordered_names_pre)) %>%
  pull(presynaptic)

left_matrix <- matrix(left_conn, ncol = length(ordered_names_pre), nrow = length(ordered_names_post_red))
dim(left_matrix)
right_matrix <- matrix(right_conn, ncol = length(ordered_names_pre), nrow = length(ordered_names_post_red))
dim(right_matrix)

# calculate pearson between between left and right matrix
pearson_left_right <- cor.test(right_matrix,
  left_matrix,
  method = "pearson"
)

# calculate pearson correlation for each left-right row pair (presyn pairs)
pearson <- matrix(nrow = ncol(left_matrix), ncol = ncol(right_matrix))

for (i in 1:(ncol(left_matrix))) {
  for (j in 1:(ncol(right_matrix))) {
    # calculate correlation, use $estimate to extract the correlation value from the cor.test output
    pearson[i, j] <- cor.test(left_matrix[i, ], right_matrix[j, ], method = "pearson")$estimate
  }
}

dim(pearson)

rownames(pearson) <- unique(left_names)
colnames(pearson) <- unique(right_names)
left_names
pearson[1:4, 1:5]

pearson_tb <- as.data.frame(pearson) %>%
  rownames_to_column(var = "left") %>%
  pivot_longer(-left, values_to = "correlation", names_to = "right") 

pearson_tb %>%
  ggplot(aes(left, right, fill = correlation)) +
  geom_tile() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = 4,
      color = rev(unlist(presyn_colour))
    ),
    axis.text.y = element_text(
      angle = 0, hjust = 1, vjust = 0.5, size = 4,
      color = rev(unlist(presyn_colour))
    ),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.ticks = element_line(linewidth = 0.1)
  ) +
  scale_fill_gradientn(
    colours = c("#0072B2", "white", "#D55E00"),
    limits = c(-1, 1),
    na.value = "white"
  ) +
  guides(size = "none") +
  labs(
    x = "left cells",
    y = "right cells",
    title = " ", fill = "correlation"
  ) +
  scale_x_discrete(
    limits =
      paste(unlist(ordered_names_pre), "left")
  ) +
  scale_y_discrete(
    limits = paste(unlist(ordered_names_pre), "right")
  )

# Saving R ggplot with R ggsave Function
ggsave("Figures/Figure4_fig_suppl3.png",
  width = 3600,
  height = 3600, limitsize = TRUE,
  units = c("px")
)
ggsave("Figures/Figure4_fig_suppl3.pdf",
  width = 4600,
  height = 4400, limitsize = TRUE,
  units = c("px")
)

# save left-right Pearson ------------------------
pearson_tb
write.table(pearson_tb, "source_data/Figure4_fig_suppl3_source_data1.txt", sep = "\t")
pearson_tb <- read.table("source_data/Figure4_fig_suppl3_source_data1.txt", sep = "\t")
