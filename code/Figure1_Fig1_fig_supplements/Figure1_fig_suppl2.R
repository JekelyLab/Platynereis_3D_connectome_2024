# R code to generate the Figure 1 figure supplement 2 of the Platynereis connectome paper
# Gaspar Jekely 2023 July

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# read neuron groups ------------------------------------------------------

SN <- nlapply(
  read.neurons.catmaid("^connectome_Sensory_neuron$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
IN <- nlapply(
  read.neurons.catmaid("^connectome_interneuron$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
MN <- nlapply(
  read.neurons.catmaid("^connectome_motorneuron$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

frag_skids <- as_tibble(catmaid_get_annotations_for_skeletons(
  "^fragmentum$", pid = 11
)) %>%
  select(skid) %>%
  unique() %>%
  pull()

#all skeletons
with_soma_skids <- as_tibble(catmaid_get_annotations_for_skeletons(
  "^with_soma$", pid = 11
)) %>%
  select(skid) %>%
  unique() %>%
  pull()

with_soma_node_count <- tibble(skid = with_soma_skids,
                          count = catmaid_get_node_count("^with_soma$", pid = 11)
)

#nodes on with_soma skeletons
total_with_soma_nodes <- with_soma_node_count %>%
  select(count) %>%
  sum()

#synapses on with_soma skeletons
with_soma_connectivity <- catmaid_query_connected(with_soma_skids, pid = 11)
num_presyn_with_soma <- sum(with_soma_connectivity$outgoing$syn.count)
num_postsyn_with_soma <- sum(with_soma_connectivity$incoming$syn.count)

#number of fragment
length(frag_skids)

frag_node_count <- tibble(skid = frag_skids,
  count = catmaid_get_node_count("^fragmentum$", pid = 11)
)

#total node count for fragments
total_frag_nodes <- frag_node_count %>%
  select(count) %>%
  sum()
total_frag_nodes

# fraction of nodes on fragments
total_frag_nodes/total_with_soma_nodes


#remove fragments of 1 node
filtered_frag_skids <- frag_node_count %>%
  filter(count>1) %>%
  pull(skid)

#single-node fragments
point_frag_skids <- frag_node_count %>%
  filter(count == 1) %>%
  pull(skid)

frag <- nlapply(
  read.neurons.catmaid(filtered_frag_skids, pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

# Prune twigs up to 2 microns long
SN_pr <- nlapply(SN, function(x) (prune_twigs(x, twig_length = 2000)))
IN_pr <- nlapply(IN, function(x) (prune_twigs(x, twig_length = 2000)))
MN_pr <- nlapply(MN, function(x) (prune_twigs(x, twig_length = 2000)))
#no pruning of fragments, many are too short

# number of presynaptic sites (total minus postsyn sites)
n_presyn_SN <- as.integer(sapply(SN, function(x) length(x$connectors$prepost) - sum(x$connectors$prepost)))
n_presyn_IN <- as.integer(sapply(IN, function(x) length(x$connectors$prepost) - sum(x$connectors$prepost)))
n_presyn_MN <- as.integer(sapply(MN, function(x) length(x$connectors$prepost) - sum(x$connectors$prepost)))
n_presyn_frag <- as.integer(sapply(frag, function(x) length(x$connectors$prepost) - sum(x$connectors$prepost)))

# number of postsynaptic sites (prepost == 1)
n_postsyn_SN <- as.integer(unlist(sapply(SN, function(x) sum(x$connectors$prepost))))
n_postsyn_IN <- as.integer(unlist(sapply(IN, function(x) sum(x$connectors$prepost))))
n_postsyn_MN <- as.integer(unlist(sapply(MN, function(x) sum(x$connectors$prepost))))
n_postsyn_frag <- as.integer(unlist(sapply(frag, function(x) sum(x$connectors$prepost))))

#synapses on fragments
point_frag_skids_partners <- catmaid_query_connected(point_frag_skids, pid = 11)

sum(n_presyn_frag)+sum(point_frag_skids_partners$outgoing$syn.count)
sum(n_postsyn_frag)+sum(point_frag_skids_partners$incoming$syn.count)

#ratio of synapses on with_soma and fragment skeletons
(sum(n_presyn_frag)+sum(point_frag_skids_partners$outgoing$syn.count))/num_presyn_with_soma
(sum(n_postsyn_frag)+sum(point_frag_skids_partners$incoming$syn.count))/num_postsyn_with_soma

# Summary statistics for neurons (e.g. cable length, number of nodes)
summ_SN <- summary(SN_pr)
summ_IN <- summary(IN_pr)
summ_MN <- summary(MN_pr)
summ_frag <- summary(frag)

# convert to tibble
summ_SN_tb <- summ_SN %>%
  rownames_to_column(var = "skid") %>%
  mutate(neuron_type = "sensory neuron") %>%
  mutate(presyn_sites = n_presyn_SN) %>%
  mutate(postsyn_sites = n_postsyn_SN)

summ_IN_tb <- summ_IN %>%
  rownames_to_column(var = "skid") %>%
  mutate(neuron_type = "interneuron") %>%
  mutate(presyn_sites = n_presyn_IN) %>%
  mutate(postsyn_sites = n_postsyn_IN)

summ_MN_tb <- summ_MN %>%
  rownames_to_column(var = "skid") %>%
  mutate(neuron_type = "motoneuron") %>%
  mutate(presyn_sites = n_presyn_MN) %>%
  mutate(postsyn_sites = n_postsyn_MN)

summ_frag_tb <- summ_frag %>%
  rownames_to_column(var = "skid") %>%
  mutate(neuron_type = "fragment") %>%
  mutate(presyn_sites = n_presyn_frag) %>%
  mutate(postsyn_sites = n_postsyn_frag)


# join tibbles
SN_IN <- full_join(summ_SN_tb, summ_IN_tb)
summ_IN_MN <- full_join(SN_IN, summ_MN_tb)
summ_All <- full_join(summ_IN_MN, summ_frag_tb)

skids <- summ_All %>%
  select(skid) %>%
  pull()

names <- as.character(catmaid_get_neuronnames(skids, pid = 11))

summ_All <- summ_All %>%
  mutate(neuron_names = names)

# save tibble
write_csv(summ_All, "source_data/Figure1_fig_suppl2_source_data_1.txt")
summ_All <- read_csv("source_data/Figure1_fig_suppl2_source_data_1.txt", show_col_types = FALSE)



# define plot themes ------------------------------------------------------

theme_h <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.key.size = unit(3, "mm")
  )

theme_p <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(3, "mm")
  )

# generate plots ----------------------------------------------------------

# plot histograms of presyn sites

h2 <- summ_All %>%
  ggplot() +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "sensory neuron"),
    aes(presyn_sites),
    fill = "#E69F00", alpha = 1, bins = 40,
    color = "grey", linewidth = 0.02
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "interneuron"),
    aes(presyn_sites),
    fill = "#CC79A7", alpha = 0.5, bins = 40,
    color = "black", linewidth = 0.02
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "motoneuron"),
    aes(presyn_sites),
    fill = "#0072B2", alpha = 0.8, bins = 40
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(1, 400)
  ) +
  labs(x = "# presynapses", y = "count") +
  theme_h

h2

# plot histograms of postsyn sites

h3 <- summ_All %>%
  ggplot() +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "sensory neuron"),
    aes(postsyn_sites),
    fill = "#E69F00", alpha = 1, bins = 40,
    color = "grey", linewidth = 0.02
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "interneuron"),
    aes(postsyn_sites),
    fill = "#CC79A7", alpha = 0.5, bins = 40,
    color = "black", linewidth = 0.02
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "motoneuron"),
    aes(postsyn_sites),
    fill = "#0072B2", alpha = 0.6, bins = 40
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(1, 400)
  ) +
  labs(x = "# postsynapses", y = "count") +
  theme_h

h3

# plot histograms of cable length sites

h4 <- summ_All %>%
  ggplot() +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "sensory neuron"),
    aes(cable.length / 1000),
    fill = "#E69F00", alpha = 1, bins = 40,
    color = "grey", linewidth = 0.1
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "interneuron"),
    aes(cable.length / 1000),
    fill = "#CC79A7", alpha = 0.5, bins = 40,
    color = "black", linewidth = 0.02
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "motoneuron"),
    aes(cable.length / 1000),
    fill = "#0072B2", alpha = 0.6, bins = 40
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(10, 2000)
  ) +
  labs(x = "cable length [µm]", y = "count") +
  theme_h

h4


h5 <- summ_All %>%
  ggplot() +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "sensory neuron"),
    aes((postsyn_sites - presyn_sites) / (postsyn_sites + presyn_sites)),
    fill = "#E69F00", alpha = 1, bins = 57,
    color = "grey", linewidth = 0.1
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "interneuron"),
    aes((postsyn_sites - presyn_sites) / (postsyn_sites + presyn_sites)),
    fill = "#CC79A7", alpha = 0.5, bins = 60,
    color = "black", linewidth = 0.02
  ) +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "motoneuron"),
    aes((postsyn_sites - presyn_sites) / (postsyn_sites + presyn_sites)),
    fill = "#0072B2", alpha = 0.6, bins = 55
  ) +
  scale_x_continuous(
    trans = "identity",
    limits = c(-1, 1)
  ) +
  labs(x = "(I-O)/(I+O)", y = "count") +
  theme_h

h5


# plot of presyn sites over cable length

p2 <- summ_All %>%
  filter(neuron_type != "fragment") %>%
  ggplot(aes(cable.length / 1000,
    presyn_sites,
    colour = neuron_type, alpha = neuron_type,
    shape = neuron_type, size = postsyn_sites
  )) +
  geom_jitter(stroke = 0, width = 0, height = 0.01) +
  scale_color_manual(values = list(
    interneuron = "#CC79A7", motoneuron = "#0072B2",
    `sensory neuron` = "#E69F00"
  )) +
  scale_alpha_manual(values = c(0.5, 0.8, 1)) +
  scale_size_area(max_size = 4) +
  labs(x = "cable length [µm]", y = "presynapses") +
  scale_x_continuous(
    trans = "log10",
    limits = c(50, 2000)
  ) + # change "log10" to "identity" to remove log scale
  theme_p +
  geom_text(
    data = summ_All %>% 
      filter(neuron_type != "fragment") %>%
      filter(presyn_sites > 100), # Filter data first
    aes(label = neuron_names), size = 2.5, alpha = 0.7, nudge_x = -0.1, check_overlap = TRUE, col = "black"
  ) + # Apply guides function
  guides(
    size = guide_legend("postsynapses"), shape = guide_legend("neuron type"),
    alpha = guide_legend("neuron type"),
    colour = guide_legend(override.aes = list(size = 4), "neuron type")
  ) +
  guides(colour = "none", shape = "none", alpha = "none")

p2

# plot of presyn sites over segment number
p3 <- summ_All %>%
  filter(neuron_type != "fragment") %>%
  ggplot(aes(cable.length / 1000,
    postsyn_sites,
    colour = neuron_type, alpha = neuron_type,
    shape = neuron_type, size = presyn_sites
  )) +
  geom_jitter(stroke = 0, width = 0, height = 0.01) +
  scale_color_manual(values = list(
    interneuron = "#CC79A7", motoneuron = "#0072B2",
    `sensory neuron` = "#E69F00"
  )) +
  scale_alpha_manual(values = c(0.5, 0.8, 1)) +
  scale_size_area(max_size = 4) +
  labs(x = "cable length [µm]", y = "postsynapses") +
  scale_x_continuous(
    trans = "log10",
    limits = c(50, 2000)
  ) + # change "log10" to "identity" to remove log scale
  theme_p +
  geom_text(
    data = summ_All %>% 
      filter(neuron_type != "fragment") %>%
      filter(presyn_sites > 100), # Filter data first
    aes(label = neuron_names), size = 2.5, alpha = 0.7, nudge_x = -0.1, check_overlap = TRUE, col = "black"
  ) + # Apply guides function
  guides(
    size = guide_legend("presynapses"), shape = guide_legend("neuron type"),
    alpha = guide_legend("neuron type"),
    colour = guide_legend(override.aes = list(size = 4), "neuron type")
  ) +
  guides(colour = "none", shape = "none", alpha = "none")

p3

# plot of cable length  over segment number

p4 <- summ_All %>%
  filter(neuron_type != "fragment") %>%
  ggplot(aes(segments,
    cable.length / 1000,
    colour = neuron_type, alpha = neuron_type,
    shape = neuron_type, size = presyn_sites
  )) +
  geom_jitter(stroke = 0, width = 0, height = 0.01) +
  scale_color_manual(values = list(
    interneuron = "#CC79A7", motoneuron = "#0072B2",
    `sensory neuron` = "#E69F00"
  )) +
  scale_alpha_manual(values = c(0.5, 0.8, 1)) +
  scale_size_area(max_size = 4) +
  labs(x = "# segments", y = "cable length [µm]") +
  scale_x_continuous(
    trans = "log10",
    limits = c(1, 200)
  ) + # change "log10" to "identity" to remove log scale
  scale_y_continuous(
    limits = c(0, 1500)
  ) +
  theme_p +
  geom_text(
    data = summ_All %>% 
      filter(neuron_type != "fragment") %>%
      filter(presyn_sites > 100), # Filter data first
    aes(label = neuron_names), size = 2.5, alpha = 0.7, nudge_x = -0.1, check_overlap = TRUE, col = "black"
  ) + # Apply guides function
  guides(
    size = guide_legend("# segments"),
    shape = guide_legend("neuron type"),
    colour = guide_legend(override.aes = list(size = 4), "neuron type")
  ) +
  guides(alpha = "none")

p4

# plot fragment statistics

h_frag <- summ_All %>%
  ggplot() +
  geom_histogram(
    data = summ_All %>% filter(neuron_type == "fragment"),
    aes(cable.length / 1000),
    fill = "grey", alpha = 1, bins = 40,
    color = "black", linewidth = 0.05
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(0.01, 1000),
    breaks = c(1, 10, 100, 1000)
  ) +
  labs(x = "cable length [µm]", y = "count") +
  theme_h

h_frag

#median cable length of fragments
summ_All %>%
  filter(neuron_type == "fragment") %>%
  select(cable.length) %>%
  pull() %>% median()/1000

# calculate correlations ------------------------------------------------------------
head(summ_All)

seg <- summ_All %>%
  select(segments) %>%
  pull()
pre <- summ_All %>%
  select(presyn_sites) %>%
  pull()
post <- summ_All %>%
  select(postsyn_sites) %>%
  pull()
cable <- summ_All %>%
  select(cable.length) %>%
  pull()


cor.test(pre, post)

pre_seg <- cor.test(pre, seg)
cor.test(post, seg)
pre_cable <- cor.test(pre, cable)
signif(as.numeric(pre_cable$estimate), 2)
post_cable <- cor.test(post, cable)

cor.test(post, cable)
seg_cable <- cor.test(seg, cable)


# assemble figure ---------------------------------------------------------

h_frag_lab <- ggdraw(h_frag) + draw_label(
  "fragments",
  x = 0.3, y = 0.98, size = 11
)

top_panels <- plot_grid(h4, h_frag_lab, h3, h2, h5,
                        rel_widths = c(1, 1, 1, 1, 1),
                        ncol = 5,
                        labels = c("A", "B", "C", "D", "E"),
                        label_size = 12, 
                        label_fontface = "plain"
                        )


p2_labelled <- ggdraw(p2) + draw_label(
  paste("Pearson's:", signif(as.numeric(pre_cable$estimate), 2)),
  x = 0.3, y = 0.98, size = 11
)


p3_labelled <- ggdraw(p3) + 
  draw_label(
    paste("Pearson's:", signif(as.numeric(post_cable$estimate), 2)),
    x = 0.3, y = 0.98, size = 11
  )

p4_labelled <- ggdraw(p4) +
  draw_label(
    paste("Pearson's:", signif(as.numeric(seg_cable$estimate), 2)),
    x = 0.3, y = 0.98, size = 11
  )



bottom_panels <- plot_grid(p2_labelled, p3_labelled, p4_labelled,
                        rel_widths = c(1, 1, 1),
                        ncol = 3,
                        labels = c("F", "G", "H"),
                        label_size = 12, 
                        label_fontface = "plain"
)

layout <- "
A
#
B
"

Fig1_fig_suppl2 <- top_panels +
  bottom_panels +
  plot_layout(design = layout, heights = c(1, 0.02, 1.2),
              widths = c(1, 1)) &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure1_fig_suppl2.png",
  limitsize = FALSE,
  units = c("px"), Fig1_fig_suppl2,
  width = 4000, height = 1800
)

ggsave("Figures/Figure1_fig_suppl2.pdf",
  limitsize = FALSE,
  units = c("px"), Fig1_fig_suppl2,
  width = 4000, height = 1800
)
