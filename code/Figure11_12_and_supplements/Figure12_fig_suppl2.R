# Code to generate Figure12-fig-suppl2 of the Platynereis 3d connectome paper
# Gaspar Jekely 2023
# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# load cell type connectivity ---------------

syn_tb <- readRDS("source_data/Figure4_source_data1.rds")
celltypes_table <- read_csv("supplements/Supplementary_Table1.csv")

syn_weights <- syn_tb %>%
  activate(edges) %>%
  pull(synapses)

syn_tb <- syn_tb %>%
  mutate(authority = centrality_authority(weights = syn_weights))

trunk_only_celltypes <- celltypes_table %>% 
  filter(!str_detect(`soma position`, "head")) %>%
  select(`name of cell type`) %>%
  pull()

# binarize and only consider >2 synapses
syn_tb_bin <- syn_tb %>%
  activate(edges) %>%
  mutate(syn_filtered = ifelse(synapses <= 3, 0, 1))

syn_tb_bin_summed <- syn_tb_bin %>%
  activate(edges) %>%
  select(from, syn_filtered) %>%
  group_by(from) %>%
  mutate(postsyn_partners = sum(syn_filtered)) %>%
  as_tibble() %>%
  ungroup() %>%
  group_by(to) %>%
  mutate(presyn_partners = sum(syn_filtered)) %>%
  ungroup()

#add names of from and to nodes
syn_tb_bin_summed_named <- syn_tb_bin_summed %>%
  mutate(presyn_name = syn_tb %>%
           activate(nodes) %>%
           select(name) %>%
           as_tibble() %>%
           slice(from) %>%
           pull()
  ) %>%
  mutate(postsyn_name = syn_tb %>%
           activate(nodes) %>%
           select(name) %>%
           as_tibble() %>%
           slice(to) %>%
           pull()
  ) %>%
  mutate(presyn_type = syn_tb %>%
           activate(nodes) %>%
           select(group) %>%
           as_tibble() %>%
           slice(from) %>%
           pull()
  ) %>%
  mutate(postsyn_type = syn_tb %>%
           activate(nodes) %>%
           select(group) %>%
           as_tibble() %>%
           slice(to) %>%
           pull()
  ) 

syn_tb_bin_summed_named_pre_trunk <- syn_tb_bin_summed_named %>%
  filter(presyn_name %in% trunk_only_celltypes)

syn_tb_bin_summed_named_post_trunk <- syn_tb_bin_summed_named %>%
  filter(postsyn_name %in% trunk_only_celltypes)


#plot theme
plot_theme <-  theme(
  text = element_text(size = 10),
  axis.text.y = element_text(angle = 90, hjust = 1, size = 10),
  axis.title = element_text(size = 12),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.length.y = unit(1, "mm"),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.key.size = unit(3, "mm"),
  legend.position = "top",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 8),
  axis.line = element_blank()
)

#filter cells to plot
presyn_ordered_names <- syn_tb_bin_summed_named_pre_trunk %>%
  select(presyn_name, presyn_type, postsyn_partners) %>%
  filter(postsyn_partners > 4) %>%
  arrange(desc(postsyn_partners)) %>%
  pull(presyn_name) %>%
  unique()

postsyn_ordered_names <- syn_tb_bin_summed_named_post_trunk %>%
  select(postsyn_name, postsyn_type, presyn_partners) %>%
  filter(presyn_partners > 4) %>%
  arrange(desc(presyn_partners)) %>%
  pull(postsyn_name) %>%
  unique()

# plot by number of postsyn cell types ----------
plot_pre_trunk <- syn_tb_bin_summed_named_pre_trunk %>%
  select(presyn_name, presyn_type, postsyn_partners) %>%
  unique() %>%
  ggplot(aes(presyn_name, postsyn_partners, fill = presyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "motoneuron" = "#0072B2",
    "interneuron" = "#CC79A7"
  )) +
  guides(fill = "none") +
  geom_col() +
  scale_x_discrete(limits = presyn_ordered_names) +
  labs(
    y = "# postsyn cell types", 
    x = "trunk neuronal cell types"
  ) +
  plot_theme +
  geom_text(aes(label = (presyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(limits = c(0, 30),
                     breaks = c(0,5,10,15,20,25))


# plot by number of presyn cell types -------

plot_post_trunk <- syn_tb_bin_summed_named_post_trunk %>%
  select(postsyn_name, postsyn_type, presyn_partners) %>%
  unique() %>%
  ggplot(aes(postsyn_name, presyn_partners, fill = postsyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "motoneuron" = "#0072B2",
    "interneuron" = "#CC79A7",
    "effector" = "#BBBBBB"
  ),
  labels = c(
    "sensory neuron" = "SN",
    "interneuron" = "IN",
    "motoneuron" = "MN",
    "effector" = "effector"
  )
  ) +
  geom_col() +
  scale_x_discrete(limits = postsyn_ordered_names) +
  labs(
    y = "# presyn neuron types", 
    x = "trunk cell types",
    fill = ""
  ) +
  plot_theme +
  geom_text(aes(label = (postsyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(limits = c(0, 30),
                     breaks = c(0,5,10,15))

# betweenness, pagerank, degree ------------

between_names <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(betweenness)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name) %>%
  pull()

pagerank_names <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(pagerank)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name) %>%
  pull()

authority_names <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(authority)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name) %>%
  pull()

weighted_degree_names <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(value)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name) %>%
  pull()

betweenness_plot <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(betweenness)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name, betweenness, group) %>%
  ggplot(aes(name, betweenness, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "motoneuron" = "#0072B2",
    "interneuron" = "#CC79A7"
  )) +
  scale_x_discrete(limits = between_names) +
  labs(
    y = "betweenness", 
    x = "trunk cell types",
    fill = ""
  ) +
  plot_theme +
  geom_text(aes(label = (between_names), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  guides(fill = "none") +
  scale_y_continuous(limits = c(0, 4500))
syn_tb
pagerank_plot <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(pagerank)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name, pagerank, group) %>%
  ggplot(aes(name, pagerank, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "motoneuron" = "#0072B2",
    "interneuron" = "#CC79A7",
    "effector" = "#BBBBBB"
  )) +
  scale_x_discrete(limits = pagerank_names) +
  labs(
    y = "pagerank", 
    x = "trunk cell types",
    fill = ""
  ) +
  plot_theme +
  geom_text(aes(label = (pagerank_names), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  guides(fill = "none") +
  scale_y_continuous(limits = c(0, 0.03))

weighted_degree_plot <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(value)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name, value, group) %>%
  ggplot(aes(name, value, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "motoneuron" = "#0072B2",
    "interneuron" = "#CC79A7"
  )) +
  scale_x_discrete(limits = weighted_degree_names) +
  labs(
    y = "weighted degree", 
    x = "trunk cell types",
    fill = ""
  ) +
  plot_theme +
  geom_text(aes(label = (weighted_degree_names), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  guides(fill = "none") +
  scale_y_continuous(limits = c(0, 110))


authority_plot <- syn_tb %>%
  filter(name %in% trunk_only_celltypes) %>%
  arrange(desc(authority)) %>%
  slice(1:28) %>%
  as_tibble() %>%
  select(name, authority, group) %>%
  ggplot(aes(name, authority, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "motoneuron" = "#0072B2",
    "interneuron" = "#CC79A7",
    "effector" = "#BBBBBB"
  )) +
  scale_x_discrete(limits = authority_names) +
  labs(
    y = "authority", 
    x = "trunk cell types",
    fill = ""
  ) +
  plot_theme +
  geom_text(aes(label = (authority_names), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  guides(fill = "none") +
  scale_y_continuous(limits = c(0, 0.0115))

authority_plot

# assemble figure ----------------

layout = "
ABC
DEF
"

Figure12_fig_suppl2 <- plot_pre_trunk + 
  plot_post_trunk +
  weighted_degree_plot + pagerank_plot + 
  betweenness_plot + authority_plot +
  plot_layout(
    design = layout
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure12_fig_suppl2.png",
       limitsize = FALSE,
       units = c("px"), Figure12_fig_suppl2,
       width = 4800, height = 2550, bg = "white"
)

ggsave("Figures/Figure12_fig_suppl2.pdf",
       limitsize = FALSE,
       units = c("px"), Figure12_fig_suppl2, 
       width = 4800, height = 2550
)

