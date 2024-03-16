# Code to generate Figure10-fig-suppl1 of the Platynereis 3d connectome paper
# Gaspar Jekely 2023
# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# load cell type connectivity ---------------

syn_tb <- readRDS("source_data/Figure4_source_data1.rds")
celltypes_table <- read_csv("supplements/Supplementary_Table1.csv")

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

plot_pre_trunk_hist <- syn_tb_bin_summed_named_pre_trunk %>%
  select(presyn_name, postsyn_partners) %>%
  unique() %>%
  ggplot(aes(postsyn_partners)) +
  geom_histogram(binwidth = 1, color = "white") +
  labs(
    x = "# postsyn cell types", 
    y = "# trunk neuronal cell types"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(size = 10, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 12),
                     breaks = c(0, 2, 4, 6, 8, 10))  +
  scale_x_continuous(breaks = c(0, 5,10,24)) 

plot_pre_trunk_hist

# plot by number of presyn cell types -------

plot_post_trunk_hist <- syn_tb_bin_summed_named_post_trunk %>%
  select(postsyn_name, presyn_partners) %>%
  unique() %>%
  ggplot(aes(presyn_partners)) +
  geom_histogram(binwidth = 1 , color = "white") +
  labs(
    x = "# presyn cell types", 
    y = "# trunk neuronal cell types"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(size = 10,  hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 35),
                     breaks = c(0,5,10, 20, 30)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15)) 

plot_post_trunk_hist


# plot histogram of number of cells per trunk cell type -----
celltypes_table
hist_cells_per_trunk_type <- celltypes_table %>% 
  filter(`soma position` != "head") %>%
  filter(`Sensory/inter/motor neuron` == "Sensory neuron" |
           `Sensory/inter/motor neuron` == "Interneuron" |
           `Sensory/inter/motor neuron` == "Motorneuron") %>%
  ggplot(aes(`number of cells`)) +
  geom_histogram(binwidth = 1, color = 'white') +
  labs(y = "# trunk neuronal cell types", x = "# cells per cell type") +
  scale_x_continuous(breaks = c(1, 5, 10, 20, 30)) +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 90, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1)
  ) +
  scale_y_continuous(limits = c(0, 21),
                     breaks = c(0, 5, 10, 15, 20))

hist_cells_per_trunk_type


# assemble figure ----------------

layout = "
ABC
"

Figure12_fig_suppl1 <- 
  hist_cells_per_trunk_type +
  plot_post_trunk_hist + plot_pre_trunk_hist +
  plot_layout(
    design = layout, widths = c(30, 15, 24)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure12_fig_suppl1.png",
       limitsize = FALSE,
       units = c("px"), Figure12_fig_suppl1,
       width = 3600, height = 1250, bg = "white"
)

ggsave("Figures/Figure12_fig_suppl1.pdf",
       limitsize = FALSE,
       units = c("px"), Figure12_fig_suppl1, 
       width = 4800, height = 1250
)

