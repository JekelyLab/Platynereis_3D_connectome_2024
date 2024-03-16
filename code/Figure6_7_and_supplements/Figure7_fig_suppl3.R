# Code to generate Figure7-figures supplement 3 of the Platynereis 3d connectome paper
# Gaspar Jekely 2023

# load packages and functions
source("code/libraries_functions_and_CATMAID_conn.R")

# load cell type connectivity ---------------

syn_tb <- readRDS("source_data/Figure4_source_data1.rds")
celltypes_table <- read_csv("supplements/Supplementary_Table1.csv")

head_only_celltypes <- celltypes_table %>% 
  filter(`soma position` == "head") %>%
  select(`name of cell type`) %>%
  pull()

# binarize and only consider >2 synapses
# add to edges the type of pre- and postsyn cell
syn_tb_bin <- syn_tb %>%
  activate(edges) %>%
  mutate(syn_filtered = ifelse(synapses <= 3, 0, 1)) %>%
  mutate(presyn_type = syn_tb %>% 
           as_tibble() %>% 
           slice(from) %>% 
           pull(group)
         ) %>%
  mutate(postsyn_type = syn_tb %>% 
           as_tibble() %>% 
           slice(to) %>% pull(group)
         )

#add names of from and to nodes
syn_tb_bin_named <- syn_tb_bin %>%
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

#add number of postsyn interneuron partners
syn_tb_bin_summed_named_pre_head <- syn_tb_bin_named %>%
  filter(presyn_name %in% head_only_celltypes) %>%
  activate(edges) %>%
  select(from, syn_filtered, presyn_type, postsyn_type, presyn_name, postsyn_name) %>%
  filter(postsyn_type == "interneuron") %>%
  group_by(from) %>%
  mutate(postsyn_IN_partners = sum(syn_filtered)) %>%
  as_tibble() %>%
  ungroup()

syn_tb_bin_summed_named_pre_head

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
presyn_ordered_names <- syn_tb_bin_summed_named_pre_head %>%
  select(presyn_name, presyn_type, postsyn_IN_partners, postsyn_type) %>%
  filter(postsyn_IN_partners > 1) %>%
  filter(presyn_type == "sensory neuron") %>%
  arrange(desc(postsyn_IN_partners)) %>%
  pull(presyn_name) %>%
  unique()
presyn_ordered_names

# plot by number of postsyn head IN cell types for head SNs ----------
plot_postIN_head <- syn_tb_bin_summed_named_pre_head %>%
  select(presyn_name, presyn_type, postsyn_IN_partners) %>%
  unique() %>%
  ggplot(aes(presyn_name, postsyn_IN_partners, fill = presyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "interneuron" = "#0072B2",
    "motoneuron" = "#CC79A7"
  )) +
  guides(fill = "none") +
  geom_col() +
  scale_x_discrete(limits = presyn_ordered_names) +
  labs(
    y = "# postsyn IN cell types", 
    x = "head SN cell types"
  ) +
  plot_theme +
  geom_text(aes(label = (presyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
    limits = c(0, max(syn_tb_bin_summed_named_pre_head %>%
                        select(postsyn_IN_partners)) +2
    )
  )
plot_postIN_head


# plot for presyn partners

#add number of postsyn interneuron partners
syn_tb_bin_named_summed_SN <- syn_tb_bin_named %>%
  filter(presyn_name %in% head_only_celltypes) %>%
  filter(postsyn_name %in% head_only_celltypes) %>%
  activate(edges) %>%
  select(to, syn_filtered, presyn_type, postsyn_type, presyn_name, postsyn_name) %>%
  filter(presyn_type == "sensory neuron") %>%
  group_by(to) %>%
  mutate(presyn_SN_partners = sum(syn_filtered)) %>%
  as_tibble() %>%
  ungroup()
syn_tb_bin_named_summed_SN

postsyn_ordered_names <- syn_tb_bin_named_summed_SN %>%
  select(postsyn_name, postsyn_type, presyn_SN_partners) %>%
  filter(presyn_SN_partners > 1) %>%
  arrange(desc(presyn_SN_partners)) %>%
  filter(postsyn_type == "interneuron" | postsyn_type == "motoneuron") %>%
  pull(postsyn_name) %>%
  unique()
postsyn_ordered_names
# plot by number of presyn head SN cell types for head INs ----------
plot_pre_SN_head <- syn_tb_bin_named_summed_SN %>%
  select(postsyn_name, postsyn_type, presyn_SN_partners) %>%
  unique() %>%
  ggplot(aes(postsyn_name, presyn_SN_partners, fill = postsyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "interneuron" = "#0072B2",
    "motoneuron" = "#CC79A7"
  )) +
  guides(fill = "none") +
  geom_col() +
  scale_x_discrete(limits = postsyn_ordered_names) +
  labs(
    y = "# presyn head SN cell types", 
    x = "head IN and MN cell types"
  ) +
  plot_theme +
  geom_text(aes(label = (postsyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 11)
                     )
  
plot_pre_SN_head

# IN to IN connections ------------------

#filter cells to plot
presyn_ordered_IN_names <- syn_tb_bin_summed_named_pre_head %>%
  select(presyn_name, presyn_type, postsyn_IN_partners, postsyn_type) %>%
  filter(postsyn_IN_partners > 1) %>%
  filter(presyn_type == "interneuron") %>%
  arrange(desc(postsyn_IN_partners)) %>%
  pull(presyn_name) %>%
  unique()
presyn_ordered_IN_names

# plot by number of postsyn head IN cell types for head SNs ----------
plot_IN_IN_post_head <- syn_tb_bin_summed_named_pre_head %>%
  select(presyn_name, presyn_type, postsyn_IN_partners) %>%
  unique() %>%
  ggplot(aes(presyn_name, postsyn_IN_partners, fill = presyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "interneuron" = "#0072B2",
    "motoneuron" = "#CC79A7"
  )) +
  guides(fill = "none") +
  geom_col() +
  scale_x_discrete(limits = presyn_ordered_IN_names) +
  labs(
    y = "# postsyn head IN cell types", 
    x = "head IN cell types"
  ) +
  plot_theme +
  geom_text(aes(label = (presyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     limits = c(0, 11
                     )
  )
plot_IN_IN_post_head


# IN to IN pre -----------------

syn_tb_bin_named_summed_IN_pre <- syn_tb_bin_named %>%
  activate(edges) %>%
  select(to, syn_filtered, presyn_type, postsyn_type, presyn_name, postsyn_name) %>%
  filter(presyn_type == "interneuron") %>%
  group_by(to) %>%
  mutate(presyn_IN_partners = sum(syn_filtered)) %>%
  as_tibble() %>%
  ungroup()
syn_tb_bin_named_summed_IN_pre

syn_tb_bin_named_summed_IN_pre_head <- syn_tb_bin_named_summed_IN_pre %>%
  filter(postsyn_name %in% head_only_celltypes)


presyn_ordered_IN_post_names <- syn_tb_bin_named_summed_IN_pre_head %>%
  select(postsyn_name, postsyn_type, presyn_IN_partners) %>%
  filter(presyn_IN_partners > 1) %>%
  filter(postsyn_type == "interneuron") %>%
  arrange(desc(presyn_IN_partners)) %>%
  pull(postsyn_name) %>%
  unique()
presyn_ordered_IN_post_names

plot_IN_IN_pre_head <- syn_tb_bin_named_summed_IN_pre_head %>%
  select(postsyn_name, postsyn_type, presyn_IN_partners) %>%
  unique() %>%
  ggplot(aes(postsyn_name, presyn_IN_partners, fill = postsyn_type)) +
  scale_fill_manual(values = c(
    "sensory neuron" = "#E69F00",
    "interneuron" = "#0072B2",
    "motoneuron" = "#CC79A7"
  )) +
  guides(fill = "none") +
  geom_col() +
  scale_x_discrete(limits = presyn_ordered_IN_post_names) +
  labs(
    y = "# presyn head IN cell types", 
    x = "head IN cell types"
  ) +
  plot_theme +
  geom_text(aes(label = (postsyn_name), angle = 90),
            cex = 3, col = "black",
            position = position_dodge2(width = 0.9), 
            vjust = 0.6, hjust = -0.1
  ) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     limits = c(0, 11)
  )
plot_IN_IN_pre_head

# anatomy ------------

#INW and presynaptic partners
{
INW_presyn <- syn_tb_bin_named_summed_IN_pre_head %>%
  filter(postsyn_name == "INW") %>%
  filter(syn_filtered == 1) %>%
  pull(presyn_name)

INW_presyn_annot <- sapply(INW_presyn, function(x){
  celltypes_table %>% 
    filter(`soma position` == "head") %>%
    filter(`name of cell type` == x) %>%
    pull(`CATMAID annotation`)
  }
)

INW <- nlapply(
  read.neurons.catmaid("^celltype117$", pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

INW_presyn_annot[1]

yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                           invertFaces = T, conn = NULL, pid = 11)
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                        function(x) smooth_neuron(x, sigma=6000))
scalebar_50um_anterior = read.neurons.catmaid("^scalebar_50um_anterior$", pid=11)

INW_pre1 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[1], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW_pre2 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[2], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW_pre3 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[3], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW_pre4 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[4], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW_pre5 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[5], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW_pre6 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[6], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW_pre7 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[7], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)
INW_pre8 <- nlapply(
  read.neurons.catmaid(INW_presyn_annot[8], pid = 11),
  function(x) smooth_neuron(x, sigma = 6000)
)

plot_background()
plot3d(
  INW, soma = T, lwd = 4,
  add = T, alpha = 1, col = "black"
)
texts3d(75000, 90000, 300, text = "INW", col = "black", cex = 3)

plot3d(
  INW_pre1, soma = T, lwd = 3,
  add = T, alpha = 1, col = bluepurple[3]
)
texts3d(56000, 40000, 300, text = INW_presyn[1], col = "black", cex = 2)

plot3d(
  INW_pre2, soma = T, lwd = 6,
  add = T, alpha = 1, col = bluepurple[5]
)
texts3d(90000, 34000, 300, text = INW_presyn[2], col = "black", cex = 2)
plot3d(
  INW_pre4, soma = T, lwd = 8,
  add = T, alpha = 1, col = bluepurple[3]
)
texts3d(28000, 65000, 3000, text = "INMB", col = "black", cex = 2)

plot3d(
  INW_pre5, soma = T, lwd = 4,
  add = T, alpha = 1, col = bluepurple[7]
)

plot3d(
  INW_pre6, soma = T, lwd = 3,
  add = T, alpha = 1, col = bluepurple[9]
)

plot3d(
  INW_pre7, soma = T, lwd = 2,
  add = T, alpha = 1, col = bluepurple[7]
)
texts3d(55000, 78000, 300, text = INW_presyn[7], col = "black", cex = 2)

plot3d(
  INW_pre8, soma = T, lwd = 5,
  add = T, alpha = 1, col = bluepurple[4]
)
texts3d(60000, 52000, 300, text = INW_presyn[8], col = "black", cex = 2)


rgl.snapshot("pictures/INW_presyn_partners.png")
close3d()

}

# INRGWa and presyn partners
{
  INRGWa_presyn <- syn_tb_bin_named_summed_SN %>%
    filter(postsyn_name == "INRGWa") %>%
    filter(syn_filtered == 1) %>%
    pull(presyn_name)
  
  INRGWa_presyn_annot <- sapply(INRGWa_presyn, function(x){
    celltypes_table %>% 
      filter(`soma position` == "head") %>%
      filter(`name of cell type` == x) %>%
      pull(`CATMAID annotation`)
  }
  )
  
  INRGWa <- nlapply(
    read.neurons.catmaid("^celltype6$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  INRGWa_presyn_annot
  
  INRGWa_pre1 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[1], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre2 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[2], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre3 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[3], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre4 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[4], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre5 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[5], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre6 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[6], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre7 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[7], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre8 <- nlapply(
    read.neurons.catmaid(INRGWa_presyn_annot[8], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
}
#plotting
{  
  plot_background()
  plot3d(
    INRGWa, soma = T, lwd = 2,
    add = T, alpha = 1, col = "black"
  )
  texts3d(57000, 53000, 300, text = "INRGWa", col = "black", cex = 3)
  
  plot3d(
    INRGWa_pre1, soma = T, lwd = 3,
    add = T, alpha = 1, col = oranges[3]
  )
  texts3d(96000, 32000, 300, text = INRGWa_presyn[1], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre2, soma = T, lwd = 1,
    add = T, alpha = 1, col = oranges[5]
  )
  texts3d(30000, 33000, 300, text = "SNnuch", col = "black", cex = 2)
  plot3d(
    INRGWa_pre4, soma = T, lwd = 8,
    add = T, alpha = 1, col = oranges[3]
  )
  texts3d(29000, 64000, 3000, text = INRGWa_presyn[4], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre5, soma = T, lwd = 4,
    add = T, alpha = 1, col = oranges[7]
  )
  texts3d(72000, 60000, 300, text = INRGWa_presyn[5], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre6, soma = T, lwd = 3,
    add = T, alpha = 1, col = oranges[9]
  )
  texts3d(62000, 20000, 300, text = INRGWa_presyn[6], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre7, soma = T, lwd = 2,
    add = T, alpha = 1, col = oranges[7]
  )
  texts3d(78000, 32000, 300, text = INRGWa_presyn[7], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre8, soma = T, lwd = 5,
    add = T, alpha = 1, col = oranges[4]
  )
  texts3d(84000, 40000, 300, text = INRGWa_presyn[8], col = "black", cex = 2)
  
  plot3d(
    scalebar_50um_anterior, lwd = 4,
    add = T, alpha = 1, col = "black"
  )
  texts3d(104000, 116000, 300, text = "50 μm", col = "black", cex = 2)
  
  rgl.snapshot("pictures/INRGWa_presyn_partners.png")
 
  
}
close3d() 


# INRGWa and postsyn partners
{
  INRGWa_postsyn <- syn_tb_bin_summed_named_pre_head %>%
    filter(presyn_name == "INRGWa") %>%
    filter(syn_filtered == 1) %>%
    pull(postsyn_name)
  
  INRGWa_postsyn_annot <- sapply(INRGWa_postsyn, function(x){
    celltypes_table %>% 
      filter(`soma position` == "head") %>%
      filter(`name of cell type` == x) %>%
      pull(`CATMAID annotation`)
  }
  )
  
  INRGWa <- nlapply(
    read.neurons.catmaid("^celltype6$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  INRGWa_postsyn_annot
  
  INRGWa_pre1 <- nlapply(
    read.neurons.catmaid(INRGWa_postsyn_annot[1], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre2 <- nlapply(
    read.neurons.catmaid(INRGWa_postsyn_annot[2], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre3 <- nlapply(
    read.neurons.catmaid(INRGWa_postsyn_annot[3], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre4 <- nlapply(
    read.neurons.catmaid(INRGWa_postsyn_annot[4], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre5 <- nlapply(
    read.neurons.catmaid(INRGWa_postsyn_annot[5], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre6 <- nlapply(
    read.neurons.catmaid(INRGWa_postsyn_annot[6], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  INRGWa_pre7 <- nlapply(
    read.neurons.catmaid(INRGWa_postsyn_annot[7], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
}  
#plotting
{  
  plot_background()
  plot3d(
    INRGWa, soma = T, lwd = 2,
    add = T, alpha = 1, col = "black"
  )
  texts3d(50000, 43000, 300, text = "INRGWa", col = "black", cex = 3)
  
  plot3d(
    INRGWa_pre2, soma = T, lwd = 1,
    add = T, alpha = 1, col = bluepurple[5]
  )
  texts3d(56000, 57000, 300, text = INRGWa_postsyn[2], col = "black", cex = 2)
  plot3d(
    INRGWa_pre3, soma = T, lwd = 5,
    add = T, alpha = 1, col = bluepurple[5]
  )
  texts3d(86000, 57000, 300, text = INRGWa_postsyn[3], col = "black", cex = 2)

  plot3d(
    INRGWa_pre4, soma = T, lwd = 8,
    add = T, alpha = 1, col = bluepurple[3]
  )
  texts3d(49000, 69000, 3000, text = INRGWa_postsyn[4], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre5, soma = T, lwd = 4,
    add = T, alpha = 1, col = bluepurple[7]
  )
  texts3d(48000, 85000, 300, text = INRGWa_postsyn[5], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre6, soma = T, lwd = 5,
    add = T, alpha = 1, col = bluepurple[9]
  )
  texts3d(94000, 40000, 300, text = INRGWa_postsyn[6], col = "black", cex = 2)
  
  plot3d(
    INRGWa_pre7, soma = T, lwd = 2,
    add = T, alpha = 1, col = bluepurple[7]
  )
  texts3d(78000, 82000, 300, text = INRGWa_postsyn[7], col = "black", cex = 2)
  
  rgl.snapshot("pictures/INRGWa_postsyn_partners.png")
}
close3d() 

# hCR partners

{
  hCR_presyn <- syn_tb_bin_summed_named_pre_head %>%
    filter(presyn_name == "hCR") %>%
    filter(syn_filtered == 1) %>%
    pull(postsyn_name)
  hCR_presyn
  hCR <- nlapply(
    read.neurons.catmaid("^celltype100$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  hCR_presyn_annot <- sapply(hCR_presyn, function(x){
    celltypes_table %>% 
      filter(`name of cell type` == x) %>%
      pull(`CATMAID annotation`)
  }
  )
  hCR_presyn_annot
  
  hCR_pre1 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[1], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  hCR_pre2 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[2], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  hCR_pre3 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[3], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  hCR_pre4 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[4], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  hCR_pre5 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[5], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  hCR_pre6 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[6], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  hCR_pre7 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[7], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  hCR_pre8 <- nlapply(
    read.neurons.catmaid(hCR_presyn_annot[8], pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )  
  
}

#plotting
{  
  plot_background_ventral_no_ac()
  par3d(zoom=0.4)
  par3d(windowRect = c(0, 0, 800, 800))
  nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
  
  plot3d(
    hCR, soma = T, lwd = 2,
    add = T, alpha = 1, col = "black"
  )
  texts3d(40000, 53000, 300, text = "hCR", col = "black", cex = 3)
  
  plot3d(
    hCR_pre1, soma = T, lwd = 3,
    add = T, alpha = 1, col = oranges[3]
  )
  texts3d(76000, 57000, 300, text = hCR_presyn[1], col = "black", cex = 2)
  
  plot3d(
    hCR_pre2, soma = T, lwd = 1,
    add = T, alpha = 1, col = oranges[5]
  )
  texts3d(62000, 63000, 300, text =  hCR_presyn[2], col = "black", cex = 2)
  plot3d(
    hCR_pre4, soma = T, lwd = 8,
    add = T, alpha = 1, col = oranges[3]
  )
  texts3d(64000, 120000, 3000, text = hCR_presyn[4], col = "black", cex = 2)
  
  plot3d(
    hCR_pre5, soma = T, lwd = 4,
    add = T, alpha = 1, col = oranges[7]
  )
  texts3d(33000, 115000, 300, text = hCR_presyn[5], col = "black", cex = 2)
  
  plot3d(
    hCR_pre6, soma = T, lwd = 3,
    add = T, alpha = 1, col = oranges[9]
  )
  texts3d(78000, 118000, 300, text = hCR_presyn[6], col = "black", cex = 2)
  
  plot3d(
    hCR_pre7, soma = T, lwd = 2,
    add = T, alpha = 1, col = oranges[7]
  )
  texts3d(117000, 110000, 300, text = hCR_presyn[7], col = "black", cex = 2)
  
  plot3d(
    hCR_pre8, soma = T, lwd = 5,
    add = T, alpha = 1, col = oranges[4]
  )
  
  rgl.snapshot("pictures/hCR_postsyn_partners.png")
  
}
close3d()  

# assemble figure -------------------

panel_hCR <- ggdraw() + draw_image(readPNG(
  "pictures/hCR_postsyn_partners.png"
  )) +
  draw_label(
    "hCR and postsynaptic partners",
    x = 0.5, y = 0.98, color = "black", size = 12
  )

panel_INRGWa_pre <- ggdraw() + draw_image(readPNG(
  "pictures/INRGWa_presyn_partners.png"
)) +
  draw_label(
    "INRGWa and presynaptic partners",
    x = 0.5, y = 0.98, color = "black", size = 12
  )

panel_INRGWa_post <- ggdraw() + draw_image(readPNG(
  "pictures/INRGWa_postsyn_partners.png"
)) +
  draw_label(
    "INRGWa and postsynaptic partners",
    x = 0.5, y = 0.98, color = "black", size = 12
  )
panel_INW_pre <- ggdraw() + draw_image(readPNG(
  "pictures/INW_presyn_partners.png"
)) +
  draw_label(
    "INW and presynaptic partners",
    x = 0.5, y = 0.98, color = "black", size = 12
  )

Fig7_fig_suppl3_top <- plot_grid(
  plot_postIN_head, plot_pre_SN_head, plot_IN_IN_post_head, plot_IN_IN_pre_head,
  labels = c("A", "B", "C", "D"), ncol = 4, 
  rel_widths = c(13, 18, 28, 25),
  label_size = 12, label_y = 1, label_x = 0,
  label_fontfamily = "sans", label_fontface = "plain"
  )

layout <- "
AAAA
####
BCDE
"

Fig7_fig_suppl3 <- Fig7_fig_suppl3_top +
  panel_hCR + panel_INRGWa_pre +  panel_INRGWa_post + panel_INW_pre +
  plot_layout(design = layout, heights = c(0.7, 0.02, 1)) +
  plot_annotation(tag_levels = list(c("", 
                                      "E", "F", "G", "H"))) &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure7_fig_suppl3.png",
       limitsize = FALSE,
       units = c("px"), Fig7_fig_suppl3,
       width = 4600, height = 2400, bg = "white"
)

ggsave("Figures/Figure7_fig_suppl3.pdf",
       limitsize = FALSE,
       units = c("px"), Fig7_fig_suppl3, 
       width = 4600, height = 2400
)
