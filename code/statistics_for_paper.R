# statistics for paper - to be sourced and inserted into the text
# Gaspar Jekely 2023

library(tidygraph)
library(dplyr)
library(tibble)
library(igraph)
library(catmaid)
library(readr)
#source("code/CATMAID_connection.R")

# read graph -------------------------
connect.tb <- readRDS("supplements/connectome_graph_tibble.rds")

n_nodes_connectome <- connect.tb %>%
  activate(nodes) %>%
  pull(skids) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")

n_fragments <- connect.tb %>%
  activate(nodes) %>%
  filter(class == "fragmentum") %>%
  pull(skids) %>%
  length() %>%
  format(scientific = FALSE, big.mark = ",")

n_edges <- connect.tb %>%
  activate(edges) %>%
  pull(from) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")

n_synapses <- connect.tb %>%
  activate(edges) %>%
  pull(weight) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")

density <- graph.density(as.igraph(connect.tb), loop = FALSE)
density <- round(density, 5)


connect.tb %>%
  filter(class == "Sensory neuron" |
           class == "interneuron" |
           class == "motorneuron") %>%
  pull(skids) %>% length()

# cell type graph -------------------
  
Platy_celltype_graph <- readRDS("source_data/Figure4_source_data1.rds")

Platy_celltype_graph_n_synapses <- Platy_celltype_graph %>%
  activate(edges) %>%
  pull(synapses) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")
Platy_celltype_graph_n_nodes <- Platy_celltype_graph %>%
  activate(nodes) %>%
  pull(name) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")
Platy_celltype_graph_n_edges <- Platy_celltype_graph %>%
  activate(edges)%>%
  pull(synapses) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")
# number of skeletons stats -----------------

frag_all_annot <- read.table("data/frag_all_annot.txt", sep = "\t")
with_soma_all_annot <- read.table("data/with_soma_all_annot.txt", sep = "\t")
with_soma_node_count <- read.table("data/with_soma_node_count.txt", sep = "\t")
fragment_node_count <-  read.table("data/fragment_node_count.txt", sep = "\t")

num_soma <- with_soma_all_annot %>%
  pull(skid) %>% unique() %>%
  length() %>%
  format(scientific = FALSE, big.mark = ",")

num_frag <- frag_all_annot %>%
  pull(skid) %>% unique() %>%
  length() %>%
  format(scientific = FALSE, big.mark = ",")

with_soma_node_count_sum <- sum(with_soma_node_count) %>%
  format(scientific = FALSE, big.mark = ",")
fragment_node_count_sum <- sum(fragment_node_count) %>%
  format(scientific = FALSE, big.mark = ",")

fraction_fragment <- round(sum(fragment_node_count)/sum(with_soma_node_count)*100, 1)

n_pigment <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "pigment cell") %>%
  pull(skid) %>% unique() %>%
  length()

n_gland <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "gland cell") %>%
  pull(skid) %>% unique() %>%
  length()

n_glia <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "glia cell") %>%
  pull(skid) %>% unique() %>%
  length()

n_follicle <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "follicle cell") %>%
  pull(skid) %>% unique() %>%
  length()

n_MUS <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "muscle") %>%
  pull(skid) %>% unique() %>%
  length()

n_ciliated <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "ciliated cell") %>%
  pull(skid) %>% unique() %>%
  length()

n_epidermis <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "epithelia_cell") %>%
  pull(skid) %>% unique() %>%
  length()
n_chaeta <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "celltype_non_neuronal22") %>%
  pull(skid) %>% unique() %>%
  length()
n_aciculae <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "celltype_non_neuronal23") %>%
  pull(skid) %>% unique() %>%
  length()

n_pnb <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "pnb") %>%
  pull(skid) %>% unique() %>%
  length()

n_antenna <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "antennae_cell") %>%
  pull(skid) %>% unique() %>%
  length()

n_palp <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "palp sensory neuron") %>%
  pull(skid) %>% unique() %>%
  length()

n_stomodeum <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "stomodeum") %>%
  pull(skid) %>% unique() %>%
  length()

n_SN_stomodeum <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "SN stomodeum") %>%
  pull(skid) %>% unique() %>%
  length()

n_dividing <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "dividing") %>%
  pull(skid) %>% unique() %>%
  length()


head_trunk_skids <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "head_trunk") %>%
  pull(skid) %>% unique() 

n_head_trunk_asc <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(skid %in% head_trunk_skids) %>%
  filter(annotation == "torso") %>%
  pull(skid) %>% unique() %>%
  length()

n_head_trunk_desc <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(skid %in% head_trunk_skids) %>%
  filter(annotation == "episphere") %>%
  pull(skid) %>% unique() %>%
  length()


neuro_celltype_skids <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "celltype") %>%
  pull(skid)

n_neuro_celltypes <- length(neuro_celltype_skids)

with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "immature neuron")

n_all_neuronal_skids <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "interneuron" |
           annotation == "motorneuron" |
           annotation == "Sensory neuron" |
           annotation == "immature sensory neuron" |
           annotation == "immature neuron") %>%
  unique() %>% pull(skid) %>% length()

n_non_neuro_celltypes <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "celltype_non_neuronal") %>%
  pull(skid) %>% length()
n_non_neuro_celltypes

num_remaining_cells <- with_soma_all_annot %>%
  pull(skid) %>% unique() %>%
  length()-n_neuro_celltypes 
num_remaining_cells <- num_remaining_cells %>%
  format(scientific = FALSE, big.mark = ",")

n_non_celltype <- with_soma_all_annot %>%
  pull(skid) %>% unique() %>%
  length()-n_non_neuro_celltypes-n_neuro_celltypes
n_non_celltype <- n_non_celltype %>%
  format(scientific = FALSE, big.mark = ",")

celltype_skids <- with_soma_all_annot %>%
  as_tibble() %>%
  filter(annotation == "celltype" | annotation == "celltype_non_neuronal")

# cell-type graph stats
Platy_celltype_graph <- readRDS("source_data/Figure4_source_data1.rds")
Platy_celltype_graph <- Platy_celltype_graph %>%
  mutate(degree = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE
  ))

SN_celltype_graph <- Platy_celltype_graph %>%
  filter(group == "sensory neuron") %>%
  pull(name) %>% length()
IN_celltype_graph <- Platy_celltype_graph %>%
  filter(group == "interneuron") %>%
  pull(name) %>% length()
MN_celltype_graph <- Platy_celltype_graph %>%
  filter(group == "motoneuron") %>%
  pull(name) %>% length()

effector_celltype_graph  <- Platy_celltype_graph %>%
  filter(group == "effector") %>%
  pull(name) %>% length()

# connectivity stats ---------------

with_soma_and_fragment_conn <- read.table("data/with_soma_and_fragment_conn.txt", sep = "\t")

with_soma_presyn_sum <- with_soma_and_fragment_conn %>%
  filter(skeleton == "with_soma") %>%
  filter(connection_type == "presynapse") %>%
  pull(syn.count) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")

with_soma_postsyn_sum <- with_soma_and_fragment_conn %>%
  filter(skeleton == "with_soma") %>%
  filter(connection_type == "postsynapse") %>%
  pull(syn.count) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")

fragment_presyn_sum <- with_soma_and_fragment_conn %>%
  filter(skeleton == "fragment") %>%
  filter(connection_type == "presynapse") %>%
  pull(syn.count) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")

fragment_postsyn_sum <- with_soma_and_fragment_conn %>%
  filter(skeleton == "fragment") %>%
  filter(connection_type == "postsynapse") %>%
  pull(syn.count) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")

percent_postsyn_on_frag <- round(with_soma_and_fragment_conn %>%
  filter(skeleton == "fragment") %>%
  filter(connection_type == "postsynapse") %>%
  pull(syn.count) %>% sum()/with_soma_and_fragment_conn %>%
  filter(skeleton == "with_soma") %>%
  filter(connection_type == "postsynapse") %>%
  pull(syn.count) %>% sum() * 100, 1)

percent_presyn_on_frag <- round(with_soma_and_fragment_conn %>%
  filter(skeleton == "fragment") %>%
  filter(connection_type == "presynapse") %>%
  pull(syn.count) %>% sum()/with_soma_and_fragment_conn %>%
  filter(skeleton == "with_soma") %>%
  filter(connection_type == "presynapse") %>%
  pull(syn.count) %>% sum() * 100, 1)

percent_postsyn_on_frag <- round(
  with_soma_and_fragment_conn %>%
  filter(skeleton == "fragment") %>%
  filter(connection_type == "postsynapse") %>%
  pull(syn.count) %>% sum()/with_soma_and_fragment_conn %>%
  filter(skeleton == "with_soma") %>%
  filter(connection_type == "postsynapse") %>%
  pull(syn.count) %>% sum() * 100, 1)

# cable length quantifications -------------

summ_All <- read_csv("source_data/Figure1_fig_suppl2_source_data_1.txt", show_col_types = FALSE)
summ_All <- summ_All %>%
  mutate(cable.length.um = cable.length/1000)

# mean cable length
median_cable_length_SN <- summ_All %>%
  filter(neuron_type == "sensory neuron") %>%
  select(cable.length.um) %>% 
  pull() %>%
  median() %>%
  format(scientific = FALSE, big.mark = ",", digits = 4)

median_cable_length_IN <- summ_All %>%
  filter(neuron_type == "interneuron") %>%
  select(cable.length.um) %>% 
  pull() %>%
  median() %>%
  format(scientific = FALSE, big.mark = ",", digits = 4)

median_cable_length_MN <- summ_All %>%
  filter(neuron_type == "motoneuron") %>%
  select(cable.length.um) %>% 
  pull() %>%
  median() %>%
  format(scientific = FALSE, big.mark = ",", digits = 4)

median_cable_length_all <-summ_All %>%
  filter(neuron_type != "fragment") %>%
  select(cable.length.um) %>% 
  pull() %>%
  median() %>%
  format(scientific = FALSE, big.mark = ",", digits = 4)

median_cable_length_frag <- summ_All %>%
  filter(neuron_type == "fragment") %>%
  select(cable.length.um) %>% 
  pull() %>%
  median() %>%
  format(scientific = FALSE, big.mark = ",", digits = 3)

summ_All %>%
  filter(neuron_type != "fragment") %>%
  arrange(desc(cable.length)) %>%
  select(neuron_names, cable.length) %>%
  slice(1:10)

summ_All %>%
  filter(neuron_type != "fragment") %>%
  arrange(desc(presyn_sites)) %>%
  select(neuron_names, presyn_sites) %>%
  slice(1:10)

summ_All %>%
  filter(neuron_type != "fragment") %>%
  arrange(desc(postsyn_sites)) %>%
  select(neuron_names, postsyn_sites) %>%
  slice(1:10)

median_presyn_sites <- summ_All %>%
  filter(neuron_type != "fragment") %>%
  select(presyn_sites) %>%
  pull() %>%
  median()

median_postsyn_sites <- summ_All %>%
  filter(neuron_type != "fragment") %>%
  select(postsyn_sites) %>%
  pull() %>%
  median()

# node properties ---------------
library(influenceR)
connect.tb %>%
  filter(node_is_cut == TRUE) %>%
  pull(names)
connect.tb %>%
  mutate(node_is_keyplayer = node_is_keyplayer(k = 10)) %>%
  filter(node_is_keyplayer == TRUE) %>%
  pull(names)

connect.tb %>%
  mutate(node_is_center = node_is_center()) %>%
  filter(node_is_center == TRUE) %>%
  pull(names)


global_efficiency(as.igraph(connect.tb), weights = NULL, directed = TRUE)


average_local_efficiency(
  as.igraph(connect.tb),
  weights = NULL,
  directed = TRUE,
  mode = c("all")
)

local_efficiency <- local_efficiency(
  as.igraph(connect.tb),
  vids = V(as.igraph(connect.tb)),
  weights = E(as.igraph(connect.tb))$weight,
  directed = TRUE,
  mode = c("all")
)

connect.tb %>%
  activate(nodes) %>%
  mutate(local_efficiency = local_efficiency) %>%
  arrange(desc(local_efficiency)) %>%
  select(names, local_efficiency) %>%
  as_tibble() %>%
  slice(1:100) %>% pull(names)

