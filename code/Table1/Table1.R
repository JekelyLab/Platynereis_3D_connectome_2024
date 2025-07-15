# code to generate Table 1 of the Platynereis connectome paper
# Gaspar Jekely 2023

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/libraries_functions_and_CATMAID_conn.R")

# load cell type graph
syn_tb <- readRDS("source_data/Figure4_source_data1.rds")

# distances plot ---------------------


syn_tb_filt <- syn_tb %>%
  activate(edges) %>%
  filter(synapses>3) %>%
  activate(nodes)

# distances from SN to effectors ------------
d_eff <- distances(syn_tb_filt, v = syn_tb_filt %>%
                     filter(group == "sensory neuron") %>% pull(name),
                   to = syn_tb_filt %>%
                     filter(group == "effector") %>% pull(name),
                   mode = "out")

# sensory-motor neurons ------------------

sensory_motor_neurons <- as.data.frame(d_eff) %>% rownames_to_column("input") %>%
  pivot_longer(
    -input, names_to = "effector", values_to = "distance"
  ) %>%
  group_by(input) %>%
  mutate(shortest_path_to_effector = min(distance)) %>%
  select(input, shortest_path_to_effector) %>%
  unique() %>%
  filter(shortest_path_to_effector == 1) %>%
  pull(input)

sensory_motor_neurons_l <- paste0(sensory_motor_neurons, collapse = ", ")

# premotor SN neurons -----------
d_MN <- distances(syn_tb_filt, v = syn_tb_filt %>%
                    filter(group == "sensory neuron") %>% pull(name),
                  to = syn_tb_filt %>%
                    filter(group == "motoneuron") %>% pull(name),
                  mode = "out")

premotor_SN <- as.data.frame(d_MN) %>% rownames_to_column("input") %>%
  pivot_longer(
    -input, names_to = "MN", values_to = "distance"
  ) %>%
  filter(distance == 1)
premotor_SNs <- premotor_SN %>%
  pull(input) %>% unique()

premotor_SNs_l <- paste0(premotor_SNs, collapse = ", ")

# SN with no path to effectors ---------------
SN_with_effector_connection <- as.data.frame(d_eff) %>% rownames_to_column("input") %>%
  pivot_longer(
    -input, names_to = "effector", values_to = "distance"
  ) %>%
  group_by(input) %>%
  mutate(shortest_path_to_effector = min(distance)) %>%
  filter(shortest_path_to_effector != Inf) %>%
  pull(input) %>% unique()

SN_with_effector_connection_l <- paste0(SN_with_effector_connection, collapse = ", ")

all_SN <- as.data.frame(d_eff) %>% rownames_to_column("input") %>%
  pivot_longer(
    -input, names_to = "effector", values_to = "distance"
  ) %>%
  pull(input) %>% unique()

SN_with_NO_effector_connection <- setdiff(all_SN, SN_with_effector_connection)

SN_with_NO_effector_connection_l <- paste0(SN_with_NO_effector_connection, collapse = ", ")

#verify that these cells have no connection to effectors
as.data.frame(d_eff) %>% rownames_to_column("input") %>%
  pivot_longer(
    -input, names_to = "effector", values_to = "distance"
  ) %>%
  group_by(input) %>%
  mutate(shortest_path_to_effector = min(distance)) %>%
  filter(input %in% SN_with_NO_effector_connection) %>%
  pull(shortest_path_to_effector) %>% unique()

# table of SN types
table_SN <- plot_ly(
  type = "table",
  columnwidth = c(10, 20, 20, 20),
  columnorder = c(0, 1, 2, 3),
  header = list(
    values = c(
      "sensory-motor", "premotor SN", 
      "SN 2-5 hops from effectors", 
      "SN with no path to effectors"
    ),
    align = c("center", "center"),
    line = list(width = 1, color = "black"),
    fill = list(color = c("#CCCCCC","#CCCCCC","#CCCCCC", "#CCCCCC")),
    font = list(family = "Arial", size = 14, color = "black")
  ),
  cells = list(
    values = rbind(
      sensory_motor_neurons_l,
      premotor_SNs_l,
      SN_with_effector_connection_l,
      SN_with_NO_effector_connection_l
    ),
    align = c("center", "center", "center", "center"),
    line = list(color = "black", width = 0.3),
    font = list(family = "Arial", size = 12, color = c("black"))
  )
)

data_for_table <- tibble(
  "sensory-motor" = sensory_motor_neurons_l, 
  "premotor SN" = premotor_SNs_l, 
  "SN 2-5 hops from effectors" = SN_with_effector_connection_l, 
  "SN with no path to effectors" = SN_with_NO_effector_connection_l
)

data_for_table

table_SN_tt <- tt(
  data_for_table
)

saveRDS(table_SN_tt, "data/Table1.RDS")

saveNetwork(table_SN, "pictures/Table1.html")
webshot2::webshot(
  url = "pictures/Table1.html",
  file = "Figures/Table1.pdf",
  vwidth = 1000, vheight = 310, # define the size of the browser window
  zoom = 1
)

webshot2::webshot(
  url = "pictures/Table1.html",
  file = "Figures/Table1.png",
  vwidth = 900, vheight = 350, # define the size of the browser window
  zoom = 10, cliprect = c(20, 52, 847, 275)
)

