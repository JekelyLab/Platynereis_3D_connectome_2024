# code for Figure4 figure supplement 4 of the 3d Platynereis connectome paper
#Gaspar Jekely 2023
# load packages, functions and catmaid connection
source("code/libraries_functions_and_CATMAID_conn.R")


Okabe_Ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )

color_scale <- c(Okabe_Ito, blues9)

# read graphs ----------------
Platy_celltype_graph <- readRDS("source_data/Figure4_source_data1.rds")
C_elegans <- read_tsv("source_data/Cook_et_al_Hermaphrodite.csv")
Ciona <- read.csv("source_data/Ryan_et_al_Ciona_matrix.txt", sep = "\t")

Platy_celltype_graph <- Platy_celltype_graph %>%
  mutate(weighted_degree = centrality_degree(
    weights = E(as.igraph(Platy_celltype_graph))$synapses,
    mode = "all",
    loops = TRUE,
    normalized = FALSE)
  )

Platy_celltype_graph <- Platy_celltype_graph %>%
  mutate(degree = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE,
    normalized = FALSE)
  )

# C elegans graph processing --------------

C_elegans_graph <- C_elegans %>%
  pivot_longer(-...1, names_to = "postsynaptic", values_to = "synapses") %>%
  rename(presynaptic = ...1) %>%
  filter(!is.na(synapses)) %>%
  as_tbl_graph() %>%
  mutate(degree = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE
  )) 

C_elegans_graph <- C_elegans_graph %>%
  mutate(weighted_degree = centrality_degree(
    weights = E(C_elegans_graph)$synapses,
    mode = "all",
    loops = TRUE,
    normalized = FALSE)
  )

#the C. elegans graph is composed of two subgraphs
# check connected components, the connectome should only have one component
cl <- components(as.igraph(C_elegans_graph))
# size of the largest subnetwork
length(which(cl$membership == 2))

# generate subgraph
C_elegans_graph_sub1 <- subgraph(
  as.igraph(C_elegans_graph), 
  which(cl$membership == 2)
  )

# clustering with Leiden algorithm, run with ModularityVertexPartition and resolution parameter
partition_C_elegans <- leiden(C_elegans_graph_sub1,
                              weights = E(C_elegans_graph_sub1)$weight,
                              partition_type = "RBConfigurationVertexPartition",
                              resolution_parameter = 2,
                              n_iterations = -1, seed = 42
)

mod_C_el <- round(modularity(C_elegans_graph_sub1, partition_C_elegans), 3)

# Ciona graph processing ---------------

Ciona_graph <- Ciona %>%
  pivot_longer(-X, names_to = "postsynaptic", values_to = "connectivity") %>%
  rename(presynaptic = X) %>%
  filter(!is.na(connectivity)) %>%
  as_tbl_graph() %>%
  mutate(degree = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE
  )) 

Ciona_graph <- Ciona_graph %>%
  mutate(weighted_degree = centrality_degree(
    weights = Ciona_graph %>% activate(edges) %>% pull(connectivity),
    mode = "all",
    loops = TRUE
  ))

# clustering with Leiden algorithm, run with ModularityVertexPartition and resolution parameter
partition_Ciona <- leiden(as.igraph(Ciona_graph),
                              weights = E(as.igraph(Ciona_graph))$connectivity,
                              partition_type = "RBConfigurationVertexPartition",
                              resolution_parameter = 2,
                              n_iterations = -1, seed = 42
)

mod_Ci <- round(modularity(Ciona_graph, partition_Ciona), 3)

# Ciona vis -----------------------
# convert to visNet form
Ciona.vis <- Ciona_graph %>%
  toVisNetworkData()

# assign sqrt of number of synapses to edge 'value'
Ciona.vis$edges$value <- sqrt(Ciona.vis$edges$connectivity)
Ciona.vis$nodes$value <- Ciona.vis$nodes$weighted_degree
Ciona.vis$nodes$color <- color_scale[partition_Ciona]

visNet_Ciona <- visNetwork(Ciona.vis$nodes, Ciona.vis$edges) %>%
  visIgraphLayout(layout = "layout_nicely",
                  physics = FALSE, randomSeed = 42) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 4, max = 30),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.8, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = Ciona.vis$nodes$color, border = "black"),
    scaling = list(min = 1, max = 50),
    font = list(color = "black", size = 0),
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 1500, height = 1500) %>%
  addFontAwesome()


# save as html
saveNetwork(visNet_Ciona, "pictures/visNet_Ciona.html",
            selfcontained = TRUE)

# save from web browser
webshot::webshot(
  url = "pictures/visNet_Ciona.html",
  file = "pictures/visNet_Ciona.png",
  vwidth = 1500, vheight = 1500, # define the size of the browser window
  cliprect = c(60, 60, 1450, 1450), zoom = 2, delay = 2
)


# C elegans vis -----------------------
# convert to visNet form
C_elegans.vis <- C_elegans_graph_sub1 %>%
  toVisNetworkData()

# assign edge and node 'value'
C_elegans.vis$edges$value <- C_elegans.vis$edges$synapses
C_elegans.vis$nodes$value <- C_elegans.vis$nodes$weighted_degree
C_elegans.vis$nodes$color <- color_scale[partition_C_elegans]

visNet_C_elegans <- visNetwork(C_elegans.vis$nodes, C_elegans.vis$edges) %>%
  visIgraphLayout(layout = "layout_nicely",
                  physics = FALSE, randomSeed = 42) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 4, max = 30),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.8, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = C_elegans.vis$nodes$color, border = "black"),
    scaling = list(min = 1, max = 50),
    font = list(color = "black", size = 0),
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 1500, height = 1500) %>%
  addFontAwesome()

# save as html
saveNetwork(visNet_C_elegans, "pictures/visNet_C_elegans.html",
            selfcontained = TRUE)

# save from web browser
webshot::webshot(
  url = "pictures/visNet_C_elegans.html",
  file = "pictures/visNet_C_elegans.png",
  vwidth = 1500, vheight = 1500, # define the size of the browser window
  cliprect = c(60, 60, 1200, 1200), zoom = 2, delay = 2
)

# Platynereis cell type network ----------------

Platy_celltype_graph.igraph <- as.igraph(Platy_celltype_graph)
conn_gexf <- rgexf::read.gexf("data/celltype_connectome_force_layout.gexf")

# clustering with Leiden algorithm, run with ModularityVertexPartition and resolution parameter
partition_Platy <- leiden(Platy_celltype_graph,
                              weights = E(Platy_celltype_graph.igraph)$weight,
                              partition_type = "RBConfigurationVertexPartition",
                              resolution_parameter = 2,
                              n_iterations = -1, seed = 42
)
max(partition_Platy)
mod_Platy_c <- round(modularity(Platy_celltype_graph, partition_Platy), 3)

# get coordinates from imported gephi file
coords <- as_tibble(x = conn_gexf$nodesVizAtt$position$x) %>%
  mutate(x = value) %>%
  mutate(y = conn_gexf$nodesVizAtt$position$y) %>%
  mutate(skid = conn_gexf$nodes$label)

# sort by skid
x_coord <- as.list(arrange(coords, desc(skid)) %>%
                     select(x))
y_coord <- as.list(arrange(coords, desc(skid)) %>%
                     select(y))

# assign coordinates from imported gephi graph to original CATMAID graph
Platy_celltype_graph <- activate(Platy_celltype_graph, nodes) %>%
  arrange(desc(name)) %>% # sort by skid to match coordinate list
  mutate(
    x = unlist(x_coord), # add coordinates
    y = unlist(y_coord)
  )

# convert to visNet form
Platy_celltype_graph.vis <- Platy_celltype_graph %>%
  toVisNetworkData()

# assign sqrt of number of synapses to edge 'value'
Platy_celltype_graph.vis$edges$value <- sqrt(Platy_celltype_graph.vis$edges$synapses)

#coordinates as matrix
coords <- matrix(c(Platy_celltype_graph.vis$nodes$x, Platy_celltype_graph.vis$nodes$y), ncol = 2)



Platy_celltype_graph.vis$nodes$value <- Platy_celltype_graph.vis$nodes$weighted_degree
Platy_celltype_graph.vis$nodes$color <- color_scale[partition_Platy]

visNet_Platy <- visNetwork(Platy_celltype_graph.vis$nodes, Platy_celltype_graph.vis$edges) %>%
  visIgraphLayout(type = "full", layoutMatrix = coords) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0.1),
    scaling = list(min = 4, max = 30),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.8, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(border = "black"),
    scaling = list(min = 1, max = 50),
    font = list(color = "black", size = 0),
  )  %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 1500, height = 700) %>%
  addFontAwesome()


# save as html
saveNetwork(visNet_Platy, "pictures/visNet_Platy.html",
            selfcontained = TRUE)

# save from web browser
webshot::webshot(
  url = "pictures/visNet_Platy.html",
  file = "pictures/visNet_Platy.png",
  vwidth = 1500, vheight = 700, # define the size of the browser window
  cliprect = c(80, 60, 1400, 640), zoom = 2, delay = 2
)

# Platy full network ------------------------


# read graph
Platy_full_connect.tb <- readRDS("supplements/connectome_graph_tibble.rds")
conn_graph.visn <- readRDS("supplements/connectome_graph.rds")
coords <- matrix(c(conn_graph.visn$nodes$x, conn_graph.visn$nodes$y), ncol = 2)

Platy_full_connect.tb <- Platy_full_connect.tb %>%
  mutate(degree = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE,
    normalized = FALSE)
  )

partitions <- Platy_full_connect.tb %>%
  select(partition) %>%
  pull() %>%
  max()

coords_rotated <- autoimage::rotate(
  coords,
  -pi / 3,
  pivot = c(0, 0)
)

full_Platy_visNet <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = coords_rotated) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 4, max = 30),
    color = list(inherit = TRUE, opacity = 0.5),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.5, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.1,
    color = list(border = "black"),
    opacity = 1,
    scaling = list(min = 1, max = 50),
    font = list(size = 0)
  ) %>%
  visOptions(
    highlightNearest = list(
      enabled = TRUE, degree = 1,
      algorithm = "hierarchical", labelOnly = FALSE
    ),
    width = 2500, height = 2500, autoResize = FALSE
  )

full_Platy_visNet
# save as html
saveNetwork(full_Platy_visNet, "pictures/Full_connectomeNo_text_modules.html", selfcontained = TRUE)

webshot::webshot(
  url = "pictures/Full_connectomeNo_text_modules.html",
  file = "pictures/Full_connectomeNo_text_modules.png",
  vwidth = 2550, vheight = 2550, # define the size of the browser window
  cliprect = c(90, 105, 2420, 2400), zoom = 2
)

# statistics comparison ------------------------


ge_C_el <- round(global_efficiency(C_elegans_graph_sub1), 3) %>%
  format(scientific = FALSE, big.mark = ",")
ge_Ci <- round(global_efficiency(as.igraph(Ciona_graph)), 3) %>%
  format(scientific = FALSE, big.mark = ",")
ge_Pdu_c <- round(global_efficiency(as.igraph(Platy_celltype_graph)), 3) %>%
  format(scientific = FALSE, big.mark = ",")
ge_Pdu_f <- round(global_efficiency(Platy_full_connect.tb), 3) %>%
  format(scientific = FALSE, big.mark = ",")

is.directed(C_elegans_graph_sub1)
is.directed(as.igraph(Ciona_graph))
is.directed(as.igraph(Platy_celltype_graph))
is.directed(Platy_full_connect.tb)

gd_C_el <- round(graph.density(C_elegans_graph_sub1), 3) %>%
  format(scientific = FALSE, big.mark = ",")
gd_Ci <- round(graph.density(Ciona_graph), 3) %>%
  format(scientific = FALSE, big.mark = ",")
gd_Pdu_c <- round(graph.density(Platy_celltype_graph), 3) %>%
  format(scientific = FALSE, big.mark = ",")
gd_Pdu_f <- round(graph.density(Platy_full_connect.tb), 3) %>%
  format(scientific = FALSE, big.mark = ",")

# Average clustering coefficient
tr_C_el <- round(transitivity(C_elegans_graph_sub1, type = "average"), 3) %>%
  format(scientific = FALSE, big.mark = ",")
tr_Ci <- round(transitivity(Ciona_graph, type = "average"), 3) %>%
  format(scientific = FALSE, big.mark = ",")
tr_Pdu_c <- round(transitivity(Platy_celltype_graph, type = "average"), 3) %>%
  format(scientific = FALSE, big.mark = ",")
tr_Pdu_f <- round(transitivity(Platy_full_connect.tb, type = "average"), 3) %>%
  format(scientific = FALSE, big.mark = ",")

#diameter
di_C_el <- round(diameter(C_elegans_graph_sub1), 3) %>%
  format(scientific = FALSE, big.mark = ",")
di_Ci <- round(diameter(Ciona_graph), 3) %>%
  format(scientific = FALSE, big.mark = ",")
di_Pdu_c <- round(diameter(Platy_celltype_graph), 3) %>%
  format(scientific = FALSE, big.mark = ",")
di_Pdu_f <- round(diameter(Platy_full_connect.tb), 3) %>%
  format(scientific = FALSE, big.mark = ",")

#mean distance
md_C_el <- round(mean_distance(C_elegans_graph_sub1), 2) %>%
  format(scientific = FALSE, big.mark = ",")
md_Ci <- round(mean_distance(Ciona_graph), 2) %>%
  format(scientific = FALSE, big.mark = ",")
md_Pdu_c <- round(mean_distance(Platy_celltype_graph), 2) %>%
  format(scientific = FALSE, big.mark = ",")
md_Pdu_f <- round(mean_distance(Platy_full_connect.tb), 2) %>%
  format(scientific = FALSE, big.mark = ",")
                  
#average degree
avgd_C_el <- round(C_elegans_graph_sub1 %>%
  as_tbl_graph() %>%
  pull(degree) %>% mean()/2, 2) %>%
  format(scientific = FALSE, big.mark = ",")
avgd_Ci <- round(Ciona_graph %>%
  pull(degree) %>% mean()/2, 2) %>%
  format(scientific = FALSE, big.mark = ",")
avgd_Pdu_c <- round(Platy_celltype_graph %>%
  pull(degree) %>% mean()/2, 2) %>%
  format(scientific = FALSE, big.mark = ",")
avgd_Pdu_f <- round(Platy_full_connect.tb %>%
  pull(degree) %>% mean()/2, 2) %>%
  format(scientific = FALSE, big.mark = ",")

# n synapses, nodes and edges
C_elegans_n_synapses <- C_elegans_graph_sub1 %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  pull(synapses) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")
C_elegans_n_edges <- C_elegans_graph_sub1 %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  pull(synapses) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")
C_elegans_n_nodes <- C_elegans_graph_sub1 %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  pull(name) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")

Ciona_n_edges <- Ciona_graph %>%
  activate(edges) %>%
  pull(connectivity) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")
Ciona_n_nodes <- Ciona_graph %>%
  activate(nodes) %>%
  pull(name) %>% length()  %>%
  format(scientific = FALSE, big.mark = ",") #no synapse count for Ciona, edge weights are not syn count

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

#modularity
mod_Platy_c <- round(modularity(as.igraph(Platy_celltype_graph), 
                          leiden(Platy_celltype_graph,
                          weights = E(as.igraph(Platy_celltype_graph))$synapses,
                          partition_type = "RBConfigurationVertexPartition",
                          resolution_parameter = 2,
                          n_iterations = -1, seed = 42)
                          ),
                     3)

Platy_full_n_synapses <- Platy_full_connect.tb %>%
  activate(edges) %>%
  pull(weight) %>% sum() %>%
  format(scientific = FALSE, big.mark = ",")
Platy_full_n_edges <- Platy_full_connect.tb %>%
  activate(edges) %>%
  pull(weight) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")
Platy_full_n_nodes <- Platy_full_connect.tb %>%
  activate(nodes) %>%
  pull(skids) %>% length() %>%
  format(scientific = FALSE, big.mark = ",")

#modularity
mod_Platy_f <- round(modularity(as.igraph(Platy_full_connect.tb), 
           leiden(as.igraph(Platy_full_connect.tb),
           weights = E(as.igraph(Platy_full_connect.tb))$weight,
           partition_type = "RBConfigurationVertexPartition",
           resolution_parameter = 2,
           n_iterations = -1, seed = 42)
           ), 3)


# Erdos graph -------------------------

Erdös_graphs_100_wg <- lapply(1:100, function(x){
  Erdös_graph <- erdos.renyi.game(
    length(V(Platy_celltype_graph)), 
    length(E(Platy_celltype_graph)),
    type = "gnm", directed = TRUE, loops = FALSE
  )
  
  E(Erdös_graph)$weight <- E(Platy_celltype_graph)$synapses
  x <- Erdös_graph
}
)


#calculate mean distance
md_Erdos <- round(mean(sapply(Erdös_graphs_100_wg, function(x) mean_distance(x))), 2)
#calculate diameter
di_Erdos <- round(mean(sapply(Erdös_graphs_100_wg, function(x) diameter(x))), 3)
Erdos_n_nodes <- length(V(Platy_celltype_graph)) 
Erdos_n_edges <- length(E(Platy_celltype_graph))
Erdos_n_synapses <- mean(sapply(Erdös_graphs_100_wg, function(x) sum(E(x)$weight)))

ge_Erdos <- round(mean(sapply(Erdös_graphs_100_wg, function(x) global_efficiency(x))), 3)

gd_Erdos <- round(mean(sapply(Erdös_graphs_100_wg, function(x) graph.density(x))), 3)
avgd_Erdos <- round(mean(sapply(Erdös_graphs_100_wg, function(x) degree(x)))/2, 2)

#transitivity (clustering coefficient) 
tr_Erdos <- round(mean(sapply(
  Erdös_graphs_100_wg, 
  function(x) transitivity(x, type = "average")
)), 3)

#modularity
mod_erdos <- round(mean(sapply(
  Erdös_graphs_100_wg,
  function(x)
  modularity(x, leiden(
    x,
    weights = E(x)$weight,
    partition_type = "RBConfigurationVertexPartition",
    resolution_parameter = 2,
    n_iterations = -1, seed = 42)
    ))
  ), 3)

# Erdos vis -----------------------

Erdös_graph <- erdos.renyi.game(
  length(V(Platy_celltype_graph)), 
  length(E(Platy_celltype_graph)),
  type = "gnm", directed = TRUE, loops = FALSE
)

E(Erdös_graph)$weight <- E(Platy_celltype_graph)$synapses

Erdös_graph <- as_tbl_graph(Erdös_graph) %>%
  mutate(weighted_degree = centrality_degree(
    weights = E(Erdös_graph)$weight,
    mode = "all",
    loops = TRUE
  ))

# convert to visNet form
erdos.vis <- Erdös_graph %>%
  toVisNetworkData()

partition_erdos <- leiden(Erdös_graph,
                          weights = E(Erdös_graph)$weight,
                          partition_type = "RBConfigurationVertexPartition",
                          resolution_parameter = 2,
                          n_iterations = -1, seed = 42
)

# assign edge and node 'value'
erdos.vis$edges$value <- erdos.vis$edges$weight
erdos.vis$nodes$value <- erdos.vis$nodes$weighted_degree
color_erdos <- c(color_scale, hcl.colors(12, palette = "viridis"))
erdos.vis$nodes$color <- color_erdos[partition_erdos]
erdos.vis$nodes$color

visNet_erdos <- visNetwork(erdos.vis$nodes, erdos.vis$edges) %>%
  visIgraphLayout(layout = "layout_nicely",
                  physics = FALSE, randomSeed = 42) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0),
    scaling = list(min = 4, max = 30),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.8, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = erdos.vis$nodes$color, border = "black"),
    scaling = list(min = 1, max = 50),
    font = list(color = "black", size = 0),
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 1500, height = 1500) %>%
  addFontAwesome()

# save as html
saveNetwork(visNet_erdos, "pictures/visNet_erdos.html",
            selfcontained = TRUE)

# save from web browser
webshot::webshot(
  url = "pictures/visNet_erdos.html",
  file = "pictures/visNet_erdos.png",
  vwidth = 1500, vheight = 1500, # define the size of the browser window
  cliprect = c(0, 0, 1500, 1500), zoom = 5, delay = 2
)

# table ---------------------

table <- plot_ly(
  type = "table",
  columnwidth = c(8, 6, 6, 6, 6, 6),
  columnorder = c(0, 1, 2, 3 , 4, 5),
  header = list(
    values = c(
      "graph property", 
      "<i>Platynereis</i> full",
      "<i>Platynereis</i> cell-type",
      "<i>C. elegans</i> full",
      "<i>Ciona</i> full",
      "random graph"),
    align = c("center", "center", "center", "center", "center", "center"),
    line = list(width = 1, color = "black"),
    fill = list(color = c("#CCCCCC", "#AAAAAA", "#AAAAAA", "#AAAAAA", "#AAAAAA", "#AAAAAA")),
    font = list(family = "Arial", size = 14, color = "black")
  ),
  cells = list(
    values = rbind(
      c("# nodes", "# edges",
        "# synapses", "diameter", "mean distance", "avg. degree", 
        "density",
        "global efficiency", "local avg. transitivity",
        "modularity"
      ),
      c(
        Platy_full_n_nodes, Platy_full_n_edges, 
        Platy_full_n_synapses, di_Pdu_f, md_Pdu_f, 
        avgd_Pdu_f, gd_Pdu_f, ge_Pdu_f, tr_Pdu_f,
        mod_Platy_f
      ),
      c(
        Platy_celltype_graph_n_nodes, Platy_celltype_graph_n_edges, 
        Platy_celltype_graph_n_synapses, di_Pdu_c, md_Pdu_c, 
        avgd_Pdu_c, gd_Pdu_c, ge_Pdu_c, tr_Pdu_c,
        mod_Platy_c
      ),
      c(
        C_elegans_n_nodes, C_elegans_n_edges,
        C_elegans_n_synapses, di_C_el, md_C_el, 
        avgd_C_el, gd_C_el, ge_C_el, tr_C_el,
        mod_C_el
      ),
      c(
        Ciona_n_nodes, Ciona_n_edges, "n.a", di_Ci, md_Ci, 
        avgd_Ci, gd_Ci, ge_Ci, tr_Ci,
        mod_Ci
      ),
      c(
        Erdos_n_nodes, Erdos_n_edges, Erdos_n_synapses, di_Erdos, md_Erdos, 
        avgd_Erdos, gd_Erdos, ge_Erdos, tr_Erdos,
        mod_erdos
      )
    ),
    align = c("center", "center","center","center","center","center"),
    line = list(color = "black", width = 0.3),
    font = list(family = "Arial", size = 12, color = c("black"))
  )
)

table
saveNetwork(table, "pictures/Fig4_fig_suppl4_network_stats_table.html")
webshot::webshot(
  url = "pictures/Fig4_fig_suppl4_network_stats_table.html",
  file = "pictures/Fig4_fig_suppl4_network_stats_table.png",
  vwidth = 600, vheight = 400, # define the size of the browser window
  cliprect = c(20, 55, 540, 325), zoom = 10
)

# assemble figure --------------------

panel_C_elegans <- ggdraw() + draw_image(readPNG("pictures/visNet_C_elegans.png")) +
  draw_label("C. elegans full connectome", x = 0.5, y = 0.97, size = 11)
panel_Ciona <- ggdraw() + draw_image(readPNG("pictures/visNet_Ciona.png")) +
  draw_label("Ciona larva full connectome", x = 0.5, y = 0.97, size = 11)

panel_erdos <- ggdraw() + draw_image(readPNG("pictures/visNet_erdos.png")) +
  draw_label("random graph", x = 0.5, y = 0.97, size = 11)

panel_Platy <- ggdraw() + draw_image(readPNG("pictures/visNet_Platy.png")) +
  draw_label("Platynereis cell-type connectome", x = 0.5, y = 0.97, size = 11)
panel_Platy_full <- ggdraw() + draw_image(readPNG("pictures/Full_connectomeNo_text_modules.png")) +
  draw_label("Platynereis full connectome", x = 0.4, y = 0.97, size = 11)

panel_table <- ggdraw() + draw_image(readPNG("pictures/Fig4_fig_suppl4_network_stats_table.png"))

layout <- "
AABB
CDEF
"

Fig <- panel_Platy_full + panel_table +
  panel_Platy + panel_C_elegans + panel_Ciona + panel_erdos +
  plot_layout(design = layout, widths = c(1, 1, 1, 1),
              heights = c(2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=12, face='plain'))

ggsave("Figures/Figure4_fig_suppl4.png", limitsize = FALSE, 
       units = c("px"), Fig, 
       width = 4200, height = 3000, bg='white'
)  

ggsave("Figures/Figure4_fig_suppl4.pdf", limitsize = FALSE, 
       units = c("px"), Fig, 
       width = 4200, height = 3000
)  

