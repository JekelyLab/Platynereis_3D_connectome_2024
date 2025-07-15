# code to generate Figure 4 fig suppl 5 of the Platynereis connectome paper
# Gaspar Jekely 2024

# load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# network plot ----------------------------------

# load network
syn_tb <- readRDS("source_data/Figure4_source_data1.rds")
conn_gexf <- rgexf::read.gexf("data/celltype_connectome_force_layout.gexf")

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


# add node in and out weight and I/O ratio --------------

syn.igraph <- as.igraph(syn_tb)
# calculate node weighted degree
weights_in <- strength(
  syn.igraph, vids = V(syn.igraph), 
  mode = "in", loops = TRUE
)
weights_out <- strength(
  syn.igraph, vids = V(syn.igraph), 
  mode = "out", loops = TRUE
)

# assign weighted degree to nodes
V(syn.igraph)$weights_in <- weights_in
V(syn.igraph)$weights_out <- weights_out

syn_tb <- syn.igraph %>%
  as_tbl_graph()

# assign coordinates from imported gephi graph to original CATMAID graph
syn_tb <- activate(syn_tb, nodes) %>%
  arrange(desc(name)) %>% # sort by skid to match coordinate list
  mutate(
    x = unlist(x_coord), # add coordinates
    y = unlist(y_coord)
  )

syn_tb_IO <- syn_tb %>%
  activate(nodes) %>%
  mutate(
    IO_rel_diff = round(((weights_in-weights_out)/
      (weights_in+weights_out)+1)*15), digits = 2)
  
level <- syn_tb_IO %>%
  as_tibble() %>%
  select(IO_rel_diff) %>%
  pull()
level

# graph visualisation -----------------------------------------------------

# convert to visNet form
syn.vis <- syn_tb_IO %>%
  toVisNetworkData()

# assign sqrt of number of synapses to edge 'value'
syn.vis$edges$value <- sqrt(syn.vis$edges$synapses)

#level	: Number. Default to undefined. When using the hierarchical layout, the level determines where the node is going to be positioned.
syn.vis$nodes$level <- level
#hierarchical layout

visNet <- visNetwork(syn.vis$nodes, syn.vis$edges) %>%
  visIgraphLayout(layout = "layout_nicely", physics = FALSE, 
                  randomSeed = 42, type="square") %>%
  visHierarchicalLayout(levelSeparation=280, 
                        nodeSpacing=70,
                        direction='LR',
                        blockShifting = TRUE,
                        sortMethod='hubsize',
                        shakeTowards='roots') %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0.2),
    scaling = list(min = 3, max = 30),
    color = list(inherit = TRUE, opacity = 0.5),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 1, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = syn.vis$nodes$color, border = "black"),
    scaling = list(min = 10, max = 50),
    font = list(color = "black", size = 40),
  ) %>%
  visGroups(
    groupname = "sensory neuron", color = "#E69F00", shape = "dot",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "interneuron", shape = "square",
    opacity = 1, color = "#CC79A7"
  ) %>%
  visGroups(
    groupname = "motoneuron", shape = "diamond",
    opacity = 1, color = "#0072B2"
  ) %>%
  visGroups(
    groupname = "effector", shape = "triangle",
    opacity = 1, color = "#CCCCCC"
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 6000, height = 3000) %>%
  addFontAwesome()

# save as html
saveNetwork(visNet, "pictures/network_with_IO_layout.html",
            selfcontained = TRUE
)
# save from web browser
webshot::webshot(
  url = "pictures/network_with_IO_layout.html",
  file = "Figures/Figure4_fig_suppl5.png",
  vwidth = 6000, vheight = 3000, # define the size of the browser window
  cliprect = c(100, 280, 5500, 2750), zoom = 1, delay = 1
)

write_rds(
  syn.vis, "source_data/Figure4_fig_suppl5_source_data1.bz2", 
  compress = "bz2"
  )
syn.vis <- read_rds("source_data/Figure4_fig_suppl5_source_data1.bz2")

