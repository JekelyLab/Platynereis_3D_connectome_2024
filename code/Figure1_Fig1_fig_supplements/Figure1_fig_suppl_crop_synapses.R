# code to generate synapses Fig supplement of the Platynereis 3d connectome paper
# Sanja Jasek & Gaspar Jekely 2024

source("code/Natverse_functions_and_conn.R")

dir.create("synapse_tiff_stacks")

# get all synapses from CATMAID
all_syn_connectors <- catmaid_fetch(
  path = "11/connectors/",
  body = list(
    relation_type = "presynaptic_to",
    relation_type = "postsynaptic_to",
    with_partners = "true"
  )
)

as_tibble(all_syn_connectors$connectors)

all_syn_connectors$connectors[[1]][[3]]

# cut out volume is not isometric, because it looks weird
half_bb_size_xy=1800
half_bb_size_z=1000
for (nconnector in seq_along(all_syn_connectors$connectors)) {
  #for (nconnector in seq(1,5)) {
  print(nconnector)
  id=all_syn_connectors$connectors[[nconnector]][[1]]
  x=all_syn_connectors$connectors[[nconnector]][[2]]
  y=all_syn_connectors$connectors[[nconnector]][[3]]
  z=all_syn_connectors$connectors[[nconnector]][[4]]
  catmaid_fetch(
    path = "11/crop/",
    body = list(
      stack_ids=28,
      min_x=x-half_bb_size_xy,
      min_y=y-half_bb_size_xy,
      min_z=z-half_bb_size_z,
      max_x=x+half_bb_size_xy,
      max_y=y+half_bb_size_xy,
      max_z=z+half_bb_size_z,
      zoom_level=0,
      single_channel='true',
      rotationcw=0
    )
  )
}

# server needs time to process this
Sys.sleep(200)

# download stacks ---------------------------------------------------------
# the exact file names are needed for download,
# to get them we need to check messages from the server
# originally I was checking which messages are new, and downloading those
# but this could lead to false positives
# so I decided that's it's probably safe to download the last n messages
# where n is the number of jobs

# It would be cool if I can put the synapse ID as file name,
# but it's not guaranteed that jobs will return in the same order

setwd("synapse_tiff_stacks")

cat_messages <- catmaid_fetch(
  path = "messages/list"
)

for (i in seq_along(all_syn_connectors$connectors)) {
  cat_message <- cat_messages[[i]]$text
  link <- regmatches(cat_message, 
                     regexpr("/crop/download/crop_.{6}.tiff", cat_message))
  filename <- regmatches(link, 
                         regexpr("crop_.{6}.tiff", link))
  full_link <- paste("https:/catmaid.ex.ac.uk", link, sep = "")
  download.file(full_link, destfile = filename)
}

# filter duplicate synapses

# make table with xyz positions of synapses
# conn_xyz <- data.frame()
# for (nconnector in seq_along(all_syn_connectors$connectors)) {
#   print(nconnector)
#   id=all_syn_connectors$connectors[[nconnector]][[1]]
#   x=all_syn_connectors$connectors[[nconnector]][[2]]
#   y=all_syn_connectors$connectors[[nconnector]][[3]]
#   z=all_syn_connectors$connectors[[nconnector]][[4]]
#   connector <- data.frame(
#     id=id,
#     x=x,
#     y=y,
#     z=z
#   )
#   conn_xyz <- rbind(conn_xyz, connector)
# }

# if a synapse has multiple postsynaptic sites, we were marking them as multiple synapses
# this is not ideal, so we try not ot have duplicates by discarding connectors that are too close to each other
# 
# min_distance <- 2000
# conn_xyz_curated <- data.frame()
# for (i in seq_along(nrow(conn_xyz))) {
#   xyz1 <- c(conn_xyz$x[[i]], conn_xyz$y[[i]], conn_xyz$z[[i]])
#   for (j in seq(i,nrow(conn_xyz))) {
#     xyz2 <- c(conn_xyz$x[[j]], conn_xyz$y[[j]], conn_xyz$z[[j]])
#     xyz <- rbind(xyz1, xyz2)
#     edistance <- dist(xyz)
#   }
# }

