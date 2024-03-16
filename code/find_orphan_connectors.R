library(catmaid)
library(tidyverse)

source("~/R/conn.R")

get_half_links <- function(connector_list) {
  for (i in (1:length(connector_list$partners))) {
    if ( length(connector_list$partners[[i]]) < 2 ) {
      # connectors with less than 2 partners
      # return treenode, not connector, because the connectors at same index don't match the connectors
      print(connector_list$partners[[i]][[1]][[2]])
    }
  }
}

get_sonnectors_with_single_node_partners <- function(connector_list) {
  for (i in (1:length(connector_list$partners))) {
    #print(i)
    if ( length(connector_list$partners[[i]]) > 1 ) {
      # connectors with 2 partners, but one of them only has a single node
      skid1 <- connector_list$partners[[i]][[1]][[3]]
      cable_length1 <- catmaid_fetch(path = paste("11/skeletons/", skid1, "/cable-length", sep = ""))[[1]]
      skid2 <- connector_list$partners[[i]][[2]][[3]]
      cable_length2 <- catmaid_fetch(path = paste("11/skeletons/", skid2, "/cable-length", sep = ""))[[1]]
      if (cable_length1 == 0 | cable_length2 == 0) {
        # only return treenode of fist partner, no matter which partner has a single node
        # because it's equally easy to find it in catmaid with either ID, but reduces code here
        print(connector_list$partners[[i]][[1]][[2]])
      }
    }
  }
}


print("Presynaptic half links")
all_pre_connectors <- catmaid_fetch(path = "11/connectors/",
                                    body = list(relation_type="presynaptic_to", 
                                                with_partners="true"))
get_half_links(all_pre_connectors)


print("Postsynaptic half links")
all_post_connectors <- catmaid_fetch(path = "11/connectors/",
                                     body = list(relation_type="postsynaptic_to", 
                                                 with_partners="true"))
get_half_links(all_post_connectors)

print("Single node skeletons with synapses")
get_sonnectors_with_single_node_partners(all_pre_connectors)
