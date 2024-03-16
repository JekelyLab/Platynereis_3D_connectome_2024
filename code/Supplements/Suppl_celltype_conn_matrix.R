#code to generate the Supplementary figure cell type connectivity matrix of the 3 day old Platynereis larva 
#described in Veraszto et al.
#Gaspar Jekely 2022 

#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

# read all skids based on cell type annotations ---------------------------

skids_celltypelist = list()
n_neuron_types = 202
n_non_neuron_types = 91 #we ignore the yolk (91) and some mesoderm cells (91)

#first we read all skids for Sensory cell types from 1-200
counter=0
for (i in c(1:n_neuron_types)){
  annotation = paste("annotation:^celltype", i, "$", sep="")
  #retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "Sensory neuron")
  print (i)
  if (length(skids)>0) {
    counter <- counter+1
    skids_celltypelist[[counter]] <- skids
    print ("Sensory neuron")
  }
}
length(skids_celltypelist)
counter

#we read all skids for interneuron cell types from 1-200
#we do not reset the counter and continue to fill the skids_celltypelist
for (i in c(1:n_neuron_types)){
  annotation = paste("annotation:^celltype", i, "$", sep="")
  #retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "interneuron")
  print (i)
  if (length(skids)>0) {
    counter <- counter+1
    skids_celltypelist[[counter]] <- skids
    print ("interneuron")
  }
}
length(skids_celltypelist)
counter

#we read all skids for motorneuron cell types from 1-200
for (i in c(1:n_neuron_types)){
  annotation = paste("annotation:^celltype", i, "$", sep="")
  #retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "motorneuron")
  print (i)
  if (length(skids)>0) {
    counter <- counter+1
    skids_celltypelist[[counter]] <- skids
    print ("motorneuron")
  }
}
length(skids_celltypelist)
counter
#the skids_celltypelist has the skids ordered by Sensory, inter, motorneuron

#we read the skids for all non-neuronal cell types from 1-90
#in the order 'ciliated cell', 'gland cell', 'pigment cell', 'muscle'
#continue to use the counter
for (i in c(1:n_non_neuron_types)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "ciliated cell")
  print (i)
  if (length(skids)>0) {
    counter <- counter+1
    skids_celltypelist[[counter]] <- skids
    print ("ciliated cell")
  }
}
length(skids_celltypelist)
counter


for (i in c(1:n_non_neuron_types)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "gland cell")
  print (i)
  if (length(skids)>0) {
    counter <- counter+1
    skids_celltypelist[[counter]] <- skids
    print ("gland cell")
  }
}
length(skids_celltypelist)
counter

for (i in c(1:n_non_neuron_types)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "pigment cell")
  print (i)
  if (length(skids)>0) {
    counter <- counter+1
    skids_celltypelist[[counter]] <- skids
    print ("pigment cell")
  }
}
length(skids_celltypelist)
counter

for (i in c(1:n_non_neuron_types)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #retrieve all skids that match the two annotations
  skids <- skids_by_2annotations(annotation, "muscle")
  print (i)
  if (length(skids)>0) {
    counter <- counter+1
    skids_celltypelist[[counter]] <- skids
    print ("muscle")
  }
}
length(skids_celltypelist)
counter

#we ignore other non-neuronal cell types

# retrieve connectivity from CATMAID --------------------------------------

#define empty synapse matrix 
synapse_matrix <- matrix(NA, nrow = (length(skids_celltypelist)), 
                         ncol = (length(skids_celltypelist)))
dim(synapse_matrix)

#iterate through cell type list, retrieve and count connections and populate connectivity matrix
#Runs till 275
for (i in 1:(length(skids_celltypelist))){
  for (j in 1:(length(skids_celltypelist))){
    cellgroup_conn <- catmaid_get_connectors_between(pre=unlist(skids_celltypelist[i]), 
         post=unlist(skids_celltypelist[j]), pid=11, conn = conn_http1)
    N_synapses <- nrow(cellgroup_conn)
    if(length(N_synapses)==0) N_synapses <- 0 #if there are no synapses, need to add 0
    synapse_matrix[i, j] <- N_synapses
  }
  print (i)
}

# get neuron names and rename matrix col and row ---------------------------------

#we use the Catmaid neuron name before the "_" to parse the generic cell type name

#get the neuron name of the first skid in the skids list
celltype_names<-list()

for(i in 1:length(skids_celltypelist)){
  name<-catmaid_get_neuronnames(skids_celltypelist[[i]][1],pid=11)
  name <- sub("_.*$", "", name )
  celltype_names[i] <- name
}

#check duplicated names
celltype_names[duplicated(celltype_names)]
#check length
length(celltype_names)

#assign column and row names
rownames(synapse_matrix) <- celltype_names
colnames(synapse_matrix) <- celltype_names


# search for annotations --------------------------------------------------

#these are the annotations to search for
annot_to_search <- c("Sensory neuron", "interneuron", "motorneuron", "effector")

#iterate through skid list, get annotations for the first skid in every cell type
#check for the presence of the annotations and add the annotation to the type_of_cell list
type_of_cell <- list()
for(i in seq_along(skids_celltypelist)){
  annot <- catmaid_get_annotations_for_skeletons(skids=skids_celltypelist[[i]][1], pid = 11)
  for (j in seq_along(annot_to_search)) {
    if(sum(annot[,'annotation'] %in% annot_to_search[j]) ==1 )
    {
      type_of_cell[i] <- annot_to_search[j]
      break()
    } else {type_of_cell[i] <- "NULL" }
  }
}


# convert to tibble and plot conn matrix ----------------------------------

syn_tb <- as_tibble(synapse_matrix) %>%
  mutate(presyn_cell = unlist(celltype_names)) %>%
  mutate(type_presyn = unlist(type_of_cell)) %>%
  pivot_longer(-c(presyn_cell, type_presyn), names_to = "postsyn_cell", 
               values_to = "synapses") %>%
  filter(synapses>0) %>%
  group_by(postsyn_cell) %>%
  mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))

write_csv(syn_tb, "supplements/Suppl_celltype_connectivity_tibble.csv")
syn_tb <- read_csv("supplements/Suppl_celltype_connectivity_tibble.csv")

#these are the neurons
celltype_names[1:200]
#plot 
syn_plot <- syn_tb %>%
  ggplot(aes(x=postsyn_cell, y=presyn_cell
             )) +
  geom_tile(colour = "grey50", size = 0.02, fill = NA) +
  geom_point(aes(size = sqrt(synapses), color = synapse_fraction), 
             shape = 15, stroke = 0) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text (angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"), 
        axis.text.y = element_text (angle = 0, hjust = 1, vjust = 0.5, size = 10, color = "black"),
        axis.title.x = element_text (size=15),
        axis.title.y = element_text (size=15)
        ) +
  #order data as in the celltype_names list
  scale_x_discrete(limits = as.character(celltype_names)) +
  scale_y_discrete(limits = as.character(rev(celltype_names[1:200]))) +
  scale_colour_gradient2(
    low = "white",
    mid = "#0072B2",
    high ="#D55E00",
    midpoint = 0.5,
    space = "Lab",
    na.value = "white",
    guide = "colourbar",
    aesthetics = "colour"
)


ggsave("supplements/Suppl_celltype_connectivity_matrix.pdf", limitsize = FALSE, 
       units = c("px"), syn_plot, width = 14000, height = 10000, bg='white')


