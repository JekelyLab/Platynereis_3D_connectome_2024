#R/natverse code to generate Figure 9 fig suppl 2 for the Platynereis 3d connectome paper
#Gaspar Jekely Feb 2022

#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

#annotations of all MB celltypes and their direct pre and postsyn partners
MB_annotations <- c("^celltype18$", "^celltype20$","^celltype187$",
                    "^celltype110$","^celltype131$","^celltype138$","^celltype133$","^celltype137$",
                    "^INMBtype1$","^INMBtype2$","^INMBtype3$","^INMBtype4$","^INMBtype5$",
                    "^INMBtype6$","^INMBtype7$","^INMBtype8$","^INMBtype9$","^INMBtype10$",
                    "^celltype183$","^celltype191$","^celltype185$","^celltype184$","^celltype58$",
                    "^celltype195$","^celltype194$","^celltype140$","^celltype192$","^MBprojMouth$",
                    "^celltype17$","^celltype167$","^celltype168$","^celltype5$","^palp sensory neuron$","^antennae_cell$",
                    "^celltype199$","^celltype153$",
                    "^celltype196$","^celltype127$","^celltype21$", "^celltype117$","^celltype116$",
                    "^celltype125$", "^celltype119$", 
                    "^celltype19$")

MB_celltypes <- c("SNhook","SNlasso","SNtorii",
                  "SN_NS1","SN_NS17","SN_NS18","SN_NS19","SN_NS27",
                  "INMBtype1","INMBtype2","INMBtype3","INMBtype4","INMBtype5",
                  "INMBtype6","INMBtype7","INMBtype8","INMBtype9","INMBtype10",
                  "INhorn","INtorii","INMBdescFMRF","INMBPDF","INrope",
                  "INMBdesc2","INMBdesc3","INbigloop","INUturnMB","MBprojMouth",
                  "SNhorn","SNstiff","SNbronto","cPRC","palp","antenna",
                  "INdecussM","INdecusshook",  
                  "INpara", "INZ","INlasso", "INW", "INfoot", 
                  "INproT2", "INhook", 
                  "MNant")

length(MB_annotations)
length(MB_celltypes)

#retrieve skids matching two annotations (left and right side of MB neurons separately)
skids_left <- lapply(MB_annotations, function(x) skids_by_2annotations(x, "left_side"))
skids_right <- lapply(MB_annotations, function(x) skids_by_2annotations(x, "right_side"))

#use sapply to iterate through the left and right skids list and get all to all synaptic connectivity
left_synapses <- sapply(skids_left, function(x) lapply(skids_left, function(y) 
  catmaid_get_connectors_between(pre_skids = unlist(x), 
                               post_skids = unlist(y), pid=11)) )

right_synapses <- sapply(skids_right, function(x) lapply(skids_right, function(y) 
  catmaid_get_connectors_between(pre_skids = unlist(x), 
                                 post_skids = unlist(y), pid=11)) )

#count the number of synapses
left_synapse_list <- lapply(left_synapses, function(x) length(x$connector_id))
right_synapse_list <- lapply(right_synapses, function(x) length(x$connector_id))

#convert synapse list into a matrix of appropriate dimensions
synapse_matrix_l = matrix(unlist(left_synapse_list), byrow=TRUE, nrow=length(MB_annotations))
rownames(synapse_matrix_l) <- MB_celltypes
colnames(synapse_matrix_l) <- MB_celltypes
synapse_matrix_l

synapse_matrix_r = matrix(unlist(right_synapse_list), byrow=TRUE, nrow=length(MB_annotations))
rownames(synapse_matrix_r) <- MB_celltypes
colnames(synapse_matrix_r) <- MB_celltypes
synapse_matrix_r

write.csv(synapse_matrix_l, 
          "source_data/Figure9_fig_suppl2_source_data3.txt"
)
write.csv(synapse_matrix_r, 
          "source_data/Figure9_fig_suppl2_source_data4.txt"
)
#plot with ggplot
{
pl1 <-   as.data.frame(synapse_matrix_l) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
    ggplot() +
    geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction), 
               stroke = 0)  + 
    theme(
      axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
      axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
      axis.title.x = element_text (size=15),
      axis.title.y = element_text (size=15)
    ) +
    labs(x="postsynaptic neuron type",y="presynaptic neuron type",title="left side") +
    scale_size_area(max_size=5) +
    guides(color = 'none', size = 'none')  + 
    #order data as in the celltype_names list
    scale_x_discrete(limits = as.character(MB_celltypes)) +
    scale_y_discrete(limits = as.character(rev(MB_celltypes))) +
    scale_colour_gradient2(
      low = "#0072B2",
      mid = "#D55E00",
      high ="#D55E00",
      midpoint = 0.5,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    )+
    #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
    theme(panel.background = element_rect(fill = "grey98", color = "white"))

pl1 
# Saving R ggplot with R ggsave Function
ggsave("pictures/MBleft_syn_matrix.png", 
         width = 2500, 
         height = 2500, limitsize = TRUE, 
         units = c("px"))
}


{
pl2 <-   as.data.frame(synapse_matrix_r) %>%
    rownames_to_column(var = "presyn_cell_group") %>%
    pivot_longer(-presyn_cell_group, names_to = "postsyn_cell_group", values_to = "synapses")%>%
    group_by(postsyn_cell_group) %>%
    mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE))%>%
    ggplot() +
    geom_point(aes(x = postsyn_cell_group, y = presyn_cell_group, size = sqrt(synapses), color = synapse_fraction), 
               stroke = 0)  + 
    theme(
      axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
      axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
      axis.title.x = element_text (size=15),
      axis.title.y = element_text (size=15)
    ) +
    labs(x="postsynaptic neuron type",y="presynaptic neuron type",title="right side") +
    scale_size_area(max_size=5) +
    guides(color = 'legend')  + 
    #order data as in the celltype_names list
    scale_x_discrete(limits = as.character(MB_celltypes)) +
    scale_y_discrete(limits = as.character(rev(MB_celltypes))) +
    scale_colour_gradient2(
      low = "#0072B2",
      mid = "#D55E00",
      high ="#D55E00",
      midpoint = 0.5,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    )+
    #        scale_fill_gradient2(low = "#74ADD1", mid= "#FDAE61", high = "#D73027") +
    theme(panel.background = element_rect(fill = "grey98", color = "white"))
pl2  
# Saving R ggplot with R ggsave Function
ggsave("pictures/MBright_syn_matrix.png", 
         width = 2900, 
         height = 2500, limitsize = TRUE, 
         units = c("px"))
}


# correlation calculations ---------------

#calculate pearson between left and right-side synapse lists/matrices
length(left_synapse_list)

cor.test(synapse_matrix_l,synapse_matrix_r, method = "pearson")
cor.test(unlist(left_synapse_list),unlist(right_synapse_list), method="pearson")

#randomised
cor.test(unlist(sample(left_synapse_list)),unlist(sample(right_synapse_list)), method="pearson")

# Sholl analysis ----------------------

{
#define empty matrices
sholl_left_mat <- matrix(data=NA, nrow=length(skids_left), ncol=100)
sholl_right_mat <- matrix(data=NA, nrow=length(skids_left), ncol=100)

for (i in 1:length(skids_left)){    #iterate through the cell type list
  #read the neuron list based on the left side skids
  left_neurons <- nlapply(read.neurons.catmaid(unlist(skids_left[i]), pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  right_neurons <- nlapply(read.neurons.catmaid(unlist(skids_right[i]), pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
  
  print(left_neurons)
  
  #do sholl analysis on left and right cells separately
  sholl_left <- sholl_analysis(
    left_neurons,
    start = soma(left_neurons),
    starting.radius = 1000,
    ending.radius = 100000,
    radius.step = 1000
  )
  
  print(right_neurons)
  
  sholl_right <- sholl_analysis(
    right_neurons,
    start = soma(right_neurons),
    starting.radius = 1000,
    ending.radius = 100000,
    radius.step = 1000
  )
  
  sholl_intersections_left <- lapply(sholl_left, function(x) (x$intersections))
  sholl_means_left <- rowMeans(as.data.frame(sholl_intersections_left[1:length(sholl_intersections_left)]))
  
  sholl_intersections_right <- lapply(sholl_right, function(x) (x$intersections))
  sholl_means_right <- rowMeans(as.data.frame(sholl_intersections_right[1:length(sholl_intersections_right)]))
  
  print (i)
  sholl_left_mat[i,] <- sholl_means_left
  sholl_right_mat[i,] <- sholl_means_right
}   


#plot all sholl plots
matplot(t(sholl_right_mat), type = "l")

#define smoothing function
lowess_fun <- function(x,f=0){
  return(lowess(x, f=0.04, iter=100))
}

lowess_sholl_left <- apply(sholl_left_mat, 1, lowess_fun,f=0.04)
lowess_sholl_right <- apply(sholl_right_mat, 1, lowess_fun,f=0.04)

#test plot
plot(lowess_sholl_left[[2]])

#create smoothed matrix
sholl_left_smoothed_mat <- sholl_left_mat
sholl_right_smoothed_mat <- sholl_right_mat

#add smoothed values
for (i in c(1:length(lowess_sholl_left))){
  sholl_left_smoothed_mat[i,] <- lowess_sholl_left[[i]]$y
  sholl_right_smoothed_mat[i,] <- lowess_sholl_right[[i]]$y
}

#plot the smoothed sholl diagrams
matplot(t(sholl_left_smoothed_mat), type = "l")
sholl_left_smoothed_mat

#assign rownames
rownames(sholl_left_smoothed_mat) <- MB_celltypes
rownames(sholl_right_smoothed_mat) <- MB_celltypes

#remove all-zero rows
sholl_left_smoothed_mat.no0 = sholl_left_smoothed_mat[ rowSums(sholl_left_smoothed_mat)!=0, ] 
sholl_right_smoothed_mat.no0 = sholl_right_smoothed_mat[ rowSums(sholl_right_smoothed_mat)!=0, ] 

#plot heatmap with heatmaply
library(heatmaply)

#plot as heatmap
hm1 <- as.data.frame((sholl_left_smoothed_mat.no0)) %>%
  rownames_to_column(var = "neuron")%>%
  pivot_longer(-neuron,  values_to = "Sholl")%>%
  mutate(distance = sub('V','',name))%>%
  ggplot() +
  geom_raster(aes(x = distance, y = neuron, fill = Sholl))+ 
  theme(
    axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
    axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
    axis.title.x = element_text (size=15),
    axis.title.y = element_text (size=15)
  ) +
  #order data as in the celltype_names list
  scale_x_discrete(limits = factor(1:100), breaks=seq(10, 100, by=10)) +
  scale_y_discrete(limits = as.character(rev(MB_celltypes))) +
  labs(x="distance from soma (µm)",y="neuron type", title="left side") +
  scale_size_area(max_size=5) +
  guides(fill = 'none') +
  scale_fill_gradientn(
    colours=c('white','#0072B2','#E69F00','#E69F00','#E69F00','#D55E00','#D55E00','#D55E00'),
    space = "Lab"
)
hm1

# Save image
ggsave("pictures/MBtypes_left_Sholl.png", 
       width = 2500, 
       height = 2500, limitsize = TRUE, 
       units = c("px"))

hm2 <- as.data.frame((sholl_right_smoothed_mat.no0)) %>%
  rownames_to_column(var = "neuron")%>%
  pivot_longer(-neuron,  values_to = "Sholl")%>%
  mutate(distance = sub('V','',name))%>%
  ggplot() +
  geom_raster(aes(x = distance, y = neuron, fill = Sholl))+ 
  theme(
    axis.text.x = element_text (angle = 90,hjust = 1, vjust = 0.5, size=10), 
    axis.text.y = element_text (angle = 0,hjust = 1, vjust = 0.5, size=10),
    axis.title.x = element_text (size=15),
    axis.title.y = element_text (size=15)
  ) +
  #order data as in the celltype_names list
  scale_x_discrete(limits = factor(1:100), breaks=seq(10, 100, by=10)) +
  scale_y_discrete(limits = as.character(rev(MB_celltypes))) +
  labs(x="distance from soma (µm)",y="neuron type",title="right side") +
  scale_size_area(max_size=5) +
  guides(color = 'legend') +
  scale_fill_gradientn(
    colours=c('white','#0072B2','#E69F00','#E69F00','#E69F00','#D55E00','#D55E00','#D55E00'),
    space = "Lab"
  )
hm2

#save source data
write.csv(sholl_left_smoothed_mat.no0, 
          "source_data/Figure9_fig_suppl2_source_data1.txt"
)
write.csv(sholl_right_smoothed_mat.no0, 
          "source_data/Figure9_fig_suppl2_source_data2.txt"
)

# Save image
ggsave("pictures/MBtypes_right_Sholl.png", 
       width = 2900, 
       height = 2500, limitsize = TRUE, 
       units = c("px"))

}

# make figure ---------------------

layout="AB
##
CD"


img1 <- readPNG("pictures/MBtypes_left_Sholl.png")
panelA <- ggdraw() + draw_image(img1, scale = 1)
img2 <- readPNG("pictures/MBtypes_right_Sholl.png")
panelB <- ggdraw() + draw_image(img2, scale = 1)
img3 <- readPNG("pictures/MBleft_syn_matrix.png")
panelC <- ggdraw() + draw_image(img3, scale = 1)
img4 <- readPNG("pictures/MBright_syn_matrix.png")
panelD <- ggdraw() + draw_image(img4, scale = 1)


Fig_MB_left_right <- panelA + panelB + panelC + panelD +
  plot_layout(design = layout, heights = c(1,0.05,1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face='plain'))


ggsave("Figures/Figure9_fig_Suppl2.pdf", limitsize = FALSE, 
       units = c("px"), Fig_MB_left_right, width = 5400, height = 4500)


ggsave("Figures/Figure9_fig_Suppl2.png", limitsize = FALSE, 
       units = c("px"), Fig_MB_left_right, width = 5400, height = 4500, bg='white')
