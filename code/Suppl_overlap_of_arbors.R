# load nat and all associated packages, incl catmaid
library(natverse)
library(nat)
library(heatmaply)
library(factoextra)

conn = source("~/R/conn.R")
#for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
#for this we configure to http/1.1
conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))

####################################

#set working directory
setwd("~/git/Veraszto_et_al_2021_connectome/")
mainDir = getwd()
dir.create(file.path(mainDir, "Figures/Suppl_Figure_predicted_synapses"), showWarnings = FALSE)
setwd('Figures/Suppl_Figure_predicted_synapses')

##########################################################
#first we demonstrate the overlap algorythm on PRC and IN1 cells
#read the neurons from Catmaid, without smoothing

PRC = read.neurons.catmaid('celltype1$', pid=11, fetch.annotations = F)
IN1 = read.neurons.catmaid('celltype2$', pid=11, fetch.annotations = F)


#Calculate overlap between neurons  
system.time(overlap_PRC_IN1 <-  overlap_score(PRC, IN1, delta = 160, progress = TRUE))

#get synapses between PRC and IN1 cells
PRC_IN1_connectors <- catmaid_get_connectors_between(pre_skids = names(PRC), post_skids = names(IN1), pid=11)
#tabulate the results
PRC_IN1_connectors_tab=table(PRC_IN1_connectors[,1:2])
# create empty connectivity matrix
PRC_IN1_conn_mat = matrix(0,nrow=length(PRC),ncol=length(IN1))
# name rows and cols for skids
colnames(PRC_IN1_conn_mat) = names(IN1)
rownames(PRC_IN1_conn_mat) = names(PRC)
# populate connectivity matrix based on np_conn_table
for (x in rownames(PRC_IN1_connectors_tab)){
  for (y in colnames(PRC_IN1_connectors_tab)){
    PRC_IN1_conn_mat[x,y] = PRC_IN1_connectors_tab[x,y]
}}

#show heatmaps
heatmaply(overlap_PRC_IN1,Rowv = NA, Colv = NA)
heatmaply(PRC_IN1_conn_mat,Rowv = NA, Colv = NA)

#calculate correlation
cor(c(PRC_IN1_conn_mat), c(overlap_PRC_IN1))

##############################################################
#this section calculates the overlap score between all cells 
#that are annotated as connectome_neuron in the episphere

#function to retrieve skids based on two annotations
skids_by_2annotations <- function(annotation1,annotation2){
  annotations_cells = list()
  annotation1 <- paste("^", annotation1, "$", sep ="")
  annotations_cells[[1]] <- catmaid_get_annotations_for_skeletons(annotation1, pid = 11, conn=conn_http1)
  #we retrieve those skeletons that are also annotated with right_side
  return(unlist(lapply(annotations_cells,function(x) x[x$annotation==annotation2,1])))
}

episphere_skidsMN <- skids_by_2annotations('motorneuron','episphere')
episphere_skidsSN <- skids_by_2annotations('Sensory neuron','episphere')
episphere_skidsIN <- skids_by_2annotations('connectome_interneuron','episphere')
episphere_skids <- c(episphere_skidsMN,episphere_skidsSN,episphere_skidsIN)

#total number of unique episphere skids of SN IN MN cells
length(unique(episphere_skids))

episphere_skeletons <- read.neurons.catmaid(unique(episphere_skids), pid=11, conn=conn_http1)

#calculate overlap between all episphere neurons
system.time(overlap_episphere <- overlap_score(episphere_skeletons, episphere_skeletons, delta=160, progress = TRUE))

#save table
write.csv(overlap_episphere, file = "overlap_episphere.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")


heatmaply(overlap_episphere, Rowv=NA,Colv=NA)

overlap_kmeans <- kmeans(overlap_episphere, centers=50, nstart = 25)

fviz_cluster(overlap_kmeans, data=overlap_episphere)

kmeans()
overlap_episphere_dist = as.dist(overlap_episphere)
overlap_episphere_dist[1:10]

overlap_episphere_clust = hclust(overlap_episphere_dist, method="ward.D2")


# 3d plotting of clustering results
nopen3d(); plot3d(overlap_episphere_clust, k=7, db=episphere_skeletons, soma=F); nview3d("frontal")

#############################################################
#calculate potential synapses for all neurons categorised as celltypes in the connectome
#############################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)
library(heatmaply)

# catmaid connection, needs username, password AND token
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info

#read the neurons from Catmaid, without smoothing with sigma
neurons_in_celltypes = read.neurons.catmaid('^celltype$', pid=11, fetch.annotations = F)

## S3 method for class 'neuron'
plot3d(neurons_in_celltypes, WithLine = TRUE, WithNodes = TRUE,PlotSubTrees = TRUE,add = TRUE, col = NULL, soma = T)
close3d()
#get synapses between PRC and IN1 cells
neurons_in_celltypes_conn <- catmaid_get_connectors_between(pre_skids = names(neurons_in_celltypes), post_skids = names(neurons_in_celltypes), pid=11)

#tabulate the results
neurons_in_celltypes_conn_tab=table(neurons_in_celltypes_conn[,1:2])

# create empty connectivity matrix
neurons_in_celltypes_conn_mat = matrix(0,nrow=length(neurons_in_celltypes),ncol=length(neurons_in_celltypes))

# name rows and cols for skids
colnames(neurons_in_celltypes_conn_mat) = names(neurons_in_celltypes)
rownames(neurons_in_celltypes_conn_mat) = names(neurons_in_celltypes)

# populate connectivity matrix based on np_conn_table
for (x in rownames(neurons_in_celltypes_conn_tab)){
  for (y in colnames(neurons_in_celltypes_conn_tab)){
    neurons_in_celltypes_conn_mat[x,y] = neurons_in_celltypes_conn_tab[x,y]
  }
}

#plot heatmap with heatmaply
heatmaply(sqrt(neurons_in_celltypes_conn_mat),
          column_text_angle = 90,row_text_angle = 0,
          fontsize_row = 3,
          fontsize_col = 3,
          Rowv=F, Colv=F,
          show_dendrogram = c(F, F),
          grid_gap = 0,
          hide_colorbar = F,
          plot_method = c("ggplot"),
          revC=F,
          col=c('grey20','cyan','orange'),
          xlab = "postsynaptic celltypes, neuronal 1-182 and non-neuronal 1-90",
          ylab = "presynaptic neuronal celltypes 1-182",
          main = "connectivity matrix of celltypes (sqrt of summed synapses)",
          
)

write.csv(neurons_in_celltypes_conn_mat, file = "Neurons_in_celltypes_connectivity_matrix.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")


###############################################
#Include the parallel library. If the next line does not work, run install.packages(“parallel”) first
#library(parallel)
# Use the detectCores() function to find the number of cores in system
#no_cores <- detectCores()
#no_cores
# Setup cluster
#clust <- makeCluster(no_cores) #This line will take time
#it is important to close the cluster at the end of execution step so that core memory is released
#stopCluster(clust)
###############################################


potential_synapses_celltypes <- potential_synapses(
  neurons_in_celltypes,
  neurons_in_celltypes,
  s=160,  #calibrated based on PRCs and IN1s 
  method = c("approx"), 
  .parallel=TRUE
)


#transpose matrix to make it compatible with the synaptic matrix
potential_synapses_celltypes <- t(potential_synapses_celltypes)

