# load nat and all associated packages, incl catmaid
library(natverse)
library(nat)
source("~/R/conn.R")
library(heatmaply)
#use multiple cores
library(doMC)
library(parallel)

####################################
#this section was used to test the potential_synapses algorhythm on PRC and IN1 cells

#read the neurons from Catmaid, without smoothing with sigma
PRC = read.neurons.catmaid('celltype1$', pid=11, fetch.annotations = F)
IN1 = read.neurons.catmaid('celltype2$', pid=11, fetch.annotations = F)

## S3 method for class 'neuron'
plot3d(PRC, WithLine = TRUE, WithNodes = TRUE,PlotSubTrees = TRUE,add = TRUE, col = NULL, soma = T)
plot3d(IN1, WithLine = TRUE, WithNodes = TRUE,PlotSubTrees = TRUE,add = TRUE, col = NULL, soma = T)


# Use the detectCores() function to find the number of cores in system
no_cores <- detectCores()
no_cores
registerDoMC(cores = no_cores)

#Calculate number of potential synapses between two neurons  
#This implements the method of Stepanyants and Chklovskii 

system.time( 
potential_synapses_PRC_IN1 <- potential_synapses(
  PRC,
  IN1,
  s=160,  #calibrated based on PRCs and IN1s 
  method = c("approx"), 
  .parallel=TRUE
)

)

system.time(PRC_d2<-dotprops(PRC,resample=1,k=5,.parallel=TRUE))


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
  }
}

heatmaply(PRC_IN1_conn_mat,
          colors = c("white","#F28D24","#AEC7E8","#2078B5","#2078B5"),
          column_text_angle = 90,row_text_angle = 0,
          Rowv=NA,Colv = NA,
          show_dendrogram = c(F, F),
          grid_gap = 1,
          hide_colorbar = F,
          plot_method = c("ggplot"),
          fontsize_row = 10,
          fontsize_col = 10,
          revC=F
)

#transpose matrix to make it compatible with synapse matrix
potential_synapses_PRC_IN1 <- t(potential_synapses_PRC_IN1)

heatmaply(potential_synapses_PRC_IN1/4,  #division by 4 scales predicted synapse numbers to real synapse numbers
          colors = c("white","#F28D24","#AEC7E8","#2078B5","#2078B5"),
          column_text_angle = 90,row_text_angle = 0,
          Rowv=NA,Colv = NA,
          show_dendrogram = c(F, F),
          grid_gap = 1,
          hide_colorbar = F,
          plot_method = c("ggplot"),
          fontsize_row = 10,
          fontsize_col = 10,
          revC=F
)

##############################################################
#this section was used to test which s value gives the highes correlation 
#calculate correlation by converting to vectors first
cor(c(PRC_IN1_conn_mat), c(potential_synapses_PRC_IN1))

corr_synapses_100_2000 <- list()
for(i in seq(100,300, by=10)){
  potential_synapses_PRC_IN1 <- potential_synapses(
  PRC,
  IN1,
  s=i,  #calibrated based on PRCs and IN1s, maximum correlation at 170 nm
  method = c("approx")
  )
  print (i)
  print(cor(c(PRC_IN1_conn_mat), c(t(potential_synapses_PRC_IN1))))
  corr_synapses_100_2000[i] <- cor(c(PRC_IN1_conn_mat), c(t(potential_synapses_PRC_IN1)))
}



#############################################################
#calculate potential synapses for all neurons categorised as celltypes in the connectome
#############################################################

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#set working directory
setwd('/Figures/Suppl_Figure_predicted_synapses')

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
library(heatmaply)

#plot heatmap
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

