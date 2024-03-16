#This code was used to generate the full connectivity matrix of the 3 day old Platynereis larva described in Veraszto et al. 2021
#Gaspar Jekely 2021 Feb

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)
library(data.table)

#set working directory
setwd('/Users/gj274/OneDrive\ -\ University\ of\ Exeter/Paper/Connectome/Supplements_source_data/Full_Connectome/')

# catmaid connection, needs username, password AND token
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info


#retrieve all skeletons that are part of the connectome
connectome_cells <- read.neurons.catmaid("^connectome$", pid=11, fetch.annotations = FALSE) 

#retrieve all annotations for the same skeletons
connectome_annotations <- catmaid_get_annotations_for_skeletons("^connectome$", pid = 11)

#retrieve all by all connectivity information
connectome_connectivity <- catmaid_get_connectors_between(pre=names(connectome_cells), 
                                                          post=names(connectome_cells), pid=11)
#retrieve all skeleton names
skeleton_names <- catmaid_get_neuronnames(names(connectome_cells), pid=11)

#list of skids
names(connectome_cells)[1]
#list of names
skeleton_names


names_to_skids <- data.table(names(connectome_cells),skeleton_names)
colnames(names_to_skids) <- c('skids','skeleton names')
names_to_skids[1:2,1:2]

# use table() to cross-tabulate number of connections between skids
connectome_conn_table = table(connectome_connectivity[,1:2])
connectome_conn_table[1:10,1:10]
dim(connectome_conn_table)


counter <- 0
col_names=list(1:length(colnames(connectome_conn_table)))
for (i in colnames(connectome_conn_table)){
  counter <- counter+1
  index <- grep(i, names(connectome_cells))
  col_names[counter] <- skeleton_names[index]
}

length(col_names)
length(colnames(connectome_conn_table))

counter <- 0
row_names=list(1:length(rownames(connectome_conn_table)))
for (i in rownames(connectome_conn_table)){
  counter <- counter+1
  index <- grep(i, names(connectome_cells))
  row_names[counter] <- skeleton_names[index]
}

length(row_names)
length(rownames(connectome_conn_table))


write.csv(connectome_conn_table, file = "connectome_conn_table_skids.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")


rownames(connectome_conn_table) <- row_names
colnames(connectome_conn_table) <- col_names
connectome_conn_table[1:10,1:10]

write.csv(connectome_conn_table, file = "connectome_conn_table_names.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")

write.csv(names_to_skids, file = "connectome_skids_to_names.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")

write.csv(connectome_annotations, file = "connectome_annotations.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")

connectome_annotations
