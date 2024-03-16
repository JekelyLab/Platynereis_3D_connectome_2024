#This code was used to generate the matrix sowing connectivity by segment and cell type
#Gaspar Jekely 2021 Jan

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
Sys.setenv('R_MAX_VSIZE'=8000000000)
options(timeout = 4000000) 

# load nat and all associated packages, incl catmaid
library(natverse)

#load RColorBrewer for color palette 
library(RColorBrewer)

#set working directory
setwd("~/git/Veraszto_et_al_2021_connectome/")
mainDir = getwd()
dir.create(file.path(mainDir, "statistics"), showWarnings = FALSE)
dir.create(file.path(mainDir, "statistics/segmental_celltype_connectivity"), showWarnings = FALSE)
setwd('statistics/segmental_celltype_connectivity')


# catmaid connection, needs username, password AND token - weird!
# can run this separate file using source function
source("~/R/conn.R")
# best practice is to store this info in your .Renviron file and R will
# automatically read and store it on start-up - you don't have to deal with it,
# and your code won't contain potentially compromising info

#this is a function to find common substrings in strings - only needed if the desired celltype name is the intersect of all names in the set 
intersect2 <- function (x, y)  
{
  y <- as.vector(y)
  y[match(as.vector(x), y, 0L)]
}

celltypelist = list()
annotation_celltypelist = list()
number_of_neuron_types=0

#first we read all celltypes from 1-182 and all annotations
for (i in c(1:182)){
  number_of_neuron_types <- number_of_neuron_types + 1
  annotation = paste("annotation:^celltype", i, "$", sep="")
    #read presyn neuron group by annotation, also retrieve all annotation
  celltypelist[[i]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = FALSE) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
  }

#we read all non-neuronal celltypes from 1-90 and all annotations
for (i in c(1:90)){
  annotation = paste("annotation:^celltype_non_neuronal", i, "$", sep="")
  #read presyn neuron group by annotation, also retrieve all annotation
  celltypelist[[i+number_of_neuron_types]] <- read.neurons.catmaid(annotation, pid=11, fetch.annotations = FALSE) #we add the next celltype to the celltype list
  #retrieve all annotations for the same neurons and create the annotations data frames
  annotation_celltypelist[[i+number_of_neuron_types]] <- catmaid_get_annotations_for_skeletons(annotation, pid = 11)
}

#define the six body regions, matching the catmaid annotations
body_regions <- c('episphere','segment_0', 'segment_1', 'segment_2', 'segment_3', 'pygidium')

#define empty synapse list with the right dimensions
synapse_list <- vector("list", length(annotation_celltypelist)*length(body_regions)*length(annotation_celltypelist)*length(body_regions))

#we iterate through the celltype lists and body regions and retrieve connectivity between pre- and postsyn sets of skids
list_position=0
for (i in c(1:6)){  #iterate through the six body regions
  body_region_i = body_regions[i]
  cycle=0
  print (body_regions[i])
  for (df1 in annotation_celltypelist){    #iterate through the celltype list  (all presyn cells per body region)
    presyn_skids <- df1[df1$annotation == body_region_i,1]
    cycle <- cycle + 1
    print (cycle)
    
    for (j in c(1:6)){  #nested iteration through the six body regions (all postsyn cells)
      body_region_j = body_regions[j] 
        for (df2 in annotation_celltypelist){  #nested iteration through the celltype list (postsyn cells)
        postsyn_skids <- df2[df2$annotation == body_region_j,1]
        assign("cellgroup_conn", NULL, envir = .GlobalEnv)  #in every iteration we empty the connectivity list
        # get connectors betwwen neurons of interest
        if (length(presyn_skids)==0 | length(postsyn_skids)==0) {
          N_synapses <- 0
      } else {
          cellgroup_conn = catmaid_get_connectors_between(pre=presyn_skids, post=postsyn_skids, pid=11)
          N_synapses=nrow(cellgroup_conn)
         if(length(cellgroup_conn) == 0) {N_synapses <- 0}
        }
        list_position <- list_position + 1
        synapse_list [[list_position]] <- N_synapses
        }
    }
    }
print(list_position)
  }

#convert synapse list into a matrix of appropriate dimensions
length(synapse_list)/(length(annotation_celltypelist)*6)
synapse_matrix = matrix(unlist(synapse_list), byrow=TRUE, nrow=length(annotation_celltypelist)*length(body_regions) )

#we make a celltype name list using the first cell in every set (not perfect) and _body_region
celltype_names=list()
for (i in 1:6){
for (df in celltypelist){
#to retrieves all the neuron names in one element of the celltype list
  neuro_names <- as.character(attr(df,"df")$name)
  #to retrieve the common characters in the neuron names (as generic celltype name)
  #print (paste(Reduce(intersect2, strsplit(neuro_names, NULL)), collapse = ''))
  #celltype_names[[length(celltype_names) + 1]] <- paste(Reduce(intersect2, strsplit(neuro_names, NULL)), collapse = '')
  celltype_names[[length(celltype_names) + 1]] <- paste(neuro_names[1], body_regions[i], sep = "_")
}
}

synapse_matrix = as.data.frame(synapse_matrix)

#assign column names to matrix
synapse_matrix=setNames(synapse_matrix, as.character(celltype_names))

#assign row names to matrix
rownames(synapse_matrix) <- as.character(celltype_names)
synapse_matrix = as.matrix(synapse_matrix)

#we remove all-zero rows and columns
synapse_matrix.no0 = synapse_matrix[ rowSums(synapse_matrix)!=0, ] 
synapse_matrix.no0 = synapse_matrix.no0[ ,colSums(synapse_matrix.no0[1:nrow(synapse_matrix.no0),1:ncol(synapse_matrix.no0)])!=0]

#check dimensions of reduced matrix
ncol(synapse_matrix.no0)
nrow(synapse_matrix.no0)

#save the matrix as pdf
pdf('Segment_sorted_connectivity_matrix.pdf', width=7.5, height=7)
heatmap((synapse_matrix.no0),   #show matrix 
        Rowv=NA, Colv=NA,
        cexRow = 0.07, cexCol = 0.07, revC=T, 
        scale = 'none',
        col=hcl.colors(300, "Oslo", alpha = 1, rev = FALSE, fixup = TRUE), #col=brewer.pal(9, 'YlOrRd'), #col=terrain.colors(500, rev=T),
        symm=F, margins= c(3,2))
dev.off()


write.csv(synapse_matrix.no0, file = "synapse_matrix_by_celltype_and_segment.csv",
          quote = FALSE,
          eol = "\n", na = "NA",
          fileEncoding = "")


#head
sum(synapse_matrix.no0[1:119,1:131])
sum(synapse_matrix.no0[120:126,1:131])
sum(synapse_matrix.no0[127:162,1:131])
sum(synapse_matrix.no0[163:203,1:131])
sum(synapse_matrix.no0[204:234,1:131])
sum(synapse_matrix.no0[235:248,1:131])


#segment_0
(synapse_matrix.no0[120:126,132:145])

sum(synapse_matrix.no0[1:119,132:145])
sum(synapse_matrix.no0[120:126,132:145])
sum(synapse_matrix.no0[127:162,132:145])
sum(synapse_matrix.no0[163:203,132:145])
sum(synapse_matrix.no0[204:234,132:145])
sum(synapse_matrix.no0[235:248,132:145])



#segment_1
(synapse_matrix.no0[127:162,146:213])

sum(synapse_matrix.no0[1:119,146:213])
sum(synapse_matrix.no0[120:126,146:213])
sum(synapse_matrix.no0[127:162,146:213])
sum(synapse_matrix.no0[163:203,146:213])
sum(synapse_matrix.no0[204:234,146:213])
sum(synapse_matrix.no0[235:248,146:213])



#segment_2
(synapse_matrix.no0[163:203,214:292])

sum(synapse_matrix.no0[1:119,214:292])
sum(synapse_matrix.no0[120:126,214:292])
sum(synapse_matrix.no0[127:162,214:292])
sum(synapse_matrix.no0[163:203,214:292])
sum(synapse_matrix.no0[204:234,214:292])
sum(synapse_matrix.no0[235:248,214:292])


#segment_3
(synapse_matrix.no0[204:234,293:363])

sum(synapse_matrix.no0[1:119,293:363])
sum(synapse_matrix.no0[120:126,293:363])
sum(synapse_matrix.no0[127:162,293:363])
sum(synapse_matrix.no0[163:203,293:363])
sum(synapse_matrix.no0[204:234,293:363])
sum(synapse_matrix.no0[235:248,293:363])

#pygidium
(synapse_matrix.no0[235:248,364:380])

sum(synapse_matrix.no0[1:119,364:380])
sum(synapse_matrix.no0[120:126,364:380])
sum(synapse_matrix.no0[127:162,364:380])
sum(synapse_matrix.no0[163:203,364:380])
sum(synapse_matrix.no0[204:234,364:380])
sum(synapse_matrix.no0[235:248,364:380])


#separating MUS from neurons for the three main segments

#head
(synapse_matrix.no0[1:114,1:114])

#pygidium
#neurons
(synapse_matrix.no0[1,364:376])

#neurons
(synapse_matrix.no0[1,1:114])

#segment_1
(synapse_matrix.no0[127:162,146:213])

#neurons:
(synapse_matrix.no0[1,146:177])
#MUS:
(synapse_matrix.no0[1,184:212])

#neurons sg1 input from body regions:
sum(synapse_matrix.no0[1:119,146:177])
sum(synapse_matrix.no0[120:126,146:177])
sum(synapse_matrix.no0[127:162,146:177])
sum(synapse_matrix.no0[163:203,146:177])
sum(synapse_matrix.no0[204:234,146:177])
sum(synapse_matrix.no0[235:248,146:177])

#MUS sg1 input from body regions:
sum(synapse_matrix.no0[1:119,184:212])
sum(synapse_matrix.no0[120:126,184:212])
sum(synapse_matrix.no0[127:162,184:212])
sum(synapse_matrix.no0[163:203,184:212])
sum(synapse_matrix.no0[204:234,184:212])
sum(synapse_matrix.no0[235:248,184:212])


#segment_2
(synapse_matrix.no0[163:203,214:292])


#neurons:
(synapse_matrix.no0[1,214:250])
#MUS:
(synapse_matrix.no0[1,258:291])

#neurons sg1 input from body regions:
sum(synapse_matrix.no0[1:119,214:250])
sum(synapse_matrix.no0[120:126,214:250])
sum(synapse_matrix.no0[127:162,214:250])
sum(synapse_matrix.no0[163:203,214:250])
sum(synapse_matrix.no0[204:234,214:250])
sum(synapse_matrix.no0[235:248,214:250])

#MUS sg1 input from body regions:
sum(synapse_matrix.no0[1:119,258:291])
sum(synapse_matrix.no0[120:126,258:291])
sum(synapse_matrix.no0[127:162,258:291])
sum(synapse_matrix.no0[163:203,258:291])
sum(synapse_matrix.no0[204:234,258:291])
sum(synapse_matrix.no0[235:248,258:291])



#segment_3
(synapse_matrix.no0[204:234,293:363])

#neurons:
(synapse_matrix.no0[1,293:320])
#MUS:
(synapse_matrix.no0[1,332:362])

#neurons sg1 input from body regions:
sum(synapse_matrix.no0[1:119,293:320])
sum(synapse_matrix.no0[120:126,293:320])
sum(synapse_matrix.no0[127:162,293:320])
sum(synapse_matrix.no0[163:203,293:320])
sum(synapse_matrix.no0[204:234,293:320])
sum(synapse_matrix.no0[235:248,293:320])

#MUS sg1 input from body regions:
sum(synapse_matrix.no0[1:119,332:362])
sum(synapse_matrix.no0[120:126,332:362])
sum(synapse_matrix.no0[127:162,332:362])
sum(synapse_matrix.no0[163:203,332:362])
sum(synapse_matrix.no0[204:234,332:362])
sum(synapse_matrix.no0[235:248,332:362])

#episphere
#neurons
(synapse_matrix.no0[1:114,1:114])

#neurons head input from body regions:
sum(synapse_matrix.no0[1:119,1:114])
sum(synapse_matrix.no0[120:126,1:114])
sum(synapse_matrix.no0[127:162,1:114])
sum(synapse_matrix.no0[163:203,1:114])
sum(synapse_matrix.no0[204:234,1:114])
sum(synapse_matrix.no0[235:248,1:114])

#pygidium
#neurons
(synapse_matrix.no0[1,364:376])

#neurons pyg input from body regions:
sum(synapse_matrix.no0[1:119,364:376])
sum(synapse_matrix.no0[120:126,364:376])
sum(synapse_matrix.no0[127:162,364:376])
sum(synapse_matrix.no0[163:203,364:376])
sum(synapse_matrix.no0[204:234,364:376])
sum(synapse_matrix.no0[235:248,364:376])

#plotting for segment 1
sg1_inputs <- vector("list", 12)


#neurons sg1 input from body regions:
sg1_inputs[1] <- sum(synapse_matrix.no0[1:119,146:177])
sg1_inputs[2] <- sum(synapse_matrix.no0[120:126,146:177])
sg1_inputs[3] <- sum(synapse_matrix.no0[127:162,146:177])
sg1_inputs[4] <- sum(synapse_matrix.no0[163:203,146:177])
sg1_inputs[5] <- sum(synapse_matrix.no0[204:234,146:177])
sg1_inputs[6] <- sum(synapse_matrix.no0[235:248,146:177])

#MUS sg1 input from body regions:
sg1_inputs[7] <- sum(synapse_matrix.no0[1:119,184:212])
sg1_inputs[8] <- sum(synapse_matrix.no0[120:126,184:212])
sg1_inputs[9] <- sum(synapse_matrix.no0[127:162,184:212])
sg1_inputs[10] <- sum(synapse_matrix.no0[163:203,184:212])
sg1_inputs[11] <- sum(synapse_matrix.no0[204:234,184:212])
sg1_inputs[12] <- sum(synapse_matrix.no0[235:248,184:212])

sg1_inputs_mt = matrix(unlist(sg1_inputs), byrow=TRUE, nrow=2)
sg1_inputs_mt[1,] = sg1_inputs_mt[1,]/sum(sg1_inputs_mt[1,])
sg1_inputs_mt[1,]
sg1_inputs_mt[2,] = sg1_inputs_mt[2,]/sum(sg1_inputs_mt[2,])
sg1_inputs_mt[2,]

barplot(sg1_inputs_mt, beside=T, ylim=c(0:1))


#plotting for segment 2
sg2_inputs <- vector("list", 12)

#neurons sg1 input from body regions:
sg2_inputs[1] <- sum(synapse_matrix.no0[1:119,214:250])
sg2_inputs[2] <- sum(synapse_matrix.no0[120:126,214:250])
sg2_inputs[3] <- sum(synapse_matrix.no0[127:162,214:250])
sg2_inputs[4] <- sum(synapse_matrix.no0[163:203,214:250])
sg2_inputs[5] <- sum(synapse_matrix.no0[204:234,214:250])
sg2_inputs[6] <- sum(synapse_matrix.no0[235:248,214:250])

#MUS sg1 input from body regions:
sg2_inputs[7] <- sum(synapse_matrix.no0[1:119,258:291])
sg2_inputs[8] <- sum(synapse_matrix.no0[120:126,258:291])
sg2_inputs[9] <- sum(synapse_matrix.no0[127:162,258:291])
sg2_inputs[10] <- sum(synapse_matrix.no0[163:203,258:291])
sg2_inputs[11] <- sum(synapse_matrix.no0[204:234,258:291])
sg2_inputs[12] <- sum(synapse_matrix.no0[235:248,258:291])


sg2_inputs_mt = matrix(unlist(sg2_inputs), byrow=TRUE, nrow=2)
sg2_inputs_mt[1,] = sg2_inputs_mt[1,]/sum(sg2_inputs_mt[1,])
sg2_inputs_mt[1,]
sg2_inputs_mt[2,] = sg2_inputs_mt[2,]/sum(sg2_inputs_mt[2,])
sg2_inputs_mt[2,]

barplot(sg2_inputs_mt, beside=T, ylim=c(0:1))


#plotting for segment 3
sg3_inputs <- vector("list", 12)

#neurons sg1 input from body regions:
sg3_inputs[1] <- sum(synapse_matrix.no0[1:119,293:320])
sg3_inputs[2] <- sum(synapse_matrix.no0[120:126,293:320])
sg3_inputs[3] <- sum(synapse_matrix.no0[127:162,293:320])
sg3_inputs[4] <- sum(synapse_matrix.no0[163:203,293:320])
sg3_inputs[5] <- sum(synapse_matrix.no0[204:234,293:320])
sg3_inputs[6] <- sum(synapse_matrix.no0[235:248,293:320])

#MUS sg1 input from body regions:
sg3_inputs[7] <- sum(synapse_matrix.no0[1:119,332:362])
sg3_inputs[8] <- sum(synapse_matrix.no0[120:126,332:362])
sg3_inputs[9] <- sum(synapse_matrix.no0[127:162,332:362])
sg3_inputs[10] <- sum(synapse_matrix.no0[163:203,332:362])
sg3_inputs[11] <- sum(synapse_matrix.no0[204:234,332:362])
sg3_inputs[12] <- sum(synapse_matrix.no0[235:248,332:362])


sg3_inputs_mt = matrix(unlist(sg3_inputs), byrow=TRUE, nrow=2)
sg3_inputs_mt[1,] = sg3_inputs_mt[1,]/sum(sg3_inputs_mt[1,])
sg3_inputs_mt[1,]
sg3_inputs_mt[2,] = sg3_inputs_mt[2,]/sum(sg3_inputs_mt[2,])
sg3_inputs_mt[2,]

#plotting for head

head_inputs <- vector("list", 12)

head_inputs[1] <- sum(synapse_matrix.no0[1:119,1:114])
head_inputs[2] <- sum(synapse_matrix.no0[120:126,1:114])
head_inputs[3] <- sum(synapse_matrix.no0[127:162,1:114])
head_inputs[4] <- sum(synapse_matrix.no0[163:203,1:114])
head_inputs[5] <- sum(synapse_matrix.no0[204:234,1:114])
head_inputs[6] <- sum(synapse_matrix.no0[235:248,1:114])


head_inputs_mt = matrix(unlist(head_inputs), byrow=TRUE, nrow=1)
head_inputs_mt[1,] = head_inputs_mt[1,]/sum(head_inputs_mt[1,])
head_inputs_mt[1,]


#plotting for pyg
pyg_inputs <- vector("list", 6)

pyg_inputs[1] <- sum(synapse_matrix.no0[1:119,364:376])
pyg_inputs[2] <- sum(synapse_matrix.no0[120:126,364:376])
pyg_inputs[3] <- sum(synapse_matrix.no0[127:162,364:376])
pyg_inputs[4] <- sum(synapse_matrix.no0[163:203,364:376])
pyg_inputs[5] <- sum(synapse_matrix.no0[204:234,364:376])
pyg_inputs[6] <- sum(synapse_matrix.no0[235:248,364:376])

pyg_inputs_mt = matrix(unlist(pyg_inputs), byrow=TRUE, nrow=1)
pyg_inputs_mt[1,] = pyg_inputs_mt[1,]/sum(pyg_inputs_mt[1,])
pyg_inputs_mt[1,]


#assign row names to matrix
colnames(sg1_inputs_mt) <- c('head','sg0', 'sg1', 'sg2','sg3', 'pyg')
pdf('Inputs_to_sg_1_neurons and MUS.pdf', width=5, height=5)
barplot(sg1_inputs_mt, beside=T, ylim=c(0:1), col =c("#D73027", "#4575B4"), ylab = "% input", space=c(0,0.5), 
        cex.axis = 1,
       cex.names = 1, main = "Inputs to segment 1")
dev.off()

colnames(sg2_inputs_mt) <- c('head','sg0', 'sg1', 'sg2','sg3', 'pyg')
pdf('Inputs_to_sg_2_neurons and MUS.pdf', width=5, height=5)
barplot(sg2_inputs_mt, beside=T, ylim=c(0:1), col =c("#D73027", "#4575B4"), ylab = "% input", space=c(0,0.5), 
        cex.axis = 1,
        cex.names = 1, main = "Inputs to segment 2")
dev.off()

colnames(sg3_inputs_mt) <- c('head','sg0', 'sg1', 'sg2','sg3', 'pyg')
pdf('Inputs_to_sg_3_neurons and MUS.pdf', width=5, height=5)
barplot(sg3_inputs_mt, beside=T, ylim=c(0:1), col =c("#D73027", "#4575B4"), ylab = "% input", space=c(0,0.5), 
        cex.axis = 1, horiz=FALSE,
        cex.names = 1, main = "Inputs to segment 3")
dev.off()

colnames(head_inputs_mt) <- c('head','sg0', 'sg1', 'sg2','sg3', 'pyg')
pdf('Inputs_to_head_neurons.pdf', width=5, height=5)
barplot(head_inputs_mt, beside=T, ylim=c(0:1), col =c("#D73027"), ylab = "% input", 
        cex.axis = 1, horiz=FALSE,
        cex.names = 1, main = "Inputs to head")
dev.off()


colnames(pyg_inputs_mt) <- c('head','sg0', 'sg1', 'sg2','sg3', 'pyg')
pdf('Inputs_to_pyg_neurons.pdf', width=5, height=5)
barplot(pyg_inputs_mt, beside=T, ylim=c(0:1), col =c("#D73027"), ylab = "% input", 
        cex.axis = 1, horiz=FALSE,
        cex.names = 1, main = "Inputs to pygidium")
dev.off()


