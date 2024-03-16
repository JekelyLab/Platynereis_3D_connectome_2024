# Supplemental table 1 of the Platynereis 3d connectome paper
# retrieve cell types connectivity matrix, list the number of cells per celltype, name, annotation
# and main pre and postsynaptic partners

# load packages, functions and catmaid connectivity
source("code/Natverse_functions_and_conn.R")

# retrieve skids
# define empty cell type skids list
celltype_skids <- list()
counts <- 0
# iterate through the neuronal cell types
for (i in c(1:202)) {
  annotation <- paste("annotation:^celltype", i, "$", sep = "")
  # use the sourced skids_by_2annotations function to retrieve all skids per celltype that are
  # also annotated with with_soma
  skids <- skids_by_2annotations(annotation, "with_soma")
  # if no skids are returned, next
  if (identical(skids, integer(0))) {
    print("this celltype is missing!")
  }
  counts <- counts + 1
  celltype_skids[[counts]] <- skids
  print(i)
  print(skids)
}

# define empty cell type non-neuronal skids list
celltype_NN_skids <- list()
counts <- 0
# iterate through the non-neuronal cell types excluding EC, yolk and mesoderm
for (i in c(1:92)) {
  annotation <- paste("annotation:^celltype_non_neuronal", i, "$", sep = "")
  # use the sourced skids_by_2annotations function to retrieve all skids per celltype that are
  # also annotated with with_soma
  skids <- skids_by_2annotations(annotation, "with_soma")
  # if no skids are returned, next
  if (identical(skids, integer(0))) {
    print("this celltype is missing!")
  }
  counts <- counts + 1
  celltype_NN_skids[[counts]] <- skids
  print(i)
  print(skids)
}

all_celltypes_skids <- c(celltype_skids, celltype_NN_skids)

# get the neuron name of the first skid in the skids list
celltype_names <- list()
for (i in 1:length(all_celltypes_skids)) {
  name <- catmaid_get_neuronnames(all_celltypes_skids[[i]][1], pid = 11)
  name <- sub("_.*$", "", name)
  celltype_names[i] <- name
}

# check duplicated names
celltype_names[duplicated(celltype_names)]

# get some parameters from the lists of cell types
number_cells <- unlist(lapply(all_celltypes_skids, function(x) length(x)))

annotations <- list()
# list of annotations cell type 1-202
annotations <- list()
for (i in c(1:202)) {
  annotations[i] <- paste("celltype", i, sep = "")
}
# list of annotations cell type non-neuronal 1-92
annotations2 <- list()
for (i in c(1:92)) {
  annotations2[i] <- paste("celltype_non_neuronal", i, sep = "")
}
# combine the two
annotations <- c(annotations, annotations2)

# retrieve all annotations for all skids per cell type
all_annot_per_celltype <- lapply(all_celltypes_skids, function(x) catmaid_get_annotations_for_skeletons(unlist(x), pid = 11))

all_annot_per_celltype
# convert list of lists to a list of vectors
all_annot_per_celltype <- lapply(all_annot_per_celltype, unlist)
length(all_annot_per_celltype)

# match the list of annotations
# sensory, motor, inter, effectors
{
  toMatch <- c("Sensory neuron")
  SN <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("interneuron")
  IN <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("motorneuron")
  MN <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("muscle")
  Mus <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("ciliated cell")
  cil <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("pigment cell")
  pig <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("gland cell")
  gland <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("glia cell")
  glia <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))

  SN_IN_MN <- c(1:length(annotations))
  SN_IN_MN <- lapply(SN_IN_MN, function(x) c("", "", "", "", "", "", "", ""))
  # write SN_IN_MN list with annotations
  for (i in 1:length(SN_IN_MN)) {
    if (SN[i] == TRUE) SN_IN_MN[[i]][1] <- "Sensory neuron"
    if (IN[i] == TRUE) SN_IN_MN[[i]][2] <- "Interneuron"
    if (MN[i] == TRUE) SN_IN_MN[[i]][3] <- "Motoneuron"
    if (Mus[i] == TRUE) SN_IN_MN[[i]][4] <- "muscle"
    if (cil[i] == TRUE) SN_IN_MN[[i]][5] <- "multiciliated cell"
    if (pig[i] == TRUE) SN_IN_MN[[i]][6] <- "pigment cell"
    if (gland[i] == TRUE) SN_IN_MN[[i]][7] <- "gland cell"
    if (glia[i] == TRUE) SN_IN_MN[[i]][8] <- "glia cell"
  }
  # collapse list into single character string
  SN_IN_MN <- lapply(SN_IN_MN, function(x) paste(unlist(x), collapse = " "))
  SN_IN_MN <- str_trim(SN_IN_MN, side = c("both"))
}

# query body position and match the list of annotations
{
  toMatch <- c("episphere")
  head <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("segment_0")
  sg0 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("segment_1")
  sg1 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("segment_2")
  sg2 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("segment_3")
  sg3 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("pygidium")
  pyg <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))

  head_tail <- c(1:length(annotations))
  head_tail <- lapply(head_tail, function(x) c("", "", "", "", "", ""))
  # write list with annotations
  for (i in 1:length(head_tail)) {
    if (head[i] == TRUE) head_tail[[i]][1] <- "head"
    if (sg0[i] == TRUE) head_tail[[i]][2] <- "sg0"
    if (sg1[i] == TRUE) head_tail[[i]][3] <- "sg1"
    if (sg2[i] == TRUE) head_tail[[i]][4] <- "sg2"
    if (sg3[i] == TRUE) head_tail[[i]][5] <- "sg3"
    if (pyg[i] == TRUE) head_tail[[i]][6] <- "pygidium"
  }
  # collapse list into single character string
  head_tail <- lapply(head_tail, function(x) paste(unlist(x), collapse = " "))
  head_tail <- str_trim(head_tail, side = c("both"))
}

# query ganglia and other anatomical terms
{
  toMatch <- c("nuchal organ")
  nuch <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Dorsal_sensory_organ")
  dso <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Dorsolateral sense organs")
  dlso <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("antennae_cell")
  ant <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Apical_organ")
  AO <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("visual system")
  vis <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("mushroom body")
  mush <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("eyespot_PRC")
  eyesp <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("ciliomotor")
  CM <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("periphery")
  periph <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("VNC_ventral")
  VNC <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Dorsolateral IN cluster")
  DLIN <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("mechanosensory_girdle")
  MGir <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("central_brain")
  CeBr <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("ventral_head")
  VeHe <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("descending_head")
  Desc <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("decussating")
  Decus <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("biaxonal cell")
  biax <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("parapodial muscle complex")
  parapod <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))

  head_glomeruli <- c(1:length(annotations))
  head_glomeruli <- lapply(head_glomeruli, function(x) c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""))
  # write list with annotations
  for (i in 1:length(head_glomeruli)) {
    if (nuch[i] == TRUE) head_glomeruli[[i]][1] <- "nuchal organ"
    if (dso[i] == TRUE) head_glomeruli[[i]][2] <- "dorsal sensory organ"
    if (dlso[i] == TRUE) head_glomeruli[[i]][3] <- "dorsolateral sensory organ"
    if (ant[i] == TRUE) head_glomeruli[[i]][4] <- "antenna"
    if (AO[i] == TRUE) head_glomeruli[[i]][5] <- "apical organ"
    if (vis[i] == TRUE) head_glomeruli[[i]][6] <- "visual system"
    if (mush[i] == TRUE) head_glomeruli[[i]][7] <- "mushroom body"
    if (eyesp[i] == TRUE) head_glomeruli[[i]][8] <- "eyespot"
    if (CM[i] == TRUE) head_glomeruli[[i]][9] <- "ciliomotor"
    if (VNC[i] == TRUE) head_glomeruli[[i]][10] <- "ventral nerve cord"
    if (periph[i] == TRUE) head_glomeruli[[i]][11] <- "peripheral NS"
    if (DLIN[i] == TRUE) head_glomeruli[[i]][12] <- "dorsolateral interneuron cluster"
    if (MGir[i] == TRUE) head_glomeruli[[i]][13] <- "mechanosensory girdle"
    if (CeBr[i] == TRUE) head_glomeruli[[i]][14] <- "central brain"
    if (VeHe[i] == TRUE) head_glomeruli[[i]][15] <- "ventral head"
    if (Desc[i] == TRUE) head_glomeruli[[i]][16] <- "descending"
    if (Decus[i] == TRUE) head_glomeruli[[i]][17] <- "decussating"
    if (biax[i] == TRUE) head_glomeruli[[i]][18] <- "biaxonal neuron"
    if (parapod[i] == TRUE) head_glomeruli[[i]][19] <- "parapodia"
  }
  # collapse list into single character string
  head_glomeruli <- lapply(head_glomeruli, function(x) paste(unlist(x), collapse = " "))
  head_glomeruli <- str_trim(head_glomeruli, side = c("both"))
}

# query transmitter phenotypes
{
  toMatch <- c("dense cored vesicles")
  DCV <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("cholinergic")
  ACh <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("serotonergic")
  Ser <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("glutamatergic")
  Glu <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("adrenergic")
  Adr <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("dopaminergic")
  Dop <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("FVa_positive_neurons")
  FVa <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("FVRIa_positive_neurons")
  FVRI <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("PDF_positive_neurons")
  PDF <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("ATO_positive_neurons")
  ATO <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("FMRFa_positive_neurons")
  FMRF <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("leucokinin_positive_neurons")
  Leu <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("MIP_positive_neuron")
  MIP <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("proenkephalin_positive_neuron")
  PROENK <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("RGWa_positive_neurons")
  RGW <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))

  transmitter <- c(1:length(annotations))
  transmitter <- lapply(transmitter, function(x) c("", "", "", "", "", "", "", "", "", "", "", "", "", "", ""))
  # write list with annotations
  for (i in 1:length(head_glomeruli)) {
    if (DCV[i] == TRUE) transmitter[[i]][15] <- "dense cored vesicles"
    if (ACh[i] == TRUE) transmitter[[i]][2] <- "cholinergic"
    if (Ser[i] == TRUE) transmitter[[i]][3] <- "serotonergic"
    if (Glu[i] == TRUE) transmitter[[i]][4] <- "glutamatergic"
    if (Adr[i] == TRUE) transmitter[[i]][5] <- "adrenergic"
    if (Dop[i] == TRUE) transmitter[[i]][6] <- "dopaminergic"
    if (FVa[i] == TRUE) transmitter[[i]][7] <- "FVamide"
    if (FVRI[i] == TRUE) transmitter[[i]][8] <- "FVRIamide"
    if (PDF[i] == TRUE) transmitter[[i]][9] <- "pigment dispersing factor"
    if (ATO[i] == TRUE) transmitter[[i]][10] <- "orexin/allatotropin"
    if (FMRF[i] == TRUE) transmitter[[i]][11] <- "FMRFamide"
    if (Leu[i] == TRUE) transmitter[[i]][12] <- "leucokinin"
    if (MIP[i] == TRUE) transmitter[[i]][13] <- "myoinhibitory peptide"
    if (PROENK[i] == TRUE) transmitter[[i]][14] <- "proenkephalin"
    if (RGW[i] == TRUE) transmitter[[i]][1] <- "RGWamide"
  }
  # collapse list into single character string
  transmitter <- lapply(transmitter, function(x) paste(unlist(x), collapse = " "))
  transmitter <- str_trim(transmitter, side = c("both"))
}

# query reference annotations
{
  toMatch <- c("Bezares-Calderon_et_al_2018")
  publ1 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  publ1
  toMatch <- c("Williams_et_al_2017")
  publ2 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Veraszto_et_al_2017")
  publ3 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Veraszto_et_al_2018")
  publ4 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Randel_et_al_2015")
  publ5 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Shahidi_et_al_2015")
  publ6 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  toMatch <- c("Jasek_et_al_2022")
  publ7 <- unlist(lapply(all_annot_per_celltype, function(x) all(toMatch %in% x)))
  reference <- c(1:length(annotations))
  reference <- lapply(reference, function(x) c("", "", "", "", "", "", ""))
  # write list with annotations
  for (i in 1:length(reference)) {
    if (publ1[i] == TRUE) reference[[i]][1] <- "Bezares-Calderon et al. 2018"
    if (publ2[i] == TRUE) reference[[i]][2] <- "Williams et al. 2017"
    if (publ3[i] == TRUE) reference[[i]][3] <- "Veraszto et al. 2017"
    if (publ4[i] == TRUE) reference[[i]][4] <- "Veraszto et al. 2018"
    if (publ5[i] == TRUE) reference[[i]][5] <- "Randel et al. 2015"
    if (publ6[i] == TRUE) reference[[i]][6] <- "Shahidi et al. 2015"
    if (publ7[i] == TRUE) reference[[i]][7] <- "Jasek et al. 2022"
  }

  # collapse list into single character string
  reference <- lapply(reference, function(x) paste(unlist(x), collapse = " "))
  reference <- str_trim(reference, side = c("both"))
  # add 'this study' to cells without a reference annotation
  reference <- lapply(reference, function(x) gsub("^$", "this study", x))
}


# convert to character and substitute c()
celltype_skids_char <- gsub("([c(])|([)])", "", as.character(all_celltypes_skids))

# retrieve synaptic connectivity and generate matrix
# define empty synapse list with the right dimensions
celltype_synapse_list <- vector("list", length(celltype_skids) * length(celltype_skids))

# retrieve all against all number of synapses
for (i in 1:length(all_celltypes_skids)) {
  for (j in 1:length(all_celltypes_skids)) {
    # get connectors between two cell lists
    presyn_skids <- unlist(all_celltypes_skids[i])
    postsyn_skids <- unlist(all_celltypes_skids[j])
    connectivity <- catmaid_get_connectors_between(
      pre = presyn_skids,
      post = postsyn_skids, pid = 11
    )
    # check the number of synapses from group1 -> group2
    N_synapses <- dim(connectivity)[1]
    counter <- ((i * length(all_celltypes_skids) - length(all_celltypes_skids)) + j)
    print(counter)
    # add value to synapse list
    print(N_synapses)
    # change NULL to 0
    if (is.null(N_synapses)) {
      N_synapses <- 0
    }
    celltype_synapse_list[[counter]] <- N_synapses
  }
}

# convert synapse list into a matrix of appropriate dimensions
celltype_synapse_matrix <- matrix(unlist(celltype_synapse_list), byrow = TRUE, nrow = length(all_celltypes_skids))

# add row and column names to matrix
rownames(celltype_synapse_matrix) <- unlist(celltype_names)
colnames(celltype_synapse_matrix) <- unlist(celltype_names)

write.csv2(celltype_synapse_matrix, file = "data/all_celltypes_synapse_matrix.csv")
# celltype_synapse_matrix <-  as.matrix(read.csv(file="data/all_celltypes_synapse_matrix.csv", sep=";"))
dim(celltype_synapse_matrix)
celltype_synapse_matrix[1:3, 1:4]


# get all postsyn partners of all cell types with a cut off
all_postsyn <- list()
for (j in c(1:length(annotations))) {
  postsyn_partners <- list()
  for (i in c(1:length(annotations))) {
    if (celltype_synapse_matrix[j, i] > 10) {
      #    print (unlist(celltype_names[i]))
      postsyn_partners <- append(postsyn_partners, unlist(celltype_names[i]))
    }
  }
  print(unlist(postsyn_partners))
  all_postsyn[[j]] <- postsyn_partners
}

# get all presyn partners of all cell types with a cut off
all_presyn <- list()
for (j in c(1:length(annotations))) {
  presyn_partners <- list()
  for (i in c(1:length(annotations))) {
    if (celltype_synapse_matrix[i, j] > 10) {
      #    print (unlist(celltype_names[i]))
      presyn_partners <- append(presyn_partners, unlist(celltype_names[i]))
    }
  }
  print(unlist(presyn_partners))
  all_presyn[[j]] <- presyn_partners
}

# convert list of lists to a list of vectors
all_presyn <- lapply(all_presyn, unlist)
all_postsyn <- lapply(all_postsyn, unlist)
celltype_names
annotations
# construct the tibble with the cell types, numbers and annotations
celltype_tibble <- tibble("name of cell type" = unlist(celltype_names)) %>%
  mutate("number of cells" = number_cells) %>%
  mutate("CATMAID annotation" = unlist(annotations)) %>%
  mutate("Sensory/inter/motor neuron" = str_trim(unlist(SN_IN_MN), side = c("both"))) %>%
  mutate("soma position" = str_trim(unlist(head_tail), side = c("both"))) %>%
  mutate("region" = str_trim(unlist(head_glomeruli), side = c("both"))) %>%
  mutate("transmitter phenotype" = str_trim(unlist(transmitter), side = c("both"))) %>%
  mutate(
    "postsynaptic partners (<10 synapses)" =
      gsub(patt = 'c\\(|\\)|\\,|\\"|NULL', repl = "", all_postsyn)
  ) %>%
  mutate(
    "presynaptic partners (<10 synapses)" =
      gsub(patt = 'c\\(|\\)|\\,|\\"|NULL', repl = "", all_presyn)
  ) %>%
  mutate("reference" = str_trim(unlist(reference), side = c("both")))

# write table
readr::write_csv(celltype_tibble, file = "supplements/Supplementary_Table1.csv", na = "", quote = "none")

