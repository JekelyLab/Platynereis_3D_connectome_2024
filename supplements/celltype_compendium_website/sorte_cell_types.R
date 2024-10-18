#filter celltypes

library(dplyr)
celltypes <- rio::import("supplements/Supplementary_Table1.csv")

#sensory cells
celltypes |>
  filter(`Sensory/inter/motor neuron` == "Sensory neuron") |>
  filter(region != "mushroom body" & 
           region !=  "visual system" & 
           region !=  "apical organ" & 
           region !=  "nuchal organ" &
           region !=  "eyespot" &
           !grepl("girdle", region)) |>
  filter(!grepl("DLSO", `name of cell type`)) |>
  select(`name of cell type`, `CATMAID annotation`)


#interneuron
celltypes |>
  filter(`Sensory/inter/motor neuron` == "Interneuron") |>
#  filter(`soma position` == "pygidium") |>
  filter(region !=  "visual system" & region !=  "apical organ" & region !=  "ciliomotor" &
           !grepl("central brain", region) &
           !grepl("girdle", region) &
           !grepl("dorsolateral", region) &
           !grepl("mushroom", region)) |>
  filter(!grepl("DLSO", `name of cell type`)) |>
  select(`name of cell type`, `CATMAID annotation`, region)



#motoneuron
celltypes |>
  filter(`Sensory/inter/motor neuron` == "Motoneuron") |>
  filter(!grepl("ciliomotor", region) &
#           !grepl("girdle", region) &
#           !grepl("dorsolateral", region) &
            !grepl("mushroom", region)) |>
  select(`name of cell type`, `CATMAID annotation`, region)

