library(catmaid)
library(tidyverse)

source("~/R/conn.R")

labels <- catmaid_get_label_stats(pid = 11)

gold <- labels %>% 
  filter(str_detect(labelName, "gold")) %>%
  select(labelName, skeletonID)

neuropeptides <- c("ATO", "FMRFa", "FVa", "FVRIa", "leucokinin", "luqin", "MIP", "PDF", "proenkephalin", "RGWa", "RYa")

skids_siGOLD <- unique(gold$skeletonID)

has_siGOLD <- function(skid, neuropeptide) {
  gold_score_threhold=3.4
  gold_on_skid <- gold %>%
    filter(skeletonID==skid) %>%
    select(labelName)
  neuropeptide_gold <- gold_on_skid %>%
    filter(str_detect(labelName, neuropeptide)) %>%
    pull()
  gold_score=0
  for (neuropep_gold in neuropeptide_gold) {
    ngold <- as.numeric(str_extract(neuropep_gold, "[0-9]+"))
    if (is.na(ngold)) {
      break
    }
    ngold <- as.numeric(ngold)
    # if gold label is on soma, give it fraction of score
    # because soma has big area, so much higher chance of non-specific labeling
    if (is.na(str_extract(neuropep_gold, "soma"))) {
      # gold is not on soma
      gold_score_to_add = ngold
    } else {
      # gold is on soma
      gold_score_to_add = ngold / 3
    }
    # adding 0.2 to score in each iteration to add more weight
    # to labeling in multiple layers, rather than same amount 
    # of labeling in one layer
    gold_score <- gold_score_to_add + gold_score + 0.2
  } 
  if ( gold_score > gold_score_threhold) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


add_siGOLD_annotation <- function(skid, neuropeptide) {
  if (has_siGOLD(skid, neuropeptide)) {
    print(paste("annotation to add:", skid, neuropeptide))
    siGOLD_annotation <- paste("expression:", neuropeptide, ":siGOLD:catmaid", sep = "")
    return(skid, neuropeptide)
    #catmaid_set_annotations_for_skeletons(skid, siGOLD_annotation, pid = 11)
  }
}

siGOLD_positive_neurons <- data.frame()
for (skid in skids_siGOLD) {
  for (neuropeptide in neuropeptides) {
    #add_siGOLD_annotation(skid, neuropeptide)
    if (has_siGOLD(skid, neuropeptide)) {
      df <- data.frame(skid=skid, neuropeptide=neuropeptide)
      siGOLD_positive_neurons <- rbind(siGOLD_positive_neurons, df)
    }
  }
}

remove_siGOLD_annotation <- function(skid, neuropeptide) {
  if (has_siGOLD(skid, neuropeptide)) {
    print(paste("skid", "neuropeptide"))
    siGOLD_annotation <- paste("expression:", neuropeptide, ":siGOLD:catmaid", sep = "")
    # catmaid_remove_annotations_for_skeletons function is seriously broken
    # TODO: use catmaid_fetch instead
  }
}

siGOLD_negative_neurons <- data.frame()
for (neuropeptide in neuropeptides) {
  siGOLD_annotation <- paste("expression:", neuropeptide, ":siGOLD:catmaid", sep = "")
  neurons_annotated <- catmaid_query_by_annotation(siGOLD_annotation, pid=11)
  if (is.null(skids_annotated)) {
    break
  }
  skids_annotated <- neurons_annotated %>%
    select(skid) %>%
    pull()
  for (skid in skids_annotated) {
    if (!has_siGOLD(skid, neuropeptide)) {
      df <- data.frame(skid=skid, neuropeptide=neuropeptide)
      siGOLD_negative_neurons <- rbind(siGOLD_negative_neurons, df)
    }
  }
}