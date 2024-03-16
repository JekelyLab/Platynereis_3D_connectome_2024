rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.

library(tidyverse)
library(tidygraph)
library(stats)

segmental_colors <- brewer.pal(6, 'Paired')


#read graph
conn.tb <- readRDS("supplements/connectome_graph_tibble.rds")
conn.tb

conn.tb %>%
  as_tibble() %>%
  select(celltype_annotation )

conn.tb <- conn.tb %>%
  as_tibble() %>%
  group_by(celltype_annotation) %>%
  mutate(count = n()) %>%
  filter(!segment == 'fragment') %>%
  filter(!celltype_annotation == 'not_celltype')
conn.tb
conn.tb %>%
  ggplot(aes(partition, skids, color = segment, size = eigen)) +
  geom_point() +
  scale_color_manual(values = segmental_colors)

# create a distance matrix between categories based on their count
conn.tb %>%
  group_by(celltype_annotation) %>%
  dist(count)

# group the data by category
grouped_data <- conn.tb %>%
  group_by(celltype_annotation)

# count the number of rows for each category
data_with_counts <- grouped_data %>%
  summarise(count = n()) %>%
  ungroup()

# create a distance matrix between categories based on their count
distance_matrix <- dist(data_with_counts$count)

  
distance_matrix
# perform hierarchical clustering on the distance matrix
hclust_object <- hclust(distance_matrix)

# reorder the categories based on the hierarchical clustering
ordered_categories <- conn.tb$celltype_annotation[hclust_object$order]
length(ordered_categories)
ordered_categories
conn.tb %>%
  group_by(celltype_annotation) %>%
  summarize(count = n())

# order the data based on the ordered categories
data_ordered <- conn.tb[match(ordered_categories, conn.tb$celltype_annotation), ]
data_ordered <- conn.tb[order(match(conn.tb$celltype_annotation, ordered_categories)), ]

data_ordered
conn.tb

conn.tb %>%
  group_by(celltype_annotation) %>%
  ggplot(aes(partition, celltype_annotation, color = segment, size = count)) +
  geom_point() +
  scale_color_manual(values = segmental_colors) +
  theme(axis.text = element_text(size = 2))


