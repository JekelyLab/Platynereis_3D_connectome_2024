# code to generate the cell type connectivity matrix of the 3 day old Platynereis larva
# Gaspar Jekely 2022-2023

source("code/libraries_functions_and_CATMAID_conn.R")

# load the cell type graph generated in the Figure4 code
syn_tb <- readRDS("source_data/Figure4_source_data1.rds")

names <- syn_tb %>%
  activate(nodes) %>%
  select(name) %>%
  pull()
length(names)
syn_plot <- syn_tb %>%
  activate(edges) %>%
  as_tibble() %>%
  mutate(from_name = names[from]) %>%
  mutate(to_name = names[to]) %>%
  ggplot(aes(x = from_name, y = to_name)) +
  geom_tile(colour = "grey60", linewidth = 0.05, fill = NA) +
  geom_point(aes(size = sqrt(synapses), color = synapses),
    shape = 15, stroke = 0
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.key.size = unit(30, "pt"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 24)
  ) +
  coord_flip() +
  # order data as in the cell type_names list
  scale_x_discrete(limits = as.character(rev(names[1:202])), name = "presynaptic") +
  scale_y_discrete(limits = as.character(names), name = "postsynaptic") +
  scale_colour_gradient2(
    low = "white",
    mid = "#0072B2",
    high = "#E69F00",
    midpoint = 0.5,
    space = "Lab",
    na.value = "white",
    guide = "colourbar",
    aesthetics = "colour"
  )

ggsave("Figures/Figure4_fig_suppl1.pdf",
  limitsize = FALSE,
  units = c("px"), syn_plot, width = 14000, height = 10000
)

ggsave("Figures/Figure4_fig_suppl1.png",
  limitsize = FALSE,
  units = c("px"), syn_plot, width = 14000, height = 10000, bg = "white"
)
