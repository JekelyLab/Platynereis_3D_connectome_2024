#code to generate Figure14 figure supplement 2 of the Platynereis 3d connectome paper
#Gaspar Jekely Feb-Dec 2022


#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

imgTEM <- readPNG("pictures/Girdle_TEM_VNC_40um.png")

panelTEM <- ggdraw() + draw_image(imgTEM, scale = 1) + 
  draw_label("sensory neuron", color = "#E69F00", size = 14,
             x = 0.25, y = 1, fontface = "bold") +
  draw_label("interneuron", color = "#56B4E9", size = 14,
             x = 0.4, y = 1, fontface = "bold") +
  draw_label("motoneuron", color = "red", size = 14,
             x = 0.55, y = 1, fontface = "bold") +
  draw_label(expression(paste("10 ", mu, "m")), x = 0.82, y = 0.08, fontfamily = "sans", fontface = "plain",
             color = "white", size = 15) + 
  draw_line(x=c(0.7, 0.95), y=0.05, color='white', size=2)

layout = "A"

Fig_VNC <- panelTEM + 
  plot_layout(design = layout) +
  theme(plot.tag = element_text(size = 12, face='plain'))

ggsave("Figures/Figure14_fig_suppl2.png", limitsize = FALSE, 
       units = c("px"), Fig_VNC, width = 4200, height = 1900, bg='white')

ggsave("Figures/Figure14_fig_suppl2.pdf", limitsize = FALSE, 
       units = c("px"), Fig_VNC, width = 4200, height = 1900)
