#code to generate Figure12 figure supplement 3 of the Platynereis 3d connectome paper
#Gaspar Jekely 2023

#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

#redefine read and plot a neuron function
read_plot_neuron_ventral <- function(annotation, color){
  neuron1 = nlapply(read.neurons.catmaid(annotation, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  plot_background_ventral_no_ac()
  plot3d(neuron1, soma=TRUE, lwd=4,
         rev = FALSE, fixup = F, add=T, forceClipregion = F, alpha=1,
         col=color)
  nview3d("ventral", extramat=rotationMatrix(0.4, 1, 0.048, 0.5))
  par3d(zoom=0.52)
}

# read cells by segment --------------
{
sg0 = nlapply(read.neurons.catmaid("^segment_0$", pid=11),
              function(x) smooth_neuron(x, sigma=6000))
sg1 = nlapply(read.neurons.catmaid("^segment_1$", pid=11),
              function(x) smooth_neuron(x, sigma=6000))
sg2 = nlapply(read.neurons.catmaid("^segment_2$", pid=11),
              function(x) smooth_neuron(x, sigma=6000))
sg3 = nlapply(read.neurons.catmaid("^segment_3$", pid=11),
              function(x) smooth_neuron(x, sigma=6000))
pyg = nlapply(read.neurons.catmaid("^pygidium$", pid=11),
              function(x) smooth_neuron(x, sigma=6000))
}

# list segment-specific cell types

{
seg <- c(
  "segment_0", "segment_1", "segment_2", 
  "segment_3", "pygidium"
  )

for (cell in 1:202){
  celltype <- paste("^celltype", cell, "$", sep = "")
  for (i in 1:5){
  present_in_seg[i] <- length(skids_by_2annotations(celltype, seg[i]))
  }
  if (sum(as.numeric(present_in_seg != 0), na.rm = TRUE) == 1)
    print (paste("celltype", cell, sep = ""))
}  


#outcome of segment-specific tally
# c("celltype36", "celltype38", "celltype47", "celltype54", "celltype59", "celltype60", "celltype67", "celltype70", "celltype72", "celltype74", "celltype77", "celltype78", "celltype82", "celltype87", "celltype90", "celltype93", "celltype98", "celltype102", "celltype108", "celltype115", "celltype142", "celltype145", "celltype147", "celltype149", "celltype150", "celltype151", "celltype157", "celltype158", "celltype167", "celltype169", "celltype170", "celltype176", "celltype177", "celltype180", "celltype182", "celltype200", "celltype202")
}

# read and plot neurons ----------

# sg0 specific ----------------

read_plot_neuron_ventral("^celltype167$", blues[c(7,9,5,6)])
plot3d(sg0, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
plot3d(scalebar_50um_ventral, lwd=3, add=T, color = 'black')
rgl.snapshot("pictures/seg_spec_SNstiff.png")
close3d()

# sg1 specific ----------------

read_plot_neuron_ventral("^celltype102$", blues[7])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_doCRunp.png")
close3d()

read_plot_neuron_ventral("^celltype38$", blues[c(4,7)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_MS5.png")
close3d()

read_plot_neuron_ventral("^celltype59$", blues[c(9,6)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_loop.png")
close3d()

read_plot_neuron_ventral("^celltype142$", bluepurple[c(6,9)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_Ser_tr1.png")
close3d()

read_plot_neuron_ventral("^celltype87$", bluepurple[c(5,8)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_MC3cover.png")
close3d()

read_plot_neuron_ventral("^celltype182$", bluepurple[c(2,3,4,5,6,7,8,9)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_MNsmile.png")
close3d()

read_plot_neuron_ventral("^celltype151$", bluepurple[c(5,9)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_MNladder.png")
close3d()

read_plot_neuron_ventral("^celltype60$", blues[c(5,9)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INCM.png")
close3d()

read_plot_neuron_ventral("^celltype77$", blues[c(5,8)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INsplitPB_RF_Ya.png")
close3d()

read_plot_neuron_ventral("^celltype78$", blues[c(9,6)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INsplitPBant.png")
close3d()

read_plot_neuron_ventral("^celltype149$", blues[c(6,9)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INsplitBronto.png")
close3d()

read_plot_neuron_ventral("^celltype150$", blues[c(5,9)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INleucoPU.png")
close3d()

read_plot_neuron_ventral("^celltype158$", blues[c(9,6)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INcomm_DownL.png")
close3d()

read_plot_neuron_ventral("^celltype200$", blues[c(9,6)])
plot3d(sg1, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INiaStiff.png")
close3d()


# sg2 specific ------------------

read_plot_neuron_ventral("^celltype47$", bluepurple[c(5,9)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_MNspinning.png")
close3d()

read_plot_neuron_ventral('celltype67', bluepurple[c(9,6)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_MNbow.png")
close3d()

read_plot_neuron_ventral("^celltype202$", blues[c(6,9)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INsplitSpin.png")
close3d()

read_plot_neuron_ventral("^celltype108$", blues[c(5,9,7)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INbiaxLeg.png")
close3d()

read_plot_neuron_ventral("^celltype157$", blues[c(5,9)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INsplitVent.png")
close3d()

read_plot_neuron_ventral("^celltype169$", blues[c(5,9,7)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INcommascFV.png")
close3d()

read_plot_neuron_ventral("^celltype180$", blues[c(5,9)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INbiaxHmid.png")
close3d()

read_plot_neuron_ventral("^celltype176$", blues[c(7,9)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INcosplit.png")
close3d()

read_plot_neuron_ventral("^celltype177$", blues[c(8,9)])
plot3d(sg2, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INcommdescFVa.png")
close3d()



# segment 3 specific ------------------

read_plot_neuron_ventral("^celltype70$", blues[c(8,9)])
plot3d(sg3, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INchaeMech.png")
close3d()

read_plot_neuron_ventral("^celltype74$", blues[c(8,9)])
plot3d(sg3, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INsplitCRATO.png")
close3d()

# pygidium specific ---------------------

read_plot_neuron_ventral("^celltype82$", bluepurple[c(5,9,7)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_cioMNcover.png")
close3d()

read_plot_neuron_ventral("^celltype170$", blues[c(5,9)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_SNpygM.png")
close3d()

read_plot_neuron_ventral("^celltype115$", blues[c(6,9)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_SNPDF-pyg.png")
close3d()

read_plot_neuron_ventral("^celltype147$", blues[c(5,9,6,7)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INasc-pyg.png")
close3d()

read_plot_neuron_ventral("^celltype145$", blues[c(9,7)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INATOpyg.png")
close3d()

read_plot_neuron_ventral("^celltype72$", blues[c(9,7)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_INFVa_pyg.png")
close3d()

read_plot_neuron_ventral("^celltype54$", blues[9])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_pygPBunp.png")
close3d()

read_plot_neuron_ventral("^celltype93$", blues[c(8,9)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_pygCirrusPUSM.png")
close3d()

read_plot_neuron_ventral("^celltype98$", blues[c(2,3,4,5,6,7,8,9,2,3,4,5,6,7,8,9,6)])
plot3d(pyg, soma=T, lwd=1,
       add=T, alpha=0.08, col='grey90') 
rgl.snapshot("pictures/seg_spec_pygCirrusPU.png")
close3d()


# assemble figure ---------------------------------------------------------

panel_SNstiff <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_SNstiff.png")
) + 
  draw_label("SNstiff", x = 0.5, y = 0.95, size = 11) +
  draw_label(expression(paste("50 ", mu, " m")), 
             x = 0.84, y = 0.07, size = 9, color = "black")  +
    geom_segment(aes(x = 0.1,
                     y = 0.9,
                     xend = 0.1,
                     yend = 0.82),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) +
    geom_segment(aes(x = 0.1,
                     y = 0.82,
                     xend = 0.1,
                     yend = 0.9),
                 arrow = arrow(type = 'closed', length = unit(0.8, "mm"))) + 
    draw_label("a", x = 0.1, y = 0.93, size = 8) +
    draw_label("p", x = 0.1, y = 0.79, size = 8) 


panel_doCRunp <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_doCRunp.png")
) + 
  draw_label("doCRunp", x = 0.5, y = 0.95, size = 11)

panel_MS5 <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_MS5.png")
) + 
  draw_label("MS5", x = 0.5, y = 0.95, size = 11)

panel_Loop <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_loop.png")
) + 
  draw_label("Loop", x = 0.5, y = 0.95, size = 11)

panel_Ser_tr1 <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_Ser_tr1.png")
) + 
  draw_label("Ser-tr1", x = 0.5, y = 0.95, size = 11)

panel_MC3cover <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_MC3cover.png")
) + 
  draw_label("MC3cover", x = 0.5, y = 0.95, size = 11)

panel_MNsmile <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_MNsmile.png")
) + 
  draw_label("MNsmile", x = 0.5, y = 0.95, size = 11)

panel_Mladder <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_MNladder.png")
) + 
  draw_label("MNladder", x = 0.5, y = 0.95, size = 11)

panel_INCM <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INCM.png")
) + 
  draw_label("INCM", x = 0.5, y = 0.95, size = 11)

panel_INsplitPB_RF_Ya <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INsplitPB_RF_Ya.png")
) + 
  draw_label("INsplitPB-RF/Ya", x = 0.5, y = 0.95, size = 11)

panel_INsplitPBant <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INsplitPBant.png")
) + 
  draw_label("INsplitPBant", x = 0.5, y = 0.95, size = 11)

panel_INsplitBronto <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INsplitBronto.png")
) + 
  draw_label("INsplitBronto", x = 0.5, y = 0.95, size = 11)

panel_INleucoPU <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INleucoPU.png")
) + 
  draw_label("INleucoPU", x = 0.5, y = 0.95, size = 11)

panel_INcomm_DownL <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INcomm_DownL.png")
) + 
  draw_label("INcomm-DownL", x = 0.5, y = 0.95, size = 11)

panel_INiaStiff <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INiaStiff.png")
) + 
  draw_label("INiaStiff", x = 0.5, y = 0.95, size = 11)



panel_MNspinning <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_MNspinning.png")
) + 
  draw_label("MNspinning", x = 0.5, y = 0.95, size = 11)

panel_MNbow <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_MNbow.png")
) + 
  draw_label("MNbow", x = 0.5, y = 0.95, size = 11)

panel_INsplitSpin <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INsplitSpin.png")
) + 
  draw_label("INsplitSpin", x = 0.5, y = 0.95, size = 11)

panel_INbiaxLeg <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INbiaxLeg.png")
) + 
  draw_label("INbiaxLeg", x = 0.5, y = 0.95, size = 11)

panel_INsplitVent <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INsplitVent.png")
) + 
  draw_label("INsplitVent", x = 0.5, y = 0.95, size = 11)

panel_INcommascFV <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INcommascFV.png")
) + 
  draw_label("INcommascFV", x = 0.5, y = 0.95, size = 11)

panel_INbiaxHmid <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INbiaxHmid.png")
) + 
  draw_label("INbiaxHmid", x = 0.5, y = 0.95, size = 11)

panel_INcommdescFVa <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INcommdescFVa.png")
) + 
  draw_label("INcommdescFVa", x = 0.5, y = 0.95, size = 11)

panel_INcosplit <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INcosplit.png")
) + 
  draw_label("INcosplit", x = 0.5, y = 0.95, size = 11)



panel_INchaeMech <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INchaeMech.png")
) + 
  draw_label("INchaeMech", x = 0.5, y = 0.95, size = 11)
panel_INsplitCRATO <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INsplitCRATO.png")
) + 
  draw_label("INsplitCRATO", x = 0.5, y = 0.95, size = 11)


panel_cioMNcover <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_cioMNcover.png")
) + 
  draw_label("cioMNcover", x = 0.5, y = 0.95, size = 11)

panel_SNpygM <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_SNpygM.png")
) + 
  draw_label("SNpygM", x = 0.5, y = 0.95, size = 11)

panel_SNPDF_pyg <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_SNPDF-pyg.png")
) + 
  draw_label("SNPDF-pyg", x = 0.5, y = 0.95, size = 11)

panel_INasc_pyg <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INasc-pyg.png")
) + 
  draw_label("INasc-pyg", x = 0.5, y = 0.95, size = 11)

panel_INATOpyg <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INATOpyg.png")
) + 
  draw_label("INATOpyg", x = 0.5, y = 0.95, size = 11)

panel_INFVa_pyg <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_INFVa_pyg.png")
) + 
  draw_label("INFVa-pyg", x = 0.5, y = 0.95, size = 11)

panel_pygCirrusPU <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_pygCirrusPU.png")
) + 
  draw_label("pygCirrusPU", x = 0.5, y = 0.95, size = 11)

panel_pygCirrusPUSM <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_pygCirrusPUSM.png")
) + 
  draw_label("pygCirrusPUSM", x = 0.6, y = 0.95, size = 11)

panel_pygPBunp <- ggdraw() + draw_image(
  readPNG("pictures/seg_spec_pygPBunp.png")
) + 
  draw_label("pygPBunp", x = 0.6, y = 0.96, size = 11)


#define layout

layout <- "
aAbBcCd
DeEfFgG
hHiIjJk
KlLmMnN
oOpPqQr
"
  
Figure12_fig_suppl3 <- panel_SNstiff + panel_doCRunp + panel_MS5 + panel_Loop + panel_Ser_tr1 + panel_MC3cover + panel_MNsmile + panel_Mladder + 
  panel_INCM + panel_INsplitPB_RF_Ya + panel_INiaStiff + panel_INsplitPBant + panel_INsplitBronto + panel_INleucoPU + panel_INcomm_DownL + panel_MNspinning + 
  panel_MNbow + panel_INsplitSpin + panel_INbiaxLeg + panel_INsplitVent + panel_INcommascFV + panel_INbiaxHmid + 
  panel_INcommdescFVa + panel_INcosplit +
  panel_INchaeMech + panel_INsplitCRATO + panel_cioMNcover + 
  panel_SNpygM + panel_SNPDF_pyg + panel_INasc_pyg + panel_INATOpyg + panel_INFVa_pyg +
  panel_pygCirrusPU + panel_pygCirrusPUSM + panel_pygPBunp +
  plot_layout(design = layout, heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = 'i') & 
  theme(plot.tag = element_text(size = 12, face='plain'))
  
ggsave("Figures/Figure12_fig_suppl3.png", limitsize = FALSE, 
         units = c("px"), Figure12_fig_suppl3, width = 4200, height = 4000, bg='white')

ggsave("Figures/Figure12_fig_suppl3.pdf", limitsize = FALSE, 
       units = c("px"), Figure12_fig_suppl3, width = 4200, height = 4000)




