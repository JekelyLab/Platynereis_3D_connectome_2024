#R/natverse code to generate Figure 5 of the Platynereis 3d connectome paper
#Gaspar Jekely 2023

#load natverse and other packages, some custom natverse functions and catmaid connectivity info
source("code/Natverse_functions_and_conn.R")

load_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(
    annotation, pid=11), 
    function(x) smooth_neuron(x, sigma=6000))
}

# load PDF neurons ---------------
cMNPDF <- load_neuron("^celltype10$")
SNnuch <- load_neuron("^celltype13$")
INMBPDF <- load_neuron("^celltype184$")
SNantlerPDF <- load_neuron("^celltype26$")
SNPDF_pyg <- load_neuron("^celltype115$")
SNPDF_dc <- load_neuron("^celltype27$")
SN_DLSO3_PDF <- load_neuron("^celltype40$")
SN_DLSO1.2_4 <- load_neuron("^celltype28$")
cMNdc <- load_neuron("^celltype12$")
INdescLuqinPDF <- load_neuron("^celltype175$")
hPU2l_asymPDF <- load_neuron("^celltype99$")

# load ATO neurons -------------

INATOpyg <- load_neuron("^celltype145$")
MNspider_ant <- load_neuron("^celltype61$")
INsplitCRATO <- load_neuron("^celltype74$")
INpreLadderATO <- load_neuron("^celltype154$")

# load leucokinin neurons ------------

INleucoPU <- load_neuron("^celltype150$")

# load FVamide neurons ------------

INFVa_pyg <- load_neuron("^celltype72$")
SNFVa_pyg <- load_neuron("^celltype143$")
INcommascFV <- load_neuron("^celltype169$")
INcommdescFVa <- load_neuron("^celltype177$")

# load cPRC circuit ------------

cPRC <- load_neuron("^celltype5$")
INRGW <- load_neuron("^celltype6$")
INNOS <- load_neuron("^celltype7$")

# load achatin, FMRFa and MIP neurons ---------------------

SN47Ach <- load_neuron("^celltype52$")
SNMIP <- load_neuron("^celltype24$")
SNMIP1 <- load_neuron("^celltype46$")
SNMIP4 <- load_neuron("^celltype45$")

INMBdesc2 <- load_neuron("^celltype195$")
INMBdescFMRF <- load_neuron("^celltype185$")
SN_IRP2_FMRF <- load_neuron("^celltype44$")


# load serotonin neurons ---------------------

Ser_h1 <- load_neuron("^celltype8$")
Ser_tr1 <- load_neuron("^celltype142$")
pygPBunp <- load_neuron("^celltype54$")


# load glutamatergic ---------------------

PRC <- load_neuron("^celltype1$")
hCR <- load_neuron("^celltype100$")
ventralpygCR <- load_neuron("^celltype101$")
doCRunp <- load_neuron("^celltype102$")


# load cholinergic neurons ---------------------

MC <- load_neuron("^celltype9$")
MN1 <- load_neuron("^celltype85$")
MN2 <- load_neuron("^celltype86$")
MN3 <- load_neuron("^celltype84$")
eyespotPRCR3 <- load_neuron("^celltype33$")
Loop <- load_neuron("^celltype59$")
MS1 <- load_neuron("^celltype35$")
MS2_4 <- load_neuron("^celltype36$")
MS3 <- load_neuron("^celltype37$")
MS5 <- load_neuron("^celltype38$")


# PDF plotting ----------
plot_background_ventral_no_ac()
plot3d(cMNPDF, soma = TRUE, lwd = 4, color = bluepurple[4])
plot3d(SNnuch, soma = TRUE, lwd = 2, color = bluepurple[8])
plot3d(INMBPDF, soma = TRUE, lwd = 3, color = bluepurple[6])
plot3d(SNantler, soma = TRUE, lwd = 6, color = bluepurple[9])
plot3d(SNPDF_pyg, soma = TRUE, lwd = 4, color = bluepurple[4])
plot3d(SNPDF_dc, soma = TRUE, lwd = 2, color = bluepurple[7])
plot3d(SN_DLSO3_PDF, soma = TRUE, lwd = 6, color = bluepurple[3])
plot3d(SN_DLSO1.2_4, soma = TRUE, lwd = 8, color = bluepurple[8])
plot3d(cMNdc, soma = TRUE, lwd = 8, color = bluepurple[4])
plot3d(INdescLuqinPDF, soma = TRUE, lwd = 5, color = bluepurple[5])
plot3d(hPU2l_asymPDF, soma = TRUE, lwd = 3, color = bluepurple[9])

texts3d(115000, 115000, 43000, text = "cMNPDF", cex = 2)
texts3d(124000, 44000, 44000, text = "SNnuch", cex = 2)
texts3d(20000, 74000, 37000, text = "INMBPDF", cex = 2)
texts3d(105000, 70000, 2000, text = "SNantlerPDF", cex = 2)
texts3d(55000, 138000, 155000, text = "SNPDF-pyg", cex = 2)
texts3d(35000, 43000, 23000, text = "SNPDF-dc", cex = 2)
texts3d(35000, 33000, 10000, text = "SN-DLSO3-PDF", cex = 2)
texts3d(30000, 28000, 20000, text = "SNDLSO1.2-4", cex = 2)
texts3d(74000, 4000, 43000, text = "cMNdc", cex = 2)
texts3d(38000, 48000, 30000, text = "INdescLuqinPDF", cex = 2)
texts3d(112000, 105000, 31000, text = "hPU2l-asymPDF", cex = 2)

plot3d(scalebar_50um_ventral, lwd = 3, color = "black")
par3d(zoom=0.53)

rgl.snapshot("pictures/PDF_neurons_ventral.png")
nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
rgl.pop()
#z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/PDF_neurons_frontal.png")
close3d()


# ATO plotting ----------
plot_background_ventral_no_ac()
plot3d(INATOpyg, soma = TRUE, lwd = 5, color = bluepurple[9])
plot3d(MNspider_ant, soma = TRUE, lwd = 6, color = bluepurple[3])
plot3d(INsplitCRATO, soma = TRUE, lwd = 3, color = bluepurple[6])
plot3d(INpreLadderATO, soma = TRUE, lwd = 12, color = blues[5])

texts3d(60000, 115000, 179000, text = "INATOpyg", cex = 2)
texts3d(41000, 144000, 84000, text = "MNspider-ant", cex = 2)
texts3d(27000, 134000, 117000, text = "INsplitCRATO", cex = 2)
texts3d(115000, 130000, 107000, text = "INpreLadderATO", cex = 2)

par3d(zoom=0.53)

rgl.snapshot("pictures/ATO_neurons_ventral.png")
close3d()



# leuco plotting ----------
plot_background_ventral_no_ac()
plot3d(INleucoPU, soma = TRUE, lwd = 5, color = bluepurple[9])
texts3d(55000, 130000, 51000, text = "INleucoPU", cex = 2)
par3d(zoom=0.53)

rgl.snapshot("pictures/Leuco_neurons_ventral.png")
close3d()

# FVa plotting ----------
plot_background_ventral_no_ac()
plot3d(INFVa_pyg, soma = TRUE, lwd = 5, color = bluepurple[9])
plot3d(SNFVa_pyg, soma = TRUE, lwd = 6, color = bluepurple[3])
plot3d(INcommascFV, soma = TRUE, lwd = 8, color = bluepurple[6])
plot3d(INcommdescFVa, soma = TRUE, lwd = 12, color = blues[5])

texts3d(60000, 115000, 173000, text = "INFVa-pyg", cex = 2)
texts3d(32000, 144000, 140000, text = "SNFVa-pyg", cex = 2)
texts3d(40000, 134000, 110000, text = "INcommascFV", cex = 2)
texts3d(115000, 130000, 87000, text = "INcommdescFVa", cex = 2)

par3d(zoom=0.53)

rgl.snapshot("pictures/FVa_neurons_ventral.png")
close3d()

# cPRC-INRGW plotting ----------
plot_background()

plot3d(INRGW, soma = TRUE, lwd = 1, color = bluepurple[9])
plot3d(INNOS, soma = TRUE, lwd = 4, color = bluepurple[4])
plot3d(cPRC, soma = TRUE, lwd = 3, color = oranges[4])

texts3d(100000, 36000, 4000, text = "cPRC", cex = 2)
texts3d(50000, 54000, 10000, text = "INRGW", cex = 2)
texts3d(71000, 32000, 2000, text = "INNOS", cex = 2)

rgl.snapshot("pictures/cPRC_et_al.png")
close3d()

# Achatin-FMRFa plotting ----------
plot_background()
plot3d(SN47Ach, soma = TRUE, lwd = 3, color =  oranges[5])

plot3d(INMBdesc2, soma = TRUE, lwd = 5, color =  bluepurple[4])
plot3d(INMBdescFMRF, soma = TRUE, lwd = 4, color =  bluepurple[9])
plot3d(SN_IRP2_FMRF, soma = TRUE, lwd = 3, color =  bluepurple[7])

texts3d(46000, 26000, 4000, text = "SN47Ach", cex = 2)
texts3d(36000, 72000, 4000, text = "INMBdesc2", cex = 2)
texts3d(39000, 56000, 4000, text = "INMBdescFMRF", cex = 2)
texts3d(76000, 43000, 4000, text = "SN-IRP2-FMRF", cex = 2)

rgl.snapshot("pictures/SN47Ach.png")
close3d()

# MIP plotting ---------------
plot_background()
plot3d(SNMIP, soma = TRUE, lwd = c(4,3), color =  oranges[5:6])
plot3d(SNMIP1, soma = TRUE, lwd = c(2,3), color =  oranges[7:8])
plot3d(SNMIP4, soma = TRUE, lwd = c(4,5), color =  oranges[4:5])

texts3d(63000, 103000, 4000, text = "SNMIP-vc", cex = 2)
texts3d(63000, 64000, 4000, text = "SNMIP1", cex = 2)
texts3d(60000, 29000, 4000, text = "SNMIP4", cex = 2)

rgl.snapshot("pictures/MIP.png")
close3d()

# FVa plotting ----------
plot_background_ventral_no_ac()
plot3d(Ser_h1, soma = TRUE, lwd = 3, color = "red")
plot3d(Ser_tr1, soma = TRUE, lwd = 2, color = bluepurple[6])
plot3d(pygPBunp, soma = TRUE, lwd = 5, color = bluepurple[8])

texts3d(40000, 52000, 24000, text = "Ser-h1", cex = 2)
texts3d(58000, 144000, 60000, text = "Ser-tr1", cex = 2)
texts3d(77000, 115000, 177000, text = "pygPBunp", cex = 2)
par3d(zoom=0.53)

rgl.snapshot("pictures/Ser_neurons.png")
close3d()

# plot cholinergic neurons ---------------------

plot_background_ventral_no_ac()

plot_ACh <- function(){
plot3d(MC, soma = TRUE, lwd = 7, color = "black")
plot3d(MN1, soma = TRUE, lwd = 6, color = bluepurple[4])
plot3d(MN2, soma = TRUE, lwd = 5, color = bluepurple[5])
plot3d(MN3, soma = TRUE, lwd = 4, color = bluepurple[6])
plot3d(Loop, soma = TRUE, lwd = 3, color = bluepurple[7])
plot3d(MS1, soma = TRUE, lwd = 2, color = bluepurple[8])
plot3d(MS2_4, soma = TRUE, lwd = 2, color = bluepurple[8])
plot3d(MS3, soma = TRUE, lwd = 2, color = bluepurple[8])
plot3d(MS5, soma = TRUE, lwd = 2, color = bluepurple[8])
plot3d(eyespotPRCR3, soma = TRUE, lwd = 5, color = "red")

texts3d(70000, 42000, 20000, text = "MC", cex = 2)
texts3d(62000, 11000, 72000, text = "MS4", cex = 2)
texts3d(50000, 22000, 76000, text = "MS5", cex = 2)
texts3d(77000, 58000, 7000, text = "MS1", cex = 2)
texts3d(70000, 72000, 10000, text = "MS2", cex = 2)
texts3d(48000, 32000, 14000, text = "MS3", cex = 2)
texts3d(58000, 138000, 60000, text = "Loop", cex = 2)
texts3d(20000, 72000, 20000, text = "eyespotPRCR3", cex = 2)
texts3d(68000, 104000, 20000, text = "vMN", cex = 2)
}

plot_ACh()
par3d(zoom=0.53)

rgl.snapshot("pictures/ACh_neurons_ventral.png")
nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
#z-axis clip
clipplanes3d(0, 0, -1, 105000)
par3d(windowRect = c(0, 0, 800, 800))
rgl.snapshot("pictures/ACh_neurons_frontal.png")
close3d()

# plot glutamatergic -----------

plot_background_ventral_no_ac()

plot3d(PRC, soma = TRUE, lwd = 3, color = bluepurple[4])
plot3d(hCR, soma = TRUE, lwd = 4, color = bluepurple[9])
plot3d(ventralpygCR, soma = TRUE, lwd = 2, color = bluepurple[7])
plot3d(doCRunp, soma = TRUE, lwd = 4, color = "red")

texts3d(22000, 72000, 20000, text = "PRC", cex = 2)
texts3d(30000, 72000, 1000, text = "hCR", cex = 2)
texts3d(70000, 32000, 77000, text = "doCRunp", cex = 2)
texts3d(40000, 138000, 65000, text = "trunkCR", cex = 2)
par3d(zoom=0.53)
rgl.snapshot("pictures/Glu_neurons.png")
close3d()

# network plot ----------------------------------

# load network
syn_tb <- readRDS("source_data/Figure4_source_data1.rds")
syn.igraph <- as.igraph(syn_tb)
conn_gexf <- rgexf::read.gexf("data/celltype_connectome_force_layout.gexf")

# get coordinates from imported gephi file
coords <- as_tibble(x = conn_gexf$nodesVizAtt$position$x) %>%
  mutate(x = value) %>%
  mutate(y = conn_gexf$nodesVizAtt$position$y) %>%
  mutate(skid = conn_gexf$nodes$label)

# sort by skid
x_coord <- as.list(arrange(coords, desc(skid)) %>%
                     select(x))
y_coord <- as.list(arrange(coords, desc(skid)) %>%
                     select(y))

# assign coordinates from imported gephi graph to original CATMAID graph
syn_tb <- activate(syn_tb, nodes) %>%
  arrange(desc(name)) %>% # sort by skid to match coordinate list
  mutate(
    x = unlist(x_coord), # add coordinates
    y = unlist(y_coord)
  )

# graph visualisation -----------------------------------------------------

# convert to visNet form
syn.vis <- syn_tb %>%
  toVisNetworkData()

# assign sqrt of number of synapses to edge 'value'
syn.vis$edges$value <- sqrt(syn.vis$edges$synapses)

#coordinates as matrix
coords <- matrix(c(syn.vis$nodes$x, syn.vis$nodes$y), ncol = 2)

PDF_neurons <- c(
  "cMNPDF-vcl1", "SNnuchPDF", "INMBPDF", "SNantlerPDF", 
  "SNPDF-pyg", "SNPDF-dc", "SN-DLSO3", "SN0DLSO1.2-4", 
  "cMNdc", "INdescLuqinPDF", "hPU2l-asymPDF",
  "SN-DLSO3-PDF")
ATO_neurons <- c("cMNATO", "INATOpyg", "MNspider-ant", "INsplitCRATO", "INpreLadderATO")
Leuco_neurons <- c("INleucoPU")
FVa_neurons <- c("INFVa-pyg", "SNFVa-pyg", "INcommascFV", "INcommdescFVa", "SNFVa")
RGW_neurons <- c("INRGWa")
sNPF_neurons <- c("INNOS")
pedal_neurons <- c("cPRC")
ach_neurons <- c("SN47Ach")
MIP_neurons <- c("SNMIP-vc", "SN-MIP1", "SN-MIP4")
FMRF_neurons <- c("INMBdescFMRF", 
                  "INMBdesc2", "SN-IRP2-FMRF")
Ser_neurons <- c("Ser-h1", "Ser-tr1", "pygPBunp")
Glu_neurons <- c("PRC", "hCR", "ventralpygCR", "doCRunp")
ACh_neurons <- c("MC", "MN1", "MN2", "MN3", "eyespot-PRCR3", "Loop", "MS1", "MS2", "MS3", "MS5", "MS4")
other_peptidergic_neurons <- c(
   "SN-ASTC", "SN-WLD",
  "SN-IRP2-burs", "MNbow", "SN-YF5cil")


transmitter <- syn_tb %>%
  mutate(
    transmitter = ifelse(name %in% PDF_neurons, "PDF+", 
                  ifelse(name %in% ATO_neurons, "ATO+", 
                  ifelse(name %in% Leuco_neurons, "leuco+", 
                  ifelse(name %in% FVa_neurons, "FVa+", 
                  ifelse(name %in% RGW_neurons, "RGWa+", 
                  ifelse(name %in% sNPF_neurons, "sNPF+", 
                  ifelse(name %in% pedal_neurons, "pedal+",
                  ifelse(name %in% FMRF_neurons, "FMRF+",
                  ifelse(name %in% ach_neurons, "achatin+", 
                  ifelse(name %in% MIP_neurons, "MIP+",
                  ifelse(name %in% other_peptidergic_neurons, "peptidergic+",
                  ifelse(name %in% Ser_neurons, "serotonergic", 
                  ifelse(name %in% Glu_neurons, "glutamatergic", 
                  ifelse(name %in% ACh_neurons, "cholinergic", "none")
    )))))))))))))) %>%
  select(transmitter) %>% 
  pull()
         
syn.vis$nodes$group <- transmitter
length(transmitter) - length(transmitter[transmitter=="none"])
53/202

visNet <- visNetwork(syn.vis$nodes, syn.vis$edges) %>%
  visIgraphLayout(type = "full", layoutMatrix = coords) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0.1),
    scaling = list(min = 0.1, max = 20),
    color = list(inherit = TRUE, opacity = 0.4),
    arrows = list(to = list(
      enabled = TRUE,
      scaleFactor = 0.8, type = "arrow"
    ))
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = syn.vis$nodes$color, border = "black"),
    scaling = list(min = 5, max = 50),
    font = list(color = "black", size = 40),
  ) %>%
  visGroups(
    groupname = "PDF+", color = bluepurple[9], shape = "square",
    opacity = 1
  ) %>%
  visGroups(
    groupname = "ATO+", shape = "square",
    opacity = 1, color = bluepurple[8]
  ) %>%
  visGroups(
    groupname = "leuco+", shape = "square",
    opacity = 1, color = bluepurple[7]
  ) %>%
  visGroups(
    groupname = "FVa+", shape = "square",
    opacity = 1, color = bluepurple[6]
  ) %>%
  visGroups(
    groupname = "RGWa+", shape = "square",
    opacity = 1, color = bluepurple[5]
  ) %>%
  visGroups(
    groupname = "sNPF+", shape = "square",
    opacity = 1, color = bluepurple[6]
  ) %>%
  visGroups(
    groupname = "pedal+", shape = "square",
    opacity = 1, color = bluepurple[7]
  ) %>%
  visGroups(
    groupname = "achatin+", shape = "square",
    opacity = 1, color = bluepurple[8]
  ) %>%
  visGroups(
    groupname = "MIP+", shape = "square",
    opacity = 1, color = bluepurple[7]
  ) %>%
  visGroups(
    groupname = "FMRF+", shape = "square",
    opacity = 1, color = bluepurple[4]
  ) %>%
  visGroups(
    groupname = "peptidergic+", shape = "square",
    opacity = 1, color = bluepurple[6]
  ) %>%
  visGroups(
    groupname = "serotonergic", shape = "triangle",
    opacity = 1, color = Okabe_Ito[1]
  ) %>%
  visGroups(
    groupname = "glutamatergic", shape = "triangle",
    opacity = 1, color = Okabe_Ito[2]
  ) %>%
  visGroups(
    groupname = "cholinergic", shape = "triangle",
    opacity = 1, color = Okabe_Ito[3]
  ) %>%
  visGroups(
    groupname = "none", shape = "triangle",
    opacity = 1, color = "#EEEEEE",
    font = list(color = "black", size = 20)
  ) %>%
  visOptions(
    highlightNearest = list(enabled = T, hover = T),
    width = 3000, height = 1500) %>%
  addFontAwesome()

visNet

# save as html
saveNetwork(visNet, "pictures/network_with_transmitters.html",
            selfcontained = TRUE
)
# save from web browser
webshot::webshot(
  url = "pictures/network_with_transmitters.html",
  file = "pictures/network_with_transmitters.png",
  vwidth = 3000, vheight = 1500, # define the size of the browser window
  cliprect = c(120, 130, 2800, 1360), zoom = 2, delay = 2
)

# assemble figure ----------------

panel_PDFv <-ggdraw() + draw_image(readPNG("pictures/PDF_neurons_ventral.png")) + 
  draw_label(
    "PDF+, ventral view", x = 0.5, y = 0.99, 
    color = "black", size = 11
  ) +
  draw_label(
    expression(paste("50 ", mu, " m")), 
    x = 0.83, y = 0.07, 
    color = "black", size = 10
    )


panel_PDFa <-ggdraw() + draw_image(readPNG("pictures/PDF_neurons_frontal.png")) + 
  draw_label(
    "PDF+, anterior view", x = 0.4, y = 0.99, 
    color = "black", size = 11
  )

panel_ATO <-ggdraw() + draw_image(readPNG("pictures/ATO_neurons_ventral.png")) + 
  draw_label(
    "allatotropin/orexin+", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_Leu <-ggdraw() + draw_image(readPNG("pictures/Leuco_neurons_ventral.png")) + 
  draw_label(
    "leucokinin+", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_FVa <-ggdraw() + draw_image(readPNG("pictures/FVa_neurons_ventral.png")) + 
  draw_label(
    "FVamide+", x = 0.3, y = 0.99, 
    color = "black", size = 11
  )

panel_cPRC <-ggdraw() + draw_image(readPNG("pictures/cPRC_et_al.png")) + 
  draw_label(
    "sNPF/RYamide+ (INNOS)\nRGWamide+ (INRGW)\npedal-peptide-2+ (cPRC)", 
    x = 0.1, y = 0.94, 
    color = "black", size = 11, hjust = 0
  )

panel_ach <-ggdraw() + draw_image(readPNG("pictures/SN47Ach.png")) + 
  draw_label(
    "achatin+, FMRFa+", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_MIP <-ggdraw() + draw_image(readPNG("pictures/MIP.png")) + 
  draw_label(
    "MIP+", x = 0.2, y = 0.99, 
    color = "black", size = 11
    )

panel_Ser <-ggdraw() + draw_image(readPNG("pictures/Ser_neurons.png")) + 
  draw_label(
    "serotonergic", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_ACh_v <-ggdraw() + draw_image(readPNG("pictures/ACh_neurons_ventral.png")) + 
  draw_label(
    "cholinergic, ventral view", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )
panel_ACh_a <-ggdraw() + draw_image(readPNG("pictures/ACh_neurons_frontal.png")) + 
  draw_label(
    "cholinergic, anterior view", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_Glu <-ggdraw() + draw_image(readPNG("pictures/Glu_neurons.png")) + 
  draw_label(
    "glutamatergic", x = 0.5, y = 0.99, 
    color = "black", size = 11
  )

panel_network <-ggdraw() + draw_image(readPNG("pictures/network_with_transmitters.png")) + 
  draw_label(
    "peptidergic", x = 0.1, y = 0.9, 
    color = bluepurple[7], size = 11, hjust = 0
  ) +
  draw_label(
    "serotonergic", x = 0.1, y = 0.85, 
    color = Okabe_Ito[1], size = 11, hjust = 0
  ) +
  draw_label(
    "glutamatergic", x = 0.1, y = 0.8, 
    color = Okabe_Ito[2], size = 11, hjust = 0
  ) +
  draw_label(
    "cholinergic", x = 0.1, y = 0.75, 
    color = Okabe_Ito[3], size = 11, hjust = 0
  )


layout <- "
AAAA
BCDE
FGHI
JKLM
"

Figure5 <- panel_network + 
  panel_PDFv + panel_PDFa + panel_Leu + panel_ATO + 
  panel_FVa + panel_cPRC + panel_ach + panel_MIP +
  panel_Ser + panel_Glu + panel_ACh_v + panel_ACh_a + 
  plot_layout(design = layout, heights = c(2, 1.3, 1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=12, face='plain'))

ggsave("Figures/Figure5.png", limitsize = FALSE, 
       units = c("px"), Figure5, 
       width = 3200, height = 4500, bg='white'
       )  

ggsave("Figures/Figure5.pdf", limitsize = FALSE, 
       units = c("px"), Figure5, 
       width = 3200, height = 4500
       )  




