# Generate the figure of germ layers of the Platynereis 3d connectome paper
# Gaspar Jekely 2023

# load nat and all associated packages, incl catmaid
source("code/Natverse_functions_and_conn.R")

# read cells --------------
{
  ectoderm <- nlapply(
    read.neurons.catmaid("^ectoderm$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  mesoderm <- nlapply(
    read.neurons.catmaid("^mesoderm$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  endoderm <- nlapply(
    read.neurons.catmaid("^endoderm$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  # these four dots are the most extreme points of the volume, adding them to the 3d view solves the problem with automatic zooming and movement of the field shown
  bounding_dots <- nlapply(
    read.neurons.catmaid("^bounding_dots$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  skids_ectoderm_left <- skids_by_2annotations("ectoderm", "left_side")
  ectoderm_left <- nlapply(
    read.neurons.catmaid(skids_ectoderm_left,
      pid = 11, conn = conn_http1,
      fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  skids_mesoderm_left <- skids_by_2annotations("mesoderm", "left_side")
  mesoderm_left <- nlapply(
    read.neurons.catmaid(skids_mesoderm_left,
      pid = 11, conn = conn_http1,
      fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  skids_endoderm_left <- skids_by_2annotations("endoderm", "left_side")
  endoderm_left <- nlapply(
    read.neurons.catmaid(skids_endoderm_left,
      pid = 11, conn = conn_http1,
      fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  dividing <- nlapply(
    read.neurons.catmaid("^dividing$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_dividing_left <- skids_by_2annotations("dividing", "left_side")
  dividing_left <- nlapply(
    read.neurons.catmaid(skids_dividing_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  macrophage <- nlapply(
    read.neurons.catmaid("^celltype_non_neuronal13$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_macrophage_left <- skids_by_2annotations("celltype_non_neuronal13", "left_side")
  macrophage_left <- nlapply(
    read.neurons.catmaid(skids_macrophage_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  
  epidermis <- nlapply(
    read.neurons.catmaid("^epithelia_cell$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_epidermis_left <- skids_by_2annotations("epithelia_cell", "left_side")
  epidermis_left <- nlapply(
    read.neurons.catmaid(skids_epidermis_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  stomodeum <- nlapply(
    read.neurons.catmaid("^stomodeum$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_stomodeum_left <- skids_by_2annotations("stomodeum", "left_side")
  stomodeum_left <- nlapply(
    read.neurons.catmaid(skids_stomodeum_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  stomodeum <- nlapply(
    read.neurons.catmaid("^stomodeum$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_stomodeum_left <- skids_by_2annotations("stomodeum", "left_side")
  stomodeum_left <- nlapply(
    read.neurons.catmaid(skids_stomodeum_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  hindgut <- nlapply(
    read.neurons.catmaid("^gut$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_hindgut_left <- skids_by_2annotations("gut", "left_side")
  hindgut_left <- nlapply(
    read.neurons.catmaid(skids_hindgut_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  
  pigment <- nlapply(
    read.neurons.catmaid("^pigment cell$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_pigment_left <- skids_by_2annotations("^pigment cell$", "left_side")
  pigment_left <- nlapply(
    read.neurons.catmaid(skids_pigment_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  gland <- nlapply(
    read.neurons.catmaid("^gland cell$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_gland_left <- skids_by_2annotations("^gland cell$", "left_side")
  gland_left <- nlapply(
    read.neurons.catmaid(skids_gland_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  glia <- nlapply(
    read.neurons.catmaid("^glia cell$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_glia_left <- skids_by_2annotations("glia", "left_side")
  glia_left <- nlapply(
    read.neurons.catmaid(skids_glia_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  ciliary_band <- nlapply(
    read.neurons.catmaid("^ciliated cell$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_ciliary_band_left <- skids_by_2annotations("^ciliated cell$", "left_side")
  ciliary_band_left <- nlapply(
    read.neurons.catmaid(skids_ciliary_band_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  coelothelium <- nlapply(
    read.neurons.catmaid("^celltype_non_neuronal36$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_coelothelium_left <- skids_by_2annotations("^celltype_non_neuronal36$", "left_side")
  coelothelium_left <- nlapply(
    read.neurons.catmaid(skids_coelothelium_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  
  muscle <- nlapply(
    read.neurons.catmaid("^muscle$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_mus_left <- skids_by_2annotations("^muscle$", "left_side")
  muscle_left <- nlapply(
    read.neurons.catmaid(skids_mus_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  follicle <- nlapply(
    read.neurons.catmaid("^follicle cell$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_follicle_left <- skids_by_2annotations("^follicle cell$", "left_side")
  follicle_left <- nlapply(
    read.neurons.catmaid(skids_follicle_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  
  neuroectoderm <- nlapply(
    read.neurons.catmaid("^pnb$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_neuroectoderm_left <- skids_by_2annotations("^pnb$", "left_side")
  neuroectoderm_left <- nlapply(
    read.neurons.catmaid(skids_neuroectoderm_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  
  yolk_cell <- nlapply(
    read.neurons.catmaid("^yolk$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_yolk_cell_left <- skids_by_2annotations("^yolk$", "left_side")
  yolk_cell_left <- nlapply(
    read.neurons.catmaid(skids_yolk_cell_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  yolk_blanket <- nlapply(
    read.neurons.catmaid("^yolk blanket$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_yolk_blanket_left <- skids_by_2annotations("^yolk blanket$", "left_side")
  yolk_blanket_left <- nlapply(
    read.neurons.catmaid(skids_yolk_blanket_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  nephridia <- nlapply(
    read.neurons.catmaid("^nephridia$", pid = 11),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  skids_nephridia_left <- skids_by_2annotations("^nephridia$", "left_side")
  nephridia_left <- nlapply(
    read.neurons.catmaid(skids_nephridia_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )

  
  skids_acicula_left <- skids_by_2annotations("^acicula", "left_side")
  acicula_left <- nlapply(
    read.neurons.catmaid(skids_acicula_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
  skids_chaeta_left <- skids_by_2annotations("^chaeta$", "left_side")
  chaeta_left <- nlapply(
    read.neurons.catmaid(skids_chaeta_left,
                         pid = 11, conn = conn_http1,
                         fetch.annotations = FALSE
    ),
    function(x) smooth_neuron(x, sigma = 6000)
  )
  
}

# define two windows, plot background -------------
{
  nopen3d() # opens a pannable 3d window
  mfrow3d(1, 2) # defines the two scenes
  par3d(windowRect = c(20, 30, 1200, 800)) # to define the size of the rgl window
  nview3d("ventral", extramat = rotationMatrix(0, 1, 0, 0))
  plot3d(bounding_dots,
    lwd = 1, add = T, alpha = 1,
    col = "white"
  )
  plot3d(scalebar_50um_ventral, add = T, alpha = 1,
         lwd = 2, col = "black"
  )
  nview3d("ventral", extramat=rotationMatrix(-0.01, 0, 0.1, 0))
  
  par3d(zoom = 0.48)
  next3d(clear = F)
  nview3d("right", extramat = rotationMatrix(-pi / 2, pi, -0.2, 0))
  plot3d(bounding_dots, lwd = 1, add = T, alpha = 1, col = "white")
  par3d(zoom = 0.48)
}

# plot ventral views and side views of left side cells only ------
{
  next3d(clear = F)
  
  plot3d(endoderm,
    soma = T, lwd = 1, add = T, alpha = 1, col = "grey20"
  )
  
  plot3d(mesoderm,
         soma = T, lwd = 1, add = T, alpha = 0.4, col = Okabe_Ito[2]
  )
  plot3d(ectoderm,
         soma = T, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[6]
  )
  
  plot3d(yolk, add = T, alpha = 0.1,
         col = "#E2E2E2"
  )
  
  par3d(zoom = 0.49)
  
  next3d(clear = F)

  plot3d(endoderm_left,
         soma = T, lwd = 1, add = T, alpha = 1, col = "grey20"
  )
  plot3d(mesoderm_left,
         soma = T, lwd = 1, add = T, alpha = 0.4, col = Okabe_Ito[2]
  )
  
  plot3d(ectoderm_left,
         soma = T, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[6]
  ) 
  plot3d(yolk,
         WithConnectors = F, WithNodes = F, soma = F, lwd = 2,
         add = T, alpha = 0.1,
         col = "#E2E2E2"
  )
  par3d(zoom = 0.49)

}

rgl.snapshot("pictures/germ_layers.png")
close3d()

# background function

background_anatomy <- function(){
  nopen3d() # opens a pannable 3d window
  mfrow3d(1, 2) # defines the two scenes
  par3d(windowRect = c(20, 30, 1200, 800)) # to define the size of the rgl window
  next3d(clear = F)
  
  nview3d("ventral", extramat=rotationMatrix(-0.01, 0, 0.1, 0))
  par3d(zoom = 0.49)
  plot3d(bounding_dots,
         lwd = 1, add = T, alpha = 1,
         col = "white"
  )
  plot3d(endoderm,
       soma = T, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[2]
  )

  plot3d(mesoderm,
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[2]
  )

  plot3d(ectoderm,
       soma = T, lwd = 1, add = T, alpha = 0.03, col = Okabe_Ito[2]
  )

#  plot3d(yolk, add = T, alpha = 0.1,
#       col = "#E2E2E2"
#  )

  next3d(clear = F)
  nview3d("right", extramat = rotationMatrix(-pi / 2, pi, -0.2, 0))
  par3d(zoom = 0.49)
  plot3d(bounding_dots,
         lwd = 1, add = T, alpha = 1,
         col = "white"
  )
  plot3d(endoderm_left,
       soma = T, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[2]
  )
  plot3d(mesoderm_left,
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[2]
  )
  plot3d(ectoderm_left,
       soma = T, lwd = 1, add = T, alpha = 0.03, col = Okabe_Ito[2]
  ) 
#  plot3d(yolk, add = T, alpha = 0.1,
#         col = "#E2E2E2"
#  )
}


# plot structures ------------

background_anatomy()
next3d(clear = F)
plot3d(glia,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
next3d(clear = F)
plot3d(glia_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_glia.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(pigment,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
next3d(clear = F)
plot3d(pigment_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_pigment.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(ciliary_band,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
next3d(clear = F)
plot3d(ciliary_band_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_cilia.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(gland,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
next3d(clear = F)
plot3d(gland_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_gland.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(nephridia,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
next3d(clear = F)
plot3d(nephridia_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_neph.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(muscle,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(oranges[3:9], 853, replace = TRUE)
)
next3d(clear = F)
plot3d(muscle_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(oranges[3:9], 425, replace = TRUE)
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_mus.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(dividing,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
next3d(clear = F)
plot3d(dividing_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_div.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(neuroectoderm,
       soma = T, lwd = 1, add = T, alpha = 1, col =  sample(bluepurple[1:9], 1600, replace = TRUE)
)
next3d(clear = F)
plot3d(neuroectoderm_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(bluepurple[1:9], length(skids_neuroectoderm_left), replace = TRUE)
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_neuroecto.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(follicle,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(bluepurple[2:9], 566, replace = TRUE)
)
next3d(clear = F)
plot3d(follicle_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(bluepurple[2:9], 286, replace = TRUE)
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_foll.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(epidermis,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(oranges[2:9], 1336, replace = TRUE)
)
next3d(clear = F)
plot3d(epidermis_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(oranges[2:9], 683, replace = TRUE)
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_epidermis.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(coelothelium,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(blues9[2:9], length(coelothelium), replace = TRUE)
)
next3d(clear = F)
plot3d(coelothelium_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(blues9[2:9], length(skids_coelothelium_left), replace = TRUE)
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_coelothelium.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(macrophage,
       soma = T, lwd = 1, add = T, alpha = 1, col = 'black'
)
next3d(clear = F)
plot3d(macrophage_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = 'black'
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_macrophage.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(stomodeum,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(bluepurple[2:9], length(stomodeum), replace = TRUE)
)
plot3d(hindgut,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(greens[2:9], length(hindgut), replace = TRUE)
)
next3d(clear = F)
plot3d(stomodeum_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(bluepurple[2:9], length(skids_stomodeum_left), replace = TRUE)
)
plot3d(hindgut_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = sample(greens[2:9], length(skids_hindgut_left), replace = TRUE)
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_gut.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(yolk_cell,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[7]
)
plot3d(yolk_blanket,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
next3d(clear = F)
plot3d(yolk_cell_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[7]
)
plot3d(yolk_blanket_left,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[8]
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_yolk.png")
close3d()

background_anatomy()
next3d(clear = F)
plot3d(acicula,
       soma = T, lwd = 4, add = T, alpha = 1, col = "grey20"
)
plot3d(chaeta,
       soma = T, lwd = 2, add = T, alpha = 1, col = "grey60"
)
next3d(clear = F)
plot3d(acicula_left,
       soma = T, lwd = 4, add = T, alpha = 1, col = "grey20"
)
plot3d(chaeta_left,
       soma = T, lwd = 2, add = T, alpha = 1, col = "grey60"
)
plot3d(yolk, add = T, alpha = 0.07,
       col = "#E2E2E2"
)
rgl.snapshot("pictures/germ_layers_chaeta.png")
close3d()



# assemble figure -----------------

p_germ_layers <- ggdraw() + draw_image(readPNG("pictures/germ_layers.png")) +
  draw_label(expression(paste("50 ", mu, " m")), x = 0.434, y = 0.09, size = 10) +
  draw_label("ventral view", x = 0.2, y = 0.99, size = 10) +
  draw_label("right view, left side", x = 0.7, y = 0.99, size = 10) +
  draw_label("ectoderm", x = 0.5, y = 0.3, size = 10, color=Okabe_Ito[6]) +
  draw_label("mesoderm", x = 0.5, y = 0.25, size = 10, color=Okabe_Ito[2]) +
  draw_label("endoderm", x = 0.5, y = 0.2, size = 10, color=Okabe_Ito[8])

p_glia <- ggdraw() + draw_image(readPNG("pictures/germ_layers_glia.png")) +
  draw_label("glia", x = 0.5, y = 0.3, size = 10)

p_pigment <- ggdraw() + draw_image(readPNG("pictures/germ_layers_pigment.png")) +
  draw_label("pigment", x = 0.5, y = 0.3, size = 10)

p_cilia <- ggdraw() + draw_image(readPNG("pictures/germ_layers_cilia.png")) +
  draw_label("ciliary bands", x = 0.5, y = 0.3, size = 10)

p_gland <- ggdraw() + draw_image(readPNG("pictures/germ_layers_gland.png")) +
  draw_label("glands", x = 0.5, y = 0.3, size = 10)

p_neph <- ggdraw() + draw_image(readPNG("pictures/germ_layers_neph.png")) +
  draw_label("nephridia", x = 0.5, y = 0.3, size = 10)

p_mus <- ggdraw() + draw_image(readPNG("pictures/germ_layers_mus.png")) +
  draw_label("muscle", x = 0.5, y = 0.3, size = 10)

p_div <- ggdraw() + draw_image(readPNG("pictures/germ_layers_div.png")) +
  draw_label("dividing cells", x = 0.5, y = 0.3, size = 10)

p_foll <- ggdraw() + draw_image(readPNG("pictures/germ_layers_foll.png")) +
  draw_label("follicle cells", x = 0.5, y = 0.3, size = 10)

p_neuroecto <- ggdraw() + draw_image(readPNG("pictures/germ_layers_neuroecto.png")) +
  draw_label("neuroectoderm", x = 0.51, y = 0.3, size = 10)

p_yolk <- ggdraw() + draw_image(readPNG("pictures/germ_layers_yolk.png")) +
  draw_label("yolk cells", x = 0.5, y = 0.3, size = 10, color=Okabe_Ito[7])+
  draw_label("yolk blanket", x = 0.5, y = 0.25, size = 10, color=Okabe_Ito[8])

p_chaeta <- ggdraw() + draw_image(readPNG("pictures/germ_layers_chaeta.png")) +
  draw_label("aciculae\nchaetae", x = 0.5, y = 0.3, size = 10)

p_epith <- ggdraw() + draw_image(readPNG("pictures/germ_layers_epidermis.png")) +
  draw_label("epithelium", x = 0.51, y = 0.3, size = 10)

p_coelothelium <- ggdraw() + draw_image(readPNG("pictures/germ_layers_coelothelium.png")) +
  draw_label("coelothelium", x = 0.5, y = 0.3, size = 10)

p_gut <- ggdraw() + draw_image(readPNG("pictures/germ_layers_gut.png")) +
  draw_label("stomodeum", x = 0.5, y = 0.3, size = 10, color = bluepurple[7]) +
  draw_label("hindgut", x = 0.5, y = 0.24, size = 10, color = greens[8])




layout <- "
ABC
DEF
GHI
JKL
MNO
"

Fig3_fig_suppl1 <- p_germ_layers + p_glia + p_pigment + 
  p_cilia + p_gland + p_mus +
  p_div + p_foll + p_neuroecto + 
  p_yolk + p_chaeta + p_neph + 
  p_epith + p_coelothelium + p_gut +
  plot_layout(design = layout, heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("Figures/Figure3_fig_suppl1.png",
       limitsize = FALSE,
       units = c("px"), Fig3_fig_suppl1, 
       width = 3600, height = 4000
)

ggsave("Figures/Figure3_fig_suppl1.pdf",
       limitsize = FALSE,
       units = c("px"), Fig3_fig_suppl1, 
       width = 3600, height = 4000
)

