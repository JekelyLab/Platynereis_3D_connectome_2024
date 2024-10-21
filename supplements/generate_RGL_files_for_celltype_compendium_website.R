

#| include: FALSE
library(dplyr)
library(rio)
library(tinytable)
library(knitr)
library(gt)
library(catmaid)
library(rgl)
library(htmltools)

# catmaid connection
#conn_http1 <- catmaid_login(
#  server="https://catmaid.jekelylab.ex.ac.uk/", 
#  authname="AnonymousUser",
#  config=httr::config(ssl_verifypeer=0, http_version=1)
#)

# catmaid connection, needs username, password AND token - weird!
{
  # can run this separate file using source function
  conn <- source("~/R/conn.R")
  #for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
  #for this we configure to http/1.1
#  conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
}


  
celltypes <- rio::import("supplements/Supplementary_Table1.csv")

# Neuronal cell types

## Sensory neurons




plot_background <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(0, 0, 200, 200)) #resize for frontal view
  plot3d(bounding_dots, add=T, col="white") 
  plot3d(yolk, add=T, alpha=0.1, col="#E2E2E2") 
  plot3d(acicula, soma=F, lwd=2, add=T, alpha=1, col="grey70")
  par3d(zoom=0.52)
  nview3d("frontal", extramat=rotationMatrix(0.2, 1, 0.1, 0.5))
  #z-axis clip
  #clipplanes3d(0, 0, -1, 65000)
  #y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
}

plot_background_ventral <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(windowRect = c(0, 0, 200, 200)) #resize for frontal view
  plot3d(bounding_dots, add=T, col="white") 
  plot3d(yolk, add=T, alpha=0.1, col="#E2E2E2") 
  plot3d(acicula, soma=F, lwd=2, add=T, alpha=1, col="grey70")
  par3d(zoom=0.6)
  nview3d('ventral', extramat=(rotationMatrix(0.35, 1, 0, 0)%*%rotationMatrix(0.05, 0, 0, 1)))
  #y-axis clip
  clipplanes3d(1, 0, 0.16, 7000)
  #x-axis clip
  clipplanes3d(0, -1, 0.16, 130000)
}


#load anatomy landmarks
bounding_dots = nlapply(read.neurons.catmaid("^bounding_dots$", pid=11),
                          function(x) smooth_neuron(x, sigma=6000))
yolk <- catmaid_get_volume(4, rval = c("mesh3d", "catmaidmesh", "raw"),
                             invertFaces = T, conn = NULL, pid = 11)
acicula <-   nlapply(read.neurons.catmaid("^acicula$", pid=11),
                       function(x) smooth_neuron(x, sigma=6000))


#function to plot rgl widget with neuron
plot_cell_rgl <- function(annotation, type_of_cell, background){
  ID <- celltypes |>
    filter(`Sensory/inter/motor neuron` == type_of_cell) |> 
    filter(`CATMAID annotation` == annotation) |> 
    pull(`CATMAID annotation`)

  neuron1 = nlapply(read.neurons.catmaid(ID, pid=11),
                    function(x) smooth_neuron(x, sigma=6000))
  if(background == "ventral"){
    plot_background_ventral()
  }
  else{
    plot_background()
  }
  plot3d(neuron1, WithConnectors = F, soma=TRUE, lwd=3,
         add=T, alpha=0.7, col="#0072B2")

  s <- scene3d()
  rglwidget(s, width = 400, height = 100)
}



### Photoreceptors

#### Rhabdomeric photoreceptors



# PRC
annot <- "celltype1"


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







# eyespot-PRCR3
annot <- "celltype33"


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






# eyespot-PRCR1
annot <- "celltype34"


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#### Ciliary photoreceptors {#sec-cPRC}




# cPRC
annot <- "celltype5"


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




### Mechanoreceptors

#### CR (collar receptor) hydrodynamic vibration-receptor neurons




# hCR
annot <- "celltype100"


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




# ventralpygCR
annot <- "celltype101"


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")

widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







# doCRunp
annot <- "celltype102"


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#### PU (penetrating uniciliated) mechanosensory neurons




# PU neurons celltype88-98
annot <- ("celltype88")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# PU neurons celltype88-98
annot <- ("celltype89")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






# PU neurons celltype88-98
#VentraltrunkPUunp
annot <- ("celltype90")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







# PU neurons celltype88-98
annot <- ("celltype91")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






# PU neurons celltype88-98
#hPUc1
annot <- ("celltype92")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# PU neurons celltype88-98
annot <- ("celltype93")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







# PU neurons celltype88-98
annot <- ("celltype94")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# PU neurons celltype88-98
annot <- ("celltype95")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# PU neurons celltype88-98
annot <- ("celltype96")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






# PU neurons celltype88-98
annot <- ("celltype97")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







# PU neurons celltype88-98
annot <- ("celltype98")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#### PB (penetrating biciliated) mechanosensory neurons  {#sec-PB}





# PB neurons celltype80 and 54
annot <- ("celltype80")

rgl <- plot_cell_rgl(annot, "Sensory neuron Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# PB neurons celltype80 and 54
#pygPBunp
annot <- ("celltype54")
rgl <- plot_cell_rgl(annot, "Sensory neuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#### Flow-sensory MS neurons





# MS neurons celltype35-38
annot <- ("celltype35")
rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# MS neurons celltype35-38
annot <- ("celltype36")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# MS neurons celltype35-38
annot <- ("celltype37")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








# MS neurons celltype35-38
annot <- ("celltype38")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#### Chaetal mechanosensory neurons





# chaeMech
annot <- ("celltype71")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#### Other putative mechanosensory neurons in the mechanosensory girdle





#SNantlerPDF
annot <- ("celltype26")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#SNblunt
annot <- ("celltype148")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#SNbronto
annot <- ("celltype168")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#SNFVa
annot <- ("celltype143")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SNPDF-pyg
annot <- ("celltype115")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SNpygM
annot <- ("celltype170")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






### Nuchal organ chemoreceptors




# PU neurons celltype13-15
annot <- ("celltype13")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)










# PU neurons celltype13-15
annot <- ("celltype14")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)









# PU neurons celltype13-15
annot <- ("celltype15")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






### Mushroom body sensory cells






#SNtrpa
annot <- ("celltype190")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SNtorii
annot <- ("celltype187")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#SNmus
annot <- ("celltype25")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#	SNMB2
annot <- ("celltype32")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)










#SNlasso
annot <- ("celltype20")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)









#SNhorn
annot <- ("celltype17")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#SNhook
annot <- ("celltype18")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SNgolden
annot <- ("celltype16")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SN-NS27
annot <- ("celltype137")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)










#SN-NS19
annot <- ("celltype133")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SN-NS18
annot <- ("celltype138")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)









#SN-NS17
annot <- ("celltype131")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)









#SN-NS1
annot <- ("celltype110")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





### Apical organ sensory cells




#	SNadNS22
annot <- ("celltype141")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN47Ach
annot <- ("celltype52")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	SN-YF5cil
annot <- ("celltype132")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#	SN47Ach
annot <- ("celltype52")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-WLD
annot <- ("celltype48")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#	SN-NS6
annot <- ("celltype135")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	SN-NS5
annot <- ("celltype130")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-NS4
annot <- ("celltype112")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-NS3
annot <- ("celltype111")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	SN-NS29
annot <- ("celltype136")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#	SN-NS22
annot <- ("celltype129")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	SN-NS20
annot <- ("celltype134")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	SN-NS16
annot <- ("celltype139")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	SN-NS15
annot <- ("celltype113")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-MIP4
annot <- ("celltype45")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#	SN-MIP1
annot <- ("celltype46")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)










#	SN-IRP2-FMRF
annot <- ("celltype44")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	SN-IRP2-burs
annot <- ("celltype43")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-ASTC
annot <- ("celltype39")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






### Dorsolateral sensory organ sensory cells




#	SN0DLSO1.2-4
annot <- ("celltype28")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)









#	SN-DLSO3-PDF
annot <- ("celltype40")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#	SN-DLSO3
annot <- ("celltype41")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-DLSO2
annot <- ("celltype31")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-DLSO1.3
annot <- ("celltype163")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-DLSO1.2-3
annot <- ("celltype29")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	SN-DLSO1.2-1
annot <- ("celltype30")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-DLSO1.1NP
annot <- ("celltype172")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	SN-DLSO1.1
annot <- ("celltype173")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	SN-DLSO1.0
annot <- ("celltype103")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






### Other sensory cells (unclear modality)




#SNMIP-vc
annot <- ("celltype24")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SNPDF-dc
annot <- ("celltype27")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)









#	SN-DSO
annot <- ("celltype49")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	SNbicil
annot <- ("celltype50")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SNasym
annot <- ("celltype51")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#SNaanc
annot <- ("celltype114")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#SNstiff
annot <- ("celltype167")


rgl <- plot_cell_rgl(annot, "Sensory neuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







## Interneurons

### Visual system interneurons




#IN1
annot <- ("celltype2")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INton
annot <- ("celltype3")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#INR1
annot <- ("celltype171")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Ciliomotor interneurons




#cMNATO
annot <- ("celltype11")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#cMNdc
annot <- ("celltype12")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#cMNPDF-vcl1
annot <- ("celltype10")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





### Apical organ interneurons





#INNOS
annot <- ("celltype7")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)








#INRGWa
annot <- ("celltype6")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()

path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






### Central brain interneurons




#INexsn
annot <- ("celltype188")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INface
annot <- ("celltype193")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INmidL
annot <- ("celltype118")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INmus
annot <- ("celltype174")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INpara
annot <- ("celltype196")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INsn
annot <- ("celltype57")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INsqNSasym
annot <- ("celltype179")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




### Central brain decussating interneurons




#INarc1
annot <- ("celltype55")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INarc2
annot <- ("celltype56")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INcross
annot <- ("celltype123")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INdecussfoot
annot <- ("celltype152")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INdecusshook
annot <- ("celltype153")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INdecussPre
annot <- ("celltype22")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INpro
annot <- ("celltype124")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INproT2
annot <- ("celltype125")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#INW
annot <- ("celltype117")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Central brain descending interneurons



#INbiax
annot <- ("celltype76")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)


#INdescLuqinPDF
annot <- ("celltype175")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INUturn
annot <- ("celltype144")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Other brain decussating interneurons




#INdecussM
annot <- ("celltype199")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)


### Mushroom body intrinsic interneurons




#INMBtype6
annot <- ("celltype120")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INMBtype5
annot <- ("celltype121")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INMBtype7
annot <- ("celltype122")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INMBtype9
annot <- ("celltype197")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INMBtype2
annot <- ("celltype198")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Mushroom body projection neurons




#MBmouth
annot <- ("celltype189")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INMBPDF
annot <- ("celltype184")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INbigloop
annot <- ("celltype140")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INhorn
annot <- ("celltype183")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INMBdescFMRF
annot <- ("celltype185")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INtorii
annot <- ("celltype191")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INUturnMB
annot <- ("celltype192")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INMBdesc3
annot <- ("celltype194")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INMBdesc2
annot <- ("celltype195")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INrope
annot <- ("celltype58")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Head dorsolateral interneuron 



#INpreSer
annot <- ("celltype4")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INlasso
annot <- ("celltype21")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#INpreMN
annot <- ("celltype23")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INdc
annot <- ("celltype104")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INfoot
annot <- ("celltype116")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INhook
annot <- ("celltype119")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INpear
annot <- ("celltype126")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INZ
annot <- ("celltype127")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	INlasso-postSN
annot <- ("celltype160")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INDLSO
annot <- ("celltype186")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#INmantis
annot <- ("celltype201")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#MC2biax-1
annot <- ("celltype42")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INCrossbow
annot <- ("celltype53")


rgl <- plot_cell_rgl(annot, "Interneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




### Head, mechanosensory girdle interneurons




#INsplitPUh
annot <- ("celltype83")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INsplitSpin
annot <- ("celltype202")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INMC3
annot <- ("celltype75")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Ventral nerve cord, mechanosensory girdle interneurons



#INCM
annot <- ("celltype60")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INsplitCR
annot <- ("celltype73")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INsplitCRATO
annot <- ("celltype74")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#INsplitPB-RF/Ya
annot <- ("celltype77")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#INsplitPBant
annot <- ("celltype78")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INsplitPB
annot <- ("celltype79")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INsplitBronto
annot <- ("celltype149")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#INsplitVent
annot <- ("celltype157")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Pygidial interneurons



#INasc-pyg
annot <- ("celltype147")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INFVa-pyg 
annot <- ("celltype72")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INATOpyg
annot <- ("celltype145")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





### Peripheral interneurons




#INbackcross
annot <- ("celltype105")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





### Other trunk interneurons




#INchaeMech
annot <- ("celltype70")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INbiaxLeg
annot <- ("celltype108")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INbiaxH
annot <- ("celltype109")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#INcomm-Upstairs
annot <- ("celltype146")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INleucoPU
annot <- ("celltype150")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INpreLadder
annot <- ("celltype154")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INcomm-Upcross
annot <- ("celltype155")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INcomm-DownL
annot <- ("celltype158")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INcommascFV
annot <- ("celltype169")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#INcosplit
annot <- ("celltype176")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INcommdescFVa
annot <- ("celltype177")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#INbiaxHmid
annot <- ("celltype180")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#INiaStiff
annot <- ("celltype200")


rgl <- plot_cell_rgl(annot, "Interneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




## Motor neurons

### Ciliomotor neurons



#MC
annot <- ("celltype9")


rgl <- plot_cell_rgl(annot, "Motoneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#MNakro
annot <- ("celltype164")


rgl <- plot_cell_rgl(annot, "Motoneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#	Ser-h1
annot <- ("celltype8")


rgl <- plot_cell_rgl(annot, "Motoneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	MNant
annot <- ("celltype19")


rgl <- plot_cell_rgl(annot, "Motoneuron", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	Ser-tr1
annot <- ("celltype142")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#Loop
annot <- ("celltype59")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	MN3
annot <- ("celltype84")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Mixed ciliomotor-musclemotor neurons



#	MN1
annot <- ("celltype85")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	MN2
annot <- ("celltype86")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)


### Other head motoneurons



#	MNheadV
annot <- ("celltype178")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Glandmotor neurons




#	MNgland-head
annot <- ("celltype166")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	MNspinning
annot <- ("celltype47")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



### Pigment-cell-motor neurons




#	cioMNcover
annot <- ("celltype82")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	MC3cover
annot <- ("celltype87")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







### Musclemotor neurons





#	MNspider-ant
annot <- ("celltype61")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#	MNspider-post
annot <- ("celltype62")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#	MNbiramous
annot <- ("celltype63")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	MNcrab
annot <- ("celltype65")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)

#	MNhose
annot <- ("celltype66")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	MNbow
annot <- ("celltype67")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	MNwave
annot <- ("celltype68")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	MNacicX
annot <- ("celltype69")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#	MNcommUpL
annot <- ("celltype106")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	MNob-ipsi
annot <- ("celltype107")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#MNpostv
annot <- ("celltype128")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#MNladder
annot <- ("celltype151")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#MNacic
annot <- ("celltype156")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#MNantacic
annot <- ("celltype161")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#MNpostacic
annot <- ("celltype162")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#MNantelope
annot <- ("celltype181")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#MNsmile
annot <- ("celltype182")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#MNarm
annot <- ("celltype159")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)


#MNchae
annot <- ("celltype165")


rgl <- plot_cell_rgl(annot, "Motoneuron", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




# Non-neuronal cell types

## Gland cells



#spinGland
annot <- ("celltype_non_neuronal7")


rgl <- plot_cell_rgl(annot, "gland cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#ciliatedGland
annot <- ("celltype_non_neuronal9")


rgl <- plot_cell_rgl(annot, "gland cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#MVGland
annot <- ("celltype_non_neuronal17")


rgl <- plot_cell_rgl(annot, "gland cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#HeadGland
annot <- ("celltype_non_neuronal29")


rgl <- plot_cell_rgl(annot, "gland cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#InterparaGland
annot <- ("celltype_non_neuronal30")


rgl <- plot_cell_rgl(annot, "gland cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#spinMicroGland
annot <- ("celltype_non_neuronal31")


rgl <- plot_cell_rgl(annot, "gland cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



## Glia cells




#flat-glia
annot <- ("celltype_non_neuronal15")


rgl <- plot_cell_rgl(annot, "glia cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#ECmidline
annot <- ("celltype_non_neuronal16")


rgl <- plot_cell_rgl(annot, "glia cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#Glia-pigmented
annot <- ("celltype_non_neuronal34")


rgl <- plot_cell_rgl(annot, "glia cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)


## Multiciliated cells



#nuchalCiliaCilia
annot <- ("celltype_non_neuronal4")

rgl <- plot_cell_rgl(annot, "multiciliated cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#akrotroch
annot <- ("celltype_non_neuronal1")
rgl <- plot_cell_rgl(annot, "multiciliated cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#crescentcell
annot <- ("celltype_non_neuronal2")
rgl <- plot_cell_rgl(annot, "multiciliated cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#prototroch
annot <- ("celltype_non_neuronal3")
rgl <- plot_cell_rgl(annot, "multiciliated cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#metatroch
annot <- ("celltype_non_neuronal5")
rgl <- plot_cell_rgl(annot, "multiciliated cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#paratroch
annot <- ("celltype_non_neuronal6")
rgl <- plot_cell_rgl(annot, "multiciliated cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




## Pigment cells




#covercell
annot <- ("celltype_non_neuronal8")


rgl <- plot_cell_rgl(annot, "pigment cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#eyespot-pigment-cell
annot <- ("celltype_non_neuronal10")


rgl <- plot_cell_rgl(annot, "pigment cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	eye-pigment-cell
annot <- ("celltype_non_neuronal11")


rgl <- plot_cell_rgl(annot, "pigment cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	yolk cover cell 
annot <- ("celltype_non_neuronal14")


rgl <- plot_cell_rgl(annot, "pigment cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	CB-pigment
annot <- ("celltype_non_neuronal32")


rgl <- plot_cell_rgl(annot, "pigment cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	vacuolar-cell-head
annot <- ("celltype_non_neuronal33")


rgl <- plot_cell_rgl(annot, "pigment cell", "anterior")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#	pygidial-pigment-cell
annot <- ("celltype_non_neuronal35")


rgl <- plot_cell_rgl(annot, "pigment cell", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)


## Chaetal and acicular cells




#	chaeta
annot <- ("celltype_non_neuronal22")


rgl <- plot_cell_rgl(annot, "chaetal complex", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#	acicula
annot <- ("celltype_non_neuronal23")


rgl <- plot_cell_rgl(annot, "chaetal complex", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#	acFC-noER
annot <- ("celltype_non_neuronal24")


rgl <- plot_cell_rgl(annot, "chaetal complex", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#	chaeFC-hemi
annot <- ("celltype_non_neuronal25")


rgl <- plot_cell_rgl(annot, "chaetal complex", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#	chaeFC-ER
annot <- ("celltype_non_neuronal26")


rgl <- plot_cell_rgl(annot, "chaetal complex", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#	chaeFC-noER
annot <- ("celltype_non_neuronal27")


rgl <- plot_cell_rgl(annot, "chaetal complex", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#	chaeFC-EC
annot <- ("celltype_non_neuronal28")


rgl <- plot_cell_rgl(annot, "chaetal complex", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




## Nephridial cells




#protonephridium
annot <- ("celltype_non_neuronal19")


rgl <- plot_cell_rgl(annot, "excretory system", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#nephridium
annot <- ("celltype_non_neuronal20")


rgl <- plot_cell_rgl(annot, "excretory system", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#nephridiumTip
annot <- ("celltype_non_neuronal21")


rgl <- plot_cell_rgl(annot, "excretory system", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



## Muscle cells



#muscle

for (i in 37:89){
annot <- paste("celltype_non_neuronal", i, sep = "")
rgl <- plot_cell_rgl(annot, "muscle", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)
}

## Other non-neuronal cell types






#coelothelium
annot <- ("celltype_non_neuronal36")


rgl <- plot_cell_rgl(annot, "coelothelium", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)





#yolk
annot <- ("celltype_non_neuronal91")


rgl <- plot_cell_rgl(annot, "yolk", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)






#meso
annot <- ("celltype_non_neuronal92")


rgl <- plot_cell_rgl(annot, "mesoderm", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#macrophage like
annot <- ("celltype_non_neuronal13")


rgl <- plot_cell_rgl(annot, "macrophage-like", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)



#microvillarCell
annot <- ("celltype_non_neuronal18")


rgl <- plot_cell_rgl(annot, "microvillar", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




#gland-pigment
annot <- ("celltype_non_neuronal12")


rgl <- plot_cell_rgl(annot, "gland", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)







#EC
annot <- ("celltype_non_neuronal90")


rgl <- plot_cell_rgl(annot, "epidermis", "ventral")
widgets <- rglwidget(scene3d(rgl))
# make snapshot
path1 <- paste("supplements/celltype_compendium_website/_site/snapshots/", annot, ".png", sep = "")
rgl.snapshot(path1)
close3d()
path2 <- paste("supplements/celltype_compendium_website/_site/celltype_RGLs/", annot, ".html", sep = "")
htmlwidgets::saveWidget(widgets, path2)




